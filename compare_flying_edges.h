
#include <vtkm/cont/DeviceAdapter.h>

#include <vtkm/cont/ArrayHandle.h>
#include <vtkm/cont/ArrayHandleCounting.h>
#include <vtkm/cont/DynamicArrayHandle.h>
#include <vtkm/cont/Timer.h>
#include <vtkm/Pair.h>
#include <vtkm/worklet/DispatcherMapField.h>
#include <vtkm/worklet/WorkletMapField.h>

#include <vtkContourFilter.h>
#include <vtkImageData.h>
#include <vtkNew.h>
#include <vtkSmartPointer.h>
#include <vtkTrivialProducer.h>

#include <vector>

namespace flyingedges
{
  static const vtkm::UInt8 Below = 0; //below isovalue
  static const vtkm::UInt8 Above = 1; //above isovalue
  static const vtkm::UInt8 LeftAbove = 1; //left vertex is above isovalue
  static const vtkm::UInt8 RightAbove = 2; //right vertex is above isovalue
  static const vtkm::UInt8 BothAbove = 3; //entire edge is above isovalue

  // Dealing with boundary situations when processing volumes.
  static const vtkm::UInt8 Interior = 0;
  static const vtkm::UInt8 MinBoundary = 1;
  static const vtkm::UInt8 MaxBoundary = 2;
}

namespace worklets
{

template <typename FieldType>
class ProcessXEdges : public vtkm::worklet::WorkletMapField
{
public:
  typedef void ControlSignature(FieldIn<IdType>  xRowId,
                                FieldOut<IdType> xsum,
                                FieldOut<IdType> xmin,
                                FieldOut<IdType> xmax);
  typedef void ExecutionSignature(_1, _2, _3, _4);
  typedef _1 InputDomain;

  typedef typename PortalTypes<FieldType>::PortalConst FieldPortalType;
  FieldPortalType pointData;

  typedef typename PortalTypes<vtkm::UInt8>::Portal UInt8PortalType;
  UInt8PortalType edgeXCases;

  float isovalue;
  int xdim, ydim, zdim;
  int cellsPerLayer;
  int pointsPerLayer;

  template<typename T, typename U>
  VTKM_CONT_EXPORT
  ProcessXEdges(const T& pointHandle,
                U& xEdgeCasesHandle,
                float iso,
                int dims[3]) :
         pointData( pointHandle.PrepareForInput( DeviceAdapter() ) ),
         edgeXCases( xEdgeCasesHandle.PrepareForOutput( (dims[0]-1) * (dims[1]-1) * (dims[2]-1), DeviceAdapter() ) ),
         isovalue( iso ),
         xdim(dims[0]),
         ydim(dims[1]),
         zdim(dims[2]),
         cellsPerLayer((xdim - 1) * (ydim - 1)),
         pointsPerLayer (xdim*ydim)
   {

   }

   //----------------------------------------------------------------------------
  // PASS 1: Process a single volume x-row (and all of the voxel edges that
  // compose the row). Determine the x-edges case classification, count the
  // number of x-edge intersections, and figure out where intersections along
  // the x-row begins and ends (i.e., gather information for computational
  // trimming).
  VTKM_EXEC_EXPORT
  void operator()(vtkm::Id  xRowId,
                  vtkm::Id& xsum,
                  vtkm::Id& xmin,
                  vtkm::Id& xmax) const
  {
  //I doubt that the first iteration of this is going to perform super well
  //since the kernel does sooo many mem fetches
  //compute the offset into pointData based on our xRowId
  const int x = xRowId % (xdim - 1);
  const int y = (xRowId / (xdim - 1)) % (ydim -1);
  const int z = xRowId / cellsPerLayer;

  //compute the write position for each x cell. This is not linear in memory
  //as that kills performance, instead it is strided based on the size
  //of the grid, so that we get coalesced writes on cuda. The performance
  //difference on 1024^3 grids is around 25 times!
  vtkm::Id writePos = xRowId + ((y * (xdim-1)) * (xRowId/(xdim-1)));
  const vtkm::Id writePosOffset = (xdim-1);


  //compute the first index for first cell in this x row
  const int i0 = x + y*xdim + z * pointsPerLayer;
  FieldType s0;
  FieldType s1= this->pointData.Get(i0);

  vtkm::Id tempSum = 0;
  vtkm::Id tempMin = xdim-1;
  vtkm::Id tempMax = 0;

  //for all cells in this row
  for (int i=0; i < xdim-1; ++i)
    {
    s0 = s1;
    s1 = this->pointData.Get(i0 + i+1);

    vtkm::UInt8 edgeCase = flyingedges::Below;
    if (s0 >= isovalue)
      {
      edgeCase = flyingedges::LeftAbove;
      }
    if( s1 >= isovalue)
      {
      edgeCase |= flyingedges::RightAbove;
      }

    //write back our edge case, this actually needs to be per cell
    //not per row
    this->edgeXCases.Set( writePos, edgeCase);
    writePos += writePosOffset;

    if ( edgeCase == flyingedges::LeftAbove ||
         edgeCase == flyingedges::RightAbove )
      {
      ++tempSum; //increment number of intersections along x-edge
      tempMin = ( i < tempMin ? i : tempMin);
      tempMax = i + 1;
      }//if contour interacts with this x-edge
    }

  xsum = tempSum; //write back the number of intersections along x-edge

  // The beginning and ending of intersections along the edge is used for
  // computational trimming.
  xmin = tempMin; //where intersections start along x edge
  xmax = tempMax; //where intersections end along x edge

  }
};

}

namespace fe {


static void doFlyingEdges( int vdims[3], //point dims
                           const vtkm::cont::ArrayHandle<vtkm::Float32>& field,
                           vtkm::cont::ArrayHandle<vtkm::Float32> scalarsArray,
                           vtkm::cont::ArrayHandle< vtkm::Vec<vtkm::Float32,3> > verticesArray,
                           int dims[3] //cell dims
                           )
{
  vtkm::cont::ArrayHandle< vtkm::UInt8 > xEdgeCase; //size will be dims^3
  vtkm::cont::ArrayHandle< vtkm::Id > xEdgeSum, xMinInt, xMaxInt; //size of dimsY * dimsZ

  // PASS 1: Traverse all x-rows building edge cases and counting number of
  // intersections (i.e., accumulate information necessary for later output
  // memory allocation, e.g., the number of output points along the x-rows
  // are counted).
  vtkm::Id numberOfXRows = dims[1] * dims[2];
  vtkm::cont::ArrayHandleCounting<vtkm::Id> cellCountImplicitArray(0, numberOfXRows);
  typedef worklets::ProcessXEdges< vtkm::Float32 > ProcessXEdgeFunctor;
  typedef vtkm::worklet::DispatcherMapField< ProcessXEdgeFunctor > XEdgeDispatcher;

  ProcessXEdgeFunctor xEdge(field, xEdgeCase, ISO_VALUE, vdims );
  XEdgeDispatcher xEDispatcher(xEdge);
  xEDispatcher.Invoke(cellCountImplicitArray, xEdgeSum, xMinInt, xMaxInt);


  // PASS 2: Traverse all voxel x-rows and process voxel y&z edges.  The
  // result is a count of the number of y- and z-intersections, as well as
  // the number of triangles generated along these voxel rows.

};


static void RunFlyingEdges(int vdims[3],
                           std::vector<vtkm::Float32>& buffer,
                           std::string device,
                           int MAX_NUM_TRIALS,
                           bool silent=false)
{
  int dims[3] = { vdims[0]-1, vdims[1]-1, vdims[2]-1 };

  vtkm::cont::ArrayHandle<vtkm::Float32> field;
  //construct the scheduler that will execute all the worklets
  for(int trial=0; trial < MAX_NUM_TRIALS; ++trial)
    {
    vtkm::cont::Timer<> timer;

    //setup the iso field to contour
    if(trial < MAX_NUM_TRIALS/2)
      {
      field = vtkm::cont::make_ArrayHandle(buffer);
      }

    vtkm::cont::ArrayHandle<vtkm::Float32> scalarsArray;
    vtkm::cont::ArrayHandle< vtkm::Vec<vtkm::Float32,3> > verticesArray;
    doFlyingEdges( vdims, field, scalarsArray, verticesArray, dims);

    double time = timer.GetElapsedTime();
    if(!silent)
      {
      std::cout << "num cells: " << (scalarsArray.GetNumberOfValues()/3)  << std::endl;
      std::cout << "vtkm," << device << "," << time << "," << trial << std::endl;
      }
    }

}

}