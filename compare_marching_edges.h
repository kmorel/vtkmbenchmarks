
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

namespace marchingEdges
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


namespace me
{
namespace worklets
{

typedef VTKM_DEFAULT_DEVICE_ADAPTER_TAG DeviceAdapter;

template <typename FieldType>
struct MarchingEdgeContData
{
  vtkm::cont::ArrayHandle< FieldType > scalarHandle;
  vtkm::cont::ArrayHandle< vtkm::UInt8 > xEdgeCaseHandle;

  vtkm::cont::ArrayHandle< vtkm::Id > xEdgeSumHandle;
  vtkm::cont::ArrayHandle< vtkm::Id > numTrisPerRowHandle;

  float isovalue;
  int dims[3];
  int inc1;
  int inc2;
  int sliceOffset;
  vtkm::Vec<vtkm::Float32,3> Origin;
  vtkm::Vec<vtkm::Float32,3> Spacing;

  VTKM_CONT_EXPORT
  MarchingEdgeContData( const vtkm::Vec<vtkm::Float32,3> origin,
                      const vtkm::Vec<vtkm::Float32,3> spacing,
                      int p_dims[3],
                      vtkm::cont::ArrayHandle< FieldType > scalarInHandle,
                      float iso):
    Origin(origin),
    Spacing(spacing)
   {
    this->scalarHandle = scalarInHandle;
    this->isovalue = iso;

    this->dims[0] = p_dims[0];
    this->dims[1] = p_dims[1];
    this->dims[2] = p_dims[2];

    this->inc1 = p_dims[0];
    this->inc2 = p_dims[0] * p_dims[1];
    this->sliceOffset = (p_dims[0]-1) * p_dims[1];
   }
};

//----------------------------------------------------------------------------
template <typename FieldType>
struct MarchingEdgePass1ExecData
{
  typedef typename ::worklets::PortalTypes<FieldType>::PortalConst FieldPortalType;
  typedef typename ::worklets::PortalTypes<vtkm::UInt8>::Portal UInt8PortalType;
  typedef typename ::worklets::PortalTypes<vtkm::Id>::Portal IdPortalType;

  float isovalue;
  vtkm::Id3 dims;
  int inc1;
  int inc2;
  int sliceOffset;
  FieldPortalType scalarData;
  UInt8PortalType edgeXCases;

  IdPortalType xsum;

  explicit
  MarchingEdgePass1ExecData( MarchingEdgeContData<FieldType>& metaData):
    isovalue(metaData.isovalue),
    dims(metaData.dims[0], metaData.dims[1], metaData.dims[2]),
    inc1(metaData.inc1),
    inc2(metaData.inc2),
    sliceOffset(metaData.sliceOffset),
    scalarData(metaData.scalarHandle.PrepareForInput(DeviceAdapter())),
    edgeXCases(metaData.xEdgeCaseHandle.PrepareForOutput( (metaData.dims[0]-1) * metaData.dims[1] * metaData.dims[2], DeviceAdapter() )),
    xsum(metaData.xEdgeSumHandle.PrepareForOutput(metaData.dims[1] * metaData.dims[2], DeviceAdapter()) )
  {
    std::cout << "MarchingEdgeExecData" << std::endl;
    std::cout << "algo.Dims[0]: " << dims[0] << std::endl;
    std::cout << "algo.Dims[1]: " << dims[1] << std::endl;
    std::cout << "algo.Dims[2]: " << dims[2] << std::endl;
    std::cout << "algo.SliceOffset: " << sliceOffset << std::endl;;

    std::cout << "algo.Inc[0]: " << 1 << std::endl;
    std::cout << "algo.Inc[1]: " << inc1 << std::endl;
    std::cout << "algo.Inc[2]: " << inc2 << std::endl;

  }
};

//----------------------------------------------------------------------------
template <typename FieldType>
struct MarchingEdgePass2ExecData
{
  typedef typename ::worklets::PortalTypes<FieldType>::PortalConst FieldPortalType;
  typedef typename ::worklets::PortalTypes<vtkm::UInt8>::PortalConst UInt8ConstPortalType;
  typedef typename ::worklets::PortalTypes<vtkm::Id>::Portal IdPortalType;
  typedef typename ::worklets::PortalTypes<vtkm::Id>::PortalConst IdConstPortalType;

  float isovalue;
  vtkm::Id3 dims;
  int inc1;
  int inc2;
  int sliceOffset;

  UInt8ConstPortalType edgeXCases;
  IdConstPortalType xsum;
  IdPortalType numTris;

  explicit
  MarchingEdgePass2ExecData( MarchingEdgeContData<FieldType>& metaData):
    isovalue(metaData.isovalue),
    dims(metaData.dims[0], metaData.dims[1], metaData.dims[2]),
    inc1(metaData.inc1),
    inc2(metaData.inc2),
    sliceOffset(metaData.sliceOffset),
    edgeXCases(metaData.xEdgeCaseHandle.PrepareForInput( DeviceAdapter() )),
    xsum(metaData.xEdgeSumHandle.PrepareForInput( DeviceAdapter()) ),
    numTris(metaData.numTrisPerRowHandle.PrepareForOutput( (metaData.dims[1]-1) * (metaData.dims[2]-1), DeviceAdapter()) )
  {

  }

  VTKM_EXEC_EXPORT
  vtkm::UInt8 numberOfPrimitives( vtkm::UInt8 ecase ) const
    { return EdgeCasesTable[ecase][0]; }

  VTKM_EXEC_EXPORT
  vtkm::UInt8 edgeUses( vtkm::UInt8 ecase, vtkm::UInt8 index ) const
    { return EdgeUsesTable[ecase][index] ; }

};


//----------------------------------------------------------------------------
template <typename FieldType>
class ProcessXEdges : public vtkm::worklet::WorkletMapField
{
public:
  typedef void ControlSignature(FieldIn<IdType>  id);
  typedef void ExecutionSignature(_1);
  typedef _1 InputDomain;

  MarchingEdgePass1ExecData<FieldType> d;

  VTKM_CONT_EXPORT
  ProcessXEdges( MarchingEdgePass1ExecData<FieldType>& data) :
         d( data )
   {
   }

  //----------------------------------------------------------------------------
  // PASS 1: Process a single volume x-row (and all of the voxel edges that
  // compose the row). Determine the x-edges case classification, count the
  // number of x-edge intersections, and figure out where intersections along
  // the x-row begins and ends (i.e., gather information for computational
  // trimming).
  VTKM_EXEC_EXPORT
  void operator()(const vtkm::Id&  yzid) const
  {
  //now instead of doing every 1 dim, we need to calculate out which
  //row and slice we are.
  const vtkm::Id row = yzid % this->d.dims[1];
  const vtkm::Id slice = yzid / this->d.dims[1];
  const vtkm::Id offset = slice * this->d.inc2 + row * this->d.inc1;

  // std::cout << offset << ", " << row << ", " << slice << std::endl;

  this->process(yzid, offset, row, slice);
  }

  VTKM_EXEC_EXPORT
  void process(vtkm::Id yzid, vtkm::Id offset, vtkm::Id row, vtkm::Id slice) const
  {
  const vtkm::Id nxcells = this->d.dims[0]-1;
  const vtkm::Id metaWritePos = yzid;  // no need to compute this
  const vtkm::Id writePos = slice * this->d.sliceOffset + row;

  FieldType s0=0;
  FieldType s1= this->d.scalarData.Get(offset);

  vtkm::Id tempXSum = 0;
  //for all cells in this row
  for (int i=0; i < nxcells; ++i)
    {
    s0 = s1;
    s1 = this->d.scalarData.Get(offset + i + 1);

    vtkm::UInt8 edgeCase = marchingEdges::Below;
    if (s0 >= this->d.isovalue)
      {
      edgeCase = marchingEdges::LeftAbove;
      }
    if( s1 >= this->d.isovalue)
      {
      edgeCase |= marchingEdges::RightAbove;
      }

    //write back our edge case, this actually needs to be per cell
    //not per row
    this->d.edgeXCases.Set(writePos + (i * nxcells) ,edgeCase);
    if ( edgeCase == marchingEdges::LeftAbove ||
         edgeCase == marchingEdges::RightAbove )
      {
      ++tempXSum; //increment number of intersections along x-edge
      }//if contour interacts with this x-edge
    }

  this->d.xsum.Set(metaWritePos,tempXSum); //write back the number of intersections along x-edge
  }
};


template <typename FieldType>
class ProcessYZEdges : public vtkm::worklet::WorkletMapField
{
public:
  typedef void ControlSignature(FieldIn<IdType>  id);
  typedef void ExecutionSignature(_1);
  typedef _1 InputDomain;

  MarchingEdgePass2ExecData<FieldType> d;

  VTKM_CONT_EXPORT
  ProcessYZEdges( MarchingEdgePass2ExecData<FieldType>& data) :
         d( data )
   {
   }

  //----------------------------------------------------------------------------
  // PASS 2: Process a single x-row of voxels. Count the number of y- and
  // z-intersections by topological reasoning from x-edge cases. Determine the
  // number of primitives (i.e., triangles) generated from this row. Use
  // computational trimming to reduce work. Note *ePtr[4] is four pointers to
  // four x-edge rows that bound the voxel x-row and which contain edge case
  // information.
  VTKM_EXEC_EXPORT
  void operator()(const vtkm::Id&  yzid) const
  {
  //now instead of doing every 1 dim, we need to calculate out which
  //row and slice we are.
  const vtkm::Id ycdim = (this->d.dims[1]-1);
  const vtkm::Id row = yzid % ycdim;
  const vtkm::Id slice = yzid / ycdim;
  this->process(yzid, row, slice);
  }

  VTKM_EXEC_EXPORT
  void process(vtkm::Id yzid, vtkm::Id row, vtkm::Id slice) const
  {
  const vtkm::Id nxcells = this->d.dims[0]-1;
  const vtkm::Id metaWritePos = yzid; // no need to compute this
  const vtkm::Id readPos = slice * this->d.sliceOffset + row;

  // Okay run along the x-voxels and count the number of y- and
  // z-intersections. Here we are just checking y,z edges that make up the
  // voxel axes. Also check the number of primitives generated.
  vtkm::UInt8 edgeCase[4] = {0, 0, 0, 0};

  vtkm::Id tempTriSum = 0;
  for (int i=0; i < nxcells; ++i)
    {
    edgeCase[0] = this->d.edgeXCases.Get( (i * nxcells) + readPos);
    edgeCase[1] = this->d.edgeXCases.Get( (i * nxcells) + readPos + 1 );
    edgeCase[2] = this->d.edgeXCases.Get( (i * nxcells) + readPos + this->d.sliceOffset);
    edgeCase[3] = this->d.edgeXCases.Get( (i * nxcells) + readPos + this->d.sliceOffset + 1 );

    const vtkm::UInt8 eCase = (edgeCase[0]    |
                               edgeCase[1]<<2 |
                               edgeCase[2]<<4 |
                               edgeCase[3]<<6 );

    //using the eCase for this voxel we need to lookup the number of primitives
    //that this voxel will generate.
    const vtkm::UInt8 numTris = this->d.numberOfPrimitives(eCase);
    tempTriSum += numTris;
    }//for all voxels along this x-edge

    this->d.numTris.Set(metaWritePos,tempTriSum); //write back the number of triangles
  }
};

} //worklets namespace

//now that the device adapter is included set a global typedef
//that is the chosen device tag
typedef VTKM_DEFAULT_DEVICE_ADAPTER_TAG DeviceAdapter;

static void doMarchingEdges( int pdims[3], //point dims
                           int cdims[3], //cell dims
                           const vtkm::cont::ArrayHandle<vtkm::Float32>& field,
                           vtkm::cont::ArrayHandle<vtkm::Float32>& scalarsArray,
                           vtkm::cont::ArrayHandle< vtkm::Vec<vtkm::Float32,3> >& verticesArray
                           )
{

  // Pass 0: Allocate the xEdgeCase. This is done here to make it easier
  //to pass to worklets
  vtkm::Vec<vtkm::Float32,3> origin(0,0,0);
  vtkm::Vec<vtkm::Float32,3> spacing(1,1,1);
  worklets::MarchingEdgeContData<vtkm::Float32> contData(origin,
                                                       spacing,
                                                       pdims,
                                                       field,
                                                       ISO_VALUE);

  // PASS 1: Traverse all x-rows building edge cases and counting number of
  // intersections (i.e., accumulate information necessary for later output
  // memory allocation, e.g., the number of output points along the x-rows
  // are counted).
  {
  worklets::MarchingEdgePass1ExecData<vtkm::Float32> execData( contData );
  //first pass is over edges so use point dims
  vtkm::cont::ArrayHandleCounting<vtkm::Id> passIds(0, pdims[1] * pdims[2]);

  typedef worklets::ProcessXEdges< vtkm::Float32 > ProcessXEdgeFunctor;
  typedef vtkm::worklet::DispatcherMapField< ProcessXEdgeFunctor > XEdgeDispatcher;

  ProcessXEdgeFunctor xEdge( execData );
  XEdgeDispatcher xEDispatcher(xEdge);
  xEDispatcher.Invoke(passIds);
  }

  // PASS 2: Traverse all voxel x-rows and process voxel y&z edges.  The
  // result is a count of the number of y- and z-intersections, as well as
  // the number of triangles generated along these voxel rows.
  {
  worklets::MarchingEdgePass2ExecData<vtkm::Float32> execData( contData );
  //first pass is over 'cells' so use cell dims
  vtkm::cont::ArrayHandleCounting<vtkm::Id> passIds(0, cdims[1] * cdims[2]);

  typedef worklets::ProcessYZEdges< vtkm::Float32 > ProcessYZEdgeFunctor;
  typedef vtkm::worklet::DispatcherMapField< ProcessYZEdgeFunctor > YZEdgeDispatcher;

  ProcessYZEdgeFunctor yzFunctor( execData );
  YZEdgeDispatcher yzEDispatcher(yzFunctor);
  yzEDispatcher.Invoke(passIds);
  }

  // PASS 3: Now allocate and generate output. We have the number of triangles
  // each section is going to write so we can map this to MarchingCubes now
  vtkm::cont::ArrayHandle< vtkm::Id > validCellIndicesArray;
  vtkm::cont::ArrayHandle< vtkm::Id > inputCellIterationNumber;
  vtkm::cont::ArrayHandle< vtkm::Id > numOutputTrisPerCell;
  const vtkm::Id numOutputTriangles =
      vtkm::cont::DeviceAdapterAlgorithm<DeviceAdapter>::ScanInclusive(contData.numTrisPerRowHandle,
                                                                       numOutputTrisPerCell);
  //terminate if no cells has triangles left
  std::cout << "numOutputTriangles: " << numOutputTriangles << std::endl;
  if(numOutputTriangles == 0)
    { return; }

  // vtkm::cont::ArrayHandleCounting<vtkm::Id> validCellCountImplicitArray(0, numOutputTriangles);
  // vtkm::cont::DeviceAdapterAlgorithm<DeviceAdapter>::UpperBounds(numOutputTrisPerCell,
  //                                                                validCellCountImplicitArray,
  //                                                                validCellIndicesArray);

  // numOutputTrisPerCell.ReleaseResourcesExecution();

  // //compute for each output triangle what iteration of the input cell
  // //generates it
  // vtkm::cont::DeviceAdapterAlgorithm<DeviceAdapter>::LowerBounds(validCellIndicesArray,
  //                                                                validCellIndicesArray,
  //                                                                inputCellIterationNumber);

  verticesArray.PrepareForOutput(numOutputTriangles*3, DeviceAdapter());
  scalarsArray.PrepareForOutput(numOutputTriangles*3, DeviceAdapter());

};

static void RunMarchingEdges(int pdims[3],
                           std::vector<vtkm::Float32>& buffer,
                           std::string device,
                           int MAX_NUM_TRIALS,
                           bool silent=false)
{
  int dims[3] = { pdims[0]-1, pdims[1]-1, pdims[2]-1 };

  vtkm::cont::ArrayHandle<vtkm::Float32> field;
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
    doMarchingEdges( pdims, dims, field, scalarsArray, verticesArray);

    double time = timer.GetElapsedTime();
    if(!silent)
      {
      std::cout << "num cells: " << (scalarsArray.GetNumberOfValues()/3)  << std::endl;
      std::cout << "vtkm," << device << "," << time << "," << trial << std::endl;
      }
    }

}

}
