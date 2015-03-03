
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

namespace y_stride
{
namespace worklets
{

typedef VTKM_DEFAULT_DEVICE_ADAPTER_TAG DeviceAdapter;

template <typename FieldType>
struct MarchingContData
{
  vtkm::cont::ArrayHandle< FieldType > scalarHandle;
  vtkm::cont::ArrayHandle< vtkm::Id > numTrisPerCell;


  float isovalue;
  int dims[3];
  int inc1;
  int inc2;
  vtkm::Vec<vtkm::Float32,3> Origin;
  vtkm::Vec<vtkm::Float32,3> Spacing;

  VTKM_CONT_EXPORT
  MarchingContData( const vtkm::Vec<vtkm::Float32,3> origin,
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
   }
};

//----------------------------------------------------------------------------
template <typename FieldType>
struct MarchingPass1ExecData
{
  typedef typename ::worklets::PortalTypes<FieldType>::PortalConst FieldPortalType;
  typedef typename ::worklets::PortalTypes<vtkm::Id>::Portal IdPortalType;

  float isovalue;
  vtkm::Id3 dims;
  int inc1;
  int inc2;
  FieldPortalType scalarData;

  IdPortalType numTriangles;

  explicit
  MarchingPass1ExecData( MarchingContData<FieldType>& metaData):
    isovalue(metaData.isovalue),
    dims(metaData.dims[0], metaData.dims[1], metaData.dims[2]),
    inc1(metaData.inc1),
    inc2(metaData.inc2),
    scalarData(metaData.scalarHandle.PrepareForInput(DeviceAdapter())),
    numTriangles(metaData.numTrisPerCell.PrepareForOutput( (metaData.dims[0]-1) * (metaData.dims[1]-1) * (metaData.dims[2]-1), DeviceAdapter() ))
  {
  }
};

//----------------------------------------------------------------------------
template <typename FieldType>
class ProcessCellsByYAxis : public vtkm::worklet::WorkletMapField
{
public:
  typedef void ControlSignature(FieldIn<IdType>  id);
  typedef void ExecutionSignature(_1);
  typedef _1 InputDomain;

  MarchingPass1ExecData<FieldType> d;

  VTKM_CONT_EXPORT
  ProcessCellsByYAxis( MarchingPass1ExecData<FieldType>& data) :
         d( data )
   {
   }

  //----------------------------------------------------------------------------
  // PASS 1: generate the number of triangles per cell in each Y stripe
  VTKM_EXEC_EXPORT
  void operator()(const vtkm::Id&  xzid) const
  {
  const vtkm::Id xcdim = (this->d.dims[0]-1);
  const vtkm::Id ycdim = (this->d.dims[1]-1);
  const vtkm::Id crow = xzid % xcdim;
  const vtkm::Id cslice = xzid / xcdim;
  const vtkm::Id coffset = cslice * (xcdim * ycdim) + crow;

  // const vtkm::Id x = crow;
  // const vtkm::Id y = 0;
  // const vtkm::Id z = cslice;
  const vtkm::Id poffset = (cslice * this->d.inc2) + crow;


  this->process(poffset, coffset);
  }

  VTKM_EXEC_EXPORT
  void process(vtkm::Id poffset, vtkm::Id coffset) const
  {
  const vtkm::Id nycells = this->d.dims[1]-1;
  const vtkm::Id writePos = coffset;

  FieldType f0 = this->d.scalarData.Get( poffset );
  FieldType f1 = this->d.scalarData.Get( poffset + 1 );
  FieldType f2 = this->d.scalarData.Get( poffset + this->d.inc1 + 1 );
  FieldType f3 = this->d.scalarData.Get( poffset + this->d.inc1 );
  FieldType f4 = this->d.scalarData.Get( poffset + this->d.inc2 );
  FieldType f5 = this->d.scalarData.Get( poffset + this->d.inc2 + 1 );
  FieldType f6 = this->d.scalarData.Get( poffset + this->d.inc2 + this->d.inc1 + 1 );
  FieldType f7 = this->d.scalarData.Get( poffset + this->d.inc2 + this->d.inc1);

  unsigned int cubeindex = (f0 > this->d.isovalue);
  cubeindex += (f1 > this->d.isovalue)*2;
  cubeindex += (f2 > this->d.isovalue)*4;
  cubeindex += (f3 > this->d.isovalue)*8;
  cubeindex += (f4 > this->d.isovalue)*16;
  cubeindex += (f5 > this->d.isovalue)*32;
  cubeindex += (f6 > this->d.isovalue)*64;
  cubeindex += (f7 > this->d.isovalue)*128;

  this->d.numTriangles.Set(writePos, numVerticesTable[cubeindex] / 3 );

  //for all cells in this row
  for (int i=1; i < nycells; ++i)
    {
    //update the bottom verts to be the old top verts
    //shape of voxel back face is:
    // 7 6
    // 4 5
    //shape of voxel front face is:
    // 3 2
    // 0 1
    f0 = f3;
    f3 = this->d.scalarData.Get( (this->d.inc1 * i) + poffset + this->d.inc1 );
    f1 = f2;
    f2 = this->d.scalarData.Get( (this->d.inc1 * i) + poffset + this->d.inc1 + 1 );
    f4 = f7;
    f7 = this->d.scalarData.Get( (this->d.inc1 * i) + poffset + this->d.inc1 + this->d.inc2  );
    f5 = f6;
    f6 = this->d.scalarData.Get( (this->d.inc1 * i) + poffset + this->d.inc1 + this->d.inc2 + 1 );

    cubeindex = (f0 > this->d.isovalue);
    cubeindex += (f1 > this->d.isovalue)*2;
    cubeindex += (f2 > this->d.isovalue)*4;
    cubeindex += (f3 > this->d.isovalue)*8;
    cubeindex += (f4 > this->d.isovalue)*16;
    cubeindex += (f5 > this->d.isovalue)*32;
    cubeindex += (f6 > this->d.isovalue)*64;
    cubeindex += (f7 > this->d.isovalue)*128;

    this->d.numTriangles.Set(writePos + ((this->d.dims[0]-1) * i),  numVerticesTable[cubeindex] / 3 );
    }
  }
};

} //worklets namespace

//now that the device adapter is included set a global typedef
//that is the chosen device tag
typedef VTKM_DEFAULT_DEVICE_ADAPTER_TAG DeviceAdapter;

static void doMarchingCubes( int pdims[3], //point dims
                           int cdims[3], //cell dims
                           const vtkm::cont::ArrayHandle<vtkm::Float32>& field,
                           vtkm::cont::ArrayHandle<vtkm::Float32>& scalarsArray,
                           vtkm::cont::ArrayHandle< vtkm::Vec<vtkm::Float32,3> >& verticesArray
                           )
{

  // Pass 0: allocate cont data helper
  vtkm::Vec<vtkm::Float32,3> origin(0,0,0);
  vtkm::Vec<vtkm::Float32,3> spacing(1,1,1);
  worklets::MarchingContData<vtkm::Float32> contData(origin,
                                                     spacing,
                                                     pdims,
                                                     field,
                                                     ISO_VALUE);

  //----------------------------------------------------------------------------
  // PASS 1: Process a single volume y-row  of cells.
  {
  worklets::MarchingPass1ExecData<vtkm::Float32> execData( contData );
  //first pass is over edges so use point dims
  vtkm::cont::ArrayHandleCounting<vtkm::Id> passIds(0, cdims[0] * cdims[2]);

  typedef worklets::ProcessCellsByYAxis< vtkm::Float32 > ProcessCellsYStride;
  typedef vtkm::worklet::DispatcherMapField< ProcessCellsYStride > YDispatcher;

  ProcessCellsYStride yStrideFunctor( execData );
  YDispatcher dispatcher(yStrideFunctor);
  dispatcher.Invoke(passIds);
  }

  //compute the number of valid input cells and those ids
  const vtkm::Id numOutputCells =
    vtkm::cont::DeviceAdapterAlgorithm<DeviceAdapter>::ScanInclusive(contData.numTrisPerCell,
                                                                     contData.numTrisPerCell);

  //terminate if no cells has triangles left
  if(numOutputCells == 0)
    { return; }

  vtkm::cont::ArrayHandle< vtkm::Id > validCellIndicesArray;
  vtkm::cont::ArrayHandle< vtkm::Id > inputCellIterationNumber;

  vtkm::cont::ArrayHandleCounting<vtkm::Id> validCellCountImplicitArray(0, numOutputCells);
  vtkm::cont::DeviceAdapterAlgorithm<DeviceAdapter>::UpperBounds(contData.numTrisPerCell,
                                                                 validCellCountImplicitArray,
                                                                 validCellIndicesArray);

  contData.numTrisPerCell.ReleaseResourcesExecution();

  //compute for each output triangle what iteration of the input cell
  //generates it
  vtkm::cont::DeviceAdapterAlgorithm<DeviceAdapter>::LowerBounds(validCellIndicesArray,
                                                                 validCellIndicesArray,
                                                                 inputCellIterationNumber);

  //iteration cell numbers now look like:
  // 0, 0, 0, 3, 4, 4
  // and valid input cells look like:
  // 1, 1, 1, 3, 6, 6
  //so in the iso surface, we need to talk current output index - iteration number
  //to get the proper iteration number

  //generate a single tri per cell
  const vtkm::Id numTotalVertices = numOutputCells * 3;
  typedef ::worklets::IsosurfaceSingleTri<vtkm::Float32, vtkm::Float32> IsoSurfaceFunctor;
  IsoSurfaceFunctor isosurface(ISO_VALUE,
                               pdims,
                               field,
                               field,
                               verticesArray.PrepareForOutput(numTotalVertices, DeviceAdapter()),
                               scalarsArray.PrepareForOutput(numTotalVertices, DeviceAdapter())
                               );

  vtkm::worklet::DispatcherMapField< IsoSurfaceFunctor > isosurfaceDispatcher(isosurface);
  isosurfaceDispatcher.Invoke(validCellIndicesArray, inputCellIterationNumber);

};

static void RunMarchingCubes(int pdims[3],
                           std::vector<vtkm::Float32>& buffer,
                           std::string device,
                           std::string writeLoc,
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
    doMarchingCubes( pdims, dims, field, scalarsArray, verticesArray);

    double time = timer.GetElapsedTime();
    if(!silent)
      {
      std::cout << "num cells: " << (scalarsArray.GetNumberOfValues()/3)  << std::endl;
      std::cout << "vtkm," << device << "," << time << "," << trial << std::endl;
      }

    if(trial == MAX_NUM_TRIALS-1 && !writeLoc.empty())
      {
      writeLoc += "mc_" + device + "_stride_y.ply";
      saveAsPly(verticesArray, writeLoc);
      }
    }

}

}

