
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
  vtkm::cont::ArrayHandle< vtkm::UInt8 > yEdgeCaseHandle;

  vtkm::cont::ArrayHandle< vtkm::Id3 > sumsHandle;

  vtkm::cont::ArrayHandle< vtkm::Id > numYTrisPerRowHandle;

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
  typedef typename ::worklets::PortalTypes<vtkm::Id3>::Portal Id3PortalType;

  float isovalue;
  vtkm::Id3 dims;
  int inc1;
  int inc2;
  int sliceOffset;
  FieldPortalType scalarData;
  UInt8PortalType edgeYCases;

  Id3PortalType sums;

  explicit
  MarchingEdgePass1ExecData( MarchingEdgeContData<FieldType>& metaData):
    isovalue(metaData.isovalue),
    dims(metaData.dims[0], metaData.dims[1], metaData.dims[2]),
    inc1(metaData.inc1),
    inc2(metaData.inc2),
    sliceOffset(metaData.sliceOffset),
    scalarData(metaData.scalarHandle.PrepareForInput(DeviceAdapter())),
    edgeYCases(metaData.yEdgeCaseHandle.PrepareForOutput( metaData.dims[0] * (metaData.dims[1]-1) * metaData.dims[2], DeviceAdapter() )),
    sums(metaData.sumsHandle.PrepareForOutput(metaData.dims[0] * metaData.dims[2], DeviceAdapter()) )
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
  typedef typename ::worklets::PortalTypes<vtkm::Id>::PortalConst IdConstPortalType;
  typedef typename ::worklets::PortalTypes<vtkm::Id>::Portal IdPortalType;
  typedef typename ::worklets::PortalTypes<vtkm::Id3>::Portal Id3PortalType;


  float isovalue;
  vtkm::Id3 dims;
  int inc1;
  int inc2;
  int sliceOffset;

  UInt8ConstPortalType edgeYCases;

  IdPortalType numYTris;
  Id3PortalType sums;

  explicit
  MarchingEdgePass2ExecData( MarchingEdgeContData<FieldType>& metaData):
    isovalue(metaData.isovalue),
    dims(metaData.dims[0], metaData.dims[1], metaData.dims[2]),
    inc1(metaData.inc1),
    inc2(metaData.inc2),
    sliceOffset(metaData.sliceOffset),
    edgeYCases(metaData.yEdgeCaseHandle.PrepareForInput( DeviceAdapter() )),
    numYTris(metaData.numYTrisPerRowHandle.PrepareForOutput( (metaData.dims[0]-1) * (metaData.dims[2]-1), DeviceAdapter()) ),
    sums(metaData.sumsHandle.PrepareForInPlace( DeviceAdapter()) )
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
struct MarchingEdgePass3ExecData
{
  typedef typename ::worklets::PortalTypes<FieldType>::PortalConst FieldPortalType;
  typedef typename ::worklets::PortalTypes<vtkm::UInt8>::PortalConst UInt8ConstPortalType;
  typedef typename ::worklets::PortalTypes<vtkm::Id>::PortalConst IdConstPortalType;

  typedef typename ::worklets::PortalTypes< FieldType >::Portal ScalarPortalType;
  typedef typename ::worklets::PortalTypes< vtkm::Vec<vtkm::Float32,3> >::Portal VertexPortalType;


  float isovalue;
  vtkm::Id3 dims;
  int inc1;
  int inc2;
  int sliceOffset;

  UInt8ConstPortalType edgeYCases;
  IdConstPortalType numYTris;
  ScalarPortalType scalars;
  VertexPortalType vertices;

  template<typename T, typename U>
  MarchingEdgePass3ExecData( MarchingEdgeContData<FieldType>& metaData,
                             const T& scalarsPortal,
                             const U& verticesPortal ):
    isovalue(metaData.isovalue),
    dims(metaData.dims[0], metaData.dims[1], metaData.dims[2]),
    inc1(metaData.inc1),
    inc2(metaData.inc2),
    sliceOffset(metaData.sliceOffset),
    edgeYCases(metaData.yEdgeCaseHandle.PrepareForInput( DeviceAdapter() )),
    numYTris(metaData.numYTrisPerRowHandle.PrepareForInput( DeviceAdapter()) ),
    scalars( scalarsPortal ),
    vertices( verticesPortal )

  {

  }

  VTKM_EXEC_EXPORT
  vtkm::UInt8 numberOfPrimitives( vtkm::UInt8 ecase ) const
    { return EdgeCasesTable[ecase][0]; }

};


//----------------------------------------------------------------------------
template <typename FieldType>
class ProcessYEdges : public vtkm::worklet::WorkletMapField
{
public:
  typedef void ControlSignature(FieldIn<IdType>  id);
  typedef void ExecutionSignature(_1);
  typedef _1 InputDomain;

  MarchingEdgePass1ExecData<FieldType> d;

  VTKM_CONT_EXPORT
  ProcessYEdges( MarchingEdgePass1ExecData<FieldType>& data) :
         d( data )
   {
   }

  //----------------------------------------------------------------------------
  // PASS 1: Process a single volume y-row (and all of the voxel edges that
  // compose the row). Determine the y-edges case classification.
  VTKM_EXEC_EXPORT
  void operator()(const vtkm::Id&  xzid) const
  {
  //now instead of doing every 1 dim, we need to calculate out which
  //row and slice we are.
  const vtkm::Id row = xzid % this->d.dims[0];
  const vtkm::Id slice = xzid / this->d.dims[0];
  const vtkm::Id offset = slice * this->d.inc2 + row;
  this->process(offset, row, slice);
  }

  VTKM_EXEC_EXPORT
  void process(vtkm::Id offset, vtkm::Id row, vtkm::Id slice) const
  {
  const vtkm::Id nycells = this->d.dims[1]-1;
  const vtkm::Id metaWritePos = slice * nycells + row;
  const vtkm::Id writePos = slice * this->d.sliceOffset + row;

  FieldType s0=0;
  FieldType s1= this->d.scalarData.Get(offset);

  //for all cells in this row
  vtkm::Id tempYSum = 0;
  for (int i=0; i < nycells; ++i)
    {
    s0 = s1;
    s1 = this->d.scalarData.Get(offset + this->d.inc1 * (i + 1) );

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
    this->d.edgeYCases.Set(writePos + (i * nycells), edgeCase);
    if ( edgeCase == flyingedges::LeftAbove ||
         edgeCase == flyingedges::RightAbove )
      {
      ++tempYSum;
      }
    }
  this->d.sums.Set(metaWritePos, vtkm::make_Vec(0,tempYSum,0) ); //write back the number of intersections along x-edge
  }
};


template <typename FieldType>
class ProcessXZEdges : public vtkm::worklet::WorkletMapField
{
public:
  typedef void ControlSignature(FieldIn<IdType>  id);
  typedef void ExecutionSignature(_1);
  typedef _1 InputDomain;

  MarchingEdgePass2ExecData<FieldType> d;

  VTKM_CONT_EXPORT
  ProcessXZEdges( MarchingEdgePass2ExecData<FieldType>& data) :
         d( data )
   {
   }

  //----------------------------------------------------------------------------
  // PASS 2: Process a single y-row of voxels. Count the number of x- and
  // z-intersections by topological reasoning from y-edge cases. Determine the
  // number of primitives (i.e., triangles) generated from this row.
  VTKM_EXEC_EXPORT
  void operator()(const vtkm::Id&  xzid) const
  {
  //now instead of doing every 1 dim, we need to calculate out which
  //row and slice we are.
  const vtkm::Id xcdim = (this->d.dims[0]-1);
  const vtkm::Id row = xzid % xcdim;
  const vtkm::Id slice = xzid / xcdim;
  this->process(xzid, row, slice);
  }

  VTKM_EXEC_EXPORT
  void process(vtkm::Id xzid, vtkm::Id row, vtkm::Id slice) const
  {
  const vtkm::Id nycells = this->d.dims[1]-1;
  const vtkm::Id metaWritePos = xzid; // no need to compute this
  const vtkm::Id readPos = slice * this->d.sliceOffset + row;

  // Okay run along the x-voxels and count the number of y- and
  // z-intersections. Here we are just checking y,z edges that make up the
  // voxel axes. Also check the number of primitives generated.
  vtkm::UInt8 edgeCase[4] = {0, 0, 0, 0};

  //don't use tempSums as a local write back location
  //so that the for loop doesn't require this memory
  //read before it can iterate
  vtkm::Id3 tempSums = this->d.sums.Get(metaWritePos);
  vtkm::Id tempTriSum = 0;
  vtkm::Id tempXSum = 0;
  vtkm::Id tempZSum = 0;
  for (int i=0; i < nycells; ++i)
    {
    edgeCase[0] = this->d.edgeYCases.Get( (i * nycells) + readPos);
    edgeCase[1] = this->d.edgeYCases.Get( (i * nycells) + readPos + 1 );
    edgeCase[2] = this->d.edgeYCases.Get( (i * nycells) + readPos + this->d.sliceOffset);
    edgeCase[3] = this->d.edgeYCases.Get( (i * nycells) + readPos + this->d.sliceOffset + 1 );

    const vtkm::UInt8 eCase = (edgeCase[0]    |
                               edgeCase[1]<<2 |
                               edgeCase[2]<<4 |
                               edgeCase[3]<<6 );

    //using the eCase for this voxel we need to lookup the number of primitives
    //that this voxel will generate.
    const vtkm::UInt8 numTris = this->d.numberOfPrimitives(eCase);
    if(numTris != 0)
      {
      tempTriSum += numTris;
      tempXSum += this->d.edgeUses(eCase, 2); //x-voxel axes edge always counted
      tempZSum += this->d.edgeUses(eCase, 8); //z-voxel axes edge always counted
      }
    }//for all voxels along this x-edge

  this->d.numYTris.Set(metaWritePos,tempTriSum); //write back the number of triangles

  //write back the number of intersections along x,y,z-edges
  tempSums[0] = tempXSum;
  tempSums[2] = tempZSum;
  this->d.sums.Set(metaWritePos,tempSums);
  }
};

template <typename FieldType>
class GenerateTriangles : public vtkm::worklet::WorkletMapField
{
public:
  typedef void ControlSignature(FieldIn<IdType>  id);
  typedef void ExecutionSignature(_1);
  typedef _1 InputDomain;

  MarchingEdgePass3ExecData<FieldType> d;

  VTKM_CONT_EXPORT
  GenerateTriangles( MarchingEdgePass3ExecData<FieldType>& data) :
         d( data )
   {
   }

  //----------------------------------------------------------------------------
  // PASS 3: Generate Triangles
  VTKM_EXEC_EXPORT
  void operator()(const vtkm::Id&  xzid) const
  {
  //now instead of doing every 1 dim, we need to calculate out which
  //row and slice we are.
  const vtkm::Id xcdim = (this->d.dims[0]-1);
  const vtkm::Id row = xzid % xcdim;
  const vtkm::Id slice = xzid / xcdim;
  this->process(xzid, row, slice);
  }

  VTKM_EXEC_EXPORT
  void process(vtkm::Id xzid, vtkm::Id row, vtkm::Id slice) const
  {
  const vtkm::Id nycells = this->d.dims[1]-1;
  const vtkm::Id metaReadPos = xzid; // no need to compute this
  const vtkm::Id readPos = slice * this->d.sliceOffset + row;

  // Okay run along the x-voxels and count the number of y- and
  // z-intersections. Here we are just checking y,z edges that make up the
  // voxel axes. Also check the number of primitives generated.
  vtkm::UInt8 edgeCase[4] = {0, 0, 0, 0};

  vtkm::Id currentTriIncScan = this->d.numYTris.Get(metaReadPos); //read the location to write out triangle
  for (int i=0; i < nycells; ++i)
    {
    edgeCase[0] = this->d.edgeYCases.Get( (i * nycells) + readPos);
    edgeCase[1] = this->d.edgeYCases.Get( (i * nycells) + readPos + 1 );
    edgeCase[2] = this->d.edgeYCases.Get( (i * nycells) + readPos + this->d.sliceOffset);
    edgeCase[3] = this->d.edgeYCases.Get( (i * nycells) + readPos + this->d.sliceOffset + 1 );

    }//for all voxels along this x-edge
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

  // Pass 0: Allocate the yEdgeCase. This is done here to make it easier
  //to pass to worklets
  vtkm::Vec<vtkm::Float32,3> origin(0,0,0);
  vtkm::Vec<vtkm::Float32,3> spacing(1,1,1);
  worklets::MarchingEdgeContData<vtkm::Float32> contData(origin,
                                                       spacing,
                                                       pdims,
                                                       field,
                                                       ISO_VALUE);

  //----------------------------------------------------------------------------
  // PASS 1: Process a single volume y-row (and all of the voxel edges that
  // compose the row). Determine the y-edges case classification.
  {
  worklets::MarchingEdgePass1ExecData<vtkm::Float32> execData( contData );
  //first pass is over edges so use point dims
  vtkm::cont::ArrayHandleCounting<vtkm::Id> passIds(0, pdims[0] * pdims[2]);

  typedef worklets::ProcessYEdges< vtkm::Float32 > ProcessYEdgeFunctor;
  typedef vtkm::worklet::DispatcherMapField< ProcessYEdgeFunctor > YEdgeDispatcher;

  ProcessYEdgeFunctor yEdge( execData );
  YEdgeDispatcher yEDispatcher(yEdge);
  yEDispatcher.Invoke(passIds);
  }

  // PASS 2: Process a single y-row of voxels. Count the number of x- and
  // z-intersections by topological reasoning from y-edge cases. Determine the
  // number of primitives (i.e., triangles) generated from this row..
  {
  worklets::MarchingEdgePass2ExecData<vtkm::Float32> execData( contData );
  //first pass is over 'cells' so use cell dims
  vtkm::cont::ArrayHandleCounting<vtkm::Id> passIds(0, cdims[0] * cdims[2]);

  typedef worklets::ProcessXZEdges< vtkm::Float32 > ProcessXZEdgeFunctor;
  typedef vtkm::worklet::DispatcherMapField< ProcessXZEdgeFunctor > XZEdgeDispatcher;

  ProcessXZEdgeFunctor xzFunctor( execData );
  XZEdgeDispatcher xzEDispatcher(xzFunctor);
  xzEDispatcher.Invoke(passIds);
  }

  // PASS 3: Now compute the number of output triangles, and
  // have each row write out the triangles it generated
  //
  {
  const vtkm::Id numOutputYTriangles =
      vtkm::cont::DeviceAdapterAlgorithm<DeviceAdapter>::ScanInclusive(contData.numYTrisPerRowHandle,
                                                                       contData.numYTrisPerRowHandle);
  //terminate if no cells has triangles left
  std::cout << "numOutputYTriangles: " << numOutputYTriangles << std::endl;

  //We convert the sumsHandle to an int array and run an exclusive scan. The
  //result of this operation is an array that represents the location that
  //each valid edge should write their point into memory
  vtkm::cont::DeviceAdapterAlgorithm<DeviceAdapter>::ScanExclusive(
    *reinterpret_cast< vtkm::cont::ArrayHandle< vtkm::Id >* >(&contData.sumsHandle),
    *reinterpret_cast< vtkm::cont::ArrayHandle< vtkm::Id >* >(&contData.sumsHandle));


  // const vtkm::Id  numOutputXTriangles =
  //     vtkm::cont::DeviceAdapterAlgorithm<DeviceAdapter>::ScanInclusive(contData.numXTrisPerRowHandle,
  //                                                                      contData.numXTrisPerRowHandle);

  // std::cout << "numOutputXTriangles: " << numOutputXTriangles << std::endl;


  if(numOutputYTriangles == 0)
    { return; }

  worklets::MarchingEdgePass3ExecData<vtkm::Float32> execData( contData,
                                                               scalarsArray.PrepareForOutput(numOutputYTriangles*3, DeviceAdapter()),
                                                               verticesArray.PrepareForOutput(numOutputYTriangles*3, DeviceAdapter())
                                                               );
  //first pass is over 'cells' so use cell dims
  vtkm::cont::ArrayHandleCounting<vtkm::Id> passIds(0, cdims[0] * cdims[2]);

  typedef worklets::GenerateTriangles< vtkm::Float32 > GenerateFunctor;
  typedef vtkm::worklet::DispatcherMapField< GenerateFunctor > GDispatcher;

  GenerateFunctor generateFunctor( execData );
  GDispatcher genDispatcher(generateFunctor);
  genDispatcher.Invoke(passIds);
  }

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
