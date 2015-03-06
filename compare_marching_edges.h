
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

  vtkm::cont::ArrayHandle< vtkm::Id > numTrisPerCell;

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

  explicit
  MarchingEdgePass1ExecData( MarchingEdgeContData<FieldType>& metaData):
    isovalue(metaData.isovalue),
    dims(metaData.dims[0], metaData.dims[1], metaData.dims[2]),
    inc1(metaData.inc1),
    inc2(metaData.inc2),
    sliceOffset(metaData.sliceOffset),
    scalarData(metaData.scalarHandle.PrepareForInput(DeviceAdapter())),
    edgeYCases(metaData.yEdgeCaseHandle.PrepareForOutput( metaData.dims[0] * (metaData.dims[1]-1) * metaData.dims[2], DeviceAdapter() ))
  {
    // std::cout << "MarchingEdgeExecData" << std::endl;
    // std::cout << "algo.Dims[0]: " << dims[0] << std::endl;
    // std::cout << "algo.Dims[1]: " << dims[1] << std::endl;
    // std::cout << "algo.Dims[2]: " << dims[2] << std::endl;
    // std::cout << "algo.SliceOffset: " << sliceOffset << std::endl;;

    // std::cout << "algo.Inc[0]: " << 1 << std::endl;
    // std::cout << "algo.Inc[1]: " << inc1 << std::endl;
    // std::cout << "algo.Inc[2]: " << inc2 << std::endl;

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

  IdPortalType numTris;

  explicit
  MarchingEdgePass2ExecData( MarchingEdgeContData<FieldType>& metaData):
    isovalue(metaData.isovalue),
    dims(metaData.dims[0], metaData.dims[1], metaData.dims[2]),
    inc1(metaData.inc1),
    inc2(metaData.inc2),
    sliceOffset(metaData.sliceOffset),
    edgeYCases(metaData.yEdgeCaseHandle.PrepareForInput( DeviceAdapter() )),
    numTris(metaData.numTrisPerCell.PrepareForOutput( (metaData.dims[0]-1) * (metaData.dims[1]-1) * (metaData.dims[2]-1), DeviceAdapter()) )
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
  const vtkm::Id writePos = slice * this->d.sliceOffset + row;

  FieldType s0=0;
  FieldType s1= this->d.scalarData.Get(offset);

  //for all cells in this row
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
    }
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
  const vtkm::Id ycdim = (this->d.dims[1]-1);

  const vtkm::Id row = xzid % xcdim;
  const vtkm::Id slice = xzid / xcdim;
  const vtkm::Id offset = slice * (xcdim * ycdim) + row;

  this->process(xzid, offset, row, slice);
  }

  VTKM_EXEC_EXPORT
  void process(vtkm::Id xzid, vtkm::Id writeOffset, vtkm::Id row, vtkm::Id slice) const
  {
  const vtkm::Id nycells = this->d.dims[1]-1;
  const vtkm::Id readPos = slice * this->d.sliceOffset + row;

  // Okay run along the x-voxels and count the number of y- and
  // z-intersections. Here we are just checking y,z edges that make up the
  // voxel axes. Also check the number of primitives generated.
  vtkm::UInt8 edgeCase[4] = {0, 0, 0, 0};

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

    const vtkm::UInt8 numTris = this->d.numberOfPrimitives(eCase);
    this->d.numTris.Set(writeOffset + (i * (this->d.dims[0]-1)), numTris);
    }//for all voxels along this y-edge

  }
};

/// \brief Compute isosurface vertices and scalars
///
template <typename FieldType, typename OutputType>
class IsoMarchingTri : public vtkm::worklet::WorkletMapField
{
public:
  typedef void ControlSignature(FieldIn<IdType> inputCellId,
                                FieldIn<IdType> inputIteration);
  typedef void ExecutionSignature(WorkIndex, _1, _2);
  typedef _1 InputDomain;

  typedef typename ::worklets::PortalTypes< FieldType >::PortalConst FieldPortalType;
  FieldPortalType field, source;

  typedef typename ::worklets::PortalTypes< OutputType >::Portal ScalarPortalType;
  ScalarPortalType scalars;

  typedef typename ::worklets::PortalTypes< vtkm::Vec<vtkm::Float32,3> >::Portal VertexPortalType;
  VertexPortalType vertices;

  const int xdim, ydim, zdim, cellsPerLayer, pointsPerLayer;
  const float isovalue, xmin, ymin, zmin, xmax, ymax, zmax;

  const int inputCellIdOffset;

  template<typename U, typename W, typename X>
  VTKM_CONT_EXPORT
  IsoMarchingTri( const float isovalue,
                       const int dims[3],
                       const U & field,
                       const U & source,
                       const W & vertices,
                       const X & scalars,
                       const int inputIdOffset=0):
  isovalue(isovalue),
  xdim(dims[0]), ydim(dims[1]), zdim(dims[2]),
  xmin(-1), ymin(-1), zmin(-1),
  xmax(1), ymax(1), zmax(1),
  field( field.PrepareForInput( DeviceAdapter() ) ),
  source( source.PrepareForInput( DeviceAdapter() ) ),
  vertices(vertices),
  scalars(scalars),
  cellsPerLayer((xdim-1) * (ydim-1)),
  pointsPerLayer (xdim*ydim),
  inputCellIdOffset(inputIdOffset)
  {

  }

  VTKM_EXEC_EXPORT
  void operator()(vtkm::Id outputCellId, vtkm::Id inputIndexId, vtkm::Id inputLowerBounds) const
  {
    // Get data for this cell
    const int verticesForEdge[] = { 0, 1, 1, 2, 3, 2, 0, 3,
                                    4, 5, 5, 6, 7, 6, 4, 7,
                                    0, 4, 1, 5, 2, 6, 3, 7 };

    // when operating on a slice of the data the inputIndexId
    // is relative to the start of the slice, so we need to
    // compute the proper global cell id.
    const int inputCellId = inputCellIdOffset + inputIndexId;

    const int x = inputCellId % (xdim - 1);
    const int y = (inputCellId / (xdim - 1)) % (ydim -1);
    const int z = inputCellId / cellsPerLayer;

    // Compute indices for the eight vertices of this cell
    const int i0 = x    + y*xdim + z * pointsPerLayer;
    const int i1 = i0   + 1;
    const int i2 = i0   + 1 + xdim;
    const int i3 = i0   + xdim;
    const int i4 = i0   + pointsPerLayer;
    const int i5 = i1   + pointsPerLayer;
    const int i6 = i2   + pointsPerLayer;
    const int i7 = i3   + pointsPerLayer;

    // Get the field values at these eight vertices
    float f[8];
    f[0] = this->field.Get(i0);
    f[1] = this->field.Get(i1);
    f[2] = this->field.Get(i2);
    f[3] = this->field.Get(i3);
    f[4] = this->field.Get(i4);
    f[5] = this->field.Get(i5);
    f[6] = this->field.Get(i6);
    f[7] = this->field.Get(i7);

    // Compute the Marching Cubes case number for this cell
    unsigned int cubeindex = 0;
    cubeindex += (f[0] > isovalue);
    cubeindex += (f[1] > isovalue)*2;
    cubeindex += (f[2] > isovalue)*4;
    cubeindex += (f[3] > isovalue)*8;
    cubeindex += (f[4] > isovalue)*16;
    cubeindex += (f[5] > isovalue)*32;
    cubeindex += (f[6] > isovalue)*64;
    cubeindex += (f[7] > isovalue)*128;

    // Compute the coordinates of the uniform regular grid at each of the cell's eight vertices
    vtkm::Vec<FieldType, 3> p[8];
    {
    //if we have offset and spacing, can we simplify this computation
    vtkm::Vec<FieldType, 3> offset = vtkm::make_Vec(xmin+(xmax-xmin),
                                                    ymin+(ymax-ymin),
                                                    zmin+(zmax-zmin) );

    vtkm::Vec<FieldType, 3> spacing = vtkm::make_Vec( 1.0 /(xdim-1),
                                                      1.0 /(ydim-1),
                                                      1.0 /(zdim-1));

    vtkm::Vec<FieldType, 3> firstPoint = offset * spacing *  vtkm::make_Vec( x, y, z );
    vtkm::Vec<FieldType, 3> secondPoint = offset * spacing * vtkm::make_Vec( x+1, y+1, z+1 );

    p[0] = vtkm::make_Vec( firstPoint[0],   firstPoint[1],   firstPoint[2]);
    p[1] = vtkm::make_Vec( secondPoint[0],  firstPoint[1],   firstPoint[2]);
    p[2] = vtkm::make_Vec( secondPoint[0],  secondPoint[1],  firstPoint[2]);
    p[3] = vtkm::make_Vec( firstPoint[0],   secondPoint[1],  firstPoint[2]);
    p[4] = vtkm::make_Vec( firstPoint[0],   firstPoint[1],   secondPoint[2]);
    p[5] = vtkm::make_Vec( secondPoint[0],  firstPoint[1],   secondPoint[2]);
    p[6] = vtkm::make_Vec( secondPoint[0],  secondPoint[1],  secondPoint[2]);
    p[7] = vtkm::make_Vec( firstPoint[0],   secondPoint[1],  secondPoint[2]);
    }

    // Get the scalar source values at the eight vertices
    float s[8];
    s[0] = this->source.Get(i0);
    s[1] = this->source.Get(i1);
    s[2] = this->source.Get(i2);
    s[3] = this->source.Get(i3);
    s[4] = this->source.Get(i4);
    s[5] = this->source.Get(i5);
    s[6] = this->source.Get(i6);
    s[7] = this->source.Get(i7);

    // Interpolate for vertex positions and associated scalar values
    const vtkm::Id inputIteration = (outputCellId - inputLowerBounds);
    const vtkm::Id outputVertId = outputCellId * 3;
    for (int v = 0; v < 3; v++)
      {
      const int edge = EdgeCasesTable[cubeindex][v + (inputIteration * 3)];
      const int v0   = verticesForEdge[2*edge];
      const int v1   = verticesForEdge[2*edge + 1];
      const float t  = (isovalue - f[v0]) / (f[v1] - f[v0]);

      this->vertices.Set(outputVertId + v, ::worklets::lerp(p[v0], p[v1], t));
      this->scalars.Set(outputVertId + v, ::worklets::lerp(s[v0], s[v1], t));
      }
  }
};

} //worklets namespace

//now that the device adapter is included set a global typedef
//that is the chosen device tag
typedef VTKM_DEFAULT_DEVICE_ADAPTER_TAG DeviceAdapter;

static void doMarchingEdges( int pdims[3], //point dims
                           int cdims[3], //cell dims
                           vtkm::cont::ArrayHandle<vtkm::Float32> field,
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

  contData.yEdgeCaseHandle.ReleaseResources();

  // PASS 3: Now compute the number of output triangles, and
  // have each row write out the triangles it generated
  //
  {
  const vtkm::Id numOutputCells =
      vtkm::cont::DeviceAdapterAlgorithm<DeviceAdapter>::ScanInclusive(contData.numTrisPerCell,
                                                                       contData.numTrisPerCell);

  //terminate if no cells has triangles left
  if(numOutputCells == 0)
    {return; }

  vtkm::cont::ArrayHandle< vtkm::Id > validCellIndicesArray;
  vtkm::cont::ArrayHandle< vtkm::Id > inputCellIterationNumber;

  vtkm::cont::ArrayHandleCounting<vtkm::Id> validCellCountImplicitArray(0, numOutputCells);
  vtkm::cont::DeviceAdapterAlgorithm<DeviceAdapter>::UpperBounds(contData.numTrisPerCell,
                                                                 validCellCountImplicitArray,
                                                                 validCellIndicesArray);

  contData.numTrisPerCell.ReleaseResources();

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
  typedef worklets::IsoMarchingTri<vtkm::Float32, vtkm::Float32> IsoMarchingFunctor;
  //uses the EdgeCasesTable and not the MC case table
  IsoMarchingFunctor isosurface(ISO_VALUE,
                               pdims,
                               field,
                               field,
                               verticesArray.PrepareForOutput(numTotalVertices, DeviceAdapter()),
                               scalarsArray.PrepareForOutput(numTotalVertices, DeviceAdapter())
                               );

  vtkm::worklet::DispatcherMapField< IsoMarchingFunctor > isosurfaceDispatcher(isosurface);
  isosurfaceDispatcher.Invoke(validCellIndicesArray, inputCellIterationNumber);
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
