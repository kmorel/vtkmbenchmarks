

#include <vtkm/Pair.h>
#include <vtkm/Types.h>
#include <vtkm/worklet/DispatcherMapField.h>
#include <vtkm/worklet/WorkletMapField.h>

namespace worklets
{
//now that the device adapter is included set a global typedef
//that is the chosen device tag
typedef VTKM_DEFAULT_DEVICE_ADAPTER_TAG DeviceAdapter;

/// Linear interpolation
template <typename T1, typename T2>
VTKM_EXEC_EXPORT
T1 lerp(T1 a, T1 b, T2 t)
{
  return a + t*(b-a);
}

template< typename FieldType>
struct PortalTypes
{
public:
  typedef vtkm::cont::ArrayHandle<FieldType> HandleType;
  typedef typename HandleType::template ExecutionTypes<DeviceAdapter> ExecutionTypes;

  typedef typename ExecutionTypes::Portal Portal;
  typedef typename ExecutionTypes::PortalConst PortalConst;
};


/// \brief Compute isosurface vertices, and scalars
///
template <typename FieldType, typename OutputType, int NumCellsToFuse>
class IsosurfaceFusedUniformGridFunctor : public vtkm::worklet::WorkletMapField
{
public:
  typedef void ControlSignature(FieldIn<IdType> inputCellId, FieldIn<IdType> outputVertId);
  typedef void ExecutionSignature(_1, _2);
  typedef _1 InputDomain;

  typedef typename PortalTypes< vtkm::Id >::PortalConst IdPortalType;
  IdPortalType triangleTable, vertexTable;

  typedef typename PortalTypes< FieldType >::PortalConst FieldPortalType;
  FieldPortalType field, source;

  typedef typename PortalTypes< OutputType >::Portal ScalarPortalType;
  ScalarPortalType scalars;

  typedef typename PortalTypes< vtkm::Vec<OutputType,3> >::Portal VertexPortalType;
  VertexPortalType vertices;

  const int xdim, ydim, zdim, cellsPerLayer, pointsPerLayer;
  const float isovalue, xmin, ymin, zmin, xmax, ymax, zmax;

  const int inputCellIdOffset;

  template<typename U, typename V, typename W, typename X>
  VTKM_CONT_EXPORT
  IsosurfaceFusedUniformGridFunctor(const float isovalue,
                               const int dims[3],
                               const U & field,
                               const U & source,
                               const V & vertexTable,
                               const V & triangleTable,
                               const W & vertices,
                               const X & scalars,
                               const int inputIdOffset=0):
  isovalue(isovalue),
  xdim(dims[0]), ydim(dims[1]), zdim(dims[2]),
  xmin(-1), ymin(-1), zmin(-1),
  xmax(1), ymax(1), zmax(1),
  field( field.PrepareForInput( DeviceAdapter() ) ),
  source( source.PrepareForInput( DeviceAdapter() ) ),
  vertexTable( vertexTable.PrepareForInput( DeviceAdapter() ) ),
  triangleTable( triangleTable.PrepareForInput( DeviceAdapter() ) ),
  vertices(vertices),
  scalars(scalars),
  cellsPerLayer((xdim-1) * (ydim-1)),
  pointsPerLayer (xdim*ydim),
  inputCellIdOffset(inputIdOffset)
  {

  }

  VTKM_EXEC_EXPORT
  void operator()(vtkm::Id inputIndexId, vtkm::Id firstOutputVertId) const
  {
    // when operating on a slice of the data the inputIndexId
    // is relative to the start of the slice, so we need to
    // compute the proper global cell id.
    const int firstInputCellId = (inputCellIdOffset + inputIndexId) * NumCellsToFuse;

    // Get data for this cell
    const int verticesForEdge[] = { 0, 1, 1, 2, 3, 2, 0, 3,
                                    4, 5, 5, 6, 7, 6, 4, 7,
                                    0, 4, 1, 5, 2, 6, 3, 7 };

    vtkm::Id outputVertId = firstOutputVertId;
    for(int i=0; i < NumCellsToFuse; i++)
      {
      const vtkm::Id inputCellId = firstInputCellId + i;


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

      const int numVertices  = this->vertexTable.Get(cubeindex);

      // Compute the coordinates of the uniform regular grid at each of the cell's eight vertices
      vtkm::Vec<FieldType, 3> p[8];
      p[0] = vtkm::make_Vec(xmin+(xmax-xmin)*(1.0*x/(xdim-1)),     ymin+(ymax-ymin)*(1.0*y/(xdim-1)),     zmin+(zmax-zmin)*(1.0*z/(xdim-1)));
      p[1] = vtkm::make_Vec(xmin+(xmax-xmin)*(1.0*(x+1)/(xdim-1)), ymin+(ymax-ymin)*(1.0*y/(xdim-1)),     zmin+(zmax-zmin)*(1.0*z/(xdim-1)));
      p[2] = vtkm::make_Vec(xmin+(xmax-xmin)*(1.0*(x+1)/(xdim-1)), ymin+(ymax-ymin)*(1.0*(y+1)/(xdim-1)), zmin+(zmax-zmin)*(1.0*z/(xdim-1)));
      p[3] = vtkm::make_Vec(xmin+(xmax-xmin)*(1.0*x/(xdim-1)),     ymin+(ymax-ymin)*(1.0*(y+1)/(xdim-1)), zmin+(zmax-zmin)*(1.0*z/(xdim-1)));
      p[4] = vtkm::make_Vec(xmin+(xmax-xmin)*(1.0*x/(xdim-1)),     ymin+(ymax-ymin)*(1.0*y/(xdim-1)),     zmin+(zmax-zmin)*(1.0*(z+1)/(xdim-1)));
      p[5] = vtkm::make_Vec(xmin+(xmax-xmin)*(1.0*(x+1)/(xdim-1)), ymin+(ymax-ymin)*(1.0*y/(xdim-1)),     zmin+(zmax-zmin)*(1.0*(z+1)/(xdim-1)));
      p[6] = vtkm::make_Vec(xmin+(xmax-xmin)*(1.0*(x+1)/(xdim-1)), ymin+(ymax-ymin)*(1.0*(y+1)/(xdim-1)), zmin+(zmax-zmin)*(1.0*(z+1)/(xdim-1)));
      p[7] = vtkm::make_Vec(xmin+(xmax-xmin)*(1.0*x/(xdim-1)),     ymin+(ymax-ymin)*(1.0*(y+1)/(xdim-1)), zmin+(zmax-zmin)*(1.0*(z+1)/(xdim-1)));

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
      for (int v = 0; v < numVertices; v++)
        {
        const int edge = this->triangleTable.Get(cubeindex*16 + v);
        const int v0   = verticesForEdge[2*edge];
        const int v1   = verticesForEdge[2*edge + 1];
        const float t  = (isovalue - f[v0]) / (f[v1] - f[v0]);

        this->vertices.Set(outputVertId + v, lerp(p[v0], p[v1], t));
        this->scalars.Set(outputVertId + v, lerp(s[v0], s[v1], t));
        }

      outputVertId += numVertices;
      }
  }
};

/// \brief Computes Marching Cubes case number for each cell, along with the number of triangles generated by that case
///
template <typename FieldType, typename CellNumType>
class ClassifyCellOutputTri : public vtkm::worklet::WorkletMapField
{
public:
  typedef void ControlSignature(FieldIn<IdType> inputCellId, FieldOut<AllTypes> numCellsOut);
  typedef _2 ExecutionSignature(_1);
  typedef _1 InputDomain;

  typedef typename PortalTypes<FieldType>::PortalConst FieldPortalType;
  FieldPortalType pointData;

  typedef typename PortalTypes<vtkm::Id>::PortalConst TablePortalType;
  TablePortalType vertTable;

  float isovalue;
  int xdim, ydim, zdim;
  int cellsPerLayer;
  int pointsPerLayer;

  template<typename T, typename U>
  VTKM_CONT_EXPORT
  ClassifyCellOutputTri(const T& pointHandle,
               const U& vertTableHandle,
               float iso,
               int dims[3]) :
         pointData( pointHandle.PrepareForInput( DeviceAdapter() ) ),
         vertTable( vertTableHandle.PrepareForInput( DeviceAdapter() ) ),
         isovalue( iso ),
         xdim(dims[0]),
         ydim(dims[1]),
         zdim(dims[2]),
         cellsPerLayer((xdim - 1) * (ydim - 1)),
         pointsPerLayer (xdim*ydim)
   {

   }

  VTKM_EXEC_EXPORT
  CellNumType operator()(vtkm::Id cellId) const
  {
    // Compute 3D indices of this cell
    const int x = cellId % (xdim - 1);
    const int y = (cellId / (xdim - 1)) % (ydim -1);
    const int z = cellId / cellsPerLayer;

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
    const float f0 = this->pointData.Get(i0);
    const float f1 = this->pointData.Get(i1);
    const float f2 = this->pointData.Get(i2);
    const float f3 = this->pointData.Get(i3);
    const float f4 = this->pointData.Get(i4);
    const float f5 = this->pointData.Get(i5);
    const float f6 = this->pointData.Get(i6);
    const float f7 = this->pointData.Get(i7);

    // Compute the Marching Cubes case number for this cell
    unsigned int cubeindex = (f0 > isovalue);
    cubeindex += (f1 > isovalue)*2;
    cubeindex += (f2 > isovalue)*4;
    cubeindex += (f3 > isovalue)*8;
    cubeindex += (f4 > isovalue)*16;
    cubeindex += (f5 > isovalue)*32;
    cubeindex += (f6 > isovalue)*64;
    cubeindex += (f7 > isovalue)*128;

    return CellNumType(this->vertTable.Get(cubeindex) / 3);
  }
};

/// \brief Compute isosurface vertices and scalars
///
template <typename FieldType, typename OutputType>
class IsosurfaceSingleTri : public vtkm::worklet::WorkletMapField
{
public:
  typedef void ControlSignature(FieldIn<IdType> inputCellId,
                                FieldIn<IdType> inputIteration);
  typedef void ExecutionSignature(WorkIndex, _1, _2);
  typedef _1 InputDomain;

  typedef typename PortalTypes< vtkm::Id >::PortalConst IdPortalType;
  IdPortalType triangleTable, vertexTable;

  typedef typename PortalTypes< FieldType >::PortalConst FieldPortalType;
  FieldPortalType field, source;

  typedef typename PortalTypes< OutputType >::Portal ScalarPortalType;
  ScalarPortalType scalars;

  typedef typename PortalTypes< vtkm::Vec<vtkm::Float32,3> >::Portal VertexPortalType;
  VertexPortalType vertices;

  const int xdim, ydim, zdim, cellsPerLayer, pointsPerLayer;
  const float isovalue, xmin, ymin, zmin, xmax, ymax, zmax;

  const int inputCellIdOffset;

  template<typename U, typename V, typename W, typename X>
  VTKM_CONT_EXPORT
  IsosurfaceSingleTri( const float isovalue,
                       const int dims[3],
                       const U & field,
                       const U & source,
                       const V & vertexTable,
                       const V & triangleTable,
                       const W & vertices,
                       const X & scalars,
                       const int inputIdOffset=0):
  isovalue(isovalue),
  xdim(dims[0]), ydim(dims[1]), zdim(dims[2]),
  xmin(-1), ymin(-1), zmin(-1),
  xmax(1), ymax(1), zmax(1),
  field( field.PrepareForInput( DeviceAdapter() ) ),
  source( source.PrepareForInput( DeviceAdapter() ) ),
  vertexTable( vertexTable.PrepareForInput( DeviceAdapter() ) ),
  triangleTable( triangleTable.PrepareForInput( DeviceAdapter() ) ),
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

    const int numVertices  = this->vertexTable.Get(cubeindex);

    // Compute the coordinates of the uniform regular grid at each of the cell's eight vertices
    vtkm::Vec<FieldType, 3> p[8];
    p[0] = vtkm::make_Vec(xmin+(xmax-xmin)*(1.0*x/(xdim-1)),     ymin+(ymax-ymin)*(1.0*y/(xdim-1)),     zmin+(zmax-zmin)*(1.0*z/(xdim-1)));
    p[1] = vtkm::make_Vec(xmin+(xmax-xmin)*(1.0*(x+1)/(xdim-1)), ymin+(ymax-ymin)*(1.0*y/(xdim-1)),     zmin+(zmax-zmin)*(1.0*z/(xdim-1)));
    p[2] = vtkm::make_Vec(xmin+(xmax-xmin)*(1.0*(x+1)/(xdim-1)), ymin+(ymax-ymin)*(1.0*(y+1)/(xdim-1)), zmin+(zmax-zmin)*(1.0*z/(xdim-1)));
    p[3] = vtkm::make_Vec(xmin+(xmax-xmin)*(1.0*x/(xdim-1)),     ymin+(ymax-ymin)*(1.0*(y+1)/(xdim-1)), zmin+(zmax-zmin)*(1.0*z/(xdim-1)));
    p[4] = vtkm::make_Vec(xmin+(xmax-xmin)*(1.0*x/(xdim-1)),     ymin+(ymax-ymin)*(1.0*y/(xdim-1)),     zmin+(zmax-zmin)*(1.0*(z+1)/(xdim-1)));
    p[5] = vtkm::make_Vec(xmin+(xmax-xmin)*(1.0*(x+1)/(xdim-1)), ymin+(ymax-ymin)*(1.0*y/(xdim-1)),     zmin+(zmax-zmin)*(1.0*(z+1)/(xdim-1)));
    p[6] = vtkm::make_Vec(xmin+(xmax-xmin)*(1.0*(x+1)/(xdim-1)), ymin+(ymax-ymin)*(1.0*(y+1)/(xdim-1)), zmin+(zmax-zmin)*(1.0*(z+1)/(xdim-1)));
    p[7] = vtkm::make_Vec(xmin+(xmax-xmin)*(1.0*x/(xdim-1)),     ymin+(ymax-ymin)*(1.0*(y+1)/(xdim-1)), zmin+(zmax-zmin)*(1.0*(z+1)/(xdim-1)));

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
    const vtkm::Id cellOffset = cubeindex*16 + (inputIteration * 3);
    for (int v = 0; v < 3; v++)
      {
      const int edge = this->triangleTable.Get( cellOffset + v );
      const int v0   = verticesForEdge[2*edge];
      const int v1   = verticesForEdge[2*edge + 1];
      const float t  = (isovalue - f[v0]) / (f[v1] - f[v0]);

      this->vertices.Set(outputVertId + v, lerp(p[v0], p[v1], t));
      this->scalars.Set(outputVertId + v, lerp(s[v0], s[v1], t));
      }
  }
};

}