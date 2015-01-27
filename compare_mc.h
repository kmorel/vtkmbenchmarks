
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
#include "isosurface.h"

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

/// Vector cross-product
template <typename T>
VTKM_EXEC_EXPORT
T cross(T a, T b)
{
  return vtkm::make_Vec(a[1]*b[2] - a[2]*b[1], a[2]*b[0] - a[0]*b[2], a[0]*b[1] - a[1]*b[0]);
}

/// Vector normalization
template <typename T>
VTKM_EXEC_EXPORT
T normalize(T v)
{
  return ((1.0f / sqrt(vtkm::dot(v,v))) * v);
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



/// \brief Computes Marching Cubes case number for each cell, along with the number of vertices generated by that case
///
template <typename FieldType>
class ClassifyCell : public vtkm::worklet::WorkletMapField
{
public:
  typedef void ControlSignature(FieldIn<IdType> inputCellId, FieldOut<AllTypes> outputCaseInfo);
  typedef _2 ExecutionSignature(_1);
  typedef _1 InputDomain;

  typedef typename PortalTypes<FieldType>::PortalConst FieldPortalType;
  FieldPortalType pointData;

  typedef typename PortalTypes<vtkm::Id>::PortalConst TablePortalType;
  TablePortalType vertexTable;

  float isovalue;
  int xdim, ydim, zdim;
  int cellsPerLayer;
  int pointsPerLayer;

  template<typename T, typename U>
  VTKM_CONT_EXPORT
  ClassifyCell(const T& pointHandle,
               const U& vertTableHandle,
               float iso,
               int dims[3]) :
         pointData( pointHandle.PrepareForInput( DeviceAdapter() ) ),
         vertexTable( vertTableHandle.PrepareForInput( DeviceAdapter() ) ),
         isovalue( iso ),
         xdim(dims[0]),
         ydim(dims[1]),
         zdim(dims[2]),
         cellsPerLayer((xdim - 1) * (ydim - 1)),
         pointsPerLayer (xdim*ydim)
   {

   }

  VTKM_EXEC_EXPORT
  vtkm::Pair<vtkm::Id, vtkm::Id> operator()(const vtkm::Id &cellId) const
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

    // Return the Marching Cubes case number and the number of vertices this case generates
    return vtkm::make_Pair(cubeindex, this->vertexTable.Get(cubeindex));
  }
};


/// \brief Return whether the cell generates geometry or not
///
class IsValidCell : public vtkm::worklet::WorkletMapField
{
public:
  typedef void ControlSignature(FieldIn<AllTypes> inputCaseInfo, FieldOut<IdType> outputIsValid);
  typedef _2 ExecutionSignature(_1);
  typedef _1 InputDomain;

  VTKM_CONT_EXPORT
  IsValidCell() { };

  VTKM_EXEC_EXPORT
  vtkm::Id operator()(const vtkm::Pair<vtkm::Id, vtkm::Id> &caseInfo) const
  {
    return caseInfo.second != 0;
  }
};

/// \brief Return a value of the specified field of pairs at the given index, similar to a permutation iterator
///
class Permute : public vtkm::worklet::WorkletMapField
{
public:
  typedef void ControlSignature(FieldIn<IdType> inputIndex, FieldOut<IdType> outputVerticesEnum);
  typedef _2 ExecutionSignature(_1);
  typedef _1 InputDomain;

  typedef PortalTypes< vtkm::Pair< vtkm::Id, vtkm::Id> >::PortalConst InputDomainPortalType;
  InputDomainPortalType inputDomain;

  template<typename T>
  VTKM_CONT_EXPORT
  Permute(const T& inputArrayHandle) :
    inputDomain(inputArrayHandle.PrepareForInput(DeviceAdapter()))
    {

    }

  VTKM_EXEC_EXPORT
  vtkm::Id operator()(const vtkm::Id &index) const
  {
    return this->inputDomain.Get(index).second;
  }
};



/// \brief Compute isosurface vertices, normals, and scalars
///
template <typename FieldType, typename OutputType>
class IsosurfaceFunctorUniformGrid : public vtkm::worklet::WorkletMapField
{
public:
  typedef void ControlSignature(FieldIn<IdType> inputCellId, FieldIn<IdType> outputVertId);
  typedef void ExecutionSignature(_1, _2);
  typedef _1 InputDomain;

  typedef typename PortalTypes< vtkm::Id >::PortalConst IdPortalType;
  IdPortalType outputVerticesEnum, triangleTable;

  typedef typename PortalTypes< vtkm::Pair<vtkm::Id, vtkm::Id> >::PortalConst IdPairPortalType;
  IdPairPortalType caseInfo;

  typedef typename PortalTypes< FieldType >::PortalConst FieldPortalType;
  FieldPortalType field, source;

  typedef typename PortalTypes< OutputType >::Portal ScalarPortalType;
  ScalarPortalType scalars;

  typedef typename PortalTypes< vtkm::Vec<OutputType,3> >::Portal VertexPortalType;
  VertexPortalType vertices;

  const int xdim, ydim, zdim, cellsPerLayer;
  const float isovalue, xmin, ymin, zmin, xmax, ymax, zmax;

  template<typename T, typename U, typename V, typename W, typename X>
  VTKM_CONT_EXPORT
  IsosurfaceFunctorUniformGrid(const float isovalue,
                               const int dims[3],
                               const T & caseInfo,
                               const U & field,
                               const U & source,
                               const V & triangleTable,
                               const W & vertices,
                               const X & scalars):
  isovalue(isovalue),
  xdim(dims[0]), ydim(dims[1]), zdim(dims[2]),
  xmin(-1), ymin(-1), zmin(-1),
  xmax(1), ymax(1), zmax(1),
  caseInfo( caseInfo.PrepareForInput( DeviceAdapter() ) ),
  field( field.PrepareForInput( DeviceAdapter() ) ),
  source( source.PrepareForInput( DeviceAdapter() ) ),
  triangleTable( triangleTable.PrepareForInput( DeviceAdapter() ) ),
  vertices(vertices),
  scalars(scalars),
  cellsPerLayer((xdim-1) * (ydim-1))
  {

  }

  VTKM_EXEC_EXPORT
  void operator()(vtkm::Id inputCellId, vtkm::Id outputVertId) const
  {
    // Get data for this cell
    const int cubeindex    = caseInfo.Get(inputCellId).first;
    const int numVertices  = caseInfo.Get(inputCellId).second;
    const int verticesForEdge[] = { 0, 1, 1, 2, 3, 2, 0, 3,
                              4, 5, 5, 6, 7, 6, 4, 7,
                              0, 4, 1, 5, 2, 6, 3, 7 };

    // Compute 3D indices of this cell
    int x = inputCellId % (xdim - 1);
    int y = (inputCellId / (xdim - 1)) % (ydim - 1);
    int z = inputCellId / cellsPerLayer;

    // Compute indices for the eight vertices of this cell
    int i[8];
    i[0] = x      + y*xdim + z * xdim * ydim;
    i[1] = i[0]   + 1;
    i[2] = i[0]   + 1 + xdim;
    i[3] = i[0]   + xdim;
    i[4] = i[0]   + xdim * ydim;
    i[5] = i[1]   + xdim * ydim;
    i[6] = i[2]   + xdim * ydim;
    i[7] = i[3]   + xdim * ydim;

    // Get the field values at these eight vertices
    float f[8];
    f[0] = this->field.Get(i[0]);
    f[1] = this->field.Get(i[1]);
    f[2] = this->field.Get(i[2]);
    f[3] = this->field.Get(i[3]);
    f[4] = this->field.Get(i[4]);
    f[5] = this->field.Get(i[5]);
    f[6] = this->field.Get(i[6]);
    f[7] = this->field.Get(i[7]);

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
    s[0] = this->source.Get(i[0]);
    s[1] = this->source.Get(i[1]);
    s[2] = this->source.Get(i[2]);
    s[3] = this->source.Get(i[3]);
    s[4] = this->source.Get(i[4]);
    s[5] = this->source.Get(i[5]);
    s[6] = this->source.Get(i[6]);
    s[7] = this->source.Get(i[7]);

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
  }
};


static const vtkm::Float32 ISO_VALUE=0.07;

static void RunvtkmMarchingCubes(int vdims[3],
                                 std::vector<vtkm::Float32>& buffer,
                                 std::string device,
                                 int MAX_NUM_TRIALS,
                                 bool silent=false)
{
  int dims[3] = { vdims[0]-1, vdims[1]-1, vdims[2]-1 };
  int dim3 = dims[0] * dims[1] * dims[2];

  //construct the scheduler that will execute all the worklets
  for(int i=0; i < MAX_NUM_TRIALS; ++i)
    {
    vtkm::cont::Timer<> timer;
    vtkm::cont::ArrayHandle<vtkm::Float32> field = vtkm::cont::make_ArrayHandle(buffer);

    // Set up the Marching Cubes tables
    vtkm::cont::ArrayHandle<vtkm::Id> vertexTableArray = vtkm::cont::make_ArrayHandle(numVerticesTable, 256);
    vtkm::cont::ArrayHandle<vtkm::Id> triangleTableArray = vtkm::cont::make_ArrayHandle(triTable, 256*16);

    // Call the ClassifyCell functor to compute the Marching Cubes case numbers for each cell, and the number of vertices to be generated
    vtkm::cont::ArrayHandleCounting<vtkm::Id> cellCountImplicitArray(0, dim3);
    vtkm::cont::ArrayHandle< vtkm::Pair<vtkm::Id, vtkm::Id> > caseInfoArray;

    typedef ClassifyCell< vtkm::Float32 > CellClassifyFunctor;
    typedef vtkm::worklet::DispatcherMapField< CellClassifyFunctor > ClassifyDispatcher;

    //classify each cell
    CellClassifyFunctor cellClassify(field, vertexTableArray, ISO_VALUE, vdims );
    ClassifyDispatcher classifyCellDispatcher(cellClassify);
    classifyCellDispatcher.Invoke(cellCountImplicitArray, caseInfoArray);

    // Determine which cells are "valid" (i.e., produce geometry), and perform
    //an inclusive scan to get running total of the number of "valid" cells
    vtkm::cont::ArrayHandle<vtkm::Id> validCellIndicesArray;
    {
    vtkm::cont::ArrayHandle<vtkm::Id> validCellEnumArray;
    vtkm::worklet::DispatcherMapField<IsValidCell> isValidCellDispatcher;
    isValidCellDispatcher.Invoke(caseInfoArray, validCellEnumArray);
    vtkm::cont::DeviceAdapterAlgorithm<DeviceAdapter>::ScanInclusive(validCellEnumArray, validCellEnumArray);


    // Return if no cells generate geometry
    unsigned int numValidCells = validCellEnumArray.GetPortalConstControl().Get(validCellEnumArray.GetNumberOfValues() - 1);
    if (numValidCells == 0) return;

    // Use UpperBounds, a "permutation", and an exclusive scan to compute the starting output vertex index for each cell

    vtkm::cont::ArrayHandleCounting<vtkm::Id> validCellCountImplicitArray(0, numValidCells);
    vtkm::cont::DeviceAdapterAlgorithm<DeviceAdapter>::UpperBounds(validCellEnumArray,
                                                                   validCellCountImplicitArray,
                                                                   validCellIndicesArray);
    }

    vtkm::cont::ArrayHandle<vtkm::Id> outputVerticesEnumArray;

    {
    vtkm::cont::ArrayHandle<vtkm::Id> validVerticesArray;
    vtkm::worklet::DispatcherMapField< Permute > permuteDispatcher( (Permute(caseInfoArray)) );
    permuteDispatcher.Invoke(validCellIndicesArray, validVerticesArray);
    vtkm::cont::DeviceAdapterAlgorithm<DeviceAdapter>::ScanExclusive(validVerticesArray,
                                                                     outputVerticesEnumArray);
    }

    // Determine the total number of output vertices
    int numTotalVertices = caseInfoArray.GetPortalConstControl().Get(validCellIndicesArray.GetPortalConstControl().Get(validCellIndicesArray.GetNumberOfValues()-1)).second +
                           outputVerticesEnumArray.GetPortalConstControl().Get(outputVerticesEnumArray.GetNumberOfValues()-1);

    vtkm::cont::ArrayHandle<vtkm::Float32> scalarsArray;
    vtkm::cont::ArrayHandle< vtkm::Vec<vtkm::Float32,3> > verticesArray;



    typedef IsosurfaceFunctorUniformGrid<vtkm::Float32, vtkm::Float32> IsoSurfaceFunctor;
    IsoSurfaceFunctor isosurface(ISO_VALUE,
                                 vdims,
                                 caseInfoArray,
                                 field,
                                 field,
                                 triangleTableArray,
                                 verticesArray.PrepareForOutput(numTotalVertices, DeviceAdapter()),
                                 scalarsArray.PrepareForOutput(numTotalVertices, DeviceAdapter())
                                 );

    vtkm::worklet::DispatcherMapField< IsoSurfaceFunctor > isosurfaceDispatcher(isosurface);

    isosurfaceDispatcher.Invoke(validCellIndicesArray, outputVerticesEnumArray);


    double time = timer.GetElapsedTime();
    if(!silent)
      std::cout << "vtkm," << device << "," << time << "," << i << std::endl;
    }
}

static void RunVTKMarchingCubes(vtkImageData* image, int MAX_NUM_TRIALS)
{
  vtkNew<vtkTrivialProducer> producer;
  producer->SetOutput(image);
  producer->Update();

  for(int i=0; i < MAX_NUM_TRIALS; ++i)
    {

    vtkNew<vtkContourFilter> marching;
    marching->SetInputConnection(producer->GetOutputPort());

    vtkm::cont::Timer<> timer;

    marching->ComputeGradientsOff();
    marching->ComputeNormalsOn();
    marching->ComputeScalarsOn();
    marching->SetNumberOfContours(1);
    marching->SetValue(0, ISO_VALUE);

    marching->Update();

    double time = timer.GetElapsedTime();
    std::cout << "VTK,Serial," << time << "," << i << std::endl;
    }
}
