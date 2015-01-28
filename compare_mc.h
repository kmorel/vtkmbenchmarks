
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
template <typename FieldType, int NumCellsToFuse>
class ClassifyCell : public vtkm::worklet::WorkletMapField
{
public:
  typedef void ControlSignature(FieldIn<IdType> inputCellId, FieldOut<AllTypes> hasOutput, FieldOut<AllTypes> numCellsOut);
  typedef _3 ExecutionSignature(_1, _2);
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
  ClassifyCell(const T& pointHandle,
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
  vtkm::Id operator()(vtkm::Id firstCellId, vtkm::Id& hasOutput) const
  {
    vtkm::Id vertCount=0;
    for(int i=0; i < NumCellsToFuse; i++)
      {
      const vtkm::Id cellId = (firstCellId * NumCellsToFuse) + i;

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

      vertCount += this->vertTable.Get(cubeindex);
      }

    // Return the number of triangles this case generates
    hasOutput = (vertCount == 0) ? 0 : 1;
    return vertCount;
  }
};


/// \brief Compute isosurface vertices, normals, and scalars
///
template <typename FieldType, typename OutputType, int NumCellsToFuse>
class IsosurfaceFunctorUniformGrid : public vtkm::worklet::WorkletMapField
{
public:
  typedef void ControlSignature(FieldIn<IdType> inputCellId, FieldIn<IdType> outputVertId);
  typedef void ExecutionSignature(_1, _2);
  typedef _1 InputDomain;

  typedef typename PortalTypes< vtkm::Id >::PortalConst IdPortalType;
  IdPortalType outputVerticesEnum, triangleTable, vertexTable;

  typedef typename PortalTypes< vtkm::Pair<vtkm::Id, vtkm::Id> >::PortalConst IdPairPortalType;
  IdPairPortalType caseInfo;

  typedef typename PortalTypes< FieldType >::PortalConst FieldPortalType;
  FieldPortalType field, source;

  typedef typename PortalTypes< OutputType >::Portal ScalarPortalType;
  ScalarPortalType scalars;

  typedef typename PortalTypes< vtkm::Vec<OutputType,3> >::Portal VertexPortalType;
  VertexPortalType vertices;

  const int xdim, ydim, zdim, cellsPerLayer, pointsPerLayer;
  const float isovalue, xmin, ymin, zmin, xmax, ymax, zmax;

  template<typename U, typename V, typename W, typename X>
  VTKM_CONT_EXPORT
  IsosurfaceFunctorUniformGrid(const float isovalue,
                               const int dims[3],
                               const U & field,
                               const U & source,
                               const V & vertexTable,
                               const V & triangleTable,
                               const W & vertices,
                               const X & scalars):
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
  pointsPerLayer (xdim*ydim)
  {

  }

  VTKM_EXEC_EXPORT
  void operator()(vtkm::Id firstInputCellId, vtkm::Id firstOutputVertId) const
  {
    // Get data for this cell
    const int verticesForEdge[] = { 0, 1, 1, 2, 3, 2, 0, 3,
                                    4, 5, 5, 6, 7, 6, 4, 7,
                                    0, 4, 1, 5, 2, 6, 3, 7 };

    vtkm::Id outputVertId = firstOutputVertId;
    for(int i=0; i < NumCellsToFuse; i++)
      {
      const vtkm::Id inputCellId = (firstInputCellId * NumCellsToFuse) + i;


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


static const vtkm::Float32 ISO_VALUE=0.07;


template< int NumCellsToFuse>
static void doMarchingCubes( int vdims[3],
                             const vtkm::cont::ArrayHandle<vtkm::Float32>& field,
                             vtkm::cont::ArrayHandle<vtkm::Float32> scalarsArray,
                             vtkm::cont::ArrayHandle< vtkm::Vec<vtkm::Float32,3> > verticesArray,
                             int count)
{
  // Set up the Marching Cubes tables
  vtkm::cont::ArrayHandle<vtkm::Id> vertexTableArray = vtkm::cont::make_ArrayHandle(numVerticesTable, 256);
  vtkm::cont::ArrayHandle<vtkm::Id> triangleTableArray = vtkm::cont::make_ArrayHandle(triTable, 256*16);

  vtkm::cont::ArrayHandle< vtkm::Id > cellHasOutput;
  vtkm::cont::ArrayHandle< vtkm::Id > numOutputVertsPerCell;

  // Call the ClassifyCell functor to compute the Marching Cubes case
  //numbers for each cell, and the number of vertices to be generated
  vtkm::cont::ArrayHandleCounting<vtkm::Id> cellCountImplicitArray(0, count);
  typedef ClassifyCell< vtkm::Float32, NumCellsToFuse > CellClassifyFunctor;
  typedef vtkm::worklet::DispatcherMapField< CellClassifyFunctor > ClassifyDispatcher;

  CellClassifyFunctor cellClassify(field, vertexTableArray, ISO_VALUE, vdims );
  ClassifyDispatcher classifyCellDispatcher(cellClassify);
  classifyCellDispatcher.Invoke(cellCountImplicitArray, cellHasOutput, numOutputVertsPerCell);

  // Determine which cells are "valid" (i.e., produce geometry), and perform
  // an inclusive scan to get running total of the number of "valid" cells
  vtkm::cont::ArrayHandle<vtkm::Id> validCellIndicesArray;
  vtkm::cont::ArrayHandle<vtkm::Id> outputVerticesLocArray;
  int numTotalVertices = 0;
  {
  const vtkm::Id numValidInputCells =
      vtkm::cont::DeviceAdapterAlgorithm<DeviceAdapter>::ScanInclusive(cellHasOutput,
                                                                       cellHasOutput);
  numTotalVertices =
      vtkm::cont::DeviceAdapterAlgorithm<DeviceAdapter>::ScanInclusive(numOutputVertsPerCell,
                                                                       outputVerticesLocArray);
  // Return if no cells generate geometry
  if (numTotalVertices == 0) return;

  //compute the original cell ids that will produce output. We do this by
  //figuring out for each output cell value what input cell generates it
  vtkm::cont::ArrayHandleCounting<vtkm::Id> validCellCountImplicitArray(0, numValidInputCells);
  vtkm::cont::DeviceAdapterAlgorithm<DeviceAdapter>::UpperBounds(cellHasOutput,
                                                                 validCellCountImplicitArray,
                                                                 validCellIndicesArray);
  }

  cellHasOutput.ReleaseResourcesExecution();
  numOutputVertsPerCell.ReleaseResourcesExecution();

  typedef IsosurfaceFunctorUniformGrid<vtkm::Float32, vtkm::Float32, NumCellsToFuse> IsoSurfaceFunctor;
  IsoSurfaceFunctor isosurface(ISO_VALUE,
                               vdims,
                               field,
                               field,
                               vertexTableArray,
                               triangleTableArray,
                               verticesArray.PrepareForOutput(numTotalVertices, DeviceAdapter()),
                               scalarsArray.PrepareForOutput(numTotalVertices, DeviceAdapter())
                               );

  vtkm::worklet::DispatcherMapField< IsoSurfaceFunctor > isosurfaceDispatcher(isosurface);

  isosurfaceDispatcher.Invoke(validCellIndicesArray, outputVerticesLocArray);
};

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
    const bool fuse4Cells = (dim3%4 == 0);
    const bool fuse3Cells = (dim3%3 == 0);

    //setup the iso field to contour
    vtkm::cont::ArrayHandle<vtkm::Float32> field = vtkm::cont::make_ArrayHandle(buffer);

    //classify each cell, and merge classification of cells based on if
    //we can fuse 3 or 4 cells at a time
    vtkm::cont::ArrayHandle<vtkm::Float32> scalarsArray;
    vtkm::cont::ArrayHandle< vtkm::Vec<vtkm::Float32,3> > verticesArray;
    if(fuse4Cells)
      {
      doMarchingCubes<4>( vdims, field, scalarsArray, verticesArray, dim3/4);
      }
    else if(fuse3Cells)
      {
      doMarchingCubes<3>( vdims, field, scalarsArray, verticesArray, dim3/3);
      }
    else
      {
      doMarchingCubes<1>( vdims, field, scalarsArray, verticesArray, dim3);
      }

    double time = timer.GetElapsedTime();
    if(!silent)
      {
      std::cout << (scalarsArray.GetNumberOfValues() / 3) << std::endl;
      std::cout << "vtkm," << device << "," << time << "," << i << std::endl;
      }
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

    std::cout << marching->GetOutput()->GetNumberOfCells() << std::endl;
    std::cout << "VTK,Serial," << time << "," << i << std::endl;
    }
}
