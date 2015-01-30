
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

namespace per_tri {

//now that the device adapter is included set a global typedef
//that is the chosen device tag
typedef VTKM_DEFAULT_DEVICE_ADAPTER_TAG DeviceAdapter;

static void doMarchingCubes( int vdims[3],
                             const vtkm::cont::ArrayHandle<vtkm::Float32>& field,
                             std::vector< vtkm::cont::ArrayHandle< vtkm::Float32> >& scalarsArray,
                             std::vector< vtkm::cont::ArrayHandle< vtkm::Vec<vtkm::Float32,3> > >& verticesArray,
                             int count)
{
  // Set up the Marching Cubes tables
  vtkm::cont::ArrayHandle<vtkm::Id> vertexTableArray = vtkm::cont::make_ArrayHandle(numVerticesTable, 256);
  vtkm::cont::ArrayHandle<vtkm::Id> triangleTableArray = vtkm::cont::make_ArrayHandle(triTable, 256*16);

  vtkm::cont::ArrayHandle< vtkm::Id > cellHasOutput;
  vtkm::cont::ArrayHandle< vtkm::UInt8 > numOutputVertsPerCell;

  // Call the ClassifyCell functor to compute the Marching Cubes case
  //numbers for each cell, and the number of vertices to be generated
  vtkm::cont::ArrayHandleCounting<vtkm::Id> cellCountImplicitArray(0, count);
  typedef worklets::ClassifyCell< vtkm::Float32, vtkm::Id, vtkm::UInt8 > CellClassifyFunctor;
  typedef vtkm::worklet::DispatcherMapField< CellClassifyFunctor > ClassifyDispatcher;

  CellClassifyFunctor cellClassify(field, vertexTableArray, ISO_VALUE, vdims );
  ClassifyDispatcher classifyCellDispatcher(cellClassify);
  classifyCellDispatcher.Invoke(cellCountImplicitArray, cellHasOutput, numOutputVertsPerCell);

  for(int i=0; i < 5; ++i)
    {
    vtkm::cont::ArrayHandle< vtkm::Id > cellOutputIncSum;
    vtkm::cont::ArrayHandle< vtkm::Id > validCellIndicesArray;

    const vtkm::Id numValidInputCells =
      vtkm::cont::DeviceAdapterAlgorithm<DeviceAdapter>::ScanInclusive(cellHasOutput,
                                                                       cellOutputIncSum);

    vtkm::cont::ArrayHandleCounting<vtkm::Id> validCellCountImplicitArray(0, numValidInputCells);
    vtkm::cont::DeviceAdapterAlgorithm<DeviceAdapter>::UpperBounds(cellOutputIncSum,
                                                                   validCellCountImplicitArray,
                                                                   validCellIndicesArray);

    //terminate if no cells has triangles left
    if(numValidInputCells == 0)
      { break; }

    //allocate
    //generate a single tri per cell
    vtkm::cont::ArrayHandle< vtkm::Float32 > scalars;
    vtkm::cont::ArrayHandle< vtkm::Vec<vtkm::Float32,3> > verts;

    const vtkm::Id numTotalVertices = numValidInputCells * 3;
    typedef worklets::IsosurfaceSingleTri<vtkm::Float32, vtkm::Float32> IsoSurfaceFunctor;
    IsoSurfaceFunctor isosurface(i,
                                 ISO_VALUE,
                                 vdims,
                                 field,
                                 field,
                                 vertexTableArray,
                                 triangleTableArray,
                                 verts.PrepareForOutput(numTotalVertices, DeviceAdapter()),
                                 scalars.PrepareForOutput(numTotalVertices, DeviceAdapter())
                                 );

    vtkm::worklet::DispatcherMapField< IsoSurfaceFunctor > isosurfaceDispatcher(isosurface);
    isosurfaceDispatcher.Invoke(validCellIndicesArray, validCellCountImplicitArray);

    scalarsArray.push_back(scalars);
    verticesArray.push_back(verts);

    //decrement the count in numOutputVertsPerCell, and cellHasOutput
    typedef worklets::DecrementCounts DecCountFunctor;
    typedef vtkm::worklet::DispatcherMapField< DecCountFunctor > DecDispatcher;

    DecCountFunctor decFunc(numValidInputCells, cellHasOutput, numOutputVertsPerCell);
    DecDispatcher decCellCountDispatcher( decFunc );
    decCellCountDispatcher.Invoke(validCellIndicesArray);
    }

};

static void RunMarchingCubes(int vdims[3],
                                 std::vector<vtkm::Float32>& buffer,
                                 std::string device,
                                 int MAX_NUM_TRIALS,
                                 bool silent=false)
{
  int dims[3] = { vdims[0]-1, vdims[1]-1, vdims[2]-1 };
  int dim3 = dims[0] * dims[1] * dims[2];


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

    std::vector< vtkm::cont::ArrayHandle< vtkm::Float32 > > scalarsArrays;
    std::vector< vtkm::cont::ArrayHandle< vtkm::Vec<vtkm::Float32,3> > > verticesArrays;
    doMarchingCubes( vdims, field, scalarsArrays, verticesArrays, dim3);

    double time = timer.GetElapsedTime();

    std::size_t numCells = 0;
    for(std::size_t layerIndex=0; layerIndex < scalarsArrays.size(); ++layerIndex)
      {
      numCells += (scalarsArrays[layerIndex].GetNumberOfValues()/3);
      }

    if(!silent)
      {
      std::cout << "num cells: " << numCells  << std::endl;
      std::cout << "vtkm," << device << "," << time << "," << trial << std::endl;
      }
    }

}
}

