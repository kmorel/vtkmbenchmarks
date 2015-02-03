
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

  vtkm::cont::ArrayHandle< vtkm::Id > numOutputTrisPerCell;
  // Call the ClassifyCell functor to compute the Marching Cubes case
  //numbers for each cell, and the number of vertices to be generated
  vtkm::cont::ArrayHandleCounting<vtkm::Id> cellCountImplicitArray(0, count);
  typedef worklets::ClassifyCellOutputTri< vtkm::Float32, vtkm::Id > CellClassifyFunctor;
  typedef vtkm::worklet::DispatcherMapField< CellClassifyFunctor > ClassifyDispatcher;

  CellClassifyFunctor cellClassify(field, vertexTableArray, ISO_VALUE, vdims );
  ClassifyDispatcher classifyCellDispatcher(cellClassify);
  classifyCellDispatcher.Invoke(cellCountImplicitArray, numOutputTrisPerCell);

  vtkm::cont::ArrayHandle< vtkm::Id > validCellIndicesArray;
  vtkm::cont::ArrayHandle< vtkm::Id > inputCellIterationNumber;

  //compute the number of valid input cells and those ids
  const vtkm::Id numOutputCells =
    vtkm::cont::DeviceAdapterAlgorithm<DeviceAdapter>::ScanInclusive(numOutputTrisPerCell,
                                                                     numOutputTrisPerCell);

  //terminate if no cells has triangles left
  if(numOutputCells == 0)
    { return; }

  vtkm::cont::ArrayHandleCounting<vtkm::Id> validCellCountImplicitArray(0, numOutputCells);
  vtkm::cont::DeviceAdapterAlgorithm<DeviceAdapter>::UpperBounds(numOutputTrisPerCell,
                                                                 validCellCountImplicitArray,
                                                                 validCellIndicesArray);

  numOutputTrisPerCell.ReleaseResourcesExecution();

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
  vtkm::cont::ArrayHandle< vtkm::Float32 > scalars;
  vtkm::cont::ArrayHandle< vtkm::Vec<vtkm::Float32,3> > verts;

  const vtkm::Id numTotalVertices = numOutputCells * 3;
  typedef worklets::IsosurfaceSingleTri<vtkm::Float32, vtkm::Float32> IsoSurfaceFunctor;
  IsoSurfaceFunctor isosurface(ISO_VALUE,
                               vdims,
                               field,
                               field,
                               vertexTableArray,
                               triangleTableArray,
                               verts.PrepareForOutput(numTotalVertices, DeviceAdapter()),
                               scalars.PrepareForOutput(numTotalVertices, DeviceAdapter())
                               );

  vtkm::worklet::DispatcherMapField< IsoSurfaceFunctor > isosurfaceDispatcher(isosurface);
  isosurfaceDispatcher.Invoke(validCellIndicesArray, inputCellIterationNumber);

  scalarsArray.push_back(scalars);
  verticesArray.push_back(verts);
}


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

