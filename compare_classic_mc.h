
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

namespace mc {

//now that the device adapter is included set a global typedef
//that is the chosen device tag
typedef VTKM_DEFAULT_DEVICE_ADAPTER_TAG DeviceAdapter;

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
  typedef worklets::FusedClassifyCell< vtkm::Float32, vtkm::Id, vtkm::Id, NumCellsToFuse > CellClassifyFunctor;
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

  typedef worklets::IsosurfaceFusedUniformGridFunctor<vtkm::Float32, vtkm::Float32, NumCellsToFuse> IsoSurfaceFunctor;
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

    vtkm::cont::ArrayHandle<vtkm::Float32> scalarsArray;
    vtkm::cont::ArrayHandle< vtkm::Vec<vtkm::Float32,3> > verticesArray;
    doMarchingCubes<1>( vdims, field, scalarsArray, verticesArray, dim3);

    double time = timer.GetElapsedTime();
    if(!silent)
      {
      std::cout << "num cells: " << (scalarsArray.GetNumberOfValues()/3)  << std::endl;
      std::cout << "vtkm," << device << "," << time << "," << trial << std::endl;
      }
    }

}

static void RunFusedMarchingCubes(int vdims[3],
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

    const bool fuse4Cells = (dims[0]%4 == 0);
    const bool fuse3Cells = (dims[0]%3 == 0);
    const bool fuse2Cells = (dims[0]%2 == 0);

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
    else if(fuse2Cells)
      {
      doMarchingCubes<2>( vdims, field, scalarsArray, verticesArray, dim3/2);
      }
    else
      {
      doMarchingCubes<1>( vdims, field, scalarsArray, verticesArray, dim3);
      }

    double time = timer.GetElapsedTime();
    if(!silent)
      {
      std::cout << "fuse4Cells: " << fuse4Cells << " fuse3Cells: " << fuse3Cells << " fuse2Cells: " << fuse2Cells << std::endl;
      std::cout << "num cells: " << (scalarsArray.GetNumberOfValues()/3)  << std::endl;
      std::cout << "vtkm," << device << "," << time << "," << trial << std::endl;
      }
    }

}
}

