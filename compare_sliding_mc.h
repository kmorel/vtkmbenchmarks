
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

namespace slide{

//now that the device adapter is included set a global typedef
//that is the chosen device tag
typedef VTKM_DEFAULT_DEVICE_ADAPTER_TAG DeviceAdapter;

static void doLayeredMarchingCubes( int vdims[3],
                             const vtkm::cont::ArrayHandle<vtkm::Float32>& field,
                             std::vector< vtkm::cont::ArrayHandle<vtkm::Float32> >& scalarsArrays,
                             std::vector< vtkm::cont::ArrayHandle< vtkm::Vec<vtkm::Float32,3> > >& verticesArrays,
                             int count,
                             int numberOfLayers)
{
  // Set up the Marching Cubes tables
  vtkm::cont::ArrayHandle<vtkm::Id> vertexTableArray = vtkm::cont::make_ArrayHandle(numVerticesTable, 256);
  vtkm::cont::ArrayHandle<vtkm::Id> triangleTableArray = vtkm::cont::make_ArrayHandle(triTable, 256*16);

  //setup the vector of outputs
  int cellsInLayer = count/numberOfLayers;
  int remainderCells = count%numberOfLayers;

  scalarsArrays.resize( numberOfLayers );
  verticesArrays.resize( numberOfLayers );
  std::cout << "numberOfLayers: " << numberOfLayers << std::endl;

  for(int i=0; i < numberOfLayers; ++i)
    {
    vtkm::cont::ArrayHandle< vtkm::Id > cellHasOutput;
    vtkm::cont::ArrayHandle< vtkm::Id > numOutputVertsPerCell;

    vtkm::Id startI = cellsInLayer * i;
    vtkm::Id size = cellsInLayer;
    if( remainderCells != 0 && (i+1) == numberOfLayers)
      { size += remainderCells; }

    // Call the ClassifyCell functor to compute the Marching Cubes case
    //numbers for each cell, and the number of vertices to be generated
    vtkm::cont::ArrayHandleCounting<vtkm::Id> cellCountImplicitArray(startI, size);
    typedef worklets::ClassifyCell< vtkm::Float32, vtkm::Id, vtkm::Id > CellClassifyFunctor;
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
    if (numTotalVertices == 0) continue;

    //compute the original cell ids that will produce output. We do this by
    //figuring out for each output cell value what input cell generates it
    vtkm::cont::ArrayHandleCounting<vtkm::Id> validCellCountImplicitArray(0, numValidInputCells);
    vtkm::cont::DeviceAdapterAlgorithm<DeviceAdapter>::UpperBounds(cellHasOutput,
                                                                   validCellCountImplicitArray,
                                                                   validCellIndicesArray);
    }

    cellHasOutput.ReleaseResourcesExecution();
    numOutputVertsPerCell.ReleaseResourcesExecution();

    typedef worklets::IsosurfaceFunctorUniformGrid<vtkm::Float32, vtkm::Float32> IsoSurfaceFunctor;

    vtkm::cont::ArrayHandle< vtkm::Float32 > scalars;
    vtkm::cont::ArrayHandle< vtkm::Vec<vtkm::Float32,3> > verts;

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

    isosurfaceDispatcher.Invoke(validCellIndicesArray, outputVerticesLocArray);

    verticesArrays[i] = verts;
    scalarsArrays[i] = scalars;
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

    const int numberOfSlices[9] = {64, 32, 25, 16, 8, 5, 4, 3, 1};
    const int widthOfEachSlice[9] = { (dims[2]/64),
                                      (dims[2]/32),
                                      (dims[2]/25),
                                      (dims[2]/16),
                                      (dims[2]/8 ),
                                      (dims[2]/5 ),
                                      (dims[2]/4 ),
                                      (dims[2]/3 ),
                                      10000000
                                     };
    //classify each cell, and merge classification of cells based on if
    //we can fuse 3 or 4 cells at a time
    std::vector< vtkm::cont::ArrayHandle<vtkm::Float32> > scalarsArrays;
    std::vector< vtkm::cont::ArrayHandle< vtkm::Vec<vtkm::Float32,3> > > verticesArrays;

    // find the best number of slices to subdivide this grid by
    for(int s=0; s < 9; ++s)
      {
      if( widthOfEachSlice[s] >= 16)
        {
        doLayeredMarchingCubes( vdims, field, scalarsArrays, verticesArrays,
                                dim3, numberOfSlices[s]);
        break;
        }
      }

    double time = timer.GetElapsedTime();
    std::size_t numCells = 0;
    for(std::size_t layerIndex=0; layerIndex < scalarsArrays.size(); ++layerIndex)
      {
      numCells += (scalarsArrays[layerIndex].GetNumberOfValues()/3);
      }

    if(!silent)
      {
      std::cout << "num cells: " << numCells << std::endl;
      std::cout << "vtkm," << device << "," << time << "," << trial << std::endl;
      }
    }

}

}