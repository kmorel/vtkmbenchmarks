
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

namespace perf {

//now that the device adapter is included set a global typedef
//that is the chosen device tag
typedef VTKM_DEFAULT_DEVICE_ADAPTER_TAG DeviceAdapter;

static void doLayeredMarchingCubes( int vdims[3],
                                   const vtkm::cont::ArrayHandle<vtkm::Float32>& field,
                                   std::vector< vtkm::cont::ArrayHandle< vtkm::Float32> >& scalarsArray,
                                   std::vector< vtkm::cont::ArrayHandle< vtkm::Vec<vtkm::Float32,3> > >& verticesArray,
                                   int count,
                                   int numberOfLayers)
{
  // Set up the Marching Cubes tables
  vtkm::cont::ArrayHandle<vtkm::Id> vertexTableArray = vtkm::cont::make_ArrayHandle(numVerticesTable, 256);
  vtkm::cont::ArrayHandle<vtkm::Id> triangleTableArray = vtkm::cont::make_ArrayHandle(triTable, 256*16);

  //setup the vector of outputs
  int cellsInLayer = count/numberOfLayers;
  int remainderCells = count%numberOfLayers;

  std::cout << "numberOfLayers: " << numberOfLayers << std::endl;
  for(int i=0; i < numberOfLayers; ++i)
    {
    vtkm::Id startI = cellsInLayer * i;
    vtkm::Id size = cellsInLayer;
    if( remainderCells != 0 && (i+1) == numberOfLayers)
      { size += remainderCells; }

    vtkm::cont::ArrayHandle< vtkm::Id > numOutputTrisPerCell;
    // Call the ClassifyCell functor to compute the Marching Cubes case
    //numbers for each cell, and the number of vertices to be generated
    vtkm::cont::ArrayHandleCounting<vtkm::Id> cellCountImplicitArray(startI, size);
    typedef worklets::ClassifyCellOutputTri< vtkm::Float32, vtkm::Id > CellClassifyFunctor;
    typedef vtkm::worklet::DispatcherMapField< CellClassifyFunctor > ClassifyDispatcher;

    CellClassifyFunctor cellClassify(field, vertexTableArray, ISO_VALUE, vdims );
    ClassifyDispatcher classifyCellDispatcher(cellClassify);
    classifyCellDispatcher.Invoke(cellCountImplicitArray, numOutputTrisPerCell);

    vtkm::cont::ArrayHandle< vtkm::Id > validCellIndicesArray;
    vtkm::cont::ArrayHandle< vtkm::Id > inputCellIterationNumber;

    // numOutputTrisPerCell to start:
    // 0, 3, 0, 1, 0, 0, 2

    //compute the number of valid input cells and those ids
    const vtkm::Id numOutputCells =
      vtkm::cont::DeviceAdapterAlgorithm<DeviceAdapter>::ScanInclusive(numOutputTrisPerCell,
                                                                       numOutputTrisPerCell);
    //go to next layer if no cells generate geometry
    if(numOutputCells == 0)
      { continue; }


    // numOutputTrisPerCell after scan inclusive:
    // 0, 3, 3, 4, 4, 4, 6

    vtkm::cont::ArrayHandleCounting<vtkm::Id> validCellCountImplicitArray(0, numOutputCells);
    vtkm::cont::DeviceAdapterAlgorithm<DeviceAdapter>::UpperBounds(numOutputTrisPerCell,
                                                                   validCellCountImplicitArray,
                                                                   validCellIndicesArray);
    //Valid Cell Indices Array
    // 1, 1, 1, 3, 6, 6

    numOutputTrisPerCell.ReleaseResourcesExecution();
    //compute for each output triangle what iteration of the input cell
    //generates it
    vtkm::cont::DeviceAdapterAlgorithm<DeviceAdapter>::LowerBounds(validCellIndicesArray,
                                                                   validCellIndicesArray,
                                                                   inputCellIterationNumber);

    //iteration cell numbers now look like:
    // 0, 0, 0, 3, 4, 4

    //so in the iso surface, we need to talk current output index - iteration number
    //to get the proper iteration number
    //to get the real input cell id we need to add StartI to each element
    //in valid cell Indices array

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
                                 scalars.PrepareForOutput(numTotalVertices, DeviceAdapter()),
                                 startI
                                 );

    vtkm::worklet::DispatcherMapField< IsoSurfaceFunctor > isosurfaceDispatcher(isosurface);
    isosurfaceDispatcher.Invoke(validCellIndicesArray, inputCellIterationNumber);

    scalarsArray.push_back(scalars);
    verticesArray.push_back(verts);
    }
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

    //because we are generating on a per output tri basis we need to have fatter
    //slices so that we can maximize the accelerator.
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

