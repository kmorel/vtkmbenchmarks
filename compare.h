
static const float ISO_VALUE=0.07;

#include "isosurface.h"
#include "worklets.h"
#include "fusedWorklets.h"

//marching cubes algorithms
#include "compare_classic_mc.h"
#include "compare_lowmem_mc.h"
#include "compare_per_tri_out_mc.h"
#include "compare_perf_mc.h"
#include "compare_sliding_mc.h"
#include "compare_vtk_mc.h"

//threshold algorithms
// #include "compare_thresh.h"

#include <vtkDataArray.h>
#include <vtkImageData.h>
#include <vtkImageResample.h>
#include <vtkNrrdReader.h>
#include <vtkPointData.h>
#include <vtkSmartPointer.h>

#include <iostream>
#include <vector>

static const int NUM_TRIALS = 8;


static vtkSmartPointer<vtkImageData>
ReadData(std::vector<vtkm::Float32> &buffer, std::string file,  double resampleSize=1.0)
{
  //make sure we are testing float benchmarks only
  assert(sizeof(float) == sizeof(vtkm::Float32));

  std::cout << "reading file: " << file << std::endl;
  vtkNew<vtkNrrdReader> reader;
  reader->SetFileName(file.c_str());
  reader->Update();

  //re-sample the dataset
  vtkNew<vtkImageResample> resample;
  resample->SetInputConnection(reader->GetOutputPort());
  resample->SetAxisMagnificationFactor(0,resampleSize);
  resample->SetAxisMagnificationFactor(1,resampleSize);
  resample->SetAxisMagnificationFactor(2,resampleSize);

  resample->Update();

  //take ref
  vtkSmartPointer<vtkImageData> image = vtkSmartPointer<vtkImageData>::New();
  vtkImageData *newImageData = vtkImageData::SafeDownCast(resample->GetOutputDataObject(0));
  image.TakeReference( newImageData );
  image->Register(NULL);

  //now set the buffer
  vtkDataArray *newData = image->GetPointData()->GetScalars();
  vtkm::Float32* rawBuffer = reinterpret_cast<vtkm::Float32*>( newData->GetVoidPointer(0) );
  buffer.resize( newData->GetNumberOfTuples() );
  std::copy(rawBuffer, rawBuffer + newData->GetNumberOfTuples(), buffer.begin() );

  return image;
}


int RunComparison(std::string device, std::string file, int pipeline, double resample_ratio)
{
  std::vector<vtkm::Float32> buffer;
  vtkSmartPointer< vtkImageData > image = ReadData(buffer, file, resample_ratio);

  //get dims of image data
  int dims[3]; image->GetDimensions(dims);
  std::cout << "data dims are: " << dims[0] << ", " << dims[1] << ", " << dims[2] << std::endl;

  //pipeline 1 is equal to threshold
  if(pipeline <= 1)
  {
    //print out header of csv
    std::cout << "Benchmarking Threshold" << std::endl;

    std::cout << "VTKM,Accelerator,Time,Trial" << std::endl;
    // RunvtkmThreshold(dims,buffer,device,NUM_TRIALS);
    if(device == "Serial")
      {
      std::cout << "VTK,Accelerator,Time,Trial" << std::endl;
      // RunVTKThreshold(image,NUM_TRIALS);
      }
  }
  else //marching cubes
  {
    std::cout << "Benchmarking Marching Cubes" << std::endl;

    std::cout << "VTKM Classic,Accelerator,Time,Trial" << std::endl;
    //Run the basic marching cubes which classifies 1 cell at a time
    //and than writes out all geometry for each input cell at a time
    try{ mc::RunMarchingCubes(dims,buffer,device,NUM_TRIALS); } catch(...) {}

    std::cout << "VTKM Fused Classic,Accelerator,Time,Trial" << std::endl;
    //Run the basic marching cubes which classifies 1/2/3/4 cells at a time
    //and than writes out all geometry for those combined cells at the same time
    try{ mc::RunFusedMarchingCubes(dims,buffer,device,NUM_TRIALS); } catch(...) {}

    std::cout << "VTKM Per Tri Output,Accelerator,Time,Trial" << std::endl;
    //Run the basic marching cubes which classifies 1 cell at a time
    //and than generate the geometry on a per output triangle basis, so we
    //have coalesced writes.
    try{ per_tri::RunMarchingCubes(dims,buffer,device,NUM_TRIALS); } catch(...) {}

    std::cout << "VTKM Sliding Window,Accelerator,Time,Trial" << std::endl;
    //Run the a sliding window marching cubes which classifies 1 cell at a time
    //the sliding window is along the Z axis, and allows us to generate a subset
    //of the triangle at a time, but requires all the input data to be uploaded.
    //The primary benifit of this approach is a reduction of memory, since the
    //number of cells we are walking is smaller
    try{ slide::RunMarchingCubes(dims,buffer,device,NUM_TRIALS); } catch(...) {}

    // std::cout << "VTKM Low Mem Inclusive Scan,Accelerator,Time,Trial" << std::endl;
    //Run the basic marching cubes which classifies 1/3/4 cells at a time
    //but use vtkm::uint8 to store the counts, and than use a modified histo
    //pyramid and inclusive scan approach to reduce the size of the lookup tables.
    // low_mem::RunMarchingCubes(dims,buffer,device,NUM_TRIALS);

    std::cout << "VTKM SuperPerf,Accelerator,Time,Trial" << std::endl;
    // Currently Combines the per tri output and the sliding window.
    // In future will need to add histo pyramid for a super fast, super low mem version
    try{  perf::RunMarchingCubes(dims,buffer,device,NUM_TRIALS); } catch(...) {}

    if(device == "Serial")
      {
      std::cout << "VTK Contour Filter,Accelerator,Time,Trial" << std::endl;
      vtk::RunContourFilter(image,NUM_TRIALS);

#ifdef VTK_HAS_FLYING_EDGES
      std::cout << "VTK Contour Filter,Accelerator,Time,Trial" << std::endl;
      vtk::RunFlyingEdges(image,NUM_TRIALS);
#endif
      }
  }

  return 0;
}
