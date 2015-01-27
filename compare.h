
#include "compare_mc.h"
// #include "compare_thresh.h"

#include <vtkDataArray.h>
#include <vtkImageData.h>
#include <vtkImageResample.h>
#include <vtkNrrdReader.h>
#include <vtkPointData.h>
#include <vtkSmartPointer.h>

#include <iostream>
#include <vector>

static const int NUM_TRIALS = 2;

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


int RunComparison(std::string device, std::string file, int pipeline)
{

  std::vector<vtkm::Float32> buffer;
  double resample_ratio = 0.3; //full data
  vtkSmartPointer< vtkImageData > image = ReadData(buffer, file, resample_ratio);

  //get dims of image data
  int dims[3]; image->GetDimensions(dims);
  std::cout << "data dims are: " << dims[0] << ", " << dims[1] << ", " << dims[2] << std::endl;

  //pipeline 1 is equal to threshold
  if(pipeline <= 1)
  {
    std::cout << "Warming up the machine" << std::endl;
    //warm up the card, whatever algorithm we run first will get a performance
    //hit for the first 10 iterations if we don't run something first
    // RunvtkmThreshold(dims,buffer,device,2,true);

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
    std::cout << "Warming up the machine" << std::endl;
    //warm up the card, whatever algorithm we run first will get a performance
    //hit for the first 10 iterations if we don't run something first
    RunvtkmMarchingCubes(dims,buffer,device,2,true);

    std::cout << "Benchmarking Marching Cubes" << std::endl;

    std::cout << "VTKM,Accelerator,Time,Trial" << std::endl;
    RunvtkmMarchingCubes(dims,buffer,device,NUM_TRIALS);

    if(device == "Serial")
      {
      RunVTKMarchingCubes(image,NUM_TRIALS);
      }
  }

  return 0;
}
