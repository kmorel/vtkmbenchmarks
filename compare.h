//=============================================================================
//
//  Copyright (c) Kitware, Inc.
//  All rights reserved.
//  See LICENSE.txt for details.
//
//  This software is distributed WITHOUT ANY WARRANTY; without even
//  the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
//  PURPOSE.  See the above copyright notice for more information.
//
//  Copyright 2012 Sandia Corporation.
//  Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
//  the U.S. Government retains certain rights in this software.
//
//=============================================================================

//marching cubes algorithms
#include "compare_vtkm_mc.h"
#include "compare_vtk_mc.h"

#include "saveAsPly.h"

#include <vtkDataArray.h>
#include <vtkImageData.h>
#include <vtkImageResample.h>
#include <vtkNrrdReader.h>
#include <vtkPointData.h>
#include <vtkSmartPointer.h>

#include <iostream>
#include <vector>

static const int NUM_TRIALS = 10;

static vtkSmartPointer<vtkImageData>
ReadData(std::vector<vtkm::Float32> &buffer, std::string file,  double resampleSize=1.0)
{
  //make sure we are testing float benchmarks only
  assert(sizeof(float) == sizeof(vtkm::Float32));

  std::cout << "loading file: " << file << std::endl;
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


int RunComparison(std::string device,
                  std::string file,
                  std::string writeLoc,
                  int targetNumCores,
                  int maxNumCores,
                  float isoValue,
                  double resampleRatio)
{
  std::vector<vtkm::Float32> buffer;
  vtkSmartPointer< vtkImageData > image = ReadData(buffer, file, resampleRatio);

  //get dims of image data
  int dims[3]; image->GetDimensions(dims);
  std::cout << "data dims are: " << dims[0] << ", " << dims[1] << ", " << dims[2] << std::endl;

/*  std::cout << "vtkSynchronizedTemplates3D,Accelerator,Cores,Time,Trial" << std::endl;
  {
  const int singleCore = 1;
  vtk::RunImageMarchingCubes(image, device,
                             singleCore, maxNumCores, isoValue, NUM_TRIALS);
  }
*/

  std::cout << "vtkmIsoSurfaceUniformGrid,Accelerator,Cores,Time,Trial" << std::endl;
  {
  vtkm::RunIsoSurfaceUniformGrid(buffer, image, device,
                                 targetNumCores, maxNumCores, isoValue, NUM_TRIALS);
  }
  return 0;
}
