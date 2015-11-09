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
#ifdef PISTON_ENABLED
#include "compare_piston_mc.h"
#endif

#include "saveAsPly.h"

#include <vtkDataArray.h>
#include <vtkImageData.h>
#include <vtkImageResample.h>
#include <vtkNew.h>
#include <vtkNrrdReader.h>
#include <vtkPointData.h>
#include <vtkSmartPointer.h>

#include <iostream>
#include <vector>

#include "ArgumentsParser.h"

static const int NUM_TRIALS = 10;

static vtkSmartPointer<vtkImageData>
ReadData(std::vector<vtkm::Float32> &buffer,
         std::string file,
         double resampleSize=1.0)
{
  //make sure we are testing float benchmarks only
  assert(sizeof(float) == sizeof(vtkm::Float32));

  vtkSmartPointer<vtkImageData> image;

  std::cout << "# Loading file: " << file << " " << resampleSize << std::endl;
  vtkNew<vtkNrrdReader> reader;
  reader->SetFileName(file.c_str());
  reader->Update();

  if (resampleSize == 1.0)
  {
    reader->Update();

    //take ref
    vtkImageData *newImageData =
        vtkImageData::SafeDownCast(reader->GetOutputDataObject(0));
    image.TakeReference( newImageData );
    image->Register(NULL);
  }
  else
  {
    //re-sample the dataset
    vtkNew<vtkImageResample> resample;
    resample->SetInputConnection(reader->GetOutputPort());
    resample->SetAxisMagnificationFactor(0,resampleSize);
    resample->SetAxisMagnificationFactor(1,resampleSize);
    resample->SetAxisMagnificationFactor(2,resampleSize);

    resample->Update();

    //take ref
    vtkImageData *newImageData =
        vtkImageData::SafeDownCast(resample->GetOutputDataObject(0));
    image.TakeReference( newImageData );
    image->Register(NULL);
  }

  //now set the buffer
  vtkDataArray *newData = image->GetPointData()->GetScalars();
  vtkm::Float32* rawBuffer = reinterpret_cast<vtkm::Float32*>( newData->GetVoidPointer(0) );
  buffer.resize( newData->GetNumberOfTuples() );
  std::copy(rawBuffer, rawBuffer + newData->GetNumberOfTuples(), buffer.begin() );

  return image;
}


int RunComparison(std::string device,
                  const vtkm::testing::ArgumentsParser &arguments)
{
  std::vector<vtkm::Float32> buffer;
  vtkSmartPointer< vtkImageData > image =
      ReadData(buffer, arguments.file(), arguments.ratio());

  std::cout << "# Filename: " << arguments.file() << std::endl;
  int dims[3]; image->GetDimensions(dims);
  std::cout << "# Data dims are: " << dims[0] << ", " << dims[1] << ", " << dims[2] << std::endl;

  vtkm::RunIsoSurfaceUniformGrid(buffer, image, device, arguments);

#ifdef PISTON_ENABLED
  piston::RunIsoSurfaceUniformGrid(buffer, image, device, arguments);
#endif

  return 0;
}
