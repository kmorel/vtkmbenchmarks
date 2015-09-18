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

#include <vtkMarchingCubes.h>
#include <vtkImageData.h>
#include <vtkNew.h>
#include <vtkSmartPointer.h>
#include <vtkTrivialProducer.h>
#include <vtkNonMergingPointLocator.h>

#include <vtkm/cont/Timer.h>

namespace vtk
{
static void RunImageMarchingCubes( vtkImageData* image,
                                   const std::string& device,
                                   int numCores,
                                   int maxNumCores,
                                   float isoValue,
                                   int MAX_NUM_TRIALS)
{
  vtkNew<vtkTrivialProducer> producer;
  producer->SetOutput(image);
  producer->Update();

  vtkm::cont::Timer<> timer;
  std::vector<double> samples;
  samples.reserve(MAX_NUM_TRIALS);

  vtkNew<vtkMarchingCubes> syncTemplates;
  syncTemplates->SetInputConnection(producer->GetOutputPort());

  vtkNonMergingPointLocator* simpleLocator = vtkNonMergingPointLocator::New();
  syncTemplates->SetLocator(simpleLocator);

  syncTemplates->ComputeGradientsOff();
  syncTemplates->ComputeNormalsOn();
  syncTemplates->ComputeScalarsOff();
  syncTemplates->SetNumberOfContours(1);

  timer.Reset();
    
  for(int i=0; i < MAX_NUM_TRIALS; ++i)
  {
    syncTemplates->SetValue(0, isoValue);
    syncTemplates->Update();

    vtkPolyData* output = syncTemplates->GetOutput();
    std::cout << (syncTemplates->GetValue(0)) << " " << output->GetNumberOfPoints() << std::endl;

    isoValue += 0.005f;
  }

  samples.push_back(timer.GetElapsedTime());

  std::sort(samples.begin(), samples.end());
  stats::Winsorize(samples, 5.0);
  std::cout << "Benchmark \'VTK Isosurface\' results:\n"
        << "\tmedian = " << stats::PercentileValue(samples, 50.0) << "s\n"
        << "\tmedian abs dev = " << stats::MedianAbsDeviation(samples) << "s\n"
        << "\tmean = " << stats::Mean(samples) << "s\n"
        << "\tstd dev = " << stats::StandardDeviation(samples) << "s\n"
        << "\tmin = " << samples.front() << "s\n"
        << "\tmax = " << samples.back() << "s\n"
        << "\t# of runs = " << samples.size() << "\n";

}

}
