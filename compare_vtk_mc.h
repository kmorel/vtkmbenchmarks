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

  for(int i=0; i < MAX_NUM_TRIALS; ++i)
    {

    vtkNew<vtkMarchingCubes> syncTemplates;
    syncTemplates->SetInputConnection(producer->GetOutputPort());

    vtkm::cont::Timer<> timer;

    syncTemplates->ComputeGradientsOff();
    syncTemplates->ComputeNormalsOn();
    syncTemplates->ComputeScalarsOn();
    syncTemplates->SetNumberOfContours(1);
    syncTemplates->SetValue(0, isoValue);

    syncTemplates->Update();

    double time = timer.GetElapsedTime();
    std::cout << "vtkMarchingCubes," << device << "," <<  numCores << "," << time << "," << i << std::endl;
    }
}

}