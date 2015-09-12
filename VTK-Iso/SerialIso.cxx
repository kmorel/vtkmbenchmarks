/*=========================================================================

  Program:   Visualization Toolkit
  Module:    ParallelIso.cxx

  Copyright (c) Ken Martin, Will Schroeder, Bill Lorensen
  All rights reserved.
  See Copyright.txt or http://www.kitware.com/Copyright.htm for details.

     This software is distributed WITHOUT ANY WARRANTY; without even
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
     PURPOSE.  See the above copyright notice for more information.

=========================================================================*/
// This example demonstrates the use of data parallelism in VTK. The
// pipeline ( vtkMPIImageReader -> vtkContourFilter -> vtkElevationFilter )
// is created in parallel and each process is assigned 1 piece to process.
// All satellite processes send the result to the first process which
// collects and renders them.

#include <string>
#include <algorithm>
#ifdef _WIN32
#include <sys/timeb.h>
#include <sys/types.h>
#else //!_WIN32
#include <limits.h>
#include <sys/time.h>
#include <unistd.h>
#endif

#include "vtkContourFilter.h"
#include "vtkDataSet.h"
#include "vtkMath.h"
#include "vtkPolyData.h"
#include "vtkTestUtilities.h"
#include "vtkImageData.h"

#include "Stats.h"

// Uncomment to display the rendered isosurface when running with -I
// TODO: I'm not sure of good values to pick here
static const float ISO_START=16;
static const float ISO_STEP=-1;
static const int ISO_NUM=15;
static const size_t DIMENSIONS = 128;


void GenerateData(vtkImageData *data){
  const int dim = DIMENSIONS;
  const int vdim = dim + 1;
  const float mins[3] = { -1, -1, -1 };
  const float maxs[3] = { 1, 1, 1 };
  data->SetExtent(0, vdim, 0, vdim, 0, vdim);
  data->AllocateScalars(VTK_FLOAT, 1);
  int* extent = data->GetExtent();
  // Port of Kewei's data generation code from his isosurface implemention
  // See isosurface.cpp:TangleField
  for (int z = extent[4]; z <= extent[5]; ++z){
    for (int y = extent[2]; y <= extent[3]; ++y){
      for (int x = extent[0]; x <= extent[1]; ++x){
        float* pixel = static_cast<float*>(data->GetScalarPointer(x, y, z));
        const float xx = 3.0 * (mins[0] + (maxs[0] - mins[0]) * (1.0 * x / dim));
        const float yy = 3.0 * (mins[1] + (maxs[1] - mins[1]) * (1.0 * y / dim));
        const float zz = 3.0 * (mins[2] + (maxs[2] - mins[2]) * (1.0 * z / dim));
        const float v = (xx*xx*xx*xx - 5.0f*xx*xx + yy*yy*yy*yy - 5.0f*yy*yy
            + zz*zz*zz*zz - 5.0f*zz*zz + 11.8f) * 0.2f + 0.5f;
        pixel[0] = v;
      }
    }
  }
}


int main()
{
  std::vector<double> samples;

  const double MAX_RUNTIME = 30;
  const size_t MAX_ITERATIONS = 1; 
  samples.reserve(MAX_ITERATIONS);
  size_t iter = 0;
  Timer timer;
  for (double el = 0.0; el < MAX_RUNTIME && iter < MAX_ITERATIONS; el += samples.back(), ++iter)
  {
    vtkImageData *data = vtkImageData::New();
    GenerateData(data);

    vtkContourFilter *iso = vtkContourFilter::New();

    iso->SetInputData(data);
    iso->SetValue(0, ISO_START);
    iso->ComputeScalarsOn();
    iso->ComputeGradientsOff();
    iso->ComputeNormalsOff();
      
    timer.Reset();
   
    for (int j = 0; j < ISO_NUM; ++j)
    {
      iso->SetValue(0, iso->GetValue(0) + ISO_STEP);
      iso->Update(); 

      vtkPolyData* output = iso->GetOutput();
      std::cout << (iso->GetValue(0)) << " " << output->GetNumberOfPoints() << std::endl;
    }
    samples.push_back(timer.GetElapsedTime());
    
    iso->Delete();
    data->Delete();
  }
    
  std::sort(samples.begin(), samples.end());
  stats::Winsorize(samples, 5.0);
  std::cout << "Benchmark \'VTK MPI Isosurface\' results:\n"
      << "\tmedian = " << stats::PercentileValue(samples, 50.0) << "s\n"
      << "\tmedian abs dev = " << stats::MedianAbsDeviation(samples) << "s\n"
      << "\tmean = " << stats::Mean(samples) << "s\n"
      << "\tstd dev = " << stats::StandardDeviation(samples) << "s\n"
      << "\tmin = " << samples.front() << "s\n"
      << "\tmax = " << samples.back() << "s\n"
      << "\t# of runs = " << samples.size() << "\n";

  return 0;
}

