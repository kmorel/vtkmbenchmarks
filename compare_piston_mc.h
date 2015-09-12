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

#include <piston/marching_cube.h>
#include <piston/image3d.h>

#include "Stats.h"

namespace piston
{

struct piston_scalar_image3d : piston::image3d<thrust::device_system_tag>
{
  typedef thrust::device_vector<vtkm::Float32> PointDataContainer;
  PointDataContainer point_data_vector;
  typedef PointDataContainer::iterator PointDataIterator;

  piston_scalar_image3d(vtkm::IdComponent xsize, vtkm::IdComponent ysize, vtkm::IdComponent zsize,
                        const std::vector<vtkm::Float32> &data)
    : piston::image3d< thrust::device_system_tag >(xsize, ysize, zsize),
      point_data_vector(data)
  {
    assert(this->NPoints == this->point_data_vector.size());
  }

  PointDataIterator point_data_begin() {
    return this->point_data_vector.begin();
  }
  PointDataIterator point_data_end() {
    return this->point_data_vector.end();
  }
};

static void RunIsoSurfaceUniformGrid(const std::vector<vtkm::Float32>& buffer,
                                     vtkImageData* image,
                                     const std::string& device,
                                     int numCores,
                                     int maxNumCores,
                                     float isoValue,
                                     int MAX_NUM_TRIALS)
{

  int dims[3];
  image->GetDimensions(dims);

  typedef piston::marching_cube< piston_scalar_image3d,
                                 piston_scalar_image3d > MC;
  
  piston_scalar_image3d pimage(dims[0],dims[1],dims[2],buffer);
  MC marching(pimage,pimage,isoValue);

  vtkm::cont::Timer<> timer;
  std::vector<double> samples;
  samples.reserve(MAX_NUM_TRIALS);
  timer.Reset();

  for(int i=0; i < MAX_NUM_TRIALS; ++i)
  {
    //piston_scalar_image3d pimage(dims[0],dims[1],dims[2],buffer);
    //MC marching(pimage,pimage,isoValue);
    marching.set_isovalue(isoValue);
    marching();

    std::cout << isoValue << " " << marching.num_total_vertices << std::endl;
    isoValue += 0.005;
  }

  samples.push_back(timer.GetElapsedTime());

  std::sort(samples.begin(), samples.end());
  stats::Winsorize(samples, 5.0);
  std::cout << "Benchmark \'Piston Isosurface\' results:\n"
        << "\tmedian = " << stats::PercentileValue(samples, 50.0) << "s\n"
        << "\tmedian abs dev = " << stats::MedianAbsDeviation(samples) << "s\n"
        << "\tmean = " << stats::Mean(samples) << "s\n"
        << "\tstd dev = " << stats::StandardDeviation(samples) << "s\n"
        << "\tmin = " << samples.front() << "s\n"
        << "\tmax = " << samples.back() << "s\n"
        << "\t# of runs = " << samples.size() << "\n";

}

}
