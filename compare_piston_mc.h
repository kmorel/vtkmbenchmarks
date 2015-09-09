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
  for(int i=0; i < MAX_NUM_TRIALS; ++i)
    {
    vtkm::cont::Timer<> timer;

    piston_scalar_image3d pimage(dims[0],dims[1],dims[2],buffer);
    MC marching(pimage,pimage,isoValue);

    marching();

    double time = timer.GetElapsedTime();
    std::cout << "piston::marching_cube," << device << "," <<  numCores << "," << time << "," << i << std::endl;
    }
}

}
