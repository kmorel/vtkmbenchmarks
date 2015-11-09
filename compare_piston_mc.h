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

#include <vtkm/cont/Timer.h>

#include "ArgumentsParser.h"

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

static void RunIsoSurfaceUniformGrid(
    const std::vector<vtkm::Float32>& buffer,
    vtkImageData* image,
    const std::string& device,
    const vtkm::testing::ArgumentsParser &arguments)
{

  int dims[3];
  image->GetDimensions(dims);

  typedef piston::marching_cube< piston_scalar_image3d,
                                 piston_scalar_image3d > MC;

  float isoValue = arguments.isovalue();

  piston_scalar_image3d pimage(dims[0],dims[1],dims[2],buffer);
  MC marching(pimage,pimage,isoValue);

  vtkm::cont::Timer<> timer;

  for (int isoIndex = 0; isoIndex < arguments.num_iso(); isoIndex++)
  {
    for (int trial = 0; trial < arguments.num_trials(); trial++)
    {
      timer.Reset();
      marching.set_isovalue(isoValue);
      marching();
      vtkm::Float64 walltime = timer.GetElapsedTime();

      std::cout
          << "- device         : " << device << std::endl
          << "  implementation : PISTON MC" << std::endl
          << "  trial          : " << trial << std::endl
          << "  seconds        : " << walltime << std::endl
          << "  isovalue       : " << isoValue << std::endl
          << "  num vertices   : " << marching.num_total_vertices << std::endl;
    }
    isoValue += arguments.isovalue_step();
  }
}

}
