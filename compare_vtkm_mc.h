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

#include <vtkm/cont/ArrayHandle.h>
#include <vtkm/cont/ArrayHandleUniformPointCoordinates.h>
#include <vtkm/cont/CellSetStructured.h>
#include <vtkm/cont/DataSet.h>
#include <vtkm/cont/Timer.h>

#include <vtkm/worklet/IsosurfaceUniformGrid.h>

#include <vtkImageData.h>

#include <vector>

#include "ArgumentsParser.h"

namespace vtkm
{
static void RunIsoSurfaceUniformGrid(
    const std::vector<vtkm::Float32>& buffer,
    vtkImageData* image,
    const std::string& device,
    const vtkm::testing::ArgumentsParser &arguments)
{

  typedef VTKM_DEFAULT_DEVICE_ADAPTER_TAG DeviceAdapter;

  int dims[3];
  image->GetDimensions(dims);

  const vtkm::Id3 pointDims(dims[0], dims[1], dims[2]);
  const vtkm::Id3 cellDims(dims[0]-1, dims[1]-1, dims[2]-1);

  vtkm::cont::ArrayHandleUniformPointCoordinates coordinates(pointDims);

  vtkm::cont::CellSetStructured<3> cellSet("cells");
  cellSet.SetPointDimensions(pointDims);

  vtkm::cont::DataSet dataSet;

  dataSet.AddCellSet(cellSet);
  dataSet.AddCoordinateSystem(
          vtkm::cont::CoordinateSystem("coordinates", 1, coordinates));


  vtkm::cont::ArrayHandle<vtkm::Float32> field = vtkm::cont::make_ArrayHandle(buffer);
  dataSet.AddField(vtkm::cont::Field("nodevar", 1, vtkm::cont::Field::ASSOC_POINTS, field));

  vtkm::cont::ArrayHandle< vtkm::Float32 > scalarsArray;
  vtkm::cont::ArrayHandle< vtkm::Vec<vtkm::Float32,3> > verticesArray;
  vtkm::cont::ArrayHandle< vtkm::Vec<vtkm::Float32,3> > normalsArray;

  vtkm::worklet::IsosurfaceFilterUniformGrid<vtkm::Float32,
                                               DeviceAdapter> isosurfaceFilter(cellDims, dataSet);

  vtkm::cont::Timer<> timer;

  float isoValue = arguments.isovalue();
  for (int isoIndex = 0; isoIndex < arguments.num_iso(); isoIndex++)
  {
    for (int trial = 0; trial < arguments.num_trials(); trial++)
    {
      timer.Reset();
      isosurfaceFilter.Run(isoValue,
                           field,
                           verticesArray,
                           normalsArray,
                           scalarsArray);
      vtkm::Float64 walltime = timer.GetElapsedTime();

      std::cout
          << "- device         : " << device << std::endl
          << "  implementation : VTK-m MC" << std::endl
          << "  trial          : " << trial << std::endl
          << "  seconds        : " << walltime << std::endl
          << "  isovalue       : " << isoValue << std::endl
          << "  num vertices   : " << verticesArray.GetNumberOfValues() << std::endl;
    }
    isoValue += arguments.isovalue_step();
  }
}

}
