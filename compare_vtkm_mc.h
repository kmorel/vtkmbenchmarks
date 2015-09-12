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

#include "Stats.h"

namespace vtkm
{
static void RunIsoSurfaceUniformGrid(const std::vector<vtkm::Float32>& buffer,
                                     vtkImageData* image,
                                     const std::string& device,
                                     int numCores,
                                     int maxNumCores,
                                     float isoValue,
                                     int MAX_NUM_TRIALS)
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
  std::vector<double> samples;
  samples.reserve(MAX_NUM_TRIALS);
  timer.Reset();

  for(int i=0; i<MAX_NUM_TRIALS; ++i)
  {
    isosurfaceFilter.Run(isoValue,
                         field,
                         verticesArray,
                         normalsArray,
                         scalarsArray);

    std::cout << isoValue << " " << verticesArray.GetNumberOfValues() << std::endl;
    isoValue += 0.005f;
  }

  samples.push_back(timer.GetElapsedTime());

  std::sort(samples.begin(), samples.end());
  stats::Winsorize(samples, 5.0);
  std::cout << "Benchmark \'VTK-m Isosurface\' results:\n"
        << "\tmedian = " << stats::PercentileValue(samples, 50.0) << "s\n"
        << "\tmedian abs dev = " << stats::MedianAbsDeviation(samples) << "s\n"
        << "\tmean = " << stats::Mean(samples) << "s\n"
        << "\tstd dev = " << stats::StandardDeviation(samples) << "s\n"
        << "\tmin = " << samples.front() << "s\n"
        << "\tmax = " << samples.back() << "s\n"
        << "\t# of runs = " << samples.size() << "\n";


}

}
