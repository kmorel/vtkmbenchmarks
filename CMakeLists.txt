cmake_minimum_required(VERSION 2.8)

project(VTKMBenchmark CXX)

set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} ${mccompare_SOURCE_DIR})

include_directories("/Users/robert/Work/vtkm/src/")
include_directories("/Users/robert/Work/vtkm/cbuild/")

find_package(Cuda)

find_package(VTK 6.2
  COMPONENTS
  vtkCommonCore
  vtkCommonDataModel
  vtkCommonExecutionModel
  vtkCommonMisc
  vtkFiltersCore
  vtkFiltersGeometry
  vtkIOLegacy
  vtkIOImage
  vtkImagingCore
  REQUIRED)

include(${VTK_USE_FILE})

set_source_files_properties(compare.h
  PROPERTIES HEADER_FILE_ONLY TRUE
  )

set(srcs
  compare.h
  compare_mc.h
  compare_thresh.h
  ArgumentsParser.cxx
  )

#-----------------------------------------------------------------------------
# setup timing libs we need to link too
set(VTKM_TIMING_LIBS)
if(UNIX AND NOT APPLE)
  set(VTKM_TIMING_LIBS rt)
endif()


add_executable(BenchmarkSerial
  ${srcs}
  mainSerial.cxx
  )

target_link_libraries(BenchmarkSerial
  vtkCommonCore
  vtkCommonDataModel
  vtkCommonExecutionModel
  vtkCommonMisc
  vtkFiltersCore
  vtkFiltersGeometry
  vtkImagingCore
  vtkIOImage
  vtkIOLegacy
  )

#Add CUDA version
cuda_add_executable(BenchmarkCuda
  ${srcs}
  mainCuda.cu
  )

target_link_libraries(BenchmarkCuda
  vtkCommonCore
  vtkCommonDataModel
  vtkCommonExecutionModel
  vtkCommonMisc
  vtkFiltersCore
  vtkFiltersGeometry
  vtkImagingCore
  vtkIOImage
  vtkIOLegacy
  )
