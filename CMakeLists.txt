cmake_minimum_required(VERSION 2.8.11)

project(VTKMBenchmark)

set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} ${mccompare_SOURCE_DIR})

find_package(VTKm REQUIRED)
include(VTKmMacros)
vtkm_configure_device(Serial)
vtkm_configure_device(Cuda)

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
option(VTK_HAS_FLYING_EDGES "Enable if you are building a vtk version that has flying edges 3d" OFF)

set(headers
  compare.h
  compare_classic_mc.h
  compare_per_tri_out_mc.h
  compare_sliding_mc.h
  compare_sliding_per_tri_out_mc.h
  compare_vtk_mc.h
  worklets.h
  )

set(srcs
  ArgumentsParser.cxx
  )

set_source_files_properties(${headers} PROPERTIES HEADER_FILE_ONLY TRUE)

add_executable(BenchmarkSerial
  ${srcs}
  ${headers}
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
  ${headers}
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

if(VTK_HAS_FLYING_EDGES)
  target_compile_definitions(BenchmarkSerial PRIVATE "VTK_HAS_FLYING_EDGES")
  target_compile_definitions(BenchmarkCuda   PRIVATE "VTK_HAS_FLYING_EDGES")
endif()
