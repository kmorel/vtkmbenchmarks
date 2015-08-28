cmake_minimum_required(VERSION 2.8.11)

project(VTKMBenchmark)

set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} ${mccompare_SOURCE_DIR})

find_package(VTKm REQUIRED)
include(VTKmMacros)
vtkm_configure_device(Serial)
vtkm_configure_device(TBB)
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

set(headers
  compare.h
  compare_vtk_mc.h
  compare_vtkm_mc.h
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


#Add TBB version
add_executable(BenchmarkTBB
  ${srcs}
  ${headers}
  mainTBB.cxx
  )

target_link_libraries(BenchmarkTBB
  vtkCommonCore
  vtkCommonDataModel
  vtkCommonExecutionModel
  vtkCommonMisc
  vtkFiltersCore
  vtkFiltersGeometry
  vtkImagingCore
  vtkIOImage
  vtkIOLegacy
  ${TBB_LIBRARIES}
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

