cmake_minimum_required(VERSION 2.8.11)

project(VTKMBenchmark)

set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} ${mccompare_SOURCE_DIR})

find_package(VTKm REQUIRED)
include(VTKmMacros)
vtkm_configure_device(Serial)

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

option(ENABLE_PISTON "Benchmark piston comparison" OFF)
if(${ENABLE_PISTON})
 find_path( PISTON_INCLUDE
    NAMES piston/piston_math.h
    DOC "Piston headers"
    )
 include_directories(${PISTON_INCLUDE})
 add_definitions("-DPISTON_ENABLED")
endif()


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

set_target_properties(BenchmarkSerial PROPERTIES COMPILE_FLAGS "-DTHRUST_DEVICE_SYSTEM=THRUST_DEVICE_SYSTEM_CPP")


#Add TBB version
option(ENABLE_TBB "Benchmark TBB Backend" ON)
if(${ENABLE_TBB})
  vtkm_configure_device(TBB)

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

  set_target_properties(BenchmarkTBB PROPERTIES COMPILE_FLAGS "-DTHRUST_DEVICE_SYSTEM=THRUST_DEVICE_SYSTEM_TBB")
endif() # TBB

#Add CUDA version
option(ENABLE_CUDA "Benchmark CUDA Backend" ON)
if(${ENABLE_CUDA})
  vtkm_configure_device(Cuda)

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

  set_target_properties(BenchmarkCuda PROPERTIES COMPILE_FLAGS "-DTHRUST_DEVICE_SYSTEM=THRUST_DEVICE_SYSTEM_CUDA")
endif() # CUDA
