cmake_minimum_required(VERSION 2.8)
 
PROJECT(VTKIsoBenchmarkProject)
 
set(VTK_DIR "/home/csewell/TGABenchmarks/VTK-Install6.1/lib/cmake/vtk-6.1/")
find_package(VTK REQUIRED)
include(${VTK_USE_FILE})

find_package(MPI)
include_directories(SYSTEM ${MPI_INCLUDE_PATH})
 
add_executable(VTKIsoBenchmark ParallelIso.cxx)
add_executable(VTKIsoSerial SerialIso.cxx)
 
if(VTK_LIBRARIES)
  target_link_libraries(VTKIsoBenchmark ${VTK_LIBRARIES} ${MPI_CPP_LIBRARIES})
  target_link_libraries(VTKIsoSerial ${VTK_LIBRARIES})
#else()
#  target_link_libraries(VTK vtkHybrid vtkWidgets)
endif()
