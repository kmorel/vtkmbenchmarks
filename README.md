## VTK-m Benchmark Suite ##

VTK-m is an open source C++  library that provides a collection of data analysis and visualization algorithms that run well on multi- and many-core processors. We currently benchmark CUDA and Serial Execution.

### How ###

VTK-m uses fine-grained concurrency for data analysis and visualization algorithms.
The basic computational unit of the VTK-m is a worklet, a function that implements the algorithm’s behavior on an element of a mesh (that is, a point, edge, face, or cell) or a small local neighborhood.

The worklet is constrained to be serial and stateless; it can access only the element passed to and from the invocation. With this constraint, the serial worklet function can be concurrently scheduled on an unlimited number of threads without the complications of threads or race conditions.

Although worklets are not allowed communication, many visualization algorithms require operations such as variable array packing and coincident topology resolution that intrinsically require significant coordination among threads. VTK-m enables such algorithms by classifying and implementing the most common and versatile communicative operations into worklet types which are managed by the dispatcher.

## Getting VTK-m ##


The VTK-m repository is located at [http://public.kitware.com/vtkm.git](http://public.kitware.com/vtkm.git)

VTK-m dependencies are:


+  [CMake 2.8.10](http://cmake.org/cmake/resources/software.html)
+  [Boost 1.50.0](http://www.boost.org) or greater
+  [Cuda Toolkit 6+](https://developer.nvidia.com/cuda-toolkit) [optional]

```
git clone http://public.kitware.com/vtkm.git vtkm
mkdir vtkm-build
cd vtkm-build
cmake-gui ../vtkm
```

## How To Use Benchmarks ##

Each program has two arguments, file and pipeline. File points
to a nrrd file, and pipeline can be 1 for threshold or 2 for marching cubes

```
./BenchmarkCuda --file=./data.nhdr --pipeline=2

```

## License ##
```
Copyright (c) Kitware, Inc.
All rights reserved.
See LICENSE.txt for details.
This software is distributed WITHOUT ANY WARRANTY; without even
the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
PURPOSE.  See the above copyright notice for more information.
Copyright 2014 Sandia Corporation.
Copyright 2014 UT-Battelle, LLC.
Copyright 2014 Los Alamos National Security.
Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
the U.S. Government retains certain rights in this software.
Under the terms of Contract DE-AC52-06NA25396 with Los Alamos National
Laboratory (LANL), the U.S. Government retains certain rights in
this software.
```