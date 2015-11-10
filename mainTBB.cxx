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
#define VTKM_DEVICE_ADAPTER VTKM_DEVICE_ADAPTER_TBB

#ifdef PISTON_ENABLED
#include <thrust/detail/config.h>
#define THRUST_DEVICE_BACKEND THRUST_DEVICE_BACKEND_TBB
#endif

#include "ArgumentsParser.h"
#include "compare.h"
#include <tbb/task_scheduler_init.h>

int main(int argc, char* argv[])
  {
  vtkm::testing::ArgumentsParser parser;
  if (!parser.parseArguments(argc, argv))
    {
    return 1;
    }

//  int maxNumCores = tbb::task_scheduler_init::default_num_threads();

  RunComparison("TBB", parser);

  return 0;
}
