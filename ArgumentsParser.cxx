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

#include "ArgumentsParser.h"

#include <vtkm/testing/OptionParser.h>
#include <iostream>
#include <sstream>
#include <string>

enum  optionIndex { UNKNOWN, HELP, FILEPATH, PIPELINE};
const vtkm::testing::option::Descriptor usage[] =
{
  {UNKNOWN,   0,"" , ""    ,      vtkm::testing::option::Arg::None, "USAGE: example [options]\n\n"
                                                                    "Options:" },
  {HELP,      0,"h" , "help",    vtkm::testing::option::Arg::None, "  --help, -h  \tPrint usage and exit." },
  {FILEPATH,      0,"", "file",      vtkm::testing::option::Arg::Optional, "  --file  \t nrrd file to read." },
  {PIPELINE,  0,"", "pipeline",  vtkm::testing::option::Arg::Optional, "  --pipeline  \t What pipeline to run ( 1 threshold, 2 marching cubes)." },
  {UNKNOWN,   0,"",  "",         vtkm::testing::option::Arg::None, "\nExample:\n"
                                                                   " example --file=./test --pipeline=1\n"},
  {0,0,0,0,0,0}
};


//-----------------------------------------------------------------------------
vtkm::testing::ArgumentsParser::ArgumentsParser():
  File(""),
  Pipeline(THRESHOLD)
{
}

//-----------------------------------------------------------------------------
vtkm::testing::ArgumentsParser::~ArgumentsParser()
{
}

//-----------------------------------------------------------------------------
bool vtkm::testing::ArgumentsParser::parseArguments(int argc, char* argv[])
{

  argc-=(argc>0);
  argv+=(argc>0); // skip program name argv[0] if present

  vtkm::testing::option::Stats  stats(usage, argc, argv);
  vtkm::testing::option::Option* options = new vtkm::testing::option::Option[stats.options_max];
  vtkm::testing::option::Option* buffer = new vtkm::testing::option::Option[stats.options_max];
  vtkm::testing::option::Parser parse(usage, argc, argv, options, buffer);

  if (parse.error())
    {
    delete[] options;
    delete[] buffer;
    return false;
    }

  if (options[HELP] || argc == 0)
    {
    vtkm::testing::option::printUsage(std::cout, usage);
    delete[] options;
    delete[] buffer;

    return false;
    }

  if ( options[FILEPATH] )
    {
    std::string sarg(options[FILEPATH].last()->arg);
    std::stringstream argstream(sarg);
    argstream >> this->File;
    }

  if ( options[PIPELINE] )
    {
    std::string sarg(options[PIPELINE].last()->arg);
    std::stringstream argstream(sarg);
    int pipelineflag = 0;
    argstream >> pipelineflag;
    if (pipelineflag == THRESHOLD)
      {
      this->Pipeline = THRESHOLD;
      }
    else if (pipelineflag == MARCHING_CUBES)
      {
      this->Pipeline = MARCHING_CUBES;
      }
    else
      {
      std::cerr << "Incorrect pipeline choice: " << pipelineflag << std::endl;
      std::cerr << "Threshold is : " << THRESHOLD  << std::endl;
      std::cerr << "Marching Cubes is : " << MARCHING_CUBES  << std::endl;
      }
    }

  delete[] options;
  delete[] buffer;
  return true;
}
