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

enum  optionIndex { UNKNOWN, HELP, FILEPATH, WRITE_LOC, ISO_VALUE, CORES, RESAMPLE_RATIO};
const vtkm::testing::option::Descriptor usage[] =
{
  {UNKNOWN,   0,"" , ""    ,      vtkm::testing::option::Arg::None, "USAGE: example [options]\n\n"
                                                                    "Options:" },
  {HELP,      0,"h" , "help",    vtkm::testing::option::Arg::None, "  --help, -h  \tPrint usage and exit." },
  {FILEPATH,      0,"", "file",      vtkm::testing::option::Arg::Optional, "  --file  \t nrrd file to read." },
  {WRITE_LOC,  0,"", "dump",  vtkm::testing::option::Arg::Optional, "  --dump  \t Folder to write dumps of the results of each algorithm." },
  {ISO_VALUE,  0,"", "isovalue",  vtkm::testing::option::Arg::Optional, "  --isovalue  \t Value to contour the dataset at." },
  {CORES,  0,"", "cores",        vtkm::testing::option::Arg::Optional, "  --cores  \t number of cores to use, 0 means all cores, -1 means test with 1 to max cores." },
  {RESAMPLE_RATIO,  0,"", "ratio",  vtkm::testing::option::Arg::Optional, "  --ratio  \t Resample ratio for the input data." },
  {UNKNOWN,   0,"",  "",         vtkm::testing::option::Arg::None, "\nExample:\n"
                                                                   " example --file=./test --pipeline=1\n"},
  {0,0,0,0,0,0}
};


//-----------------------------------------------------------------------------
vtkm::testing::ArgumentsParser::ArgumentsParser():
  File(""),
  WriteLocation(""),
  IsoValue(0.0f),
  Ratio(1.0),
  Cores(0)
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

  if ( options[WRITE_LOC] )
    {
    std::string sarg(options[WRITE_LOC].last()->arg);
    std::stringstream argstream(sarg);
    argstream >> this->WriteLocation;
    }

  if ( options[ISO_VALUE] )
    {
    std::string sarg(options[ISO_VALUE].last()->arg);
    std::stringstream argstream(sarg);
    argstream >> this->IsoValue;
    }

  if ( options[CORES] )
    {
    std::string sarg(options[CORES].last()->arg);
    std::stringstream argstream(sarg);
    argstream >> this->Cores;
    }

  if ( options[RESAMPLE_RATIO] )
    {
    std::string sarg(options[RESAMPLE_RATIO].last()->arg);
    std::stringstream argstream(sarg);
    argstream >> this->Ratio;
    }

  delete[] options;
  delete[] buffer;
  return true;
}