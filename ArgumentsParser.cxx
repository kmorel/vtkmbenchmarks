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

enum  optionIndex { UNKNOWN, HELP, FILEPATH, WRITE_LOC, ISO_VALUE, ISO_VALUE_STEP, NUM_ISO, NUM_TRIALS, RESAMPLE_RATIO};
const vtkm::testing::option::Descriptor usage[] =
{
  {UNKNOWN,        0, "" , ""    ,      vtkm::testing::option::Arg::None,     "USAGE: example [options]\n\n"
                                                                             "Options:" },
  {HELP,           0, "h", "help",      vtkm::testing::option::Arg::None,     "  --help, -h  Print usage and exit." },
  {FILEPATH,       0, "",  "file",      vtkm::testing::option::Arg::Optional, "  --file      nrrd file to read." },
  {WRITE_LOC,      0, "",  "dump",      vtkm::testing::option::Arg::Optional, "  --dump      Folder to write dumps of the results of each algorithm." },
  {ISO_VALUE,      0, "",  "isovalue",  vtkm::testing::option::Arg::Optional, "  --isovalue  Value to contour the dataset at." },
  {ISO_VALUE_STEP, 0, "",  "isostep",   vtkm::testing::option::Arg::Optional, "  --isostep   IsoValue increment."},
  {NUM_ISO,        0, "",  "numiso",    vtkm::testing::option::Arg::Optional, "  --numiso    Number of isovalues to use." },
  {NUM_TRIALS,     0, "",  "numtrials", vtkm::testing::option::Arg::Optional, "  --numtrials Number of trials per isovalue." },
  {RESAMPLE_RATIO, 0, "",  "ratio",     vtkm::testing::option::Arg::Optional, "  --ratio     Resample ratio for the input data." },
  {UNKNOWN,        0, "",  "",          vtkm::testing::option::Arg::None,     "\nExample:\n"
                                                                              " example --file=./normal_1349.nhdr\n"},
  {0,0,0,0,0,0}
};


//-----------------------------------------------------------------------------
vtkm::testing::ArgumentsParser::ArgumentsParser():
  File(""),
  WriteLocation(""),
  IsoValue(0.0f),
  IsoValueStep(0.005f),
  NumIso(10),
  NumTrials(11),
  Ratio(1.0)
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

  if ( options[ISO_VALUE_STEP] )
    {
    std::string sarg(options[ISO_VALUE_STEP].last()->arg);
    std::stringstream argstream(sarg);
    argstream >> this->IsoValueStep;
    }

  if ( options[NUM_ISO] )
    {
    std::string sarg(options[NUM_ISO].last()->arg);
    std::stringstream argstream(sarg);
    argstream >> this->NumIso;
    }

  if ( options[NUM_TRIALS] )
    {
    std::string sarg(options[NUM_TRIALS].last()->arg);
    std::stringstream argstream(sarg);
    argstream >> this->NumTrials;
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
