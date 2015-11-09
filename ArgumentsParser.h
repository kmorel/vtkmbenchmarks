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
#ifndef __argumentsParser_h
#define __argumentsParser_h

#include <string>

namespace vtkm { namespace testing {

class ArgumentsParser
{
public:
  ArgumentsParser();
  virtual ~ArgumentsParser();

  bool parseArguments(int argc, char* argv[]);

  std::string file() const
    { return this->File; }

  float isovalue() const
    { return this->IsoValue; }

  float isovalue_step() const
  { return this->IsoValueStep; }

  int num_iso() const
  { return this->NumIso; }

  int num_trials() const
  { return this->NumTrials; }

  double ratio() const
    { return this->Ratio; }

  std::string writeLocation() const
    { return this->WriteLocation; }

private:
  std::string File;
  std::string WriteLocation;
  float IsoValue;
  float IsoValueStep;
  int NumIso;
  int NumTrials;
  double Ratio;
};

}}
#endif
