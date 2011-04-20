/***************************************************************************
  Copyright (C) 2005,2007,2008 by Michael S. Kelley
  <msk@astro.umd.edu>

  ***************************************************************************/

#if !defined(__LOGFILE)
#define __LOGFILE 1

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <fstream>
#include "rundynamics.h"

using namespace std;

/** The logFile file provides runtime status. */
class logFile {
 public:
  logFile(const string);
  ~logFile();
  void flush();
  void add(const string);
  void add(const long);

 private:
  ofstream _logFile;
};

#endif
