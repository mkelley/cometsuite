/***************************************************************************

  Implements rundynamics file logs.

  Copyright (C) 2005,2006,2007,2008 by Michael S. Kelley
  <msk@astro.umd.edu>

 ***************************************************************************/

#if HAVE_CONFIG_H
#  include <config.h>
#endif /* HAVE_CONFIG_H */

#include <iostream>
#include <string>
#include <cstdlib>
#include <cstdio>
#include <ctime>
#include <fstream>
#include <sys/types.h>
#include <dirent.h>
#include "logFile.h"

using namespace std;

/** Opens the log file for run time status. */
logFile::logFile(const string basename) {
//   string filename = BASELOGNAME + nextLog() + ".txt";
  string filename = basename + ".log";
  cerr << "Run time output being sent to " << filename << endl;
  _logFile.open(filename.c_str());

  time_t t;
  time(&t);
  _logFile << "# " << filename << "  " << asctime(localtime(&t));
}

/** Closes the log file. */
logFile::~logFile() { _logFile.close(); }

/** Flushes the log file buffer. */
void logFile::flush() { _logFile.flush(); }

/** Adds a string to the log file.
 *  Newlines are not automatically added. */
void logFile::add(const string str) {
  _logFile << str;
}

/** Adds a long int to the log file. */
void logFile::add(const long l) {
  _logFile << l << " ";
}
