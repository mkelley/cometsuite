/***************************************************************************

  Tests the normal distribution code.

  Copyright (C) 2010 by Michael S. Kelley <msk@astro.umd.edu>

 ***************************************************************************/

#if HAVE_CONFIG_H
#  include <config.h>
#endif /* HAVE_CONFIG_H */

#include <iostream>
#include <fstream>
#include <string>
#include <cstdlib>
#include "Distribution.h"
#include "rundynamics.h"

#define SUBPROJECT "randtest"
#define OUTBASENAME "randtest-"
#define N 1000

using namespace std;

int main(int argc, char *argv[])
{
  ofstream outf;
  Distribution d;

  string f = string(OUTBASENAME) + "1.dat";
  outf.open(f.c_str(), ios::out | ios::trunc);
  outf << "# Test 1: normal variates.\n";
  outf << "# column  mu  sigma  minVal\n";
  outf << "#      1   0      1       0\n";
  outf << "#      2  10      2       0\n";
  outf << "#      3  10      2       5\n";

  for (int i=0; i<N; i++) {
    outf << d.dn_dx__normal(0, 1, 0) << " " <<
      d.dn_dx__normal(10, 2, 0) << " " <<
      d.dn_dx__normal(10, 2, 5) << "\n";
  }

  outf.close();

  return EXIT_SUCCESS;
}

