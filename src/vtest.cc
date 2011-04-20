/***************************************************************************

  Runs various ejection velocity methods and writes out the results.

  Copyright (C) 2008-2010 by Michael S. Kelley <msk@astro.umd.edu>

 ***************************************************************************/

#if HAVE_CONFIG_H
#  include <config.h>
#endif /* HAVE_CONFIG_H */

#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <cmath>
#include "Vej.h"
#include "Vector.h"
#include "paramSet.h"
#include "pfunctions.h"
#include "rundynamics.h"

#define SUBPROJECT "vtest"
#define OUTBASENAME "vtest-"
#define N 1000

using namespace std;

int main(int argc, char *argv[])
{
  Vej vej;
  ofstream outf;
  Vector r(_AU, 0, 0), v;
  particle p;
  paramSet parameters;
  pfunctions pf;

  stringstream str;
  str << "KERNEL: encke.bsp";
  parameters.loadParameters(str);
  p.parameters(parameters);

  /**********************************************************************/
  cout << "Test 1: isotropic outflow.\n";

  string generator = string("logbeta -5 -1; age 0 0; velocity iso 1.0; q_d cos; suncone 0 90");
  parameters.pFunc(generator);
  pf.setup(parameters, p);
  p.next();  /* shoot one off to compute the rh */

  string f = string(OUTBASENAME) + "1.dat";
  outf.open(f.c_str(), ios::out | ios::trunc);
  outf << "# Test 1: isotropic outflow.\n";
  outf << "# " << generator << "\n";
  outf << "# rh = " << p.istate().r.length() / _AU << "\n";
  outf << "# beta v_x v_y v_z\n";

  for (int i=0; i<N; i++) {
    p.next();
    outf << p.beta() << " " << p.vej() << "\n";
  }

  outf.close();

  /*********************************************************************/
  cout << "Test 2: Q_d \\propto cos(sun-zenith angle).\n";

  generator = string("logbeta -5 -1; age 0 0; velocity cos 1.0; q_d iso; suncone 0 90");
  parameters.pFunc(generator);
  pf.setup(parameters, p);

  f = string(OUTBASENAME) + "2.dat";
  outf.open(f.c_str(), ios::out | ios::trunc);
  outf << "# Test 2: Q_c \\propto cos(sun-zenith angle).\n";
  outf << "# " << generator << "\n";
  outf << "# rh = " << p.istate().r.length() / _AU << "\n";
  outf << "# beta v_x v_y v_z r_x r_y r_z\n";

  for (int i=0; i<N; i++) {
    p.next();
    outf << p.beta() << " " << p.vej() << " " << p.istate().r << "\n";
  }

  outf.close();

  /*********************************************************************/
  cout << "Test 3: qv_cos_coma.\n";

  generator = string("qv_cos_coma 0.5 4000 -5 -1");
  parameters.pFunc(generator);
  pf.setup(parameters, p);

  f = string(OUTBASENAME) + "3.dat";
  outf.open(f.c_str(), ios::out | ios::trunc);
  outf << "# Test 3: qv_cos_coma.\n";
  outf << "# " << generator << "\n";
  outf << "# beta v_x v_y v_z r_x r_y r_z\n";

  for (int i=0; i<N; i++) {
    p.next();
    outf << p.beta() << " " << p.vej() << " " << p.istate().r << "\n";
  }

  outf.close();

  return EXIT_SUCCESS;
}

