/***************************************************************************

  Runs various ejection velocity methods and writes out the results.

  Copyright (C) 2008-2010,2012 by Michael S. Kelley <msk@astro.umd.edu>

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

#define SUBPROJECT "vejtest"
#define OUTBASENAME "vejtest-"
#define N 10000

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

  string f = string(OUTBASENAME) + "iso.dat";
  outf.open(f.c_str(), ios::out | ios::trunc);
  outf.precision(3);
  outf << "# Test 1: isotropic outflow.\n";
  outf << "# The sunward hemisphere has v_ej of 1 km/s.\n";
  outf << "# The opposite hemisphere has v_ej of 0.25 km/s\n";
  outf << "# Output is relative to the comet, units of km/s\n";

  parameters.pFunc("age 0 0; radius 1 1; velocity iso 1.0; q_d iso; suncone 0 90;");
  pf.setup(parameters, p);

  // get the rh of the comet
  p.next();
  outf << "# The sun is located at " <<
    p.istate().r[0] / _AU << " " <<
    p.istate().r[1] / _AU << " " <<
    p.istate().r[2] / _AU << " AU\n";
  outf << "# v_x v_y v_z\n";

  for (int i=0; i<N; i++) {
    p.next();
    outf << p.vej() << "\n";
  }

  parameters.pFunc("velocity iso 0.25; suncone 90 180;");
  pf.setup(parameters, p);
  for (int i=0; i<N; i++) {
    p.next();
    outf << p.vej() << "\n";
  }

  outf.close();

  /**********************************************************************/
  cout << "Test 2: Q_d \\propto cos(sun-zenith angle).\n";

  f = string(OUTBASENAME) + "iso-cos.dat";
  outf.open(f.c_str(), ios::out | ios::trunc);
  outf.precision(3);
  outf << "# Test 2: cos(z_sun) outflow.\n";
  outf << "# The sunward hemisphere has v_ej of 1 km/s.\n";
  outf << "# Output is relative to the comet, units of km/s\n";
  outf << "# The sun is located at " <<
    p.istate().r[0] / _AU << " " <<
    p.istate().r[1] / _AU << " " <<
    p.istate().r[2] / _AU << " AU\n";
  outf << "# v_x v_y v_z\n";

  parameters.pFunc("velocity iso 1.0; q_d cos; suncone 0 90;");
  pf.setup(parameters, p);
  for (int i=0; i<N; i++) {
    p.next();
    outf << p.vej() << "\n";
  }

  outf.close();

  /**********************************************************************/
  cout << "Test 3: Q_d and v_ej \\propto cos(z_sun).\n";

  f = string(OUTBASENAME) + "cos-cos.dat";
  outf.open(f.c_str(), ios::out | ios::trunc);
  outf.precision(3);
  outf << "# Test 2: cos(z_sun) outflow.\n";
  outf << "# The sunward hemisphere has v_ej of 1 km/s\n";
  outf << "# at the subsolar point, and is attenuated by\n";
  outf << "# cos of z_sun.  Dust production rate\n";
  outf << "# is similarly attenuated.\n";
  outf << "# Output is relative to the comet, units of km/s\n";
  outf << "# The sun is located at " <<
    p.istate().r[0] / _AU << " " <<
    p.istate().r[1] / _AU << " " <<
    p.istate().r[2] / _AU << " AU\n";
  outf << "# v_x v_y v_z\n";

  parameters.pFunc("velocity cos 1.0;;;");
  pf.setup(parameters, p);
  for (int i=0; i<N; i++) {
    p.next();
    outf << p.vej() << "\n";
  }

  outf.close();

  /**********************************************************************/
  cout << "Test 4: Q_d and v_ej \\propto temp \\propto cos^0.25(z_sun).\n";

  f = string(OUTBASENAME) + "temp-temp.dat";
  outf.open(f.c_str(), ios::out | ios::trunc);
  outf.precision(3);
  outf << "# Test 4: cos^0.25(z_sun) outflow.\n";
  outf << "# The sunward hemisphere has v_ej of 1 km/s\n";
  outf << "# at the subsolar point, and is attenuated by\n";
  outf << "# cos^0.25 of z_sun.  Dust production rate\n";
  outf << "# is similarly attenuated.\n";
  outf << "# Output is relative to the comet, units of km/s\n";
  outf << "# The sun is located at " <<
    p.istate().r[0] / _AU << " " <<
    p.istate().r[1] / _AU << " " <<
    p.istate().r[2] / _AU << " AU\n";
  outf << "# v_x v_y v_z\n";

  parameters.pFunc("velocity temp 1.0; q_d temp");
  pf.setup(parameters, p);
  for (int i=0; i<N; i++) {
    p.next();
    outf << p.vej() << "\n";
  }

  outf.close();

  /**********************************************************************/
  cout << "Test 5: isotropic outflow, normally distributed.\n";

  f = string(OUTBASENAME) + "iso-norm.dat";
  outf.open(f.c_str(), ios::out | ios::trunc);
  outf.precision(3);
  outf << "# Test 5: isotropic outflow, normally distributed.\n";
  outf << "# The surface has v_ej of v0 = 1 km/s, sigma = 20%.\n";
  outf << "# Output is relative to the comet, units of km/s\n";
  outf << "# The sun is located at " <<
    p.istate().r[0] / _AU << " " <<
    p.istate().r[1] / _AU << " " <<
    p.istate().r[2] / _AU << " AU\n";
  outf << "# v_x v_y v_z\n";

  parameters.pFunc("velocity iso normal 0.3 0.2; q_d iso; suncone 0 180;");
  pf.setup(parameters, p);
  for (int i=0; i<N; i++) {
    p.next();
    outf << p.vej() << "\n";
  }

  outf.close();

  /**********************************************************************/
  cout << "Test 6: jet.\n";

  f = string(OUTBASENAME) + "jet.dat";
  outf.open(f.c_str(), ios::out | ios::trunc);
  outf.precision(3);
  outf << "# Test 6: jet.\n";
  outf << "# The jet is located at (lon, lat) = (0, 45).\n";
  outf << "# The opening angle is 20 deg.\n";
  outf << "# The velocity is isotropic, v_0 = 1 km/s.\n";  
  outf << "# The pole is located at (lam, bet) = (0, 90).\n";
  outf << "# The period is 8 hr and the simulation is 4 hr long.\n";
  outf << "# Output is relative to the comet, units of km/s\n";
  outf << "# The sun is located at " <<
    p.istate().r[0] / _AU << " " <<
    p.istate().r[1] / _AU << " " <<
    p.istate().r[2] / _AU << " AU\n";
  outf << "# v_x v_y v_z\n";

  parameters.pFunc("jet 0 45 45 8; pole 0 90; velocity iso 1.0");
  pf.setup(parameters, p);
  for (int i=0; i<N; i++) {
    p.next();
    outf << p.vej() << "\n";
  }

  outf.close();

  return EXIT_SUCCESS;
}

