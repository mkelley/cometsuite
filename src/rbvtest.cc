/***************************************************************************

  Generates grains to test the radius-beta-v_ej relationships.

  Copyright (C) 2009 by Michael S. Kelley <msk@astro.umd.edu>

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

#define SUBPROJECT "rbvtest"
#define OUTBASENAME "rbvtest-"
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
  cout << "Geometric grains...\n";

  string f = string(OUTBASENAME) + "g.dat";
  outf.open(f.c_str(), ios::out | ios::trunc);
  outf << "# Geometric grains\n";
  outf << "# beta = 0.57 / a / rho\n";
  outf << "# a = 0.1 to 1000 microns\n";
  outf << "# rho = 1 g/cm3\n";
  outf << "# v_ej = 0.3 / sqrt(a) [km/s]\n";
  outf << "# radius beta v_ej\n";

  parameters.pFunc("simple_coma 0.3 0 -1 3");
  pf.setup(parameters, p);
  for (int i=0; i<N; i++) {
    p.next();
    outf << p.radius() << " " << p.beta() << " " << p.vej().length() << "\n";
  }

  outf.close();

  /**********************************************************************/
  cout << "Amorphous carbon grains...\n";

  f = string(OUTBASENAME) + "ac.dat";
  outf.open(f.c_str(), ios::out | ios::trunc);
  outf << "# Amorphous carbon grains\n";
  outf << "# beta = beta(Mie theory)\n";
  outf << "# a = 0.1 to 1000 microns\n";
  outf << "# v_ej = 0.3 / sqrt(a) [km/s]\n";
  outf << "# radius beta v_ej\n";

  parameters.pFunc("simple_coma 0.3 0 -1 3; composition ac");
  pf.setup(parameters, p);
  for (int i=0; i<N; i++) {
    p.next();
    outf << p.radius() << " " << p.beta() << " " << p.vej().length() << "\n";
  }

  outf.close();

  /**********************************************************************/
  cout << "Amorphous olivine grains...\n";

  f = string(OUTBASENAME) + "ol50.dat";
  outf.open(f.c_str(), ios::out | ios::trunc);
  outf << "# Amorphous olivine 50 grains\n";
  outf << "# beta = beta(Mie theory)\n";
  outf << "# a = 0.1 to 1000 microns\n";
  outf << "# v_ej = 0.3 / sqrt(a) [km/s]\n";
  outf << "# radius beta v_ej\n";

  parameters.pFunc("simple_coma 0.3 0 -1 3; composition ol50");
  pf.setup(parameters, p);
  for (int i=0; i<N; i++) {
    p.next();
    outf << p.radius() << " " << p.beta() << " " << p.vej().length() << "\n";
  }

  outf.close();

  /**********************************************************************/
  cout << "Amorphous olivine 50 grains (rho = 1 g/cm3)...\n";

  f = string(OUTBASENAME) + "ol50-rho1.dat";
  outf.open(f.c_str(), ios::out | ios::trunc);
  outf << "# Amorphous olivine 50 grains\n";
  outf << "# beta = beta(Mie theory) * 3.3 / rho\n";
  outf << "# rho = 1 g/cm3\n";
  outf << "# a = 0.1 to 1000 microns\n";
  outf << "# v_ej = 0.3 / sqrt(a) [km/s]\n";
  outf << "# radius beta v_ej\n";

  parameters.pFunc("simple_coma 0.3 0 -1 3; composition ol50; density 1");
  pf.setup(parameters, p);
  for (int i=0; i<N; i++) {
    p.next();
    outf << p.radius() << " " << p.beta() << " " << p.vej().length() << "\n";
  }

  outf.close();

  return EXIT_SUCCESS;
}

