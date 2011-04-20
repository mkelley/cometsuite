/***************************************************************************

  Reads in accutest[123].xyz and checks the accuracy.

  Copyright (C) 2005,2008-2010 by Michael S. Kelley <msk@astro.umd.edu>

 ***************************************************************************/

#if HAVE_CONFIG_H
#  include <config.h>
#endif /* HAVE_CONFIG_H */

#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <algorithm>
#include <cstdlib>
#include <cstring>
#include "particle.h"
#include "xyzstream.h"
#include "Vector.h"
#include "state.h"
#include "getxyz.h"

#define _INFILE1 "accutest1.xyz"
#define _OUTFILE1 "accutest1.dat"
#define _COMET1 "schwassmann-wachmann 1"
#define _KERNEL1 "schwassmannwachmann1.bsp"

#define _INFILE2 "accutest2.xyz"
#define _OUTFILE2 "accutest2.dat"
#define _COMET2 "neujmin 1"
#define _KERNEL2 "neujmin1.bsp"

#define _INFILE3 "accutest3.xyz"
#define _OUTFILE3 "accutest3.dat"
#define _COMET3 "neujmin 1"
#define _KERNEL3 "neujmin1.bsp"

using namespace std;

//extern "C" int get_comet_xyz(const char *comet, const char *cometSPK, int npts,
//			 double *jd, double *r, double *v);

int main() {
  // Comet 1
  xyzstream xyzfile;
  xyzfile.xyzopen(_INFILE1, xyzstream::READ);
  if (xyzfile.fail()) {
    cerr << "Error opening " << _INFILE1 << endl;
    return EXIT_FAILURE;
  }

  // get the parameters
  paramSet parameters;
  string header = xyzfile.readHeader();
  stringstream str;
  str << header;
  parameters.loadParameters(str);
  xyzfile.initData(parameters);

  // get the number of particles
  long max;
  if (parameters.isSyndynes()) {
    max = static_cast<long>(parameters.beta().size() * parameters.steps());
  } else {
    max = parameters.nParticles();
  }

  // get the comet's position
  double jd = parameters.obsDate();
  double r[3], v[3];
  char kfile[255];
  strcpy(kfile, _KERNEL1);
  findkernel(_COMET1, kfile);
  Vector cometR, cometV;
  get_comet_xyz(_COMET1, kfile, 1, &jd, r, v);
  cometR = Vector(r);
  cometV = Vector(v);

  // try and read all the particles
  long n;
  double deltaR, deltaV;
  ofstream outf(_OUTFILE1, ios::out | ios::trunc);
  cout << _COMET1 << "\n";
  if (max <= 21) cout << "age(yr) delta_r(km) delta_v(km/s)\n";
  outf << "# age(yr) delta_r(km) delta_v(km/s)\n";
  particle p;
  for (n=0; n<max; n++) {
    xyzfile.readParticle(p);
    if (p.error) break;

    if (max <= 21) cout << n << " " << p.age() / 86400 / 365.25 <<
		     " " << (p.fstate().r - cometR).length() << " " <<
		     (p.fstate().v - cometV).length() << "\n";
    outf << n << " " << p.age() / 86400 / 365.25 <<
      " " << (p.fstate().r - cometR).length() << " " <<
      (p.fstate().v - cometV).length() << "\n";
  }

  outf.close();
  cout << "output written to " << _OUTFILE1 << "\n";
  xyzfile.close();

  // Comet 2
  xyzfile.xyzopen(_INFILE2, xyzstream::READ);
  if (xyzfile.fail()) {
    cerr << "Error opening " << _INFILE2 << endl;
    return EXIT_FAILURE;
  }

  // get the parameters
  header = xyzfile.readHeader();
  str.str("");
  str << header;
  parameters.loadParameters(str);
  xyzfile.initData(parameters);

  // get the number of particles
  if (parameters.isSyndynes()) {
    max = static_cast<long>(parameters.beta().size() * parameters.steps());
  } else {
    max = parameters.nParticles();
  }

  // get the comet's position
  jd = parameters.obsDate();
  kfile[0] = '\0';
  strcpy(kfile, _KERNEL2);
  findkernel(_COMET2, kfile);
  get_comet_xyz(_COMET2, kfile, 1, &jd, r, v);
  cometR = Vector(r);
  cometV = Vector(v);

  // try and read all the particles
  outf.open(_OUTFILE2, ios::out | ios::trunc);
  cout << _COMET2 << "\n";
  if (max <= 21) cout << "particle age(yr) delta_r(km) delta_v(km/s)\n";
  outf << "# particle age(yr) delta_r(km) delta_v(km/s)\n";
  for (n=0; n<max; n++) {
    xyzfile.readParticle(p);
    if (p.error) break;

    if (max <= 21) cout << n << " " << p.age() / 86400 / 365.25 <<
		     " " << (p.fstate().r - cometR).length() << " " <<
		     (p.fstate().v - cometV).length() << "\n";
    outf << n << " " << p.age() / 86400 / 365.25 <<
      " " << (p.fstate().r - cometR).length() << " " <<
      (p.fstate().v - cometV).length() << "\n";
  }

  outf.close();
  cout << "output written to " << _OUTFILE2 << "\n";
  xyzfile.close();

  // Comet 3
  xyzfile.xyzopen(_INFILE3, xyzstream::READ);
  if (xyzfile.fail()) {
    cerr << "Error opening " << _INFILE3 << endl;
    return EXIT_FAILURE;
  }

  // get the parameters
  header = xyzfile.readHeader();
  str.str("");
  str << header;
  parameters.loadParameters(str);
  xyzfile.initData(parameters);

  // get the number of particles
  if (parameters.isSyndynes()) {
    max = static_cast<long>(parameters.beta().size() * parameters.steps());
  } else {
    max = parameters.nParticles();
  }

  // get the comet's position
  jd = parameters.obsDate();
  kfile[0] = '\0';
  strcpy(kfile, _KERNEL3);
  findkernel(_COMET3, kfile);
  get_comet_xyz(_COMET3, kfile, 1, &jd, r, v);
  cometR = Vector(r);
  cometV = Vector(v);

  // try and read all the particles
  outf.open(_OUTFILE3, ios::out | ios::trunc);
  cout << _COMET3 << "\n";
  if (max <= 21) cout << "particle age(yr) delta_r(km) delta_v(km/s)\n";
  outf << "# particle age(yr) delta_r(km) delta_v(km/s)\n";
  for (n=0; n<max; n++) {
    xyzfile.readParticle(p);
    if (p.error) break;

    if (max <= 21) cout << n << " " << p.age() / 86400 / 365.25 <<
		     " " << (p.fstate().r - cometR).length() << " " <<
		     (p.fstate().v - cometV).length() << "\n";
    outf << n << " " << p.age() / 86400 / 365.25 <<
      " " << (p.fstate().r - cometR).length() << " " <<
      (p.fstate().v - cometV).length() << "\n";
  }

  outf.close();
  cout << "output written to " << _OUTFILE3 << "\n";
  xyzfile.close();

  return EXIT_SUCCESS;
}
