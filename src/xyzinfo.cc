/***************************************************************************

  Provides basic information on a CometSuite data file.

  Copyright (C) 2005-2010,2012 by Michael S. Kelley <msk@astro.umd.edu>

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
#include <cmath>
#include "xyzstream.h"
#include "particle.h"
#include "rundynamics.h"

#define SUBPROJECT "xyzinfo"

using namespace std;

int main(int argc, char *argv[])
{
  xyzstream xyzfile;
  if (argc > 1) {
    xyzfile.xyzopen(argv[1], xyzstream::READ);
    if (xyzfile.fail()) {
      cerr << "Error opening " << argv[1] << endl;
      return EXIT_FAILURE;
    }
  } else {
    cerr << PACKAGE_STRING << "\n";
    cerr << SUBPROJECT << " Copyright (C) 2005-2010,2012 Michael S Kelley\n";
    cerr << "Usage:\n";
    cerr << "  " << argv[0] << " xyzfile.xyz\n\n";
    return EXIT_FAILURE;
  }

  // get the parameters
  paramSet parameters;
  stringstream str;
  string header;

  if (DEBUG) cerr  << "Setting parameters from file (xyzinfo)\n";
  header = xyzfile.readHeader();
  str << header;
  parameters.loadParameters(str);
  xyzfile.initData(parameters);
  if (DEBUG) cerr  << "done\n";

  // output header information
  cout << header << endl;

  // get number of particles for syndynes or make comet / integrate xyz
  long max;
  if (parameters.isSyndynes()) {
    max = static_cast<long>(parameters.beta().size() * parameters.steps());
  } else {
    max = parameters.nParticles();
  }

  // try and read all the particles
  if (DEBUG) cerr << "Attempting to read all particles\n";
  long n;
  double lbeta_avg = 0, beta_avg = 0, radius_avg = 0,
    age_avg = 0, rh_avg = 0, rho_avg = 0;
  particle p;

  // read in grains until the expected number is reached or an error
  // occurs; update the cumulative sums after every grain
  for (n=0; n<max; n++) {
    p.radius(0);  // default condition
    xyzfile.readParticle(p);
    if (DEBUG) cerr << n << " " << p << endl;

    if (p.error) {
      cout << "Particle read error... exiting.\n\n";
      break;
    }

    lbeta_avg += log10(p.beta());
    beta_avg += p.beta();

    // Syndynes and pre-v0.7.3 xyz files may not have a correctly
    // computed radius.  In this case, compute radius from $\beta$.
    if (p.radius() <= 0) {
      p.updateRadius();
      p.updateGrainDensity();
    }

    radius_avg += p.radius();
    rho_avg += p.graindensity();
    age_avg += p.age();
    rh_avg += p.fstate().r.length() / _AU;
  }

  // Convert cumulative sums to averages
  double dn = static_cast<double>(n);
  lbeta_avg = pow(10, lbeta_avg / dn);
  beta_avg /= dn;
  radius_avg /= dn;
  rho_avg /= dn;
  age_avg /= dn;
  rh_avg /= dn;

  cout << "Average log(beta): " << lbeta_avg << endl;
  cout << "Average beta: " << beta_avg << endl;
  cout << "Average radius (micron): " << radius_avg << endl;
  cout << "Average density (g/cm3): " << rho_avg << endl;
  cout << "Average age (days): " << age_avg / 86400 << endl;
  cout << "Average rh (AU): " << rh_avg << endl;
  cout.precision(12);
  cout << "\nSuccessfully read " << n << " out of " << max <<
    " particles\n\n";

  if (n < max) return EXIT_FAILURE;
  return EXIT_SUCCESS;
}

