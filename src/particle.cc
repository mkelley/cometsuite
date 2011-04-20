/***************************************************************************

  A class to contain the total description of a particle.

  Copyright (C) 2004-2010 by Michael S. Kelley <msk@astro.umd.edu>

 ***************************************************************************/

#if HAVE_CONFIG_H
#  include <config.h>
#endif /* HAVE_CONFIG_H */

#include <iostream>
#include <cstring>
#include <string>
#include <cmath>
#include "paramSet.h"
#include "particle.h"
#include "Vector.h"
#include "state.h"
#include "rundynamics.h"

using namespace std;

/** Create a blank particle.  All parameters are set to zero or similar. */
particle::particle() {
  // Physical attributes
  beta(0);
  radius(0);
  graindensity(1);
  composition(GEOMETRIC);
  fractalDim(3.0);
  minRadius(0.1);

  // Dynamical attributes
  longlat zeroLL = {0, 0};
  Vector zeroV = Vector(0, 0, 0);
  state zeroS = {zeroV, zeroV, 0};
  istate(zeroS);
  fstate(zeroS);
  origin(zeroLL);
  vejGen.vDist.min(0);
  vejGen.vDist.max(0);
  vejGen.vDist.mu(1);
  vejGen.vDist.sigma(0.1);
  vejGen.range();
  rhlimit(-1);

  // Other attributes
  label(string(""));

  error = false;
}

/** Create a particle with this particle's parameters. */
particle::particle(const particle& p) { *this = p; }

/** Set the particle's string label. */
void particle::label(const string s) { _label = s; }

/** Set the particle's string label. */
void particle::label(char *s) { _label = s; }

/** Return the particle's string label. */
string particle::label() { return _label; }

/** Generate the particle's string label.  Only the first 16
    characters will be saved by xyzstream. */
void particle::generateLabel(const long n) {
  char l[255];
  sprintf(l, parameters().labelFormat().c_str(), n);
  _label = l;
}

/** Create a new particle, but a unique label won't be generated. */
void particle::next() { next(0); }

/** Create a new particle.
    -# Generate a new age and radius/beta/grain density:
      - Repeat until the grain is within the rh limit.
    -# Get the new ejection velocity.
    -# Set the origin on the nucleus.
    -# Update the initial state vector.
    -# Generate a new label.
 */
void particle::next(const long n) {
  double jd;
  double rc[3], vc[3];
  longlat o;
  state st;

  if (DEBUG) cerr << "particle::next()\n";

  age(ageDist.next());
  radius(radiusDist.next());
  updateGrainDensity();
  updateBeta();

  jd = parameters().obsDate() - age() / 86400;
  get_comet_xyz(parameters().comet().c_str(),
		parameters().spkKernel().c_str(),
		1, &jd, rc, vc);

  st.r = rc;
  st.v = vc;
  st.t = istate().t;  // istate().t already set above via age()

  // get an ejection velocity
  vej(vejGen.next(radius(), st.r, st.t));

  // update the initial state vector
  st.v += vej();

  // compute the origin on the nucleus
  if (vej().length() > 0) {
    o.lambda = vej().unit() * poleX() * 180.0 / M_PI + 180.0;
    o.beta = vej().unit() * poleV();
    o.beta = 90.0 - acos(o.beta) * 180.0 / M_PI;
    origin(o);
  }

  // update the particle
  istate(st);

  // generate a label
  generateLabel(n);

  if (DEBUG) cerr << "  - new particle created\n";
  if (DEBUG) cerr << "  - age beta radius= " << age() / 86400 << " " << beta() << " " << radius() << "\n";
}

/** Return the current parameter set. */
paramSet particle::parameters() { return _parameters; }

/** Set the current parameters to integrate by. */
void particle::parameters(paramSet p) { _parameters = p; }

/** Return the grain's temperature, computed via interpolation of T(D,
    a, rh) FITS files. */
double particle::temperature() {
  double rh = fstate().r.length();

  if (composition().name == GEOMETRIC)
    return 278. / sqrt(rh);
}

/** << for stream output. */
ostream& operator<<(ostream& os, particle p) {
  os << "beta: " << p.beta() << " radius: " << p.radius() << endl;
  os << "grain density: " << p.graindensity() << endl;
  os << "age: " << p.age() / 86400 << endl;
  os << "v_ej: " << p.vej().length()  << endl;
  os << "v_ej: " << p.vej()           << endl;
  os << "initial r: " << p.istate().r << endl;
  os << "initial v: " << p.istate().v << endl;
  os << "initial t: " << p.istate().t << endl;
  os << "final r: " << p.fstate().r << endl;
  os << "final v: " << p.fstate().v << endl;
  os << "final t: " << p.fstate().t << endl;
  os << "label: " << p.label() << endl;

  return os;
}
