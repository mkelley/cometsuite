/***************************************************************************

  Projects rectangular coordinates onto an observer's Celestial
  Sphere.

  Copyright (C) 2005,2006,2008,2009 by Michael S. Kelley
  <msk@astro.umd.edu>

 ***************************************************************************/

#if HAVE_CONFIG_H
#  include <config.h>
#endif /* HAVE_CONFIG_H */

#include <string>
#include <algorithm>
#include <cstdlib>
#include <fstream>
#include <valarray>
#include <sstream>
#include "Vector.h"
#include "CoordTrans.h"
#include "projection.h"
#include "longlat.h"
#include "rundynamics.h"
#include "StringConv.h"

extern "C" int get_planet_xyz(double jd, int planets, double *r);
extern "C" int get_spitzer_xyz(double jd, double *r);

using namespace std;

/** Creates a projection object with the observer at the origin. */
projection::projection() {
  _ecliptic = false;
  observer[0] = observer[1] = observer[2] = 0;
  _offset.lambda = 0;
  _offset.beta = 0;
}

/** Creates a projection object with the specified heliocentric
    rectangular coordinates as the observer. */
projection::projection(Vector o) {
  _ecliptic = false;
  observer = o;
  _offset.lambda = 0;
  _offset.beta = 0;
}

/** Copies the passed projection object. */
projection::projection(const projection& p) {
  _ecliptic = p._ecliptic;
  observer[0] = p.observer[0];
  observer[1] = p.observer[1];
  observer[2] = p.observer[2];
  _offset = p._offset;
}

/** Creates a projection object for the specified observer and time
    for the default equatorial coordinate system.  Currently, Earth,
    Spitzer, and a 3-element position are supported. */
projection::projection(string name, double jd) {
  projection(name, jd, false);
}

/** Creates a projection object for the specified observer and time
    for the equatorial or ecliptic coordinate systems.  Currently,
    Earth, Spitzer, and a 3-element position are supported. */
projection::projection(string name, double jd, bool ec) {
  string obs = name;
  transform(obs.begin(), obs.end(), obs.begin(), (int(*)(int))toupper);

  Vector r;
  if (obs.find("EARTH") != string::npos) {
    double RHp[27];
    get_planet_xyz(jd, 4, RHp);
    r = Vector(RHp+6);
  } else if (obs.find("SPITZER") != string::npos) {
    get_spitzer_xyz(jd, r.dblarr());
  } else {
    valarray<double> rr = StringConv(name).toValarray<double>();
    if (rr.size() != 3) {
      r[0] = r[1] = r[2] = 0;
    } else {
      r[0] = rr[0] * _AU;
      r[1] = rr[1] * _AU;
      r[2] = rr[2] * _AU;
    }
  }

  observer[0] = r[0];
  observer[1] = r[1];
  observer[2] = r[2];

  _offset.lambda = 0;
  _offset.beta = 0;

  _ecliptic = ec;
}

/** Returns the RA and Dec (or ecliptic lambda, beta) of an object for
    the observer. */
longlat projection::observe(const Vector target) {
  longlat coord;

  if (_ecliptic) {
    coord = getEcliptic(observer, target);
  } else {
    coord = eclipticToEqJ2000(getEcliptic(observer, target));
  }

  coord.lambda += _offset.lambda;
  coord.beta += _offset.beta;

  return coord;
}

/** Return the observer's heliocentric ecliptic coordinates. */
Vector projection::r() { return observer; }

/** Setup absolute positional offsets in units of degrees. */
void projection::offset(const valarray<float> off) {
  _offset.lambda = off[0];
  _offset.beta = off[1];
}

/** Setup absolute positional offsets in units of degrees. */
void projection::offset(const longlat off) {
  _offset = off;
}
