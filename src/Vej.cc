/***************************************************************************

  Implements velocity ejection descriptions.

  Copyright (C) 2007-2010 by Michael S. Kelley <msk@astro.umd.edu>

 ***************************************************************************/

#if HAVE_CONFIG_H
#  include <config.h>
#endif /* HAVE_CONFIG_H */

#include <cmath>
#include "Vej.h"
#include "Physical.h"
#include "Distribution.h"
#include "Vector.h"
#include "longlat.h"
#include "state.h"
#include "rundynamics.h"

using namespace std;

Vej::Vej() {
  _v0min = 0.0; // km/s
  _v0 = 1.0; // km/s
  _axis = Vector(1, 0, 0);  // AU
  _Speed = RANGE;
  _Dir = SUNWARD;
  _u1 = 0.5;
  _u2 = 0.5;

  // arbitrary values
  thetaDist.min(0);
  thetaDist.max(90);
  thetaDist.distribution(Distribution::ISOSPHERE);
  phiDist.min(0);
  phiDist.max(360);
  phiDist.distribution(Distribution::LINEAR);
  vDist.min(0);
  vDist.max(1.0);
  vDist.mu(1.0);
  vDist.sigma(0.1);
  vDist.distribution(Distribution::LINEAR);
}

/** Set the v0 parameter (km/s). */
void Vej::v0(const double vv) { _v0 = vv; }
/** Return the v0 parameter (km/s). */
double Vej::v0() { return _v0; }

/** Set the v0min parameter (km/s), used in simpleActivtyRange. */
void Vej::v0min(const double vv) { _v0min = vv; }
/** Return the v0min parameter (km/s), used in simpleActivityRange. */
double Vej::v0min() { return _v0min; }

/** Set the mu parameter (km/s). */
void Vej::mu(const double mu) { vDist.mu(mu); }
/** Return the mu parameter (km/s). */
double Vej::mu() { return vDist.mu(); }

/** Set the sigma parameter (fraction of mu). */
void Vej::sigma(const double ss) { vDist.sigma(ss); }
/** Return the sigma parameter (fraction of mu). */
double Vej::sigma() { return vDist.sigma(); }

/** Set the u1 parameter ($v_{ej} \propto a^{-u1}$). */
void Vej::u1(const double u) { _u1 = u; }
/** Return the u1 parameter. */
double Vej::u1() { return _u1; }

/** Set the u2 parameter ($v_{ej} \propto rh^{-u2}$). */
void Vej::u2(const double u) { _u2 = u; }
/** Return the u2 parameter. */
double Vej::u2() { return _u2; }

/** Set the axis of symmetry for grain ejection.  This is typically
    the direction of the sun or the north pole of the nucleus. */
void Vej::axis(const Vector a) { _axis = a; _axis = _axis.unit(); }
/** Set the axis of symmetry for grain ejection.  This is typically
    the direction of the sun or the north pole of the nucleus. */
void Vej::axis(const longlat ll) {
  _axis = longlatToVector(ll);
  _axis = _axis.unit();
}
/** Return the axis of symmetry for grain ejection. */
Vector Vej::axis() { return _axis; }

/** Set the rotation period for jets in units of seconds.  Internally,
    this also sets _angfreq (radians per s). */
void Vej::period(const double p) { _period = p; _angfreq = 2.0 * M_PI / p; }
/** Return the rotation period. */
double Vej::period() { return _period; }

/** The ejection velocity is independent of grain parameters and
    ranged from vMin() to vMax(). */
void Vej::range() {
  _Speed = RANGE;
  vDist.distribution(Distribution::LINEAR);
}

/** The ejection velocity is independent of grain parameters and
    ranged from vMin() to vMax(), but attenuated by the cos of the
    sun-zenith angle. */
void Vej::cosRange() {
  _Speed = COSRANGE;
  vDist.distribution(Distribution::LINEAR);
}

/** The ejection velocity is independent of grain parameters and
    ranged from vMin() to vMax(), but attenuated by the STM/NEATM
    temperature of the surface. */
void Vej::tempRange() {
  _Speed = TEMPRANGE;
  vDist.distribution(Distribution::LINEAR);
}

/** The ejection velocity is normally distributed about v0. */
void Vej::normal() {
  _Speed = SIMPLEACTIVITYNORMAL;
  vDist.distribution(Distribution::NORMAL);
}

/** The ejection velocity is normally distributed about v0, attenuated
    by the cos of the sun-zenith angle. */
void Vej::cosNormal() {
  _Speed = COSNORMAL;
  vDist.distribution(Distribution::NORMAL);
}

/** The ejection velocity is normally distributed about v0, attenuated
    by the STM/NEATM temperature of the surface. */
void Vej::tempNormal() {
  _Speed = TEMPNORMAL;
  vDist.distribution(Distribution::NORMAL);
}

/** The ejection velocity depends on simple coma activity. */
void Vej::simpleActivity() {
  _Speed = SIMPLEACTIVITY;
  vDist.distribution(Distribution::LINEAR);
}

/** The ejection velocity depends on simple coma activity, attenuated
    by te cos of the sun-zenith angle. */
void Vej::cosSimpleActivity() {
  _Speed = COSSIMPLEACTIVITY;
  vDist.distribution(Distribution::LINEAR);
}

/** The ejection velocity depends on simple coma activity, attenuated
    by the STM/NEATM temperature of the surface. */
void Vej::tempSimpleActivity() {
  _Speed = TEMPSIMPLEACTIVITY;
  vDist.distribution(Distribution::LINEAR);
}

/** The ejection velocity ranges from v0min to v0max. */
void Vej::simpleActivityRange() {
  _Speed = SIMPLEACTIVITYRANGE;
  vDist.distribution(Distribution::LINEAR);
}

/** The ejection velocity ranges from v0min to v0max, attenuated by
    the cos of the sun-zenith angle. */
void Vej::cosSimpleActivityRange() {
  _Speed = COSSIMPLEACTIVITYRANGE;
  vDist.distribution(Distribution::LINEAR);
}

/** The ejection velocity ranges from v0min to v0max, attenuated by
    the STM/NEATM temperature of the surface. */
void Vej::tempSimpleActivityRange() {
  _Speed = TEMPSIMPLEACTIVITYRANGE;
  vDist.distribution(Distribution::LINEAR);
}

/** The ejection velocity is normally distributed about v0. */
void Vej::simpleActivityNormal() {
  _Speed = SIMPLEACTIVITYNORMAL;
  vDist.distribution(Distribution::NORMAL);
}

/** The ejection velocity is normally distributed about v0, attenuated
    by the cos of the sun-zenith angle. */
void Vej::cosSimpleActivityNormal() {
  _Speed = COSSIMPLEACTIVITYNORMAL;
  vDist.distribution(Distribution::NORMAL);
}

/** The ejection velocity is normally distributed about v0, attenuated
    by the STM/NEATM temperature of the surface. */
void Vej::tempSimpleActivityNormal() {
  _Speed = TEMPSIMPLEACTIVITYNORMAL;
  vDist.distribution(Distribution::NORMAL);
}

/** Isotropic dust emission per solid angle. */
void Vej::isoSphere() {
  _Dir = ISOSPHERE;
  thetaDist.distribution(Distribution::ISOSPHERE);
}

/** Dust emission proportional to .... */
void Vej::cosSphere() {
  _Dir = COSSPHERE;
  thetaDist.distribution(Distribution::COSSPHERE);
}

/** Dust emission proportional to .... */
void Vej::tempSphere() {
  _Dir = TEMPSPHERE;
  thetaDist.distribution(Distribution::TEMPSPHERE);
}

/** Jet.  Dust is uniform over the jet's opening angle. */
void Vej::jet() {
  _Dir = JET;
  thetaDist.distribution(Distribution::ISOSPHERE);
}

/** Returns the magnitude of the velocity of ejection, independent of
    the direction. */
double Vej::getSpeed(const double radius, Vector r) {
  switch (_Speed) {
  case RANGE:
  case COSRANGE:
  case TEMPRANGE: break;  // min/max should already be set
  case NORMAL:
  case COSNORMAL:
  case TEMPNORMAL: break; // min/mu/sigma should already be set
  case SIMPLEACTIVITY:
  case COSSIMPLEACTIVITY:
  case TEMPSIMPLEACTIVITY:
    if (radius <= 0) {
      vDist.min(0);
      vDist.max(0);
    } else {
      vDist.min(_v0 * exp(-_u1 * log(radius)) * exp(-_u2 * log(r.length() / _AU)));
      vDist.max(vDist.min());
    }
    break;
  case SIMPLEACTIVITYRANGE:
  case COSSIMPLEACTIVITYRANGE:
  case TEMPSIMPLEACTIVITYRANGE:
    if (radius <= 0) {
      vDist.min(0);
      vDist.max(0);
    } else {
      double a = exp(-_u1 * log(radius)) * exp(-_u2 * log(r.length() / _AU));
      vDist.min(_v0min * a);
      vDist.max(   _v0 * a);
    }
    break;
  case SIMPLEACTIVITYNORMAL:
  case COSSIMPLEACTIVITYNORMAL:
  case TEMPSIMPLEACTIVITYNORMAL:
    if (radius <= 0) {
      vDist.min(0);
      vDist.mu(0);
    } else {
      vDist.min(0);
      vDist.mu(_v0 * exp(-_u1 * log(radius)) * exp(-_u2 * log(r.length() / _AU)));
    }
    break;
  }
  return vDist.next();
}

/** Returns the direction of the velocity of ejection.  r is the
    position vector, and t is the time in seconds.  If the axis of
    symmetry is the Sun, the initial state vector (istate) should
    already be set and will be used to get the sunward direction.
    Otherwise, axis() will be used as is. */
Vector Vej::getDirection(Vector r, const double t) {
  double phi;
  /* 1) rotate axis th radians around the perp vector
     2) rotate the result ph radians around the axis of symmetry
     3) make sure it is a unit vector */

  // Determine the axis of symmetry, if needed
  switch(_Dir) {
  case SUNWARD:
  case ISOSPHERE:
  case COSSPHERE:
  case TEMPSPHERE:
    axis(getSunward(r));
    phi = phiDist.next();
    break;
  case JET:
    phi = phiDist.next() + _angfreq * t;
    break;
  }

  Vector result;
  result = _axis.rotate(getPerp(_axis), thetaDist.next());
  result = result.rotate(_axis, phi);
  result = result.unit();
  return result;
}

/** Generate a new ejection velocity (calls both getSpeed() and
    getDirection()).  a is the grain size, r is the position vector,
    and t is time in seconds.  Attenuate getSpeed() if needed.  */
Vector Vej::next(const double a, Vector r, const double t) {
  Vector v = getDirection(r, t);
  double s = getSpeed(a, r);
  double w = 1;

  // attenuate if needed
  switch(_Speed) {
  case COSRANGE:
  case COSNORMAL:
  case COSSIMPLEACTIVITY:
  case COSSIMPLEACTIVITYRANGE:
  case COSSIMPLEACTIVITYNORMAL:
    w = v.cosangle(getSunward(r));
    break;
  case TEMPRANGE:
  case TEMPNORMAL:
  case TEMPSIMPLEACTIVITY:
  case TEMPSIMPLEACTIVITYRANGE:
  case TEMPSIMPLEACTIVITYNORMAL:
    w = exp(0.25 * log(v.cosangle(getSunward(r))));
    break;
  }

  return v * s * w;
}

/** Return a vector perpendicular to the input vector. */
Vector Vej::getPerp(Vector n) {
  if ((n[1] == n[2]) && (n[2] == 0)) {
    // if the input vector is the x-axis, use the z-axis
    return Vector(0, 0, 1);
  } else {
    // use the cross-product of the vector and the x-axis
    return Vector(0, -n[2], n[1]).unit();
  }
}

/** Compute the sunward unit vector given the passed state. */
Vector Vej::getSunward(Vector r) {
  return (r * -1.0).unit();
}
