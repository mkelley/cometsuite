/***************************************************************************

  Everything related to the description of a particle's dynamics.

  Copyright (C) 2008,2010,2012 by Michael S. Kelley <msk@astro.umd.edu>

 ***************************************************************************/

#if HAVE_CONFIG_H
#  include <config.h>
#endif /* HAVE_CONFIG_H */

#include "Dynamical.h"
#include "Vector.h"
#include "state.h"
#include "longlat.h"

using namespace std;

/** Constructor. */
Dynamical::Dynamical() {
  longlat zeroO = {0, 0};
  pole(zeroO);
  _origin = zeroO;

  Vector zeroV(0, 0, 0);
  _vej = zeroV;
  _istate.r = zeroV;
  _istate.v = zeroV;
  _istate.t = 0;
  _fstate.r = zeroV;
  _fstate.v = zeroV;
  _fstate.t = 0;
}

/** Return the age in seconds (at the time of observation).  This is a
    convenience function since the age is stored in istate(). */
double Dynamical::age() { return -_istate.t; }
/** Change the age (at the time of observation) to the value of a in
    seconds.  This is a convenience function since the age is stored
    in istate(). */
void Dynamical::age(const double a) { _istate.t = -a; }

/** Return the initial state vector. */
state Dynamical::istate() { return _istate; }
/** Set the initial state vector. */
void Dynamical::istate(state s) { _istate = s; }

/** Return the final state vector. */
state Dynamical::fstate() { return _fstate; }
/** Set the final state vector. */
void Dynamical::fstate(state s) { _fstate = s; }

/** Return the pole solution of the nucleus (ecliptic coordinates). */
longlat Dynamical::pole() { return _pole; }
/** Set the pole solution of the nucleus (ecliptic coordinates). */
void Dynamical::pole(longlat ll) {
  _pole = ll;
  _poleV = longlatToVector(ll).unit();

  // Define the vector at zero longitude on the equator
  if ((_poleV[1] == _poleV[2]) && (_poleV[2] == 0)) {
    // if the input vector is the x-axis, use the z-axis
    _poleX = Vector(0, 0, 1);
  } else {
    Vector X(1, 0, 0);
    _poleX = (X - _poleV * (X * _poleV)).unit();
  }

  _poleY = _poleV % _poleX;
}

/** Return the pole solution of the nucleus as a vector (ecliptic
    coordinates). */
Vector Dynamical::poleV() { return _poleV; }
/** Set the pole solution of the nucleus given a vector (ecliptic
    coordinates). */
void Dynamical::poleV(Vector v) {
  _poleV = v.unit();

  // Define the vector at zero longitude on the equator
  if ((_poleV[1] == _poleV[2]) && (_poleV[2] == 0)) {
    // if the input vector is the x-axis, use the z-axis
    _poleX = Vector(0, 0, 1);
  } else {
    Vector X(1, 0, 0);
    _poleX = (X - _poleV * (X * _poleV)).unit();
  }

  _poleY = _poleV % _poleX;
  _pole = getEcliptic(Vector(0, 0, 0), v);
}

/** Return the vector at zero longitude on the equator (ecliptic
    coordinates). */
Vector Dynamical::poleX() { return _poleX; }

/** Return the vector at 90 deg longitude on the equator (ecliptic
    coordinates). */
Vector Dynamical::poleY() { return _poleY; }

/** Return the particle's origin on the nucleus. */
longlat Dynamical::origin() { return _origin; }
/** Set the particle's origin on the nucleus. */
void Dynamical::origin(longlat ll) { _origin = ll; }

/** Return the particle's ejection velocity vector. */
Vector Dynamical::vej() { return _vej; }
/** Set the particle's ejection velocity vector. */
void Dynamical::vej(const Vector v) { _vej = v; }

/** Return the rh limit. */
double Dynamical::rhlimit() { return _rhlimit; }
/** Set the rh limit. */
void Dynamical::rhlimit(const double l) { _rhlimit = l; }
