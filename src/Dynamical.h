/***************************************************************************
  Copyright (C) 2008,2012 by Michael S. Kelley <msk@astro.umd.edu>

  ***************************************************************************/

#if !defined(__DYNAMICAL)
#define __DYNAMICAL 1

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include "state.h"
#include "longlat.h"
#include "Vej.h"
#include "Distribution.h"
#include "CoordTrans.h"

using namespace std;

/** Handles all particle dynamical parameters, including methods to
    generate them.
    - the particle's initial and final states: position, velocity, and time
    - the particle's origin (longitude, latitude) on the nucleus
    - the particle's ejection velocity
    - the heliocentric distance limit of comet activity
*/
class Dynamical : CoordTrans {
 public:
  Dynamical();

  // dynamical parameter I/O
  double  age();
  void    age(const double);
  state   istate();              // initial state
  void    istate(const state);
  state   fstate();              // final state
  void    fstate(const state);
  longlat pole();
  void    pole(const longlat);
  Vector  poleV();
  void    poleV(const Vector);
  Vector  poleX();
  Vector  poleY();
  longlat origin();
  void    origin(const longlat);
  Vector  vej();
  void    vej(const Vector);
  double  rhlimit();
  void    rhlimit(const double);

  // dynamical parameter generation
  Distribution ageDist;
  Vej vejGen;

 private:
  double _rhlimit;
  longlat _pole, _origin;
  state _istate, _fstate;
  Vector _vej, _poleV, _poleX, _poleY;
};

#endif
