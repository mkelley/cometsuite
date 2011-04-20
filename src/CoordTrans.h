/***************************************************************************
  Copyright (C) 2005,2008 by Michael S. Kelley <msk@astro.umd.edu>

  ***************************************************************************/

#if !defined(__COORDTRANS)
#define __COORDTRANS 1

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <fstream>
#include "Vector.h"
#include "longlat.h"

using namespace std;

/** Astronomical coordinate transformation routines. */
class CoordTrans {
 public:
  CoordTrans();

  // methods
  longlat getEcliptic(const Vector, const Vector);
  longlat eclipticToEqJ2000(const longlat);
  Vector  longlatToVector(const longlat);
  Vector  longlatToVector(const double, const double);
};

#endif
