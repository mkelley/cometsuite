/***************************************************************************
  Copyright (C) 2005,2006,2009,2010 by Michael S. Kelley <msk@astro.umd.edu>

 ***************************************************************************/

#if !defined(__PROJECTION)
#define __PROJECTION 1

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <string>
#include <fstream>
#include <valarray>
#include "Vector.h"
#include "CoordTrans.h"
#include "longlat.h"

using namespace std;

/** Transforms a position into RA and Dec for a particular
    observer. */
class projection : public CoordTrans {
 public:
  projection();
  projection(Vector);
  projection(const projection&);
  projection(string, double);
  projection(string, double, bool);
  projection(string, double, bool, valarray<float>);

  // methods
  void offset(const valarray<float>);
  void offset(const longlat);
  longlat observe(const Vector);
  Vector r();

 private:
  bool _ecliptic;
  Vector observer;
  longlat _offset;
};

#endif
