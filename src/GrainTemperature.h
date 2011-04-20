/***************************************************************************
  Copyright (C) 2009 by Michael S. Kelley <msk@astro.umd.edu>

  ***************************************************************************/

#if !defined(__TEMPERATUREFIT)
#define __TEMPERATUREFIT 1

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <valarray>
#include "Composition.h"

using namespace std;

/** Provides a manner for which to compute grain temperatures with
    stored T(D, a, rh) matrices. */
class GrainTemperature {
 public:
  GrainTemperature();

  void load(const materials);

  double interp(const double, const double, const double);

 private:
  // for interpolation
  int* nearest(valarray<float>, const double);

  Composition _composition;
  valarray<float> _radius, _D, _rh;
  valarray<double> _T;
  int naxis[3];
  bool dataLoaded;
};

#endif
