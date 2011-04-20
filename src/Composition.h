/***************************************************************************
  Copyright (C) 2009 by Michael S. Kelley <msk@astro.umd.edu>

  ***************************************************************************/

#if !defined(__COMPOSITION)
#define __COMPOSITION 1

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

using namespace std;

/// Material list
enum materials { GEOMETRIC, AM_CARBON, AM_OLIVINE50 };

/** Provides a framework to implement different materials for
    rundynamics.  Composition parameters should be defined in
    Physical::Physical(). */
struct Composition {
  materials name;
  double bulkdensity;

  // for radius-beta conversion
  double scale1000;
  double slope1000;
  double fit[11];
};

#endif
