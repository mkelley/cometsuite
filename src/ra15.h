/***************************************************************************
  Copyright (C) 2004,2005,2006,2007,2008 by Michael S. Kelley
  <msk@astro.umd.edu>

 ***************************************************************************/

#if !defined(__RA15)
#define __RA15 1

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include "Integrator.h"
#include "particle.h"
#include "paramSet.h"
#include "getxyz.h"

using namespace std;

/* r, v, tFinal, jd, beta, planets, tol, minStep,
   number of function calls, number of time steps. */
extern "C" int everhart(double*, double*, double, double, double, const int,
			double, double*, double*, double*, const int,
			const int);

/** Provides ra15::integrate() to Integrator. */
class ra15 : public Integrator {
 public:
  void integrate();
};

#endif
