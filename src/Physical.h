/***************************************************************************
  Copyright (C) 2008,2009 by Michael S. Kelley <msk@astro.umd.edu>

 ***************************************************************************/

#if !defined(__PHYSICAL)
#define __PHYSICAL 1

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include "Distribution.h"
#include "Composition.h"
#include "GrainTemperature.h"

using namespace std;

/** Handles all particle physical parameters, including any methods to
    generate them.  The inheriting function should consider whether a
    given parameter set makes sense (e.g., can $\beta$ = 0 and radius
    = 1).  updateBeta(), updateRadius(), and updateGrainDensity() can
    be used for this purpose.  */
class Physical {
 public:
  Physical();

  // physical parameter I/O
  double      beta();
  void        beta(const double);
  double      radius();
  void        radius(const double);
  Composition composition();
  void        composition(const int);
  double      graindensity();
  void        graindensity(const double);
  double      fractalDim();
  void        fractalDim(const double);
  double      minRadius();
  void        minRadius(const double);
  GrainTemperature grainT();
  void        grainT(const materials);

  void        updateBeta();
  void        updateBeta(const double);
  void        updateRadius();
  void        updateRadius(const double);
  void        updateGrainDensity();
  double      fit2beta(const double);

  // for random size generation
  Distribution radiusDist;

 private:
  Composition _composition;
  double _beta, _radius, _graindensity, _D, _a0;
  GrainTemperature _grainT;

  // grain library
  Composition geometric,
    am_carbon_d3000, am_carbon_d2857, am_carbon_d2727,
    am_carbon_d2609, am_carbon_d2500,
    am_olivine50_d3000, am_olivine50_d2857, am_olivine50_d2727,
    am_olivine50_d2609, am_olivine50_d2500;
};

#endif
