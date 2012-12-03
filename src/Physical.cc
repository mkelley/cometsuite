/***************************************************************************

  Implements the physical description of CometSuite grains.

  Copyright (C) 2008,2009,2012 by Michael S. Kelley <msk@astro.umd.edu>

 ***************************************************************************/

#if HAVE_CONFIG_H
#  include <config.h>
#endif /* HAVE_CONFIG_H */

#include <iostream>
#include <cmath>
#include "Physical.h"
#include "GrainTemperature.h"

using namespace std;

/** Constructor. */
Physical::Physical() {
  geometric.name = GEOMETRIC;
  geometric.bulkdensity = 1.0;
  // slope1000, scale1000 and fit should never be used for geometric
  geometric.scale1000 = 0;
  geometric.slope1000 = 0;
  geometric.fit[ 0] = 0;
  geometric.fit[ 1] = 0;
  geometric.fit[ 2] = 0;
  geometric.fit[ 3] = 0;
  geometric.fit[ 4] = 0;
  geometric.fit[ 5] = 0;
  geometric.fit[ 6] = 0;
  geometric.fit[ 7] = 0;
  geometric.fit[ 8] = 0;
  geometric.fit[ 9] = 0;
  geometric.fit[10] = 0;

  am_carbon_d3000.name = AM_CARBON;
  am_carbon_d3000.bulkdensity = 2.5;
  am_carbon_d3000.scale1000 = 0.225666;
  am_carbon_d3000.slope1000 = -1.003665;
  am_carbon_d3000.fit[ 0] = -5.06527316e-07;
  am_carbon_d3000.fit[ 1] =  7.70285780e-06;
  am_carbon_d3000.fit[ 2] = -9.59153552e-06;
  am_carbon_d3000.fit[ 3] = -3.49458578e-04;
  am_carbon_d3000.fit[ 4] =  1.42119554e-03;
  am_carbon_d3000.fit[ 5] =  3.87261875e-03;
  am_carbon_d3000.fit[ 6] = -2.88931148e-02;
  am_carbon_d3000.fit[ 7] =  2.72092622e-02;
  am_carbon_d3000.fit[ 8] =  9.48081031e-02;
  am_carbon_d3000.fit[ 9] = -1.28394360e+00;
  am_carbon_d3000.fit[10] = -1.14805226e+00;

  am_carbon_d2857.name = AM_CARBON;
  am_carbon_d2857.bulkdensity = 2.5;
  am_carbon_d2857.scale1000 = 0.308570;
  am_carbon_d2857.slope1000 = -0.858520;
  am_carbon_d2857.fit[ 0] = -4.48449385e-07;
  am_carbon_d2857.fit[ 1] =  6.68803667e-06;
  am_carbon_d2857.fit[ 2] = -6.49449185e-06;
  am_carbon_d2857.fit[ 3] = -3.14459683e-04;
  am_carbon_d2857.fit[ 4] =  1.20833040e-03;
  am_carbon_d2857.fit[ 5] =  3.69275828e-03;
  am_carbon_d2857.fit[ 6] = -2.58547676e-02;
  am_carbon_d2857.fit[ 7] =  2.45987095e-02;
  am_carbon_d2857.fit[ 8] =  8.36780673e-02;
  am_carbon_d2857.fit[ 9] = -1.12885733e+00;
  am_carbon_d2857.fit[10] = -8.21246599e-01;

  am_carbon_d2727.name = AM_CARBON;
  am_carbon_d2727.bulkdensity = 2.5;
  am_carbon_d2727.scale1000 = 0.419969;
  am_carbon_d2727.slope1000 = -0.726903;
  am_carbon_d2727.fit[ 0] = -4.05392711e-07;
  am_carbon_d2727.fit[ 1] =  5.84342739e-06;
  am_carbon_d2727.fit[ 2] = -2.78814574e-06;
  am_carbon_d2727.fit[ 3] = -2.91213810e-04;
  am_carbon_d2727.fit[ 4] =  9.99270004e-04;
  am_carbon_d2727.fit[ 5] =  3.80060169e-03;
  am_carbon_d2727.fit[ 6] = -2.33328448e-02;
  am_carbon_d2727.fit[ 7] =  1.96509102e-02;
  am_carbon_d2727.fit[ 8] =  7.77117690e-02;
  am_carbon_d2727.fit[ 9] = -9.76768785e-01;
  am_carbon_d2727.fit[10] = -5.29556032e-01;

  am_carbon_d2609.name = AM_CARBON;
  am_carbon_d2609.bulkdensity = 2.5;
  am_carbon_d2609.scale1000 = 0.562876;
  am_carbon_d2609.slope1000 = -0.610334;
  am_carbon_d2609.fit[ 0] = -3.84651390e-07;
  am_carbon_d2609.fit[ 1] =  5.29694435e-06;
  am_carbon_d2609.fit[ 2] =  1.16052167e-06;
  am_carbon_d2609.fit[ 3] = -2.84384116e-04;
  am_carbon_d2609.fit[ 4] =  8.20303384e-04;
  am_carbon_d2609.fit[ 5] =  4.22371128e-03;
  am_carbon_d2609.fit[ 6] = -2.16433669e-02;
  am_carbon_d2609.fit[ 7] =  1.24099712e-02;
  am_carbon_d2609.fit[ 8] =  7.77452561e-02;
  am_carbon_d2609.fit[ 9] = -8.30069464e-01;
  am_carbon_d2609.fit[10] = -2.75885881e-01;

  am_carbon_d2500.name = AM_CARBON;
  am_carbon_d2500.bulkdensity = 2.5;
  am_carbon_d2500.scale1000 = 0.732281;
  am_carbon_d2500.slope1000 = -0.502358;
  am_carbon_d2500.fit[ 0] = -3.71360110e-07;
  am_carbon_d2500.fit[ 1] =  4.89821904e-06;
  am_carbon_d2500.fit[ 2] =  4.55471144e-06;
  am_carbon_d2500.fit[ 3] = -2.83416679e-04;
  am_carbon_d2500.fit[ 4] =  6.83619812e-04;
  am_carbon_d2500.fit[ 5] =  4.69166303e-03;
  am_carbon_d2500.fit[ 6] = -2.07887200e-02;
  am_carbon_d2500.fit[ 7] =  5.76073758e-03;
  am_carbon_d2500.fit[ 8] =  8.39535774e-02;
  am_carbon_d2500.fit[ 9] = -6.99093103e-01;
  am_carbon_d2500.fit[10] = -6.49558715e-02;

  am_olivine50_d3000.name = AM_OLIVINE50;
  am_olivine50_d3000.bulkdensity = 3.3;
  am_olivine50_d3000.scale1000 = 0.170046;
  am_olivine50_d3000.slope1000 = -1.003555;
  am_olivine50_d3000.fit[ 0] = -6.54386456e-07;
  am_olivine50_d3000.fit[ 1] =  1.05882600e-05;
  am_olivine50_d3000.fit[ 2] = -2.25920388e-05;
  am_olivine50_d3000.fit[ 3] = -4.26746530e-04;
  am_olivine50_d3000.fit[ 4] =  2.19374014e-03;
  am_olivine50_d3000.fit[ 5] =  2.98361374e-03;
  am_olivine50_d3000.fit[ 6] = -3.79309617e-02;
  am_olivine50_d3000.fit[ 7] =  6.08255788e-02;
  am_olivine50_d3000.fit[ 8] =  5.29361872e-02;
  am_olivine50_d3000.fit[ 9] = -1.24941487e+00;
  am_olivine50_d3000.fit[10] = -1.47638834e+00;

  am_olivine50_d2857.name = AM_OLIVINE50;
  am_olivine50_d2857.bulkdensity = 3.3;
  am_olivine50_d2857.scale1000 = 0.234124;
  am_olivine50_d2857.slope1000 = -0.857340;
  am_olivine50_d2857.fit[ 0] = -6.39907671e-07;
  am_olivine50_d2857.fit[ 1] =  9.90858660e-06;
  am_olivine50_d2857.fit[ 2] = -1.55671815e-05;
  am_olivine50_d2857.fit[ 3] = -4.25865801e-04;
  am_olivine50_d2857.fit[ 4] =  1.89115631e-03;
  am_olivine50_d2857.fit[ 5] =  3.86564976e-03;
  am_olivine50_d2857.fit[ 6] = -3.46816974e-02;
  am_olivine50_d2857.fit[ 7] =  4.51559466e-02;
  am_olivine50_d2857.fit[ 8] =  5.59055807e-02;
  am_olivine50_d2857.fit[ 9] = -1.05008045e+00;
  am_olivine50_d2857.fit[10] = -1.22296183e+00;

  am_olivine50_d2727.name = AM_OLIVINE50;
  am_olivine50_d2727.bulkdensity = 3.3;
  am_olivine50_d2727.scale1000 = 0.306098;
  am_olivine50_d2727.slope1000 = -0.720557;
  am_olivine50_d2727.fit[ 0] = -6.23582425e-07;
  am_olivine50_d2727.fit[ 1] =  9.42018662e-06;
  am_olivine50_d2727.fit[ 2] = -1.11898344e-05;
  am_olivine50_d2727.fit[ 3] = -4.27991885e-04;
  am_olivine50_d2727.fit[ 4] =  1.72410167e-03;
  am_olivine50_d2727.fit[ 5] =  4.61817342e-03;
  am_olivine50_d2727.fit[ 6] = -3.40735115e-02;
  am_olivine50_d2727.fit[ 7] =  3.32849724e-02;
  am_olivine50_d2727.fit[ 8] =  7.91621514e-02;
  am_olivine50_d2727.fit[ 9] = -8.85918268e-01;
  am_olivine50_d2727.fit[10] = -1.06030697e+00;

  am_olivine50_d2609.name = AM_OLIVINE50;
  am_olivine50_d2609.bulkdensity = 3.3;
  am_olivine50_d2609.scale1000 = 0.379768;
  am_olivine50_d2609.slope1000 = -0.595100;
  am_olivine50_d2609.fit[ 0] = -5.93317247e-07;
  am_olivine50_d2609.fit[ 1] =  8.80454396e-06;
  am_olivine50_d2609.fit[ 2] = -8.01056146e-06;
  am_olivine50_d2609.fit[ 3] = -4.17591614e-04;
  am_olivine50_d2609.fit[ 4] =  1.59101641e-03;
  am_olivine50_d2609.fit[ 5] =  4.93813322e-03;
  am_olivine50_d2609.fit[ 6] = -3.39803595e-02;
  am_olivine50_d2609.fit[ 7] =  2.84211347e-02;
  am_olivine50_d2609.fit[ 8] =  1.07359016e-01;
  am_olivine50_d2609.fit[ 9] = -7.78007415e-01;
  am_olivine50_d2609.fit[10] = -9.80069105e-01;

  am_olivine50_d2500.name = AM_OLIVINE50;
  am_olivine50_d2500.bulkdensity = 3.3;
  am_olivine50_d2500.scale1000 = 0.399005;
  am_olivine50_d2500.slope1000 = -0.460211;
  am_olivine50_d2500.fit[ 0] = -5.51723976e-07;
  am_olivine50_d2500.fit[ 1] =  8.09093833e-06;
  am_olivine50_d2500.fit[ 2] = -5.51625696e-06;
  am_olivine50_d2500.fit[ 3] = -3.99746585e-04;
  am_olivine50_d2500.fit[ 4] =  1.46695958e-03;
  am_olivine50_d2500.fit[ 5] =  5.03652584e-03;
  am_olivine50_d2500.fit[ 6] = -3.37239006e-02;
  am_olivine50_d2500.fit[ 7] =  2.73624624e-02;
  am_olivine50_d2500.fit[ 8] =  1.31606954e-01;
  am_olivine50_d2500.fit[ 9] = -7.15506324e-01;
  am_olivine50_d2500.fit[10] = -9.54534540e-01;

  _beta = 0;
  _radius = 1;
  _graindensity = 1.0;
  _D = 3;
  _a0 = 0.1;
  _composition = geometric;
}

/** Return the $\beta$-value. */
double Physical::beta() { return _beta; }
/** Set the $\beta$-value.  Does not recompute radius. */
void Physical::beta(const double b) { _beta = b; }

/** Return the radius ($\mu$m). */
double Physical::radius() { return _radius; }
/** Set the radius ($\mu$m). */
void Physical::radius(const double r) { _radius = r; }

/** Return the composition. */
Composition Physical::composition() { return _composition; }
/** Set the composition. */
void Physical::composition(const int c) {
  switch (c) {
  case GEOMETRIC: _composition = geometric; break;

  case AM_CARBON:
    switch (static_cast<int>(floor(_D*1000))) {
    case 2857: _composition = am_carbon_d2857; break;
    case 2727: _composition = am_carbon_d2727; break;
    case 2609: _composition = am_carbon_d2609; break;
    case 2500: _composition = am_carbon_d2500; break;
    default:
      cerr << "Fractal dimension, D = " << _D << ", not understood, using 3.0\n";
    case 3000: _composition = am_carbon_d3000; break;
    }
    break;

  case AM_OLIVINE50:
    switch (static_cast<int>(floor(_D*1000))) {
    case 2857: _composition = am_olivine50_d2857; break;
    case 2727: _composition = am_olivine50_d2727; break;
    case 2609: _composition = am_olivine50_d2609; break;
    case 2500: _composition = am_olivine50_d2500; break;
    default:
      cerr << "Fractal dimension, D = " << _D << ", not understood, using 3.0\n";
    case 3000: _composition = am_olivine50_d3000; break;
    }
    break;
  }
  _graindensity = _composition.bulkdensity;
}
/** Set the composition. */
void Physical::composition(Composition c) { _composition = c; }

/** Return the grain density (g/cm$^3$). */
double Physical::graindensity() { return _graindensity; }
/** Set the grain density (g/cm$^3$). */
void Physical::graindensity(const double d) { _graindensity = d; }

/** Return the fractal dimension. */
double Physical::fractalDim() { return _D; }
/** Set the fractal dimension.  Update the composition, if
    necessary. */
void Physical::fractalDim(const double d) {
  _D = d;
  composition(_composition.name);
}

/** Return the minimum grain radius ($\mu$m). */
double Physical::minRadius() { return _a0; }
/** Set the minimum grain radius ($\mu$m). */
void Physical::minRadius(const double a0) { _a0 = a0; }

/** Return the grain temperature calculator. */
GrainTemperature Physical::grainT() { return _grainT; }
/** Load a composition into the grain temperature calculator. */
void Physical::grainT(const materials m) { _grainT.load(m); }

/** Update the $\beta$-value, using the current grain radius and grain
    density (see updateBeta(rho)). */
void Physical::updateBeta() { updateBeta(_graindensity); }
/** Update the $\beta$-value, using the current grain radius and
    provided grain density (ignored if the composition is not
    GEOMETRIC).  If radius <= 0, beta will be 0. */
void Physical::updateBeta(const double rho) {
  if (_radius <= 0) {
    _beta = 0;
  } else {
    if (_composition.name == GEOMETRIC)
      _beta = 0.57 / _radius / rho;
    else
      _beta = fit2beta(_radius);
  }
}

/** Update the radius, using the current grain $\beta$-value and
    grain density (see updateRadius(rho)). */
void Physical::updateRadius() { updateRadius(_graindensity); }
/** Update the radius, using the current grain $\beta$-value and
    provided grain density.  Composition and porosity are ignored,
    Q_pr = 1.  If $\beta$ <= 0, radius is set to 1 m.  This function
    is intended for syndynes and pre-v0.7.3 xyz files which did not
    include radius. */
void Physical::updateRadius(const double rho) {
  if (_beta <= 0) {
    _radius = 1e6;
  } else {
    // use the grain density and beta, Q_pr = 1
    _radius = 0.57 / rho / _beta;
  }
}

/** Update the grain density, using the current grain radius,
    composition, and fractal dimension. */
void Physical::updateGrainDensity() {
  if (_D == 3)
    _graindensity = _composition.bulkdensity;
  else
    _graindensity = _composition.bulkdensity * exp((_D - 3.0) * log(_radius / _a0));
}

/** Take a log-space polynomial fit to a $\beta$-radius relationship
    and return the $\beta$-value for the given radius.  The fits are
    defined in Physical().
*/
double Physical::fit2beta(const double a) {
  if (a > 1000.0) {
    // For large grains
    return _composition.scale1000 * exp(_composition.slope1000 * log(a));
  } else {
    double loga = log(a);
    double p = 0;
    for (int i=0; i<=10; i++)
      p = p * loga + _composition.fit[i];
    return exp(p);
  }
}

