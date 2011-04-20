/***************************************************************************

  Random variable generation.

  Copyright (C) 2008-2010 by Michael S. Kelley <msk@astro.umd.edu>

  \todo Need a random number generator that can pick numbers over a
  large dynamic range.  rand() returns an integer and the quantization
  is apparent at very small ages.

 ***************************************************************************/

#if HAVE_CONFIG_H
#  include <config.h>
#endif /* HAVE_CONFIG_H */

#include <iostream>
#include <fstream>
#include <cstdlib>
#include <ctime>
#include <cmath>
#include <queue>
#include "Distribution.h"
#include "rundynamics.h"

using namespace std;

const double Distribution::logZero = -999;
const double Distribution::emptySequence = -1e99;

Distribution::Distribution() {
  // prepare for random number generation
  // Grab a seed from /dev/urandom
  ifstream devrand("/dev/urandom", ios::in);
  unsigned int seed;
  if (devrand.fail()) {
    // instead, seed with the time
    seed = static_cast<unsigned int>(time(NULL));
  } else {
    devrand.read(reinterpret_cast<char *>(&seed), sizeof(int));
    devrand.close();
  }
  srand(seed);  // seed the random number generator

  // default min, max, mu, sigma, grid steps
  _min = 0;
  _max = 1;
  _mu = 1.0;
  _sigma = 1.0;
  _nGridSteps = 10;

  // the default sequence
  for (int i = 0; i<_nGridSteps; i++)
    _seq.push((_min - _max) / static_cast<double>(i) + _min);

  // default distribution for new()
  distribution(LINEAR);
}

/** Sets the minimum random variate value. */
void Distribution::min(const double m) {
  _min = m;
  if (m == 0) {
    _log10Min = logZero;
  } else {
    _log10Min = log10(m);
  }
}
/** Returns the minimum random variate value. */
double Distribution::min() { return _min; }

/** Sets the maximum random variate value. */
void Distribution::max(const double m) {
  _max = m;
  if (m == 0) {
    _log10Max = logZero;
  } else {
    _log10Max = log10(m);
  }
}
/** Returns the maximum random variate value. */
double Distribution::max() { return _max; }

/** Sets the minimum random variate value from log10(min). */
void Distribution::log10Min(const double m) {
  _log10Min = m;
  if (m == logZero) {
    _min = 0;
  } else {
    _min = pow(10, m);
  }
}
/** Returns log10(minimum random variate value). */
double Distribution::log10Min() { return _log10Min; }

/** Sets the maximum random variate value from log10(max). */
void Distribution::log10Max(const double m) {
  _log10Max = m;
  if (m == logZero) {
    _max = 0;
  } else {
    _max = pow(10, m);
  }
}
/** Returns log10(maximum random variate value). */
double Distribution::log10Max() { return _log10Max; }

/** Sets the mean of the normal distribution. */
void Distribution::mu(const double m) {
  _mu = m;
}
/** Returns the mean of the normal distribution. */
double Distribution::mu() { return _mu; }

/** Sets the width (sigma) of the normal distribution. */
void Distribution::sigma(const double s) {
  _sigma = s;
}
/** Returns the width (sigma) of the normal distribution. */
double Distribution::sigma() { return _sigma; }

/** Sets the variate grid steps. */
void Distribution::nGridSteps(const unsigned int gs) { _nGridSteps = gs; }
/** Returns the number of grid steps. */
unsigned int Distribution::nGridSteps() { return _nGridSteps; }

/** Sets the default sequence. */
void Distribution::setSequence(queue<double> q) { _seq = q; }
/** Returns the default sequence. */
queue<double> Distribution::getSequence() { return _seq; }

/** Returns a random value from dn/dlog(value) ~ 1 between
    log(valueMin) and log(valueMax) using the default limits. */
double Distribution::dn_dlogx__1() { return dn_dlogx__1(_log10Min, _log10Max); }

/** Returns a random value from dn/dlog(value) ~ 1 between minVal and
    maxVal.  Set log10MinVal or log10MaxVal to Distribution::logZero to
    return 0.
*/
double Distribution::dn_dlogx__1(const double log10MinVal,
				 const double log10MaxVal) {
  // if the zero flag is set, return 0
  if (log10MaxVal == logZero || log10MinVal == logZero) return 0;

  // if they are equal we are done
  if (log10MinVal == log10MaxVal) return exp(log10MinVal * log(10.0));

  // pick a value: dn/dlog(value) ~ 1; this is done by picking
  // value from a power law distribution, dn/dvalue ~ value^-1
  double u = static_cast<double>(rand()) / static_cast<double>(RAND_MAX);
  return exp((u * (log10MaxVal - log10MinVal) + log10MinVal) * log(10.0));
}

/** Returns a random value between minVal and maxVal, linearly
    distributed, using the default limits. */
double Distribution::dn_dx__1() { return dn_dx__1(_min, _max); }

/** Returns a random value between minVal and maxVal, linearly
    distributed. */
double Distribution::dn_dx__1(const double minVal, const double maxVal) {
  return (maxVal - minVal) * static_cast<double>(rand()) / 
    static_cast<double>(RAND_MAX) + minVal;
}

/** Returns a normally distributed random value centered on mu, with
    width sigma (as a fraction of of mu), and greater than minVal,
    using the default limits. */
double Distribution::dn_dx__normal() {
  return dn_dx__normal(_mu, _sigma, _min);
}

/** Returns a normally distributed random value centered on mu, with
    width sigma (as a fraction of mu), and greater than minVal.  Uses
    the formula of Box & Mueller (1958). */
double Distribution::dn_dx__normal(const double mu, const double sigma,
			      const double minVal) {
  double z = minVal - 1;
  double r1, r2;
  int i = 0;
  while (z < minVal) {
    r1 = static_cast<double>(rand()) / static_cast<double>(RAND_MAX);
    r2 = static_cast<double>(rand()) / static_cast<double>(RAND_MAX);
    z = sqrt(-2 * log(r1)) * cos(6.2831853071795862 * r2) * sigma * mu + mu;
    i++;
    if (i == 100) {
      cerr << "Distribution::dn_dx_normal()\n - Too many iterations while attempting to generate a random variate\n - (mu, sigma, minVal) = (" <<
	mu << ", " << sigma << ", " << minVal << ")\n";
      i = 0;
    }
  }
  return z;
}

/** Returns a value picked from a regular grid between minVal and
    maxVal using the default limits and grid spacing. */
double Distribution::grid() { return grid(_nGridSteps, _min, _max); }

/** Returns a value picked from a regular grid between minVal and
    maxVal.  The grid is nGrids steps across. */
double Distribution::grid(const double nGrids, const double
			     minVal, const double maxVal) {
  double u = static_cast<double>(rand()) / static_cast<double>(RAND_MAX);
  u = int(u * nGrids) / (nGrids - 1);
  return (maxVal - minVal) * u + minVal;
}

/** Returns a value picked from a regular grid spaced in log space
    between minVal and maxVal using the default limits and grid
    spacing. */
double Distribution::log_grid() { return log_grid(_nGridSteps, _log10Min, _log10Max); }

/** Returns a value picked from a regular grid spaced in log space
    between minVal and maxVal.  The grid is nGrids steps across.
    dn/dlog(val) ~ val^-1. */
double Distribution::log_grid(const double nGrids, const double
			      log10MinVal, const double log10MaxVal) {
  // if the zero flag is set, return 0
  if (log10MaxVal == logZero || log10MinVal == logZero) return 0;

  if (log10MinVal == log10MaxVal) return exp(log10MinVal * log(10.0));

  double u = static_cast<double>(rand()) / static_cast<double>(RAND_MAX);
  u = int(u * nGrids) / (nGrids - 1);
  return exp((u * (log10MaxVal - log10MinVal) + log10MinVal) * log(10.0));
}

/** Returns the next value picked from the default sequence of values.
    If the sequence is exhausted, return emptySequence. */
double Distribution::sequence() {
  return sequence(_seq);
}

/** Returns the next value picked from seq.  If the sequence is
    exhausted, return emptySequence. */
double Distribution::sequence(queue<double>& seq) {
  if (DEBUG) cerr << "Distribution::sequence()\n";
  if (seq.empty())
    return emptySequence;

  double val = seq.front();
  if (DEBUG) cerr << val << " ";
  seq.pop();
  if (DEBUG) cerr << seq.front() << endl;
  return val;
}

/** Returns theta() using default values. */
double Distribution::theta() { return theta(_min, _max); }

/** Returns a polar angle picked from a distribution that is uniformly
    distributed in solid angle on a sphere.  minVal and maxVal are in
    radians measured from the input axis. */
double Distribution::theta(const double minVal, const double maxVal) {
  // Pick a random polar angle (theta) that is part of a uniform
  // distribution in solid angle on a sphere: th =
  // acos((1-u*)cos(th_min)+u*cos(th_max)) where u is a random number
  // uniformly-distributed between 0 and 1.
  double u = static_cast<double>(rand()) / static_cast<double>(RAND_MAX);
  //  return acos(1 - u + u * cos(maxVal));  // for min = 0, max = maxVal
  return acos((1 - u) * cos(minVal) + u * cos(maxVal));
}

/** Returns theta_cos() using default values. */
double Distribution::theta_cos() { return theta_cos(_min, _max); }

/** Returns a polar angle picked from a distribution where the flux
    through a spherical solid angle is proportional to
    cos(polar_angle).  minVal and maxVal are in radians.  Typically,
    polar_angle <= 90. */
double Distribution::theta_cos(const double minVal, const double maxVal) {
  // th = acos(sqrt(1-u)) where u is uniformly-distributed between 0
  // and 1.
  double u = static_cast<double>(rand()) / static_cast<double>(RAND_MAX);
  //  return acos(sqrt(1 - u));  // min = 0, max = pi/2
  return acos(sqrt((1 - u) * cos(minVal) * cos(minVal) 
		   + u * cos(maxVal) * cos(maxVal)));
}

/** Returns theta_temp() using default values. */
double Distribution::theta_temp() { return theta_temp(_min, _max); }

/** Returns a polar angle picked from a distribution where the flux
    through a spherical solid angle is proportional to
    cos(polar_angle)^0.25 (i.e., the temperature distribution on a
    spherical nucleus).  minVal and maxVal are in radians measured.
    Typically, polar_angle <= 90. */
double Distribution::theta_temp(const double minVal, const double maxVal) {
  // th = (acos((1-u)^(4/5)) where u is uniformly-distributed between
  // 0 and 1.
  double u = static_cast<double>(rand()) / static_cast<double>(RAND_MAX);
  //  return acos(pow(1 - u, 0.8));  // min = 0, max = pi/2
  return acos(pow((1 - u) * pow(cos(minVal), 1.25) + 
		  u * pow(cos(maxVal), 1.25), 0.8));
}

/** Returns phi() using default values. */
double Distribution::phi() { return phi(_min, _max); }

/** Returns a random angle picked from uniform distribution.  minVal
    and maxVal are in radians measured from an arbitrary axis.  The
    angle may be the azimuthal angle in a distribution uniform in
    solid angle, or the theta_cos() or theta_temp()
    distributions. */
double Distribution::phi(const double minVal, const double maxVal) {
  double v = static_cast<double>(rand()) / static_cast<double>(RAND_MAX);
  return v * (maxVal - minVal) + minVal;
}

/** Set the current distribution to one of distMethod (see
    Distribution.h). */
void Distribution::distribution(const unsigned int d) {
  _distribution = d;
}

/** Return the current distribution method. */
unsigned int Distribution::distribution() {
  return _distribution;
}

/** Generate a new random variate from the current distribution. */
double Distribution::next() {
  switch (_distribution) {
  case LOG: { return dn_dlogx__1(); }
  case LINEAR: { return dn_dx__1(); }
  case NORMAL: { return dn_dx__normal(); }
  case GRID: { return grid(); }
  case LOGGRID: { return log_grid(); }
  case SEQUENCE: { return sequence(); }
  case ISOSPHERE: { return theta(); }
  case COSSPHERE: { return theta_cos(); }
  case TEMPSPHERE: { return theta_temp(); }
  }
  return -1;
}
