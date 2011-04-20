/***************************************************************************
  Copyright (C) 2008-2010 by Michael S. Kelley <msk@astro.umd.edu>

  ***************************************************************************/

#if !defined(__DISTRIBUTION)
#define __DISTRIBUTION 1

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <queue>
#include "rundynamics.h"

using namespace std;

/** Returns random variates from defined distributions. */
class Distribution {
public:
  Distribution();
  // Distribution(const double, const double);
  // Distribution(const double, const double, const unsigned int);

  /** Flag to return 0 with logarithmic distributions since log(0) does
      not exist. */
  static const double logZero;

  /** Flag that signifies our sequence is empty. */
  static const double emptySequence;

  // methods to set or return default variate limits or grid parameters
  void min(const double);
  double min();
  void max(const double);
  double max();
  void log10Min(const double);
  double log10Min();
  void log10Max(const double);
  double log10Max();
  void mu(const double);
  double mu();
  void sigma(const double);
  double sigma();
  void nGridSteps(const unsigned int);
  unsigned int nGridSteps();
  void setSequence(queue<double>);
  queue<double> getSequence();

  // methods to return random variates
  double dn_dlogx__1();
  double dn_dlogx__1(const double, const double);
  double dn_dx__1();
  double dn_dx__1(const double, const double);
  double dn_dx__normal();
  double dn_dx__normal(const double, const double, const double);
  double grid();
  double grid(const double, const double, const double);
  double log_grid();
  double log_grid(const double, const double, const double);
  double theta();
  double theta(const double, const double);
  double theta_cos();
  double theta_cos(const double, const double);
  double theta_temp();
  double theta_temp(const double, const double);
  double phi();
  double phi(const double, const double);

  // methods to return one from a list of values
  double sequence();
  double sequence(queue<double>&);

  // alternatively, we can set a default distribution and just call
  // next()
  enum distMethod { LOG, LINEAR, NORMAL, GRID, LOGGRID, SEQUENCE,
		    ISOSPHERE, COSSPHERE, TEMPSPHERE };
  void distribution(const unsigned int);
  unsigned int distribution();
  double next();

private:
  queue<double> _seq;
  double _min, _max, _log10Min, _log10Max, _mu, _sigma;
  unsigned int _nGridSteps;
  unsigned int _distribution;
};

#endif
