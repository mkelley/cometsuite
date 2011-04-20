/***************************************************************************
  Copyright (C) 2005-2010 by Michael S. Kelley <msk@astro.umd.edu>

 ***************************************************************************/

#if !defined(__PARTICLEWEIGHTS)
#define __PARTICLEWEIGHTS 1

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <string>
#include "particle.h"
#include "particleWeights.h"
#include "GrainTemperature.h"
#include "rundynamics.h"

using namespace std;

/** Handles the weighting of cometary particles.  For example, use
    these routines to determine how much flux a particle contributes
    to an image, or how a particle contributes to a particle size
    distribution given the dust production when the particle was
    released. */
class particleWeights {
 public:
  particleWeights();

  // parameter setting
  double afrhoSlope();
  void afrhoSlope(const double);
  double wavelength();
  void wavelength(const double);
  bool scatteringMode();
  void scatteringMode(const bool);
  bool thermalMode();
  void thermalMode(const bool);
  string nuclearPsd();
  void nuclearPsd(const string);

  // methods for the weights
  double afrhoWeight(const double);
  double thermalWeight(const double, const double, const double, const double);
  double scatteringWeight(const double);
  double nuclearPsdWeight(const double);
  GrainTemperature grainT();

  // utility functions
  double planck(const double, const double);

  // file scaling
  void fileScales(const vector<double>);
  vector<double> fileScales();
  double fileScale();
  void nextFileScale();

 private:
  // particle weighting parameters
  bool _scatteringMode, _thermalMode;
  int _scaleIndex;
  double _afrhoSlope, _wavelength;
  vector<double> _fileScales;
  string _nuclearPsd;
  GrainTemperature _grainT;
};

#endif
