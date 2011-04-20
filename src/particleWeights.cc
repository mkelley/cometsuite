/***************************************************************************

  Weights grains according to various schemes (e.g., thermal
  emission).

  Copyright (C) 2005-2010 by Michael S. Kelley <msk@astro.umd.edu>

 ***************************************************************************/

#if HAVE_CONFIG_H
#  include <config.h>
#endif /* HAVE_CONFIG_H */

#include <string>
#include <sstream>
#include <algorithm>
#include <cstdlib>
#include <cmath>
#include <cfloat>
#include "GrainTemperature.h"
#include "particle.h"
#include "particleWeights.h"
#include "rundynamics.h"
#include "state.h"
#include "StringConv.h"

using namespace std;

/** Default parameters. */
particleWeights::particleWeights() {
  _afrhoSlope = -2;
  _wavelength = 24;
  _scatteringMode = false;
  _thermalMode = true;
  _nuclearPsd = "none";
  _fileScales.push_back(1.0);
  _scaleIndex = -1;
}

/** Return the dust production's logarithmic slope. */
double particleWeights::afrhoSlope() { return _afrhoSlope; }
/** Set the dust production's logarithmic slope. */
void particleWeights::afrhoSlope(const double s) { _afrhoSlope = s; }

/** Return the wavelength of simulated thermal emission and light
    scattering. */
double particleWeights::wavelength() { return _wavelength; }
/** Set the wavelength for simulated thermal emission and light
    scattering. */
void particleWeights::wavelength(const double w) { _wavelength = w; }

/** Return true if we are in scattering mode. */
bool particleWeights::scatteringMode() { return _scatteringMode; }
/** Enable or disable scattering mode. */
void particleWeights::scatteringMode(const bool m) { _scatteringMode = m; }

/** Return true if we are in thermal emission mode. */
bool particleWeights::thermalMode() { return _thermalMode; }
/** Enable or disable thermal emission mode. */
void particleWeights::thermalMode(const bool m) { _thermalMode = m; }

/** Return the particle size distribution function name. */
string particleWeights::nuclearPsd() { return _nuclearPsd; }
/** Set the particle size distribution function. */
void particleWeights::nuclearPsd(const string p) {
  _nuclearPsd = p;
  transform(_nuclearPsd.begin(), _nuclearPsd.end(), _nuclearPsd.begin(),
	    (int(*)(int))tolower);
}

/** Returns the weight for a grain determined by the comet's dust
    production at t_i.  rh is in km. */
double particleWeights::afrhoWeight(const double rh) {
  return pow(rh / _AU, _afrhoSlope);
}

/** Returns the weight for a grain determined by the grains's thermal
    emission.

    Qem = 2 pi a / lambda  for a <= lambda / 2 / pi
        = 1                for a >  lambda / 2 / pi

    where a is the grain radius.  rh and Delta are in km, T is in K,
    and the returned weight is in units of Jy.
*/
double particleWeights::thermalWeight(const double radius,
				      const double rh,
				      const double Delta,
				      const double T) {
  double wt = 1;

  if ((_wavelength > 0) && _thermalMode) {
    double a = radius;
    double Qem, Delta_m;

    if (a < _wavelength / 2 / M_PI)
      Qem = a * 2 * M_PI / _wavelength;
    else
      Qem = 1.0;

    a = a * 1e-6;  // to meters
    Delta_m = Delta * 1e3;  // to meters
    wt = M_PI * a * a * Qem * 
      planck(_wavelength, T) / 
      Delta_m / Delta_m * 1e26; // Jy
  }

  return wt;
}

/** Returns the particle weight for light scattering using a rough
    approximation.

    Qsca \propto a^4 / lambda^4  for a <  lambda / 2 / pi
         \propto a^2             for a <= lambda / 2 / pi

    \todo Make this absolute flux (Jy), rather than relative flux.
 */
double particleWeights::scatteringWeight(const double radius) {
  double wt = 1;
  if ((_wavelength > 0)  && _scatteringMode) {
    double a = radius;

    if (a < _wavelength / 2 / M_PI) {
      wt = a * 2 * M_PI / _wavelength;
      wt = wt * wt * wt * wt * a * a;
    } else {
      wt = a * a;
    }
  }

  return wt;
}

/** Return the weight for a particle given the current particle size
    distribution assuming the input psd is dn/dlog(a) ~ 1.  Possible
    distributions are: none, ISM (dn/da = a^-3.5), a^x, or Hanner.
    dn/da(a = 1 micron) = 1.0 for the power-law distributions.

    \todo Don't assume the input PSD is dn/dlog(a) for syndynes.
*/
double particleWeights::nuclearPsdWeight(const double radius) {
  if (_nuclearPsd.find("none") != string::npos) return 1;

  double a = radius;

  /* dn/da = a^-3.5
     -- For radius --
     the simulation picks dn/dlog(a) = constant
     weight by N_ism/N_sim:
       dn/da|ism = a^-3.5  dn/da|sim = a^-1
       dn = da * dn/da
       -> dn_ism/dn_sim = da * a^-3.5 / (da * a^-1)

       for small intervals dn_ism/dn_sim = n_ism/n_sim = a^-2.5

     -- For beta (old method) --
     want dn/da = a^-3.5 or dn/dbeta = beta^1.5
     weight by N_ism/N_sim:
       dn/dbeta|ism = beta^1.5  dn/dbeta|sim = beta^-1
       dn = dbeta * dn/dbeta
       -> dn_ism/dn_sim = dbeta * beta^1.5 / (dbeta * beta^-1)
       for small intervals dn_ism/dn_sim = n_ism/n_sim = beta^2.5

  */
  if (_nuclearPsd.find("ism") != string::npos) return pow(a, -2.5);

  /* dn/da = a^x
     -- For radius --
     the simulation picks dn/dlog(a) = constant
     weight by N_ism/N_sim:
       dn/da|ism = a^x  dn/da|sim = a^-1
       dn = da * dn/da
       -> dn_ism/dn_sim = da * a^x / (da * a^-1)

       for small intervals dn_ism/dn_sim = n_ism/n_sim = a^(x+1)

     -- For beta (old method) --
     the simulation picks dn/dlog(a) = constant
     want dn/dbeta = beta^x
     weight by N_new/N_sim:
       dn/dbeta|new = beta^x  dn/dbeta|sim = beta^-1
       dn = dbeta * dn/dbeta
       -> dn_new/dn_sim = dbeta * beta^x / (dbeta * beta^-1)
       for small intervals dn_new/dn_sim = n_new/n_sim = beta^(x+1) */
  int c;
  if ((c = _nuclearPsd.find("a^")) != string::npos) {
    double x = StringConv(_nuclearPsd.substr(c + 2,
      _nuclearPsd.size() - 2 - c)).toScalar<double>();
    return pow(a, x + 1);
  }

  /* Hanner grain size distribution (Hanner 1983)
     dn/da = (1 - a0 / a)^M * (a0 / a)^N
       N = large particle slope
       M determines the peak grain size:
       ap = a0 = (M + N) / N
       a0 = minimum grain radius

    -> n_hanner/n_sim = (1 - a0 / a)^M * (a0/a)^N * a

     Usage: hanner a0 M N
  */
  if ((c = _nuclearPsd.find("hanner")) != string::npos) {
    vector<double> hanner(StringConv(_nuclearPsd.substr(c + 6,
      _nuclearPsd.size() - 6 - c)).toVector<double>());
    return pow(1 - hanner[0] / a, hanner[1]) * 
      pow(hanner[0] / a, hanner[2]) * a;
  }

  return 1;
}

/** Set the user defined scales for all xyz files. */
void particleWeights::fileScales(const vector<double> scales) {
  _fileScales = scales;
}

/** Return a vector of all file scales. */
vector<double> particleWeights::fileScales() { return _fileScales; }

/** Return the user defined scale for the current xyz file. */
double particleWeights::fileScale() { return _fileScales[_scaleIndex]; }

/** Set the current user defined file scale to the next available file
    scale. */
void particleWeights::nextFileScale() {
  _scaleIndex++;
  if (_scaleIndex >= _fileScales.size()) {
    _scaleIndex = _fileScales.size() - 1;
  }
}

/** Compute the Planck function: wave is in microns, temp is in
    Kelvin [units: W/m2/Hz/sr].  */
double particleWeights::planck(const double wave, const double temp) {
  const double c1 = 3.972894922e-25;  // 2 h c
  const double c2 = 1.438768660e-2;   // h c / k
  double wave_m, a;

  wave_m = wave * 1e-6;  // wavelength in meters
  a = c2 / wave_m / temp;

  if (a < DBL_MAX_EXP) {
    return c1 / (wave_m * wave_m * wave_m * (exp(a) - 1.0));
  } else {
    return 0;
  }
}
