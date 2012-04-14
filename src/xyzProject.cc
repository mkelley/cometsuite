/***************************************************************************

  Projects grains onto the Celestial Sphere.

  Copyright (C) 2005-2010,2012 by Michael S. Kelley <msk@astro.umd.edu>

 ***************************************************************************/

#if HAVE_CONFIG_H
#  include <config.h>
#endif /* HAVE_CONFIG_H */

#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <vector>
#include <cmath>
#include "xyzstream.h"
#include "Vej.h"
#include "Vector.h"
#include "projection.h"
#include "xyzProject.h"
#include "rundynamics.h"
#include "state.h"
#include "particle.h"
#include "longlat.h"
#include "paramSet.h"
#include "pfunctions.h"

using namespace std;

/** Default parameters. */
xyzProject::xyzProject() {
  _verbose = false;
  _currentFile = 0;
  _graindensity = 0;
  _ecliptic = false;
  _observerName = "Earth";
  _max = -1;
  _ageInvert = false;
  _betaInvert = false;
  _latInvert = false;
  _lonInvert = false;
  _radInvert = false;
  _sunInvert = false;
  _Coffset.lambda = 0;
  _Coffset.beta = 0;
  _npole[0] = 0.0;
  _npole[1] = 90.0;
  _rotRate = 0.0;
  _rotPhase = 0.0;

  _jet = false;
  _jetAngle = 0.0;
  _jetLoc[0] = 0.0;
  _jetLoc[1] = 0.0;

  _rhlimit = -1;
  _vejGen.v0(-1);
  _vejGen.simpleActivity();
}

/** Closes the input file. */
xyzProject::~xyzProject() { xyzfile.close(); }

/** Loads the next input file and prepares for reading particles.
    Returns true on success. */
bool xyzProject::loadXyzfile() {
  if (_currentFile >= _xyzfileNames.size()) return false;
  // reset the error indicator
  _p.error = false;

  // close the file stream if already open
  if (xyzfile.is_open()) {
    xyzfile.close();
    xyzfile.clear();
  }

  cout << "Opening " << _xyzfileNames[_currentFile] << " for input." << endl;

  // open the next file
  xyzfile.xyzopen(_xyzfileNames[_currentFile], xyzstream::READ);

  if (xyzfile.fail()) {
    cerr << "Error opening file: " << _xyzfileNames[_currentFile] << "\n";
    return false;
  }

  // get the parameters
  stringstream str;
  string header;
  header = xyzfile.readHeader();
  if (verbose()) cerr << header << endl;
  str << header;
  _parameters.loadParameters(str);
  xyzfile.initData(parameters());

  // warn if we are dealing with a Syndynes file
  if (parameters().isSyndynes()) {
    cerr << "Warning: Input file created with Syndynes.\n";
    long nParticles = static_cast<long>((parameters().beta().size() + 
					 ((parameters().orbit()>0)?1:0)) *
					parameters().steps());
    _parameters.nParticles(nParticles);
  }

  _currentFile++;
  return true;
}

/** Return the verbosity. */
bool xyzProject::verbose() { return _verbose; }
/** Set the verbosity. */
void xyzProject::verbose(const bool v) { _verbose = v; }

/** Return the input file names. */
vector<string> xyzProject::xyzfileNames() { return _xyzfileNames; }
/** Add a file to the input file names. */
void xyzProject::xyzfileNames(const string fn) { _xyzfileNames.push_back(fn); }

/** Return the current input file name. */
string xyzProject::currentFileName() { return _xyzfileNames[_currentFile-1]; }

/** Return the current input file number. */
int xyzProject::currentFileNumber() { return _currentFile; }

/** Return the default grain density. */
double xyzProject::graindensity() { return _graindensity; }
/** Set the default grain density (g/cm^3). */
void xyzProject::graindensity(const double d) { _graindensity = d; }

/** Return the ecliptic flag. */
bool xyzProject::ecliptic() { return _ecliptic; }
/** Set the ecliptic flag. */
void xyzProject::ecliptic(const bool ec) { _ecliptic = ec; }

/** Return the observer's name. */
string xyzProject::observerName() { return _observerName; }
/** Set the observer's name. */
void xyzProject::observerName(const string obs) { _observerName = obs; }

/** Return the rotation period in hours. */
double xyzProject::rotPeriod() {
  if (_rotRate == 0) {
    return 0;
  } else {
    return 0.1 / _rotRate;
  }
}
/** Set the rotation period in hours. */
void xyzProject::rotPeriod(const double P) {
  if (P == 0) {
    _rotRate = 0;
  } else {
    _rotRate = 0.1 / P;
  }
}

/** Return the rotation phase in degrees. */
double xyzProject::rotPhase() { return _rotPhase; }
/** Set the rotation phase in degrees. */
void xyzProject::rotPhase(const double phi) { _rotPhase = phi; }

/** Return the absolute positional offsets in units of degrees. */
longlat xyzProject::offset() { return _Coffset; }
/** Setup absolute positional offsets in units of degrees. */
void xyzProject::offset(const valarray<float> off) {
  _Coffset.lambda = off[0];
  _Coffset.beta = off[1];
}
/** Setup absolute positional offsets in units of degrees. */
void xyzProject::offset(const longlat off) {
  _Coffset = off;
}

/** Return the north pole. */
valarray<float> xyzProject::npole() {
  valarray<float> np(2);
  np[0] = _npole[0];
  np[1] = _npole[1];
  return np;
}
/** Set the north pole. */
void xyzProject::npole(const float* np) {
  _npole[0] = np[0];
  _npole[1] = np[1];
}
/** Set the north pole. */
void xyzProject::npole(const valarray<float> np) {
  _npole[0] = np[0];
  _npole[1] = np[1];
}

/** Return the comet's RA and Dec. */
longlat xyzProject::cometRaDec() { return _cometRaDec; }
/** Set the comet's RA and Dec. */
void xyzProject::cometRaDec(const longlat crd) { _cometRaDec = crd; }

/** Return the origin's RA and Dec. */
longlat xyzProject::originRaDec() { return _originRaDec; }
/** Set the origin's RA and Dec. */
void xyzProject::originRaDec(const longlat ord) { _originRaDec = ord; }

/** Return the origin's offset in arcseconds. */
longlat xyzProject::originOffset() { return _originOffset; }
/** Set the origin's offset in arcseconds. */
void xyzProject::originOffset(const longlat oxy) { _originOffset = oxy; }

/** Return the maximum number of particles to process. */
long xyzProject::max() { return _max; }
/** Set the maximum number of particles to process. */
void xyzProject::max(const long m) { _max = m; }

/** Return the age lower and upper limits (units: seconds). */
valarray<float> xyzProject::ageRange() { return _ageRange; }
/** Set the age lower and upper limits. */
void xyzProject::ageRange(const valarray<float> ar) {
  if (_ageRange.size() != 2)
    _ageRange.resize(2);
  _ageRange = ar;
}

/** Return the beta lower and upper limits. */
valarray<float> xyzProject::betaRange() { return _betaRange; }
/** Set the beta lower and upper limits. */
void xyzProject::betaRange(const valarray<float> br) {
  if (_betaRange.size() != 2)
    _betaRange.resize(2);
  _betaRange = br;
}

/** Return the latitude lower and upper limits. */
valarray<float> xyzProject::latRange() { return _latRange; }
/** Set the latitude lower and upper limits. */
void xyzProject::latRange(const valarray<float> lr) {
  if (_latRange.size() != 2)
    _latRange.resize(2);
  _latRange = lr;
}

/** Return the longitude lower and upper limits. */
valarray<float> xyzProject::lonRange() { return _lonRange; }
/** Set the longitude lower and upper limits. */
void xyzProject::lonRange(const valarray<float> lr) {
  if (_lonRange.size() != 2)
    _lonRange.resize(2);
  _lonRange = lr;

  // Make sure they are between 0 and 360 deg
  for (int i=0; i<2; i++) {
    while (_lonRange[i] < 0)
      _lonRange[i] += 360;
    _lonRange[i] = fmod(_lonRange[i], (float)360.0);
  }
}

/** Return the radius lower and upper limits. */
valarray<float> xyzProject::radRange() { return _radRange; }
/** Set the radius lower and upper limits. */
void xyzProject::radRange(const valarray<float> rr) {
  if (_radRange.size() != 2)
    _radRange.resize(2);
  _radRange = rr;
}

/** Return the Sun-zenith angle lower and upper limits. */
valarray<float> xyzProject::sunRange() { return _sunRange; }
/** Set the Sun-zenith angle lower and upper limits. */
void xyzProject::sunRange(const valarray<float> sr) {
  if (_sunRange.size() != 2)
    _sunRange.resize(2);
  _sunRange = sr;
}

/** True if the jet is enabled. */
bool xyzProject::jetOn() { return _jet; }
/** Set the jet status. */
void xyzProject::setJet(const bool j) { _jet = j; }

/** Return the jet's location (ecliptic lambda, beta). */
longlat xyzProject::jet() { return _jetLB; }
/** Set the jet's location (ecliptic lambda, beta). */
void xyzProject::jet(const valarray<float> jet) {
  _jetLB.lambda = jet[0];
  _jetLB.beta = jet[1];
  _jetV = longlatToVector(_jetLB);
}
/** Set the jet's location (ecliptic lambda, beta). */
void xyzProject::jet(longlat jet) {
  _jetLB.lambda = jet.lambda;
  _jetLB.beta = jet.beta;
  _jetV = longlatToVector(_jetLB);
}

/** Return the jet location (ecliptic rectangular coods). */
Vector xyzProject::jetV() { return _jet; }
/** Set the jet status. */
void xyzProject::jetV(const Vector j) {
  _jetV = j;
  _jetLB = getEcliptic(Vector(0, 0, 0), j);
}

/** Return the jet opening angle, full width (degrees). */
float xyzProject::jetAngle() { return _jetAngle; }
/** Set the jet status. */
void xyzProject::jetAngle(const float a) { _jetAngle = a; }

/** Return the simple activity ejection velocity parameter. */
double xyzProject::vLimit() { return _vejGen.v0(); }
/** Set the simple activity ejection velocity parameter.  If set to
    less than zero, the ejection velocities are unlimited.  If set to
    >=0, then the ejection velocities are limited to the range
    [0, vLimit()*sqrt(beta/r_h)]. */
void xyzProject::vLimit(const double v) { _vejGen.v0(v); }

/** Return the maximum heliocentric distance for dust production in
    AU. */
double xyzProject::rhLimit() { return _rhlimit / _AU; }
/** Set the maximum heliocentric distance for dust production in AU.
    Set to -1 for no limit. */
void xyzProject::rhLimit(const double rh) { _rhlimit = rh * _AU; }

/** Return true if the age limit logic is inverted. */
bool xyzProject::ageInvert() { return _ageInvert; }
/** Set the status of the age limit logic. */
void xyzProject::ageInvert(const bool ai) { _ageInvert = ai; }

/** Return true if the beta limit logic is inverted. */
bool xyzProject::betaInvert() { return _betaInvert; }
/** Set the status of the beta limit logic. */
void xyzProject::betaInvert(const bool bi) { _betaInvert = bi; }

/** Return true if the latitude limit logic is inverted. */
bool xyzProject::latInvert() { return _latInvert; }
/** Set the status of the latitude limit logic. */
void xyzProject::latInvert(const bool li) { _latInvert = li; }

/** Return true if the longitude limit logic is inverted. */
bool xyzProject::lonInvert() { return _lonInvert; }
/** Set the status of the longitude limit logic. */
void xyzProject::lonInvert(const bool li) { _lonInvert = li; }

/** Return true if the radius limit logic is inverted. */
bool xyzProject::radInvert() { return _radInvert; }
/** Set the status of the radius limit logic. */
void xyzProject::radInvert(const bool ri) { _radInvert = ri; }

/** Return true if the Sun-zenith angle limit logic is inverted. */
bool xyzProject::sunInvert() { return _sunInvert; }
/** Set the status of the Sun-zenith angle limit logic. */
void xyzProject::sunInvert(const bool si) { _sunInvert = si; }

/** Load the next particle in the file. */
void xyzProject::nextParticle() {
  static bool oneTimeSetup = false;
  static bool newFileSetup = false;
  static long nParticles;
  static long current = -1;
  Vector vejhat;

  if (!oneTimeSetup) {
    // load the first file on the first nextParticle call
    loadXyzfile();

    // set up the north pole unit vector and the vector along which
    // the equator crosses the Prime Meridian as defined by the Vernal
    // Equinox. [rectangular ecliptic coordinates]
    if (_npole[0] > -999) {
      // convert to radians
      float nPoleLambda = _npole[0] * M_PI / 180.0;
      float nPoleBeta   = _npole[1] * M_PI / 180.0;
      nPole[0] = cos(nPoleBeta) * cos(nPoleLambda);
      nPole[1] = cos(nPoleBeta) * sin(nPoleLambda);
      nPole[2] = sin(nPoleBeta);
      nPole = nPole.unit();

      // Project the Vernal Equinox onto the north pole, the rejection
      // is where vector along which the equator crosses the Prime
      // Meridian.  If the pole is the VE, then use ecliptic north.
      if ((_npole[0] + _npole[1]) == 0) {
	eq0 = Vector(0, 0, 1);
      } else {
	Vector X(1, 0, 0);
	eq0 = (X - nPole * (X * nPole)).unit();
      }
      eq90 = (nPole % eq0).unit();
    }

    oneTimeSetup = true;

    // we just loaded a new file
    throw(newFile);
  }

  if (!newFileSetup) {
    double jd = parameters().obsDate();
    _observer = projection(observerName(), jd, ecliptic());
    _observer.offset(_Coffset);

    // find the location of the comet
    Vector RHc, VHc;
    get_comet_xyz(parameters().comet().c_str(),
		  parameters().spkKernel().c_str(),
		  1, &jd, RHc.dblarr(), VHc.dblarr());
    cometRaDec(_observer.observe(RHc));
    // and the origin
    originRaDec(_observer.observe(Vector(0,0,0)));
    originOffset(coordOffset(originRaDec()));

    // if max is set, make sure it is less than or equal to the number
    // of particles in the file, otherwise set it to the number of
    // particles in the file
    if (_max >= 0) {
      nParticles = (_max > parameters().nParticles()) ? 
	parameters().nParticles() : _max;
    } else {
      nParticles = parameters().nParticles();
    }

    // get some needed info from the parameter set into the particle,
    // e.g., composition
    _p.parameters(parameters());
    if (parameters().pFunc().length()) {
      pfunctions pfunc;
      pfunc.setup(parameters(), _p);
      _p.grainT(_p.composition().name);
    }

    newFileSetup = true;
  }

  // about to read in the next particle
  current++;

  // Progress report every 10000 particles
  if ((current % 10000) == 0) {
    cout << "\r" << current << " completed, " << nParticles - current <<
      " remain.          ";
    cout.flush();
  }

  if (current >= nParticles) {
    cout << "\r" << current << " completed, " << nParticles - current <<
      " remain.          \n";
    cout.flush();

    // check for another file
    if (loadXyzfile()) {
      // there is...
      current = -1;
      newFileSetup = false;
      throw(newFile);
    }
    // otherwise, no more particles
    throw(nParticlesExceeded);
  }

  _p.radius(0);  // for radius test below
  xyzfile.readParticle(_p);
  if (!_p.error) {
    // Syndynes and pre-v0.7.3 xyz files may not have a correctly
    // computed radius.  In this case, compute radius from $\beta$.
    if (_p.radius() <= 0) {
      if (graindensity() > 0) {
	// use xyzProject's graindensity
	_p.updateRadius(graindensity());
      } else {
	if (_p.graindensity() > 0) {
	  // use the particle's grain density
	  _p.updateRadius();
	} else  {
	  // use the bulk density of the material
	  _p.updateRadius(_p.composition().bulkdensity);
	}
      }
    }
  } else {
    cout << "\r" << current << " completed, 0 remain.          \n";
    cout.flush();

    // check for another file
    if (loadXyzfile()) {
      // there is...
      current = -1;
      newFileSetup = false;
      throw(newFile);
    }
    // otherwise, stop with a read error
    throw(readError);
  }

  // needed for longtiude and latitude tests
  vejhat = _p.vej().unit();

  // (re)calculate this particle's ejected latitude and longitude from
  // the nucleus
  if (_npole[0] >= -999) {
    longlat origin = _p.origin();
    double x, y, offset;

    // Project vej onto the pole
    origin.beta = vejhat * nPole;
    origin.beta = 90.0 - acos(origin.beta) * 180.0 / M_PI;

    // The rejection is on the equator
    Vector rej = _p.vej() - nPole * (vejhat * nPole);
    x = rej * eq0;
    y = rej * eq90;

    origin.lambda = atan2(y, x) * 180.0 / M_PI + 360;
    origin.lambda = fmod(origin.lambda, 360); // branch cut at 0

    // Rotate the nucleus
    origin.lambda += _rotPhase;
    origin.lambda -= fmod(_rotRate * _p.age(), 360.0);

    // Branch cut at 0 deg
    while (origin.lambda < 0)
      origin.lambda += 360;
    origin.lambda = fmod(origin.lambda, 360.0);

    _p.origin(origin);
  }

  // if ranges are set, keep/throw away this particle

  // age limit
  if (_ageRange.size()) {
    bool keep = false;
    if ((_ageRange[0] <= _p.age()) && (_p.age() <= _ageRange[1]))
      keep = true;
    if (_ageInvert) keep = !keep;
    if (!keep) throw(ageLimit);
  }

  // beta limit
  if (_betaRange.size()) {
    bool keep = false;
    if ((_betaRange[0] <= _p.beta()) && (_p.beta() <= _betaRange[1]))
      keep = true;
    if (_betaInvert) keep = !keep;
    if (!keep) throw(betaLimit);
  }

  // latitude limit
  if (_latRange.size()) {
    bool keep = false;
    if ((_latRange[0] <= _p.origin().beta) &&
	(_p.origin().beta <= _latRange[1]))
      keep = true;
    if (_latInvert) keep = !keep;
    if (!keep) throw(latLimit);
  }

  // longitude limit
  if (_lonRange.size()) {
    bool keep = false;
    // This test requires origin.lambda be between 0 and 360.  Take
    // care around longitude = 0.
    if ((_lonRange[0] <= _lonRange[1])) {
      if ((_lonRange[0] <= _p.origin().lambda) &&
	  (_p.origin().lambda <= _lonRange[1]))
	keep = true;
    } else {
      if ((_lonRange[0] <= _p.origin().lambda) ||
	  (_p.origin().lambda <= _lonRange[1]))
	keep = true;
    }

    if (_lonInvert) keep = !keep;
    if (!keep) throw(lonLimit);
  }

  // radius limit
  if (_radRange.size()) {
    bool keep = false;
    if ((_radRange[0] <= _p.radius()) && (_p.radius() <= _radRange[1]))
      keep = true;
    if (_radInvert) keep = !keep;
    if (!keep) throw(radLimit);
  }

  // sun angle limit
  if (_sunRange.size()) {
    bool keep = false;
    float z_sun = acos(vejhat * (_p.istate().r.unit() * -1.0)) * 180.0 / M_PI;
    if ((_sunRange[0] <= z_sun) && (z_sun <= _sunRange[1])) keep = true;
    if (_latInvert) keep = !keep;
    if (!keep) throw(sunLimit);
  }

  // velocity limit
  if (_vejGen.v0() >= 0) {
    if (_p.vej().length() > 0) {
      // delete?
      // _vejGen.getSpeed(_p.beta(), _p.istate().r);
      throw(velocLimit);
    }
  }

  // rh limit
  if (_rhlimit > 0) {
    if (_p.istate().r.length() > _rhlimit) {
      throw(rhlimit);
    }
  }
}

/** Project the current particle onto the sky for the observer. */
longlat xyzProject::projectParticle() {
  return _observer.observe(_p.fstate().r);
}

/** Determine the current particle offset from the comet in
    arcseconds.  This is the gnomonic (TAN) projection.  See
    Calabretta and Greisen 2002. */
longlat xyzProject::particleOffset() {
  return coordOffset(projectParticle());
}

/** Determine the coordinate offset from the comet in arcseconds.  This
    is the gnomonic (TAN) projection.  See Calabretta and Greisen
    2002.

    \todo The projection doesn't agree with synover.pro or syn2ds9 for
    very large images.
*/
longlat xyzProject::coordOffset(longlat target) {
  double l1 = cometRaDec().lambda * M_PI / 180;
  double l2 = target.lambda * M_PI / 180;
  double b1 = cometRaDec().beta * M_PI / 180;
  double b2 = target.beta * M_PI / 180;

  // distance and position angle from the center of the projection
  double th, ph;
  double a, b, c; // temp variables

  a = cos(b1) * sin(b2) - sin(b1) * cos(b2) * cos(l2 - l1);
  b = cos(b2) * sin(l2 - l1);
  c = sin(b1) * sin(b2) + cos(b1) * cos(b2) * cos(l2 - l1);

  // the center of the projection is theta = 90 degr
  th = atan2(sqrt(a*a + b*b), c) + M_PI_2;
  ph = atan2(sin(l1 - l2), cos(b2) * tan(b1) - sin(b2) * cos(l1 - l2));

  /* x = 180 * cot(th) / pi * sin(ph) * 3600
     y = -180 * cot(th) / pi * cos(ph) * 3600 */
  longlat offset;
  offset.lambda = sin(ph) / tan(th) * 206264.8;
  offset.beta = cos(ph) / tan(th) * 206264.8;
  //  cerr << (th - M_PI_2) * 206265 << " " << ph * 206265 << " " << offset.lambda << " " << offset.beta << "\n";
  return offset;
}

/** Determine the particle offset from the comet in arcseconds.  This
    is the Sanson-Flamsteed projeciton.  See Calabretta and Greisen
    2002.
longlat xyzProject::particleOffset(particle p) {
  longlat radec = projectParticle(p);

  longlat offset;
  offset.lambda = (radec.lambda - cometRaDec().lambda) * 
    cos(radec.beta * M_PI / 180) * 3600;
  offset.beta = (radec.beta - cometRaDec().beta) * 3600;

  return offset;
} */

/** Return the current parameter set. */
paramSet xyzProject::parameters() { return _parameters; }

/** Set the current parameters to integrate by. */
void xyzProject::parameters(paramSet p) { _parameters = p; }

/** Return the current particle. */
particle xyzProject::p() { return _p; }

/** Set the current paraticle. */
void xyzProject::p(particle p) { _p = p; }
