/***************************************************************************

  Produces histograms from CometSuite data files.

  Copyright (C) 2005-2009 by Michael S. Kelley <msk@astro.umd.edu>

 ***************************************************************************/

#if HAVE_CONFIG_H
#  include <config.h>
#endif /* HAVE_CONFIG_H */

#include <string>
#include <vector>
#include <algorithm>
#include <cstdlib>
#include <cmath>
#include "xyzHist.h"
#include "rundynamics.h"
#include "particle.h"
#include "longlat.h"

using namespace std;

/** Default parameters. */
xyzHist::xyzHist() {
  _aperture = 3;
  _bin = sizeBin;
  _binSize = -1;
}

/** Return the output file name. */
string xyzHist::outfileName() { return _outfileName; }
/** Set the output file name. */
void xyzHist::outfileName(const string fn) { _outfileName = fn; }

/** Return the observer's aperture size. */
float xyzHist::aperture() { return _aperture; }
/** Set the observer's aperture size. */
void xyzHist::aperture(const float a) { _aperture = a; }

/** Return the particle bin type. */
int xyzHist::bin() { return _bin; }
/** Set the particle bin type. */
void xyzHist::bin(const int b) { _bin = b; }
/** Set the particle bin type. */
void xyzHist::bin(const string b) {
  string c = b;
  transform(c.begin(), c.end(), c.begin(), (int(*)(int))tolower);
  if (c.find("size") != string::npos)
    _bin = sizeBin;
  else if (c.find("beta") != string::npos)
    _bin = betaBin;
  else
    _bin = ageBin;
}

/** Return the bin size. */
float xyzHist::binSize() { return _binSize; }
/** Set the bin size. */
void xyzHist::binSize(const float bs) { _binSize = bs; }

/** Return the particle x distribution. */
vector<float> xyzHist::dn_dx() { return _dn_dx; }
/** Set the particle x distribution. */
void xyzHist::dn_dx(const vector<float> dd) { _dn_dx = dd; }

/** Return the particle x distribution error. */
vector<float> xyzHist::dn_dx_err() { return _dn_dx_err; }
/** Set the particle x distribution error. */
void xyzHist::dn_dx_err(const vector<float> dde) { _dn_dx_err = dde; }

/** Return the particle x distribution x values. */
vector<float> xyzHist::xBin() { return _xBin; }
/** Set the particle x distribution x values. */
void xyzHist::xBin(const vector<float> xb) { _xBin = xb; }

/** Load particles, project onto the image plane, remove particles
    outside the observer's aperture, sort into a histogram. */
void xyzHist::createPSD() {
  // initial bins set up
  _xBin.push_back(0);
  _dn_dx.push_back(0);
  _dn_dx_err.push_back(0);

  if (DEBUG) cerr << "Begin reading particles (xyzHist)\n";
  bool finished = false;
  while (!finished) {
    try {
      nextParticle();

      // observe the particle and get the offset in arcseconds; check if
      // the particle is in the observer's aperture
      longlat radec = projectParticle();
      longlat offset = particleOffset();
      
      double ra1 = radec.lambda * M_PI / 180,
	ra2 = cometRaDec().lambda * M_PI / 180,
	dec1 = radec.beta * M_PI / 180,
	dec2 = cometRaDec().beta * M_PI / 180;
      double rx = cos(dec1) * sin(dec2) - sin(dec1) * cos(dec2) * 
	cos(ra2 - ra1);
      double ry = cos(dec2) * sin(ra2 - ra1);
      double rz = sin(dec1) * sin(dec2) + cos(dec1) * cos(dec2) * 
	cos(ra2 - ra1);

      double d = atan2(sqrt(rx*rx + ry*ry), rz) * 206265;

      if ((d <= _aperture / 2) || (_aperture < 0)) {
	// place the particle in a bin
	double x;
	switch (_bin) {
	case sizeBin: x = _p.radius(); break;
	case betaBin: x = _p.beta(); break;
	case ageBin: x = _p.age() / 86400; break;
	}

	// weight the particle if necessary
	double rh = _p.fstate().r.length();
	double Delta = (_p.fstate().r - _observer.r()).length();

	double w = afrhoWeight(_p.istate().r.length());
	w *= thermalWeight(_p.radius(), rh, Delta,
			   _p.grainT().interp(_p.fractalDim(), _p.radius(), rh));
	w *= scatteringWeight(_p.radius());
	w *= nuclearPsdWeight(_p.radius()); 
	w *= fileScale();

	if (w > 0) {
	  if (_binSize < 0) {
	    _xBin.push_back(x);
	    _dn_dx.push_back(w);
	    _dn_dx_err.push_back(w * w);
	  } else {
	    long b = static_cast<int>(x / _binSize);
	    if (b >= _xBin.size()) {
	      // then add more elements
	      for (int i=_xBin.size(); i<=b; i++) {
		_xBin.push_back(_xBin[_xBin.size()-1] + _binSize);
		_dn_dx.push_back(0);
		_dn_dx_err.push_back(0);
	      }
	    }
	    _dn_dx[b] += w;
	    // need to think about the error a little bit more
	    _dn_dx_err[b] += w * w;
	  }
	}
      }
    } catch (xyzProjectFlags flag) {
      switch (flag) {
      case newFile:
	nextFileScale();
	break;
      case nParticlesExceeded:
      case readError:
	finished = true;
	break;
      case error:
	cerr << "Particle error!\n";
	finished = true;
	break;
      }
    }
  }
  if (DEBUG) cerr << "End reading particles (xyzHist)\n";

  // the error is the sqrt(sum(w*w)), i.e. sqrt(count in bin)
  transform(_dn_dx_err.begin(), _dn_dx_err.end(), _dn_dx_err.begin(),
	    (float(*)(float))sqrt);
}

/** Write the particle distribution output file. */
bool xyzHist::writePSD() {
  ofstream outfileStream(_outfileName.c_str(), ios::out);
  outfileStream.precision(8);

  if (!outfileStream.good()) {
    cerr << "Error opening file: " << _outfileName << "\n";
    return false;
  }

  outfileStream << "# binning ";
  switch (_bin) {
  case sizeBin: outfileStream << "size/radius (microns)\n"; break;
  case betaBin: outfileStream << "beta\n"; break;
  case ageBin: outfileStream << "age (seconds)\n"; break;
  }

  outfileStream << "# bin width = " << _binSize << "\n";
  outfileStream << "# bin\tnumber\terror\n";

  for (int i=0; i<_xBin.size(); i++)
    outfileStream << _xBin[i] << "\t" << _dn_dx[i] << "\t" <<
      _dn_dx_err[i] << "\n";
  return true;
}
