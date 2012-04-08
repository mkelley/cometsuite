/***************************************************************************

  Histogram particles onto the sky and save to a FITS file.

  Copyright (C) 2005-2010 by Michael S. Kelley <msk@astro.umd.edu>

 ***************************************************************************/

#if HAVE_CONFIG_H
#  include <config.h>
#endif /* HAVE_CONFIG_H */

#include <string>
#include <sstream>
#include <vector>
#include <valarray>
#include <CCfits/CCfits>
#include "xyzImage.h"
#include "round.h"
#include "longlat.h"
#include "rundynamics.h"

using namespace std;

/** Default image and input parameters. */
xyzImage::xyzImage() {
  _writeFits = true;
  _outfileName = "out.fits";
  _size[0] = _size[1] = 512;
  _platescale[0] = -1;
  _platescale[1] = 1;

  // declare auto-pointer to FITS. Ensures no resources
  // leaked if something fails in dynamic allocation.
  pFits = auto_ptr<CCfits::FITS>(0);
}

/** Return fits file writing status. */
bool xyzImage::writeFits() { return _writeFits; }
/** Set fits file writing status. */
void xyzImage::writeFits(const bool wf) { _writeFits = wf; }

/** Return the output file name. */
string xyzImage::outfileName() { return _outfileName; }
/** Set the output file name. */
void xyzImage::outfileName(const string fn) { _outfileName = fn; }

/** Return the image size. */
valarray<long> xyzImage::size() {
  valarray<long> sz(2);
  sz[0] = _size[0];
  sz[1] = _size[1];
  return sz;
}
/** Set the image size. */
void xyzImage::size(const long sz) {
  _size[0] = sz;
  _size[1] = sz;
}
/** Set the image size. */
void xyzImage::size(const long* sz) {
  _size[0] = sz[0];
  _size[1] = sz[1];
}
/** Set the image size. */
void xyzImage::size(const valarray<long> sz) {
  if (sz.size() > 1) {
    _size[0] = sz[0];
    _size[1] = sz[1];
  } else {
    _size[0] = sz[0];
    _size[1] = sz[0];
  }    
}

/** Return the platescale. */
valarray<float> xyzImage::platescale() {
  valarray<float> ps(2);
  ps[0] = _platescale[0];
  ps[1] = _platescale[1];
  return ps;
}
/** Set the platescale. */
void xyzImage::platescale(const float* ps) {
  _platescale[0] = ps[0];
  _platescale[1] = ps[1];
}
/** Set the platescale. */
void xyzImage::platescale(const float ps) {
  _platescale[0] = -ps;
  _platescale[1] = ps;
}
/** Set the platescale. */
void xyzImage::platescale(const valarray<float> ps) {
  if (ps.size() > 1) {
    _platescale[0] = ps[0];
    _platescale[1] = ps[1];
  } else {
    _platescale[0] = -ps[0];
    _platescale[1] = ps[0];
  }
}

/** Adds a comment to the primary image header. */
void xyzImage::addFitsComment(const string comment) {
  pFits->pHDU().writeComment(comment);
}

/** Prepare a fits image, load particles, project onto the image
    plane, sort into image arrays.  There are six arrays: simulated
    emission, number, average beta, average radius, average age, and
    average ejection latitude.  The images are fully weighted (thermal,
    scattering, afrho, and psd) except for "number".  */
void xyzImage::createImages() {
  if (_writeFits) {
    try {
      // Create a new FITS object, specifying the data type and axes for
      // the primary image. Simultaneously create (and overwrite) the
      // corresponding file.
      pFits.reset(new CCfits::FITS(string("!") + _outfileName, DOUBLE_IMG,
				   2, _size));
    } catch (CCfits::FITS::CantCreate) {
      // ... or not, as the case may be.
      throw(cantCreateFits);
    }

    // add the date to the fits header
    pFits->pHDU().writeDate();

    // add a version stamp
    pFits->pHDU().addKey(string("RDVERSN"), VERSION, "RunDynamics / xyz2fits version");
  }

  // set up the image arrays
  _center[0] = static_cast<float>(_size[0]) / 2;
  _center[1] = static_cast<float>(_size[1]) / 2;

  _weighted.resize(_size[0] * _size[1], 0);
  _number.resize(_size[0] * _size[1], 0);
  _beta.resize(_size[0] * _size[1], 0);
  _radius.resize(_size[0] * _size[1], 0);
  _age.resize(_size[0] * _size[1], 0);
  _latitude.resize(_size[0] * _size[1], 0);

  bool finished = false;
  while (!finished) {
    try {
      nextParticle();
      longlat offset = particleOffset();

      /* Determine the location in the image.  The -0.5 shifts the
	 coordinate system from the lower left corner to the center of
	 the pixel (FITS standard).  Then add this particle in; betas
	 are averaged such that they yield the correct particle size
	 for beta ~ 1/a. */
      int x = static_cast<int>(round(offset.lambda / _platescale[0] +
				     _center[0] - 0.5));
      int y = static_cast<int>(round(offset.beta / _platescale[1] +
				     _center[1] - 0.5));
      int j = index(x, y);
      if (j >= 0) {
	double rh = _p.fstate().r.length();
	double Delta = (_p.fstate().r - _observer.r()).length();

	double w = afrhoWeight(_p.istate().r.length());
	w *= thermalWeight(_p.radius(), rh, Delta,
			   _p.grainT().interp(_p.fractalDim(), _p.radius(), rh)) ;	
	w *= scatteringWeight(_p.radius());
	w *= nuclearPsdWeight(_p.radius()); 
	w *= fileScale();

	if (w > 0) {
	  _number[j]++;
	  _beta[j] += w / _p.beta();
	  _radius[j] += w * _p.radius();
	  _age[j] += w * _p.age() / 86400;
	  _latitude[j] += w * _p.origin().beta;
	  _weighted[j] += w;
	}
      }
    } catch (xyzProjectFlags flag) {
      if (flag == newFile) {
	if (_writeFits) addXyzHeader();
	nextFileScale();
      }
      if (flag == nParticlesExceeded) finished = true;
      if (flag == readError) finished = true;
    }
  }

  // Finish the averages.  If thermal or scattered is on, then convert
  // weighted to surface brightness (MJy/sr)
  for (int i=0; i<_size[0]*_size[1]; i++) {
    if (_number[i] > 0) {
      _beta[i] = 1 / (_beta[i] / _weighted[i]);
      _radius[i] = _radius[i] / _weighted[i];
      _age[i] = _age[i] / _weighted[i];
      _latitude[i] = _latitude[i] / _weighted[i];
    }
  }
  if (thermalMode() || scatteringMode()) {
    _weighted *= 2.350439e-5 * abs(_platescale[0] * _platescale[1]);
  }
}

/** Writes the fits image.  I am using the example, cookbook.cxx, from
    the CCfits documentation.  Throws xyzImageFlags on error. */
void xyzImage::writeImages() {
  if (!_writeFits) {
    cerr << "Fits writing is diabled!\n";
    throw(fitsWritingDisabled);
  }

  long nElements = _size[0] * _size[1];

  // add the primary data
  pFits->pHDU().write(1, nElements, _weighted);
  pFits->pHDU().writeComment(string("------------------------------------------"));
  pFits->pHDU().writeComment(string("Simulated dust observation, weighted by   "));
  pFits->pHDU().writeComment(string("dust production, thermal emission, light  "));
  pFits->pHDU().writeComment(string("scattering, and/or particle size          "));
  pFits->pHDU().writeComment(string("distribution.  Thermal emission is in     "));
  pFits->pHDU().writeComment(string("units of MJy/sr.  Dust production is      "));
  pFits->pHDU().writeComment(string("normalized to 1.0 at 1 AU.                "));
  pFits->pHDU().writeComment(string("Divide by 2.350439e-5 * platescale**2 to  "));
  pFits->pHDU().writeComment(string("get image weights.                        "));
  pFits->pHDU().writeComment(string("------------------------------------------"));

  // add astrometry information to the primary image
  addWCS(&(pFits->pHDU()));
  pFits->pHDU().addKey(string("BUNIT"), "MJy/sr", "Only thermal emission is calibrated");

  // add number, beta, radius, age, and latitude images as fits extensions
  vector<long int> extAx(2);
  extAx[0] = _size[0];
  extAx[1] = _size[1];
  string newName("number");
  CCfits::ExtHDU* imageExt = pFits->addImage(newName, USHORT_IMG, extAx);
  imageExt->write(1, nElements, _number);
  imageExt->writeComment(string("------------------------------------------"));
  imageExt->writeComment(string("Number of particles per pixel, unweighted."));
  imageExt->writeComment(string("------------------------------------------"));
  addWCS(imageExt); // add astrometry information
  imageExt->addKey(string("BUNIT"), "number", "");

  newName = "beta";
  imageExt = pFits->addImage(newName, FLOAT_IMG, extAx);
  imageExt->write(1, nElements, _beta);
  imageExt->writeComment(string("------------------------------------------"));
  imageExt->writeComment(string("Average beta, weighted by dust production,"));
  imageExt->writeComment(string("thermal emission, scattered light, and/or "));
  imageExt->writeComment(string("particle size distribution.  Beta is      "));
  imageExt->writeComment(string("averaged such that a ~ 1/beta(pixel)      "));
  imageExt->writeComment(string("yields the correct average particle size. "));
  imageExt->writeComment(string("------------------------------------------"));
  addWCS(imageExt); // add astrometry information
  imageExt->addKey(string("BUNIT"), "none", "");

  newName = "radius";
  imageExt = pFits->addImage(newName, FLOAT_IMG, extAx);
  imageExt->write(1, nElements, _radius);
  imageExt->writeComment(string("--------------------------------------------"));
  imageExt->writeComment(string("Average radius, weighted by dust production,"));
  imageExt->writeComment(string("thermal emission, scattered light, and/or   "));
  imageExt->writeComment(string("particle size distribution [units: microns]."));
  imageExt->writeComment(string("--------------------------------------------"));
  addWCS(imageExt); // add astrometry information
  imageExt->addKey(string("BUNIT"), "microns", "");

  newName = "age";
  imageExt = pFits->addImage(newName, FLOAT_IMG, extAx);
  imageExt->write(1, nElements, _age);
  imageExt->writeComment(string("-----------------------------------------"));
  imageExt->writeComment(string("Average age, weighted by dust production,"));
  imageExt->writeComment(string("thermal emission, scattered light, and/or"));
  imageExt->writeComment(string("particle size distribution [units: days]."));
  imageExt->writeComment(string("-----------------------------------------"));
  addWCS(imageExt); // add astrometry information
  imageExt->addKey(string("BUNIT"), "days", "");

  newName = "latitude";
  imageExt = pFits->addImage(newName, FLOAT_IMG, extAx);
  imageExt->write(1, nElements, _latitude);
  imageExt->writeComment(string("-----------------------------------"));
  imageExt->writeComment(string("Average latitude, weighted by dust "));
  imageExt->writeComment(string("production, thermal emission,      "));
  imageExt->writeComment(string("scattered light, and/or particle   "));
  imageExt->writeComment(string("size distribution [units: degrees]."));
  imageExt->writeComment(string("-----------------------------------"));
  addWCS(imageExt); // add astrometry information
  imageExt->addKey(string("BUNIT"), "degrees", "");
}

/** Returns the array index of pixel x, y. */
long xyzImage::index(const int x, const int y) {
  int j = static_cast<long>(y * _size[0] + x);

  if (x >= _size[0]) return -1;
  if (x <  0)        return -1;
  if (y >= _size[1]) return -1;
  if (y <  0)        return -1;

  return j;
}

/** Adds the world coordinate system and the origin to the fits header. */
void xyzImage::addWCS(CCfits::HDU* hdu) {
  hdu->writeComment(string("----------------------"));
  hdu->writeComment(string("Astrometry Information"));
  hdu->writeComment(string("----------------------"));

  if (ecliptic()) {
    hdu->addKey("CRVAL1", cometRaDec().lambda, "Ecliptic longitude at reference point");
    hdu->addKey("CRVAL2", cometRaDec().beta, "Ecliptic latitude at reference point");
  } else {
    hdu->addKey("CRVAL1", cometRaDec().lambda, "RA at reference point");
    hdu->addKey("CRVAL2", cometRaDec().beta, "Dec at reference point");
  }

  hdu->addKey("CRPIX1", _center[0]+1, "Pixel coordinate of reference point");
  hdu->addKey("CRPIX2", _center[1]+1, "Pixel coordinate of reference point");
  hdu->addKey("CDELT1", _platescale[0] / 3600, "Platescale [deg/pix]");
  hdu->addKey("CDELT2", _platescale[1] / 3600, "Platescale [deg/pix]");

  if (ecliptic()) {
    hdu->addKey("CTYPE1", string("ELON-TAN"), "Gnomonic projection");
    hdu->addKey("CTYPE2", string("ELAT-TAN"), "Gnomonic projection");
  } else {
    hdu->addKey("CTYPE1", string("RA---TAN"), "Gnomonic projection");
    hdu->addKey("CTYPE2", string("DEC--TAN"), "Gnomonic projection");
  }

  /*
  if (ecliptic()) {
    hdu->addKey("CTYPE1", string("ELON-SFL"), "Sanson-Flamsteed projection");
    hdu->addKey("CTYPE2", string("ELAT-SFL"), "Sanson-Flamsteed projection"); 
  } else {
    hdu->addKey("CTYPE1", string("RA---SFL"), "Sanson-Flamsteed projection");
    hdu->addKey("CTYPE2", string("DEC--SFL"), "Sanson-Flamsteed projection");
  }
  */

  if (ecliptic()) {
    hdu->addKey("ORIGINR", originRaDec().lambda, "Origin ecliptic longitude");
    hdu->addKey("ORIGIND", originRaDec().beta, "Origin ecliptic latitude");
  } else {
    hdu->addKey("ORIGINR", originRaDec().lambda, "Origin RA");
    hdu->addKey("ORIGIND", originRaDec().beta, "Origin Dec");
  }

  hdu->addKey("ORIGINX", round(originOffset().lambda / _platescale[0] +
			       _center[0] - 0.5), "Origin Y pixel coordinate");
  hdu->addKey("ORIGINY", round(originOffset().beta / _platescale[1] +
			       _center[1] - 0.5), "Origin Y pixel coordinate");
}

/** Adds the xyzfile header to the fits header. */
void xyzImage::addXyzHeader() {
  stringstream s;
  s << currentFileNumber();
  string n = s.str();

  // RunDynamics parameters
  pFits->pHDU().writeComment(string("------------------------"));
  pFits->pHDU().writeComment(string("RunDynamics Parameters"));
  pFits->pHDU().writeComment(string("------------------------"));
  pFits->pHDU().addKey(string("RDFILE") + n,
		       currentFileName(), "RunDynamics file name");
  pFits->pHDU().addKey(string("RDPROG") + n,
		       parameters().program(), "RunDynamics program name");
  pFits->pHDU().addKey(string("RDCOMET") + n,
		       parameters().comet(), "Comet name");
  pFits->pHDU().addKey(string("RDKERN") + n,
		       parameters().spkKernel(), "Comet kernel");
  pFits->pHDU().addKey(string("RDJD") + n,
		       parameters().obsDate(), "Julian date of observation");
  pFits->pHDU().addKey(string("RDXYZFL") + n,
		       parameters().outFile(), "RunDynamics output file");
  pFits->pHDU().addKey(string("RDPFUNC") + n,
		       parameters().pFunc(), "Particle generator");
  pFits->pHDU().addKey(string("RDTOL") + n,
		       parameters().tolerance(), "RA15 tolerance");
  pFits->pHDU().addKey(string("RDPLANT") + n,
		       parameters().planets(), "Planets used");
  pFits->pHDU().addKey(string("RDPLU") + n,
		       parameters().planetLookUp(), "Planet look up enabled");
  pFits->pHDU().addKey(string("RDCA") + n,
		       parameters().closeApproaches(), "Special handling of close approaches");
  pFits->pHDU().addKey(string("RDBOX") + n,
    parameters().box(), "Integration box size");
  pFits->pHDU().addKey(string("RDLTT") + n,
    parameters().lightTravelTime(), "Light travel time correction");
  stringstream b;
  for (int i = 0; i < parameters().beta().size(); i++) {
    b << " " << parameters().beta()[i];
  }
  pFits->pHDU().addKey(string("RDBETA") + n,
		       b.str(), "Syndyne beta list");
  pFits->pHDU().addKey(string("RDNDAYS") + n,
		       parameters().nDays(), "Maximum syndyne age");
  pFits->pHDU().addKey(string("RDSTEPS") + n,
		       parameters().steps(), "Syndyne steps");
  pFits->pHDU().addKey(string("RDORBIT") + n,
		       parameters().orbit(), "Size of orbit for syndynes");
  pFits->pHDU().addKey(string("RDNP") + n,
		       parameters().nParticles(), "Number of particles for Make Comet");
}
