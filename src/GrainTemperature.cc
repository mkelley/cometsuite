/***************************************************************************

  Implements grain temperatures.

  Copyright (C) 2009-2010 by Michael S. Kelley <msk@astro.umd.edu>

 ***************************************************************************/

#if HAVE_CONFIG_H
#  include <config.h>
#endif /* HAVE_CONFIG_H */

#include <valarray>
#include <cmath>
#include <CCfits/CCfits>
#include <unistd.h>
#include "GrainTemperature.h"
#include "rundynamics.h"

#ifndef _DATADIR
#define _DATADIR "./"
#endif

using namespace std;

GrainTemperature::GrainTemperature() {
  dataLoaded = false;  // set to true when a data file is loaded
}

/** Open a FITS file and populate the internal arrays. */
void GrainTemperature::load(const materials c) {
  string filename;

  // See if the FITS composition files are in the current directory
  switch (c) {
  case GEOMETRIC: filename = string("./geometric.fits"); break;
  case AM_CARBON: filename = string("./am-carbon.fits"); break;
  case AM_OLIVINE50: filename = string("./am-olivine50.fits"); break;
  }
  if (access(filename.c_str(), F_OK)) {
    // No?  OK, see if the FITS composition files are in _DATADIR
    switch (c) {
    case GEOMETRIC: filename = string(_DATADIR) + "/geometric.fits"; break;
    case AM_CARBON: filename = string(_DATADIR) + "/am-carbon.fits"; break;
    case AM_OLIVINE50: filename = string(_DATADIR) + "/am-olivine50.fits"; break;
    }

    if (access(filename.c_str(), F_OK)) {
      // Again, no?  We might be in the distribution directory.  Check
      // ../data/
      switch (c) {
      case GEOMETRIC: filename = string("../data/geometric.fits"); break;
      case AM_CARBON: filename = string("../data/am-carbon.fits"); break;
      case AM_OLIVINE50: filename = string("../data/am-olivine50.fits"); break;
      }

      // If it still doesn't work, the program will continue
      if (access(filename.c_str(), F_OK)) {
	cerr << "Cannot find FITS composition file:" << filename << "\n";
	cerr << "I searched: \"./\", \"" << _DATADIR
	     << "\", and \"../data\"" << "\n";
	cerr << "Defaulting to blackbody temperatures.\n";
      }
    }
  }

  try {
    auto_ptr<CCfits::FITS> pInfile(new CCfits::FITS(filename, CCfits::Read, true));
    CCfits::PHDU& image = pInfile->pHDU();
    image.read(_T);

    CCfits::ExtHDU& table1 = pInfile->extension(1);
    table1.column(1).read(_radius, 1, table1.column(1).rows());  // size

    CCfits::ExtHDU& table2 = pInfile->extension(2);
    table2.column(1).read(_rh, 1, table2.column(1).rows());  // rh

    CCfits::ExtHDU& table3 = pInfile->extension(3);
    table3.column(1).read(_D, 1, table3.column(1).rows());  // D

    dataLoaded = true;
    naxis[0] = _D.size();
    naxis[1] = _rh.size();
    naxis[2] = _radius.size();
  } catch (CCfits::FitsException&) {
    dataLoaded = false;    
  }
}

/** Return the temperature of the grain, interpolating to the given D,
    a (microns), and rh (km).

    The function selects the 2x2x2 cube closest in D-, radius-,
    rh-space to the desired D, radius, and rh.  Each corner of the
    cube is given a weight based on the distance to the desired point.

    \todo Reduce the amount of code.
 */
double GrainTemperature::interp(const double D,
				const double a,
				const double rh) {
  // if we don't have a solution, return the blackbody temperature
  if (!dataLoaded) return 278. / sqrt(rh / _AU);

  // find the nearest neighbors for each dimension
  int nn[3][2];  // indices for D, a, rh
  double w[3][2];  // weights
  double d;
  int* temp;

  // get the values for the nearest neighbors, and compute their
  // interpolation weights
  temp = nearest(_D, D);
  nn[0][0] = temp[0];
  nn[0][1] = temp[1];
  d = abs(_D[temp[1]] - _D[temp[0]]);  // distance between boundaries
  if (d == 0) {
    w[0][0] = 1;
    w[0][1] = 1;
  } else {
    w[0][0] = abs(D - _D[temp[1]]) / d;  // normalized distance to right boundary
    w[0][1] = abs(_D[temp[0]] - D) / d;  // normalized distance to left boundary
  }
  delete temp;

  temp = nearest(_radius, a);
  nn[1][0] = temp[0];
  nn[1][1] = temp[1];
  d = abs(_radius[temp[1]] - _radius[temp[0]]);  // distance between boundaries
  if (d == 0) {
    w[1][0] = 0.5;
    w[1][1] = 0.5;
  } else {
    w[1][0] = abs(a - _radius[temp[1]]) / d;  // normalized distance to left boundary
    w[1][1] = abs(_radius[temp[0]] - a) / d;  // normalized distance to right boundary
  }
  delete temp;

  temp = nearest(_rh, rh / _AU );
  nn[2][0] = temp[0];
  nn[2][1] = temp[1];
  d = abs(_rh[temp[1]] - _rh[temp[0]]);  // distance between boundaries
  if (d == 0) {
    w[2][0] = 0.5;
    w[2][1] = 0.5;
  } else {
    w[2][0] = abs(rh - _rh[temp[1]]) / d;  // normalized distance to left boundary
    w[2][1] = abs(_rh[temp[0]] - rh) / d;  // normalized distance to right boundary
  }
  delete temp;

  // weighted average of the temperatures
  double sum_t = 0, sum_w = 0;
  for (int i=0; i<2; i++) {
    for (int j=0; j<2; j++) {
      for (int k=0; k<2; k++) {
	int index = (nn[0][i] * naxis[2] + nn[1][j]) * naxis[1] + nn[2][k];
	sum_t += _T[index] * w[0][i] * w[1][j] * w[2][k];
	sum_w += w[0][i] * w[1][j] * w[2][k];
      }
    }
  }
  return sum_t / sum_w;
}

/** Return the indices of the two values of v nearest to vv. */
int* GrainTemperature::nearest(valarray<float> v, const double vv) {
  int *nn = new int[2];

  // is vv before i=0?
  if (((v[0] < v[1]) && (vv <= v[0])) ||
      ((v[0] > v[1]) && (vv >= v[0]))) {
    nn[0] = 0;
    nn[1] = 0;
    return nn;  // done!
  }

  for (int i=1; i<v.size(); i++) {
    // is vv between i-1 and i?
    if (((v[i-1] <= vv) && (vv <= v[i])) || ((v[i] <= vv) && (vv <= v[i-1]))) {
      // between i-1 and i
      nn[0] = i-1;
      nn[1] = i;
      return nn;  // done!
    }

    // are we on the last index?
    if (i == (v.size() - 1)) {
      // then the vv is beyond v[i]
      nn[0] = v.size() - 1;
      nn[1] = v.size() - 1;
      return nn;  // done!
    }
  }

  // should never get to this point
  nn[0] = 0;
  nn[1] = 0;
  return nn;
}
