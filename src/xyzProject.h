/***************************************************************************

  Copyright (C) 2005-2010 by Michael S. Kelley <msk@astro.umd.edu>

 ***************************************************************************/

#if !defined(__XYZPROJECT)
#define __XYZPROJECT 1

#if HAVE_CONFIG_H
#  include <config.h>
#endif /* HAVE_CONFIG_H */

#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <valarray>
#include "xyzstream.h"
#include "Vej.h"
#include "projection.h"
#include "Vector.h"
#include "rundynamics.h"
#include "longlat.h"
#include "paramSet.h"
#include "particle.h"

using namespace std;

enum xyzProjectFlags { newFile, lastParticle, noParticle, nParticlesExceeded,
		       betaLimit, ageLimit, latLimit, radLimit, sunLimit,
		       velocLimit, rhlimit, readError, error };

extern "C" int get_comet_xyz(const char *comet, const char *cometSPK, int npts,
                        double *jd, double *r, double *v);

/** Handles the loading and sky projection of xyzfiles. */
class xyzProject {
 public:
  xyzProject();
  ~xyzProject();

  // nice friend
  friend class xyzHist;
  friend class xyzImage;

  // parameter setting and file loading
  bool verbose();
  void verbose(const bool);
  vector<string> xyzfileNames();
  void xyzfileNames(const string);
  string currentFileName();
  int currentFileNumber();
  paramSet parameters();
  void parameters(paramSet);
  bool loadXyzfile();
  double graindensity();
  void graindensity(const double);
  bool ecliptic();
  void ecliptic(const bool);
  string observerName();
  valarray<float> npole();
  void npole(const float*);
  void npole(const valarray<float>);
  void observerName(const string);
  longlat offset();
  void offset(const valarray<float>);
  void offset(const longlat);
  longlat cometRaDec();
  void cometRaDec(const longlat);
  longlat originRaDec();
  void originRaDec(const longlat);
  longlat originOffset();
  void originOffset(const longlat);
  long max();
  void max(const long);
  valarray<float> ageRange();
  void ageRange(const valarray<float>);
  valarray<float> betaRange();
  void betaRange(const valarray<float>);
  valarray<float> latRange();
  void latRange(const valarray<float>);
  valarray<float> radRange();
  void radRange(const valarray<float>);
  valarray<float> sunRange();
  void sunRange(const valarray<float>);
  double vLimit();
  void vLimit(const double);
  double rhLimit();
  void rhLimit(const double);
  bool ageInvert();
  void ageInvert(const bool);
  bool betaInvert();
  void betaInvert(const bool);
  bool latInvert();
  void latInvert(const bool);
  bool radInvert();
  void radInvert(const bool);
  bool sunInvert();
  void sunInvert(const bool);

  // interfaces to particles, parameters, and positions
  particle p();
  void p(particle);
  void nextParticle();
  longlat projectParticle();
  longlat particleOffset();
  longlat coordOffset(longlat);

 private:
  bool _verbose;

  paramSet _parameters;
  particle _p;

  // file input/output parameters
  vector<string> _xyzfileNames;
  xyzstream xyzfile;
  int _currentFile;

  // default grain density
  double _graindensity;

  // projection parameters
  bool _ecliptic;
  string _observerName;
  longlat _Coffset;
  projection _observer;
  longlat _cometRaDec;
  longlat _originRaDec;
  longlat _originOffset;

  // nucleus parameters
  float _npole[2];  // lambda, beta
  Vector nPole;     // unit vector

  // particle inclusion/exclusion parameters
  long _max;
  /// \todo Change the ageRange, etc, vectors to double* (speed up?)
  valarray<float> _ageRange, _betaRange, _latRange, _radRange, _sunRange;
  double _rhlimit;
  bool _ageInvert, _betaInvert, _latInvert, _radInvert, _sunInvert;
  Vej _vejGen;
};

#endif
