/***************************************************************************
  Copyright (C) 2005 by Michael S. Kelley <msk@astro.umd.edu>

 ***************************************************************************/

#if !defined(__XYZHIST)
#define __XYZHIST 1

#if HAVE_CONFIG_H
#  include <config.h>
#endif /* HAVE_CONFIG_H */

#include <string>
#include <vector>
#include "xyzProject.h"
#include "particleWeights.h"

using namespace std;

enum binMethods { sizeBin, betaBin, ageBin };

/** Observer particle size/age/etc distributions on the sky. */
class xyzHist : public xyzProject, public particleWeights {
 public:
  xyzHist();

  // parameter setting
  string outfileName();
  void outfileName(const string);
  float aperture();
  void aperture(const float);
  int bin();
  void bin(const int);
  void bin(const string);
  float binSize();
  void binSize(const float);
  vector<float> dn_dx();
  void dn_dx(const vector<float>);
  vector<float> dn_dx_err();
  void dn_dx_err(const vector<float>);
  vector<float> xBin();
  void xBin(vector <float>);

  // processing methods
  void createPSD();
  bool writePSD();

 private:
  string _outfileName;

  // observed particle distributions and bin sizes, and bin locationss
  double _aperture;
  int _bin;
  float _binSize;
  vector<float> _dn_dx, _dn_dx_err, _xBin;
};

#endif
