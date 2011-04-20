/***************************************************************************
  Copyright (C) 2005-2010 by Michael S. Kelley <msk@astro.umd.edu>

 ***************************************************************************/

#if !defined(__XYZIMAGE)
#define __XYZIMAGE 1

#if HAVE_CONFIG_H
#  include <config.h>
#endif /* HAVE_CONFIG_H */

#include <string>
#include <valarray>
#include <CCfits/CCfits>
#include "xyzProject.h"
#include "particleWeights.h"

using namespace std;

enum xyzImageFlags { fitsWritingDisabled, cantCreateFits };

/** Creates images from xyzfiles. */
class xyzImage : public xyzProject, public particleWeights {
 public:
  xyzImage();

  // parameter setting
  bool writeFits();
  void writeFits(const bool);
  string outfileName();
  void outfileName(const string);
  valarray<long> size();
  void size(const long);
  void size(const long*);
  void size(const valarray<long>);
  valarray<float> platescale();
  void platescale(const float*);
  void platescale(const float);
  void platescale(const valarray<float>);

  // processing methods
  template <typename T>
    void addFitsKeyword(const string, const T, const string);
  void addFitsComment(const string);
  void createImages();
  void writeImages();

 private:
  bool _writeFits;
  string _outfileName;

  // image parameters
  long _size[2];
  float _platescale[2], _center[2];

  // image arrays
  valarray<double> _weighted;
  valarray<unsigned int> _number;
  valarray<float> _beta, _radius, _age, _latitude;

  // fits handling
  auto_ptr<CCfits::FITS> pFits;

  long index(const int, const int);
  void addWCS(CCfits::HDU*);
  void addXyzHeader();
};

/** Adds a keyword to the primary image header. */
template <typename T>
void xyzImage::addFitsKeyword(const string key, const T val, 
			      const string comment) {
  pFits->pHDU().addKey(key, val, comment);
}
#endif
