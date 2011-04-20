/***************************************************************************

  Convert a string to scalar, vector or valarray.

  Copyright (C) 2010 by Michael S. Kelley <msk@astro.umd.edu>

 ***************************************************************************/

#if !defined(__STRINGCONV)
#define __STRINGCONV 1

#if HAVE_CONFIG_H
#  include <config.h>
#endif /* HAVE_CONFIG_H */

#include <iostream>
#include <string>
#include <sstream>
#include <vector>
#include <valarray>
#include <cstdlib>

using namespace std;

class StringConv {
 public:
  StringConv();
  StringConv(const string);
  StringConv(const string, const string);

  void str(const string);
  void delimiter(const string);

  bool toBool();
  template<typename T> T toScalar();
  template<typename T> vector<T> toVector();
  template<typename T> valarray<T> toValarray();

 private:
  string _d, _s;
  stringstream _ss;

  void setup();
};

#endif
