/***************************************************************************

  Convert a string to a scalar, vector, or valarray.

  Copyright (C) 2010 by Michael S. Kelley <msk@astro.umd.edu>

 ***************************************************************************/

#if HAVE_CONFIG_H
#include <config.h>
#endif /* HAVE_CONFIG_H */

#include <iostream>
#include <string>
#include <sstream>
#include <valarray>
#include <vector>
#include <cstdlib>
#include "StringConv.h"

using namespace std;

/** Setup string conversion with blank strings. */
StringConv::StringConv() {
  _s = "";
  _d = "";
}

/** Setup string conversion using the default delimiters. */
StringConv::StringConv(const string s) {
  _s = s;
  _d = ",/\\:;";
}

/** Setup string conversion using the specifed delimiters. */
StringConv::StringConv(const string s, const string d) {
  _s = s;
  _d = d;
}

/** Set the string delimiter. */
void StringConv::delimiter(string d) { _d = d; }

/** Set the string. */
void StringConv::str(string s) { _s = s; }

/** Set up the stringstream for parsing. */
void StringConv::setup() {
  string str = _s;
  int c;
  for (int i = 0; i<_d.size(); i++) {
    while ((c = str.find(_d[i])) != string::npos)
      str.replace(c, 1, " ");
  }
  _ss << str;
}

/** Parse the string and return a bool. */
bool StringConv::toBool() {
  string str;
  setup();
  _ss >> str;
  transform(str.begin(), str.end(), str.begin(), (int(*)(int))toupper);
  if (str == "YES") return true;
  if (str == "NO") return false;
  if (str == "TRUE") return true;
  if (str == "FALSE") return false;
  if (str == "ON") return true;
  if (str == "OFF") return false;
  if (atoi(str.c_str()) == 0) return false;
  if (atoi(str.c_str()) != 0) return true;
  cerr << "I did not understand " << _s << ".  Assuming 'false'";
  return false;
}

/** Parses the string and returns the first scalar. */
template<typename T>
T StringConv::toScalar() {
  T v;
  setup();
  _ss >> v;
  return v;
}

/** Parses the string and returns a value array. */
template<typename T>
vector<T> StringConv::toVector() {
  vector<T> v;
  T s;
  setup();
  while (!_ss.eof()) {
    _ss >> s;
    v.push_back(s);
  }
  return v;
}

/** Parses the string and returns a value array. */
template<typename T>
valarray<T> StringConv::toValarray() {
  vector<T> v = toVector<T>();
  valarray<T> va(v.size());
  for (int i=0; i<v.size(); i++) {
    va[i] = v[i];
  }
  return va;
}

template int StringConv::toScalar<int>();
template long StringConv::toScalar<long>();
template float StringConv::toScalar<float>();
template double StringConv::toScalar<double>();
template string StringConv::toScalar<string>();

template vector<int> StringConv::toVector<int>();
template vector<float> StringConv::toVector<float>();
template vector<double> StringConv::toVector<double>();
template vector<string> StringConv::toVector<string>();

template valarray<int> StringConv::toValarray<int>();
template valarray<long> StringConv::toValarray<long>();
template valarray<float> StringConv::toValarray<float>();
template valarray<double> StringConv::toValarray<double>();
