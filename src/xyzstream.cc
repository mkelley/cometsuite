/***************************************************************************

  Data file I/O for CometSuite.

  Copyright (C) 2008-2010 by Michael S. Kelley <msk@astro.umd.edu>

 ***************************************************************************/

#if HAVE_CONFIG_H
#  include <config.h>
#endif /* HAVE_CONFIG_H */

#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <algorithm>
#include <cstdlib>
#include "xyzstream.h"
#include "particle.h"
#include "paramSet.h"
#include "rundynamics.h"
#include "Vector.h"
#include "state.h"
#include "longlat.h"

using namespace std;

xyzstream::xyzstream() {
  _initData = false;
}

/** Open an xyz file with the default methods.

    \todo Consider automatically calling read/writeHeader() and
    initData() here.

    \todo Make a definitive description of the xyz files.
*/
void xyzstream::xyzopen(string filename, const unsigned int mode) {
  if (mode == WRITE) {
    if (DEBUG) cerr << "xyzstream::open() file for writing\n";
    open(filename.c_str(), ios::binary | ios::out | ios::trunc);
  } else {
    if (DEBUG) cerr << "xyzstream::open() file for reading\n";
    open(filename.c_str(), ios::binary | ios::in);
  }
}

/** Write a time stamp to the output file. */
void xyzstream::timeStamp() {
  time_t t;
  time(&t);
  *this << "# " << asctime(localtime(&t));
}

/** Write the parameter set and saved variables description to the
    output file. */
void xyzstream::writeHeader(paramSet parameters) {
  if (DEBUG) cerr << "xyzstream::writeHeader()\n";
  *this << "# " << PACKAGE_STRING << endl;
  timeStamp();
  parameters.writeParameters(*this);

  *this << "# data file description\n";
  *this << "UNITS: " << parameters.savedDataUnits() << "\n";
  *this << "DATA: " << parameters.savedDataDescription() << "\n";
  if (DEBUG) cerr << "xyzstream::writeXyzHeader() (savedData)\n";
}

/** Returns a string from the input stream's parameter set and saved
    variables description.  Only call readHeader() or readParamSet()
    once. */
string xyzstream::readHeader() {
  string header;

  while (!eof()) {
    char buf[255];
    getline(buf, 255);
    header += buf;
    header += "\n";
    // Data is always last
    if (header.find("DATA:") != string::npos) return header;
  }

  return header;
}

/** Reads in the file header and returns a parameter set. */
paramSet xyzstream::readParameters() {
  stringstream str;
  string header;
  paramSet parameters;

  header = readHeader();
  str << header;
  parameters.loadParameters(str);
  initData(parameters);
  cout << "\n";
  parameters.writeParameters(cout);
  cout << endl;
  return parameters;
}

/** Initialize the data description so that particles may be correctly
    read and written.

    \todo This should be automatic when reading from files.
*/
void xyzstream::initData(paramSet parameters) {
  int i, j;

  _dataDescription = parameters.savedData();
  transform(_dataDescription.begin(),
	    _dataDescription.end(),
	    _dataDescription.begin(),
	    (int(*)(int))toupper);
  if (DEBUG) cerr << "xyzstream::initData(): " << _dataDescription << "\n";

  string desc = parameters.savedDataDescription();
  i = desc.find_first_not_of(" ") - 1;
  _recSize = 0;

  while (i != string::npos) {
    i++;
    int rep = 1;
    if (desc[i+1] == '[')
      rep = atoi(desc.substr(i+2, desc.find_first_of("]")-1).c_str());

    if (desc[i] == 'd') {
      _recSize += sizeof(double) * rep;
    } else if (desc[i] == 'c') {
      _recSize += sizeof(char) * rep;
    } else if (desc[i] == 'i') {
      _recSize += sizeof(int) * rep;
    }
    i = desc.find_first_of(" ", i + 1);
  }

  i = 0;
  stringstream str;
  str << _dataDescription;
  while (!str.eof()) {
    string name;
    str >> name;

    if (name == "BETA") {
      _parameterOrder[i] = BETA;
    } else if (name == "RADIUS") {
      _parameterOrder[i] = RADIUS;
    } else if (name == "COMPOSITION") {
      _parameterOrder[i] = COMPOSITION;
    } else if ((name == "GRAINDENSITY") || (name == "DENSITY")) {
      _parameterOrder[i] = GRAINDENSITY;
    } else if (name == "MASS") {
      _parameterOrder[i] = MASS;
    } else if (name == "POROSITY") {
      _parameterOrder[i] = POROSITY;
    } else if (name == "FRACTALDIM") {
      _parameterOrder[i] = FRACTALDIM;
    } else if (name == "AGE") {
      _parameterOrder[i] = AGE;
    } else if (name == "ORIGIN") {
      _parameterOrder[i] = ORIGIN;
    } else if (name == "V_EJ") {
      _parameterOrder[i] = V_EJ;
    } else if (name == "R_I") {
      _parameterOrder[i] = R_I;
    } else if (name == "V_I") {
      _parameterOrder[i] = V_I;
    } else if (name == "T_I") {
      _parameterOrder[i] = T_I;
    } else if (name == "R_F") {
      _parameterOrder[i] = R_F;
    } else if (name == "V_F") {
      _parameterOrder[i] = V_F;
    } else if (name == "T_F") {
      _parameterOrder[i] = T_F;
    } else if (name == "LABEL") {
      _parameterOrder[i] = LABEL;
    }
    i++;
  }
  _nSavedParameters = i;

  i = _dataDescription.find_first_not_of(" ");
  j = _dataDescription.find_last_not_of(" ");
  if (_dataDescription.substr(i, j - i + 1) == "RADIUS GRAINDENSITY BETA AGE ORIGIN V_EJ R_I V_I T_I R_F V_F T_F LABEL") {
    _IOMethod = ALLPARAMETERS;
  } else if (_dataDescription.substr(i, j - i + 1) == "RADIUS GRAINDENSITY BETA AGE ORIGIN V_EJ R_I V_I T_I R_F V_F T_F") {
    _IOMethod = ALLBUTLABEL;
  } else if (_dataDescription.substr(i, j - i + 1) == "RADIUS GRAINDENSITY BETA AGE ORIGIN R_I V_EJ R_F") {
    _IOMethod = DEFAULT073;
  } else if (_dataDescription.substr(i, j - i + 1) == "BETA AGE ORIGIN R_I V_EJ R_F") {
    _IOMethod = DEFAULT050;
  } else if (_dataDescription.substr(i, j - i + 1) == "BETA RADIUS COMPOSITION DENSITY MASS POROSITY AGE ORIGIN V_EJ R_I V_I T_I R_F V_F T_F") {
    _IOMethod = ORIGINALFLAVOR;
  } else {
    _IOMethod = ONEATATIME;
  }

  _initData = true;
}

/** Write the current particle to the output file.

    /todo Can this be improved?
*/
void xyzstream::writeParticle(particle p) {
  int i, byte;
  byte = 0;
  if ((_IOMethod == ALLPARAMETERS) || (_IOMethod == ALLBUTLABEL)) {
    // the whole parameter list
    // RADIUS GRAINDENSITY BETA AGE ORIGIN V_EJ R_I V_I T_I R_F V_F T_F [LABEL]
    writeD(p.radius(), byte);
    writeD(p.graindensity(), byte);
    writeD(p.beta(), byte);
    writeD(p.age(), byte);
    writeLL(p.origin(), byte);
    writeVector(p.vej(), byte);
    writeState(p.istate(), byte);
    writeState(p.fstate(), byte);
    if (_IOMethod == ALLPARAMETERS)
      writeCArray(const_cast<char *>(p.label().c_str()), 16, byte);
  } else if (_IOMethod == DEFAULT073) {
    // the default for v0.7.3+
    // RADIUS GRAINDENSITY BETA AGE ORIGIN R_I V_EJ R_F
    writeD(p.radius(), byte);
    writeD(p.graindensity(), byte);
    writeD(p.beta(), byte);
    writeD(p.age(), byte);
    writeLL(p.origin(), byte);
    writeVector(p.istate().r, byte);
    writeVector(p.vej(), byte);
    writeVector(p.fstate().r, byte);
  } else if (_IOMethod == DEFAULT050) {
    // the default for v0.5.0 - v0.7.2
    // BETA AGE ORIGIN R_I V_EJ R_F
    writeD(p.beta(), byte);
    writeD(p.age(), byte);
    writeLL(p.origin(), byte);
    writeVector(p.istate().r, byte);
    writeVector(p.vej(), byte);
    writeVector(p.fstate().r, byte);
  } else {
    // could there be a better way?
    for (i=0; i<_nSavedParameters; i++) {
      switch (_parameterOrder[i]) {
      case BETA:
	writeD(p.beta(), byte);
	break;
      case RADIUS:
	writeD(p.radius(), byte);
	break;
      case GRAINDENSITY:
	writeD(p.graindensity(), byte);
	break;
      case AGE:
	writeD(p.age(), byte);
	break;
      case ORIGIN:
	writeLL(p.origin(), byte);
	break;
      case V_EJ:
	writeVector(p.vej(), byte);
	break;
      case R_I:
	writeVector(p.istate().r, byte);
	break;
      case V_I:
	writeVector(p.istate().v, byte);
	break;
      case T_I:
	writeD(p.istate().t, byte);
	break;
      case R_F:
	writeVector(p.fstate().r, byte);
	break;
      case V_F:
	writeVector(p.fstate().v, byte);
	break;
      case T_F:
	writeD(p.fstate().t, byte);
	break;
      case LABEL:
	writeCArray(const_cast<char *>(p.label().c_str()), 16, byte);
	break;
      }
    }
  }
  write(_record, _recSize);
  if (fail()) p.error = true;
}

/** Reads a particle from an input stream.  Note that this only sets
    the particle's Physical and Dynamical variables and does not
    intialize the particle's parameter set, which was presumably set
    before this function call.

    \todo Can this be improved?
*/
void xyzstream::readParticle(particle &p) {
  int i, byte;
  state st;

  read(_record, _recSize);
  if (fail()) p.error = true;
  byte = 0;

  // if the data string is one of the defaults, then go for it,
  // otherwise, parse the data description (slower)
  if (_IOMethod == ORIGINALFLAVOR) {
    // the old default, note that COMPOSITION, MASS, POROSITY, and
    // FRACTALDIM don't exist anymore, and that DENSITY is now
    // GRAINDENSITY
    // BETA RADIUS COMPOSITION DENSITY MASS POROSITY AGE ORIGIN V_EJ R_I V_I T_I R_F V_F T_F
    p.beta(readDArray(1, byte)[0]);
    p.radius(readDArray(1, byte)[0]);
    byte += sizeof(int);  // skip composition
    p.graindensity(readDArray(1, byte)[0]);
    byte += sizeof(double);  // skip mass
    byte += sizeof(double);  // skip porosity
    p.age(readDArray(1, byte)[0]);
    p.origin(readLL(byte));
    p.vej(readVector(byte));
    p.istate(readState(byte));
    p.fstate(readState(byte));
  } else if ((_IOMethod == ALLPARAMETERS) || (_IOMethod == ALLBUTLABEL)) {
    // all paramters, with or without label
    // RADIUS GRAINDENSITY BETA AGE ORIGIN V_EJ R_I V_I T_I R_F V_F T_F [LABEL]
    p.radius(readDArray(1, byte)[0]);
    p.graindensity(readDArray(1, byte)[0]);
    p.beta(readDArray(1, byte)[0]);
    p.age(readDArray(1, byte)[0]);
    p.origin(readLL(byte));
    p.vej(readVector(byte));
    p.istate(readState(byte));
    p.fstate(readState(byte));
    if (_IOMethod == ALLPARAMETERS)
      p.label(readCArray(16, byte));
  } else if (_IOMethod == DEFAULT073) {
    // the default for v0.7.3+
    // RADIUS GRAINDENSITY BETA AGE ORIGIN R_I V_EJ R_F
    p.radius(readDArray(1, byte)[0]);
    p.graindensity(readDArray(1, byte)[0]);
    p.beta(readDArray(1, byte)[0]);
    p.age(readDArray(1, byte)[0]);
    p.origin(readLL(byte));
    st = p.istate();
    st.r = readVector(byte);
    p.istate(st);
    p.vej(readVector(byte));
    st = p.fstate();
    st.r = readVector(byte);
    p.fstate(st);
  } else if (_IOMethod == DEFAULT050) {
    // the default for v0.5.0 - v0.7.2
    // BETA AGE ORIGIN R_I V_EJ R_F
    p.beta(readDArray(1, byte)[0]);
    p.age(readDArray(1, byte)[0]);
    p.origin(readLL(byte));
    st = p.istate();
    st.r = readVector(byte);
    p.istate(st);
    p.vej(readVector(byte));
    st = p.fstate();
    st.r = readVector(byte);
    p.fstate(st);
  } else {
    // could there be a better way?
    for (i=0; i<_nSavedParameters; i++) {
      switch (_parameterOrder[i]) {
      case BETA:
	p.beta(readDArray(1, byte)[0]);
	break;
      case RADIUS:
	p.radius(readDArray(1, byte)[0]);
	break;
      case COMPOSITION:
	byte += sizeof(int);  // skip composition
	break;
      case GRAINDENSITY:
	p.graindensity(readDArray(1, byte)[0]);
	break;
      case MASS:
	byte += sizeof(double);  // skip mass
	break;
      case POROSITY:
	byte += sizeof(double);  // skip porosity
	break;
      case FRACTALDIM:
	byte += sizeof(double);  // skip fractaldim
	break;
      case AGE:
	p.age(readDArray(1, byte)[0]);
	break;
      case ORIGIN:
	p.origin(readLL(byte));
	break;
      case V_EJ:
	p.vej(readVector(byte));
	break;
      case R_I:
	st = p.istate();
	st.r = readVector(byte);
	p.istate(st);
	break;
      case V_I:
	st = p.istate();
	st.v = readVector(byte);
	p.istate(st);
	break;
      case T_I:
	st = p.istate();
	st.t = readDArray(1, byte)[0];
	p.istate(st);
	break;
      case R_F:
	st = p.fstate();
	st.r = readVector(byte);
	p.fstate(st);
	break;
      case V_F:
	st = p.fstate();
	st.v = readVector(byte);
	p.fstate(st);
	break;
      case T_F:
	st = p.fstate();
	st.t = readDArray(1, byte)[0];
	p.fstate(st);
	break;
      case LABEL:
	p.label(readCArray(16, byte));
	break;
      }
    }
  }
}

/** Returns a character array from the input record at position i.  i
    is incremented to the end of the array. */
char* xyzstream::readCArray(int n, int &i) {
  char *a = reinterpret_cast<char *>(_record+i);
  i += sizeof(char) * n;
  return a;
}

/** Returns an integer array from the input record at position i.  i is
    incremented to the end of the array. */
int* xyzstream::readIArray(int n, int &i) {
  int *a = reinterpret_cast<int *>(_record+i);
  i += sizeof(int) * n;
  return a;
}

/** Returns a double array from the input record at position i.  i is
    incremented to the end of the array. */
double* xyzstream::readDArray(int n, int &i) {
  double *a = reinterpret_cast<double *>(_record+i);
  i += sizeof(double) * n;
  return a;
}

/** Returns a Vector from the input record at position i.  i is
    incremented to the end of the Vector. */
Vector xyzstream::readVector(int &i) {
  return Vector(readDArray(3, i));
}

/** Returns a state the input record at position i (r, v, then t).  i
    is incremented to the end of the state. */
state xyzstream::readState(int &i) {
  state st;
  st.r = readVector(i);
  st.v = readVector(i);
  st.t = readDArray(1, i)[0];
  return st;
}

/** Returns 2 values from the input record at position i as a longlat.
    i is incremented to the end of the longlat. */
longlat xyzstream::readLL(int &i) {
  longlat ll;
  double *a;
  a = readDArray(2, i);
  ll.lambda = a[0];
  ll.beta = a[1];
  return ll;
}

/** Writes a character array to the output record at position i.  i
    is incremented to the end of the array. */
void xyzstream::writeCArray(char *a, int n, int &i) {
  for (int j=0; j<n; j++)
    _record[i+j] = a[j];
  i += n;
}

/** Writes an integer array to the output record at position i.  i is
    incremented to the end of the array. */
void xyzstream::writeIArray(int *a, int n, int &i) {
  char *b = reinterpret_cast<char *>(a);
  for (int j=0; j<sizeof(int)*n; j++)
    _record[i+j] = b[j];
  i += sizeof(int) * n;
}

/** Writes a double array to the output record at position i.  i is
    incremented to the end of the array. */
void xyzstream::writeDArray(double *a, int n, int &i) {
  char *b = reinterpret_cast<char *>(a);
  for (int j=0; j<sizeof(double)*n; j++)
    _record[i+j] = b[j];
  i += sizeof(double) * n;
}

/** Writes a double value to the output record at position i.  i is
    incremented to by the size of 1 double. */
void xyzstream::writeD(double a, int &i) {
  double *p;
  p = &a;
  for (int j=0; j<sizeof(double); j++)
    _record[i+j] = reinterpret_cast<char*>(p)[j];
  i += sizeof(double);
}

/** Writes a Vector to the output record at position i.  i is
    incremented to the end of the Vector. */
void xyzstream::writeVector(Vector v, int &i) {
  writeDArray(v.dblarr(), 3, i);
}

/** Writes a state to the output record at position i (r, v, then t).
    i is incremented to the end of the state. */
void xyzstream::writeState(state st, int &i) {
  writeVector(st.r, i);
  writeVector(st.v, i);
  writeDArray(&st.t, 1, i);
}

/** Writes 2 values to the output record at position i as a longlat.
    i is incremented to the end of the longlat. */
void xyzstream::writeLL(longlat ll, int &i) {
  writeDArray(&ll.lambda, 1, i);
  writeDArray(&ll.beta, 1, i);
}

