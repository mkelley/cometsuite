/***************************************************************************

  Copyright (C) 2008,2010 by Michael S. Kelley <msk@astro.umd.edu>

 ***************************************************************************/

#if !defined(__XYZSTREAM)
#define __XYZSTREAM 1

#if HAVE_CONFIG_H
#  include <config.h>
#endif /* HAVE_CONFIG_H */

#include <string>
#include <fstream>
#include "particle.h"
#include "paramSet.h"
#include "state.h"
#include "longlat.h"
#include "Vector.h"

using namespace std;

/** \todo gzip/gunzip data on the fly. */
class xyzstream : public fstream {
 public:
  xyzstream();

  // methods
  void xyzopen(string filename, const unsigned int);
  void timeStamp();
  void writeHeader(paramSet);
  string readHeader();
  paramSet readParameters();
  void initData(paramSet);
  void writeParticle(particle);
  void readParticle(particle&);

  // overloaded operators
  void operator=(xyzstream);

  // are we opening a file for reading or writing?
  enum openMethod { READ, WRITE };

  enum parameterNames { BETA, RADIUS, COMPOSITION, GRAINDENSITY, MASS,
			POROSITY, FRACTALDIM, AGE, ORIGIN, V_EJ, R_I,
			V_I, T_I, R_F, V_F, T_F, LABEL };

 private:
  bool _initData;
  int _recSize;
  string _dataDescription;
  // update this if particle descriptions get too big
  char _record[1024];
  // update this if there are ever more than 32 parameters
  int _parameterOrder[32];
  int _nSavedParameters, _IOMethod;

  enum particleIOMethod { ALLPARAMETERS, ALLBUTLABEL, DEFAULT073, DEFAULT050,
			  ORIGINALFLAVOR, ONEATATIME };

  char* readCArray(int, int&);
  int* readIArray(int, int&);
  double* readDArray(int, int&);
  Vector readVector(int&);
  state readState(int&);
  longlat readLL(int&);

  void writeCArray(char*, int, int&);
  void writeIArray(int*, int, int&);
  void writeDArray(double*, int, int&);
  void writeD(double, int&);
  void writeVector(Vector, int&);
  void writeState(state, int&);
  void writeLL(longlat, int&);
};

#endif
