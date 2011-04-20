/***************************************************************************
  Copyright (C) 2004-2006,2008-2010 by Michael S. Kelley
  <msk@astro.umd.edu>

 ***************************************************************************/

#if !defined(__PARTICLE)
#define __PARTICLE 1

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <iostream>
#include <string>
#include "Vector.h"
#include "Physical.h"
#include "Dynamical.h"
#include "Distribution.h"
#include "paramSet.h"

using namespace std;

// comet, kernel, n points, julian dates, r, v
extern "C" int get_comet_xyz(const char*, const char*, int, double*, double*,
			     double*);
extern "C" double jd2et(double);

/** Everything we need to know about a comet particle.  Creating new
    instances of this class uses a significant amount of overhead in
    the rundynamics programs.  Reuse particle instances as much as
    possible. */
class particle : public Physical, public Dynamical {
 public:
  // constructors
  particle();
  particle(const particle&);

  // methods
  paramSet parameters();
  void parameters(paramSet);
  string label();
  void label(const string);
  void label(char*);
  void generateLabel(const long n);
  void next();
  void next(const long n);
  double temperature();

  // overloaded operators
  friend ostream& operator<<(ostream&, particle);

  // vars
  bool error;

 private:
  paramSet _parameters;
  string _label;
};

#endif
