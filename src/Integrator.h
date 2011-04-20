/***************************************************************************
  Copyright (C) 2008 by Michael S. Kelley <msk@astro.umd.edu>

  ***************************************************************************/


#if !defined(__INTEGRATOR)
#define __INTEGRATOR 1

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <string>
#include "particle.h"
#include "paramSet.h"
#include "xyzstream.h"

using namespace std;

extern "C" int planet_lookup_init(double, double, int);

extern "C" int get_comet_xyz(const char*, const char*, int, double*, double*,
			     double*);

extern "C" double jd2et(double);

enum IntegrationFlags { outsideBox, outsideRhLimit };

/** The heart of rundynamics.  This class generates particles and
    integrates their positions.  However, no actual integration code
    is contained here---it is inherited from other classes.

    \todo Incorporate a Keplerian solution.
*/
class Integrator {
 public:
  // constructors
  Integrator();
  Integrator(const paramSet);
  Integrator(istream&);

  // methods
  paramSet parameters();
  void parameters(paramSet);
  void parameters(istream&);
  particle p();
  void p(particle);
  bool setup(const paramSet);
  bool setup(istream&);
  bool setup();
  int calculateOne();
  bool writeOrbit(long int);

  // integrate() will be replaced by the actual integration code, 
  // e.g., ra15::integrate()
  virtual void integrate() {};

  // vars
  string status;  // info on the last integration
  bool error;  // constructor error var

  // file input/ouput via
  xyzstream xyzfile;

  // nice friend
  friend class ra15;

 private:
  paramSet _parameters;
  particle _particle;
  long _nParticles, _currentParticle;
  bool _setup;
};

#endif
