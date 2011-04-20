/***************************************************************************
  Copyright (C) 2007,2008,2010 by Michael S. Kelley <msk@astro.umd.edu>

 ***************************************************************************/

#if !defined(__VEJ)
#define __VEJ 1

#if HAVE_CONFIG_H
#  include <config.h>
#endif /* HAVE_CONFIG_H */

#include "CoordTrans.h"
#include "Distribution.h"
#include "Vector.h"
#include "longlat.h"
#include "state.h"
#include "rundynamics.h"

using namespace std;

/** Designed to handle all things related to ejection velocities. */
class Vej : CoordTrans {
public:
  Vej();

  // ejection distributions
  Distribution thetaDist, phiDist, vDist;

  // methods to set v_ej parameters
  void v0(const double);
  double v0();
  void v0min(const double);
  double v0min();
  void mu(const double);
  double mu();
  void sigma(const double);
  double sigma();
  void u1(const double);
  double u1();
  void u2(const double);
  double u2();
  void axis(const Vector);
  void axis(const longlat);
  Vector axis();
  void period(const double);
  double period();

  void range();
  void cosRange();
  void tempRange();

  void normal();
  void cosNormal();
  void tempNormal();

  void simpleActivity();
  void cosSimpleActivity();
  void tempSimpleActivity();

  void simpleActivityRange();
  void cosSimpleActivityRange();
  void tempSimpleActivityRange();

  void simpleActivityNormal();
  void cosSimpleActivityNormal();
  void tempSimpleActivityNormal();

  void sunward();
  void isoSphere();
  void cosSphere();
  void tempSphere();
  void jet();

  // methods to generate a new v_ej
  double getSpeed(const double, Vector);
  Vector getDirection(Vector, const double);
  Vector next(const double, Vector, const double);

private:
  Vector getPerp(Vector);
  Vector getSunward(Vector);

  enum Speed { RANGE, COSRANGE, TEMPRANGE,
	       NORMAL, COSNORMAL, TEMPNORMAL,
	       SIMPLEACTIVITY, COSSIMPLEACTIVITY, TEMPSIMPLEACTIVITY,
	       SIMPLEACTIVITYRANGE, COSSIMPLEACTIVITYRANGE,
	       TEMPSIMPLEACTIVITYRANGE,
	       SIMPLEACTIVITYNORMAL, COSSIMPLEACTIVITYNORMAL,
	       TEMPSIMPLEACTIVITYNORMAL,
  } _Speed;
  enum Direction { SUNWARD, ISOSPHERE, COSSPHERE, TEMPSPHERE, JET } _Dir;

  Vector _axis;
  double _v0, _v0min, _u1, _u2, _period, _angfreq;
};

#endif
