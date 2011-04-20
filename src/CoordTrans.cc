/***************************************************************************

  Transformations between coordinate systems.

  Copyright (C) 2005,2008 by Michael S. Kelley <msk@astro.umd.edu>

 ***************************************************************************/

#if HAVE_CONFIG_H
#  include <config.h>
#endif /* HAVE_CONFIG_H */

#include <fstream>
#include <cmath>
#include "CoordTrans.h"
#include "Vector.h"
#include "longlat.h"

using namespace std;

/** An empty constructor. */
CoordTrans::CoordTrans() {
  return;
}

/** Transform heliocentric rectangular coordinates to ecliptic
    coordinates with respect to the observer. */
longlat CoordTrans::getEcliptic(const Vector RHobs, const Vector RHtarget) {
  Vector ROtarget = RHtarget;
  ROtarget = ROtarget - RHobs;

  longlat c;
  c.lambda = atan2(ROtarget[1], ROtarget[0]) * 180.0 / M_PI;
  c.beta = atan2(ROtarget[2], sqrt(ROtarget[0]*ROtarget[0] +
				   ROtarget[1]*ROtarget[1])) * 180.0 / M_PI;
  return c;
}

/** Transform ecliptic coordinates to RA-Dec in equinox J2000.
    Adapted from euler.pro W. Landsman and Daryl Yentis.  Obliquity of
    the ecliptic eps = 23.4392911111 degrees (Astronomical Almanac
    2008)
*/
longlat CoordTrans::eclipticToEqJ2000(const longlat eclip) {
  long double a = eclip.lambda * M_PI / 180.0;
  long double b = eclip.beta * M_PI / 180.0;
  long double c = 0.39777715593 * cos(b) * sin(a) + 0.91748206207 * sin(b);
  long double d = atan2(0.91748206207 * cos(b) * sin(a) - 0.39777715593 * sin(b),
		   cos(b) * cos(a));
  longlat eq;
  eq.beta = static_cast<double>(asin((c <= 1.0)?c:1.0) * 180.0 / M_PI);
  eq.lambda = static_cast<double>(fmod(d + 4.0 * M_PI, 
	        static_cast<long double>(2.0 * M_PI)) * 180.0 / M_PI);

  return eq;
}

/** Converts a set of angles (degrees) to a Vector. */
Vector CoordTrans::longlatToVector(const double l, const double b) {
  longlat angles = {l, b};
  return longlatToVector(angles);
}

/** Converts a set of angles (degrees) to a Vector. */
Vector CoordTrans::longlatToVector(const longlat a) {
  Vector axis;
  axis[0] = cos(a.beta * M_PI / 180) * cos(a.lambda * M_PI / 180);
  axis[1] = cos(a.beta * M_PI / 180) * sin(a.lambda * M_PI / 180);
  axis[2] = sin(a.beta * M_PI / 180);
  return axis.unit();
}
