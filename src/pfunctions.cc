/***************************************************************************

  Implements the particle generation functions.

  Copyright (C) 2005-2010 by Michael S. Kelley <msk@astro.umd.edu>

 ***************************************************************************/

#if HAVE_CONFIG_H
#  include <config.h>
#endif /* HAVE_CONFIG_H */

#include <iostream>
#include <cstdlib>
#include <ctime>
#include <sstream>
#include <cmath>
#include <fstream>
#include <vector>
#include "pfunctions.h"
#include "particle.h"
#include "Distribution.h"
#include "Vector.h"
#include "longlat.h"
#include "CoordTrans.h"
#include "rundynamics.h"
#include "StringConv.h"

using namespace std;

/** Initializes setup, which will keep track of whether or not the
    particle function parameter list has been parsed. */
pfunctions::pfunctions() {
  _setup = false;
}

/** Set up the given particle for particle generation.  Returns false
    if there are no errors.  The PFUNC string has the format: "pfunc
    parameters; pfunc parameters; ...", where pfunc is a module or
    template name, and parameters is the module or template parameter
    list.  The following pfuncs are allowed, or are on the to-do list.
    Logarithms are base 10, units are AU, km/s, days, and degrees.  In
    order to create beta = 0 grains (corresponding to radius = inf),
    use radius = -999:
      - Modules
        - RADIUS min max [steps]\n
	  The allowed range for grain radii.
        - LOGRADIUS log(min) log(max) [steps]\n
	  The allowed range for radii in terms of log10.
        - AGE min max [steps]\n
	  The allowed range of ages.
	- RADIUSLAW u1\n
	  v_ej is proportional to radius^-u1 [default: 0.5].
	- RHLAW u2\n
	  v_ej is proportional to rh^-u2 [default: 0.5].
        - VELOCITY ("iso"|"cos"|"temp") v0\n
	  The ejection velocity (actually speed) is proportional to v0
	  * radius^-u1 * rh^-u2, and is independent of the origin on
	  the surface (iso), proprotional to the cosine of the
	  sun-zenith angle (cos), or proportional to the surface
	  temperature (STM/NEATM) of a spherical nucleus (temp).
	- VELOCITY ("iso"|"cos"|"temp") "normal" v0 sigma [absolute]\n
	  Same as above, but here v is picked from a normal
	  distribution.  The default is centered on v' = v0 *
	  radius^-u1 * rh^-u2, with a width of sigma * v' (i.e., sigma
	  is specified as a fraction of v').  If the "absolute"
	  parameter is specified, the distribution is centered on v0
	  with a width of sigma * v0.  The velocity will be >= 0.
        - VELOCITY ("iso"|"cos"|"temp") "limit" v0min v0max [steps]\n
	  Same as above, but here v0 is picked from a range of values.
        - VELOCITY ("iso"|"cos"|"temp") "range" min max [steps]\n
	  Same as above, but here the ejection velocity is picked from
	  the absolute range min to max (independent of radius and rh).
        - Q_D ("iso"|"cos"|"temp")\n
	  The dust production is constant per solid angle (iso),
	  proportional to the cosine of the sun-zenith angle (cos), or
	  proportional to the surface temperature (STM/NEATM) of a
	  spherical nucleus (temp).
        - SUNCONE min max\n
	  The dust production is limited to a cone min to max degrees
	  from the sub-solar point.
        - RHLIMIT max\n
	  The dust production is limited to heliocentric distances
	  inside of max AU from the sun.
	- POLE nPoleLambda nPoleBeta\n
	  Set the ecliptic longitude and latitude of the nucleus north
	  pole.
	- LATITUDE thetaMin thetaMax
	  Set the ejection latitude to range from thetaMin to
	  thetaMax.
	- JET longitude latitude angle period\n
	  Eject grains in a jet.  Specify the planetocentric longitude
	  and latitude (degrees) of the jet at t=0 (the time of
	  observation), the opening angle, the rotation period of the
	  jet (i.e., nucleus) in hours.
        - COMPOSITION name\n
        - BULKDENSITY rho\n
        - FRACTALDIM D\n
      - Templates
        - BASIC_EJECTION_VELOCITY v0\n
	  Equivalent to "VELOCITY iso v0; SUNCONE 0 0".
        - ONE_JET v0 [minAge] maxAge log(radiusMin) log(radiusMax) nPoleLambda nPoleBeta jetLat alpha period\n
	  (not implemented)
        - ONE_SIMPLE_JET v0 [minAge] maxAge log(radiusMin) log(radiusMax) nPoleLambda nPoleBeta jetLat alpha period\n
	  (not implemented)
	- SIMPLE_COMA v0 [minAge] maxAge log(radiusMin) log(radiusMax)\n
	  Equivalent to "VELOCITY iso v0; AGE minAge maxAge; LOGRADIUS
	  log(radiusMin) log(radiusMax); Q_D iso; suncone 0 90".
	- ISO_COMA v0 [minAge] maxAge log(radiusMin) log(radiusMax)\n
	  Equivalent to "VELOCITY iso v0; AGE minAge maxAge; LOGRADIS
          log(radiusMin) log(radiusMax); Q_D iso; SUNCONE 0 180".
	- QV_COS_COMA v0 [minAge] maxAge log(radiusMin) log(radiusMax)\n
	  Equivalent to "VELOCITY cos v0; AGE minAge maxAge; LOGRADIUS
	  log(radiusMin) log(radiusMax); Q_D temp; SUNCONE 0 90".
	- QV_TEMP_COMA v0 [minAge] maxAge log(radiusMin) log(radiusMax)\n
	  Equivalent to "VELOCITY temp v0; AGE minAge maxAge; LOGRADIUS
	  log(radiusMin) log(radiusMax); Q_D temp; SUNCONE 0 90".
        - RTV log(radiusMin) log(radiusMax) minAge maxAge min_v0 max_v0\n
	  Equivalent to "LOGRADIUS log(radiusMin) log(radiusMax); AGE minAge
	  maxAge; VELOCITY iso min_v0 max_v0; Q_D iso; Q_D iso;
	  SUNCONE 0 180".
        - RTV_GRID steps log(radiusMin) log(radiusMax) minAge maxAge min_v0 max_v0\n
	  Equivalent to "LOGRADIUS log(radiusMin) log(radiusMax) steps; AGE
	  minAge maxAge steps; VELOCITY iso range min_v0 max_v0; Q_D
	  iso; SUNCONE 0 180".
	- RTVL log(radiusMin) log(radiusMax) minAge maxAge min_v0 max_v0 axisLam axisBet minLat maxLat
	  (not implemented)
	- DEEPIMPACT maxAge dt vmin vmax thetaMin thetaMax

    Particle fucntions are processed in the order they are found in
    the string, i.e., the last pfunc overrides any previously used
    pfuncs.  Logarithm functions are base 10.  Units are days, km/s,
    AU, and degrees.  Returns false for no error.

    \todo Use lex and yacc.
*/
bool pfunctions::setup(paramSet parameters, particle& p) {
  string pfunc, function, list;

  // Make sure the function list is terminated with a ;
  pfunc = parameters.pFunc() + ';';

  if (DEBUG) cerr << "pfunctions::setup() with " << pfunc << "\n";

  while (parse_parameters(pfunc, function, list)) {
    if (DEBUG) cerr << "  - " << function << " " << list << "\n";

    if (function == "radius") {
      radius(list, p);
      continue;
    }

    if (function == "logradius") {
      logradius(list, p);
      continue;
    }

    if (function == "age") {
      age(list, p);
      continue;
    }

    if (function == "rhlaw") {
      rhlaw(list, p);
      continue;
    }

    if (function == "radiuslaw") {
      radiuslaw(list, p);
      continue;
    }

    if (function == "velocity") {
      velocity(list, p);
      continue;
    }

    if (function == "q_d") {
      q_d(list, p);
      continue;
    }

    if (function == "rhlimit") {
      rhlimit(list, p);
      continue;
    }

    if (function == "suncone") {
      suncone(list, p);
      continue;
    }

    if (function == "pole") {
      pole(list, p);
      continue;
    }

    if (function == "latitude") {
      latitude(list, p);
      continue;
    }

    if (function == "jet") {
      jet(list, p);
      continue;
    }

    /// COMPOSITION name
    if (function == "composition") {
      composition(list, p);
      continue;
    }

    /// BULKDENSITY rho
    if (function == "bulkdensity") {
      bulkdensity(list, p);
      continue;
    }

    // FRACTALDIM D
    if (function == "fractaldim") {
      fractalDim(list, p);
      continue;
    }

    // Templates and old generators

    /* BASIC_EJECTION_VELOCITY v0
       = VELOCITY iso v0; SUNCONE 0 0
       Primarily intended for syndynes. */
    if (function == "basic_ejection_velocity") {
      // all grains are ejected directly towards the sun
      stringstream paramList;
      paramList << "iso " << list;
      velocity(paramList.str(), p);
      suncone("0 0", p);
      continue;
    }

    /** \todo Implement:
	ONE_JET v0 [minAge] maxAge log(radiusMin) log(radiusMax) nPoleLambda
	  nPoleBeta jetLat alpha period
    */

    /** \todo Implement:
        ONE_SIMPLE_JET v0 [minAge] maxAge log(radiusMin) log(radiusMax)
          nPoleLambda nPoleBeta jetLat alpha period
    */

    /* SIMPLE_COMA v0 [minAge] maxAge log(radiusMin) log(radiusMax)
       = VELOCITY iso v0; AGE minAge maxAge; LOGRADIUS log(radiusMin) log(radiusMax);
	 Q_D iso; suncone 0 90 */
    if (function == "simple_coma") {
      stringstream paramList;
      vector<string> l = StringConv(list).toVector<string>();

      paramList << "iso " << l[0];
      velocity(paramList.str(), p);
      q_d("iso", p);
      suncone("0 90", p);

      if (l.size() > 4) {
	paramList.str("");
	paramList << l[1] << " " << l[2];
	age(paramList.str(), p);

	paramList.str("");
	paramList << l[3] << " " << l[4];
	logradius(paramList.str(), p);
      } else {
	paramList.str("");
	paramList << "0 " << l[1];
	age(paramList.str(), p);

	paramList.str("");
	paramList << l[2] << " " << l[3];
	logradius(paramList.str(), p);
      }	
      continue;
    }

    /* ISO_COMA v0 [minAge] maxAge log(radiusMin) log(radiusMax)
       = VELOCITY iso v0;
         AGE minAge maxAge;
         LOGRADIUS log(radiusMin) log(radiusMax);
         Q_D iso; SUNCONE 0 180
    */
    if (function == "iso_coma") {
      stringstream paramList;
      vector<string> l = StringConv(list).toVector<string>();

      paramList << "iso " << l[0];
      velocity(paramList.str(), p);
      q_d("iso", p);
      suncone("0 180", p);

      if (l.size() > 4) {
	paramList.str("");
	paramList << l[1] << " " << l[2];
	age(paramList.str(), p);

	paramList.str("");
	paramList << l[3] << " " << l[4];
	logradius(paramList.str(), p);
      } else {
	paramList.str("");
	paramList << "0 " << l[1];
	age(paramList.str(), p);

	paramList.str("");
	paramList << l[2] << " " << l[3];
	logradius(paramList.str(), p);
      }	
      continue;
    }

    /* QV_COS_COMA v0 [minAge] maxAge log(radiusMin) log(radiusMax)
       = VELOCITY cos v0; AGE minAge maxAge; LOGRADIUS log(radiusMin) log(radiusMax);
	 Q_D temp; SUNCONE 0 90 */
    if (function == "qv_cos_coma") {
      stringstream paramList;
      vector<string> l = StringConv(list).toVector<string>();

      paramList << "cos " << l[0];
      velocity(paramList.str(), p);
      q_d("cos", p);
      suncone("0 90", p);

      if (l.size() > 4) {
	paramList.str("");
	paramList << l[1] << " " << l[2];
	age(paramList.str(), p);

	paramList.str("");
	paramList << l[3] << " " << l[4];
	logradius(paramList.str(), p);
      } else {
	paramList.str("");
	paramList << "0 " << l[1];
	age(paramList.str(), p);

	paramList.str("");
	paramList << l[2] << " " << l[3];
	logradius(paramList.str(), p);
      }	
      continue;
    }

    /* QV_TEMP_COMA v0 [minAge] maxAge log(radiusMin) log(radiusMax)
       = VELOCITY temp v0; AGE minAge maxAge; LOGRADIUS log(radiusMin) log(radiusMax);
	 Q_D temp; SUNCONE 0 90 */
    if (function == "qv_temp_coma") {
      stringstream paramList;
      vector<string> l = StringConv(list).toVector<string>();

      paramList << "temp " << l[0];
      velocity(paramList.str(), p);
      q_d("temp", p);
      suncone("0 90", p);

      if (l.size() > 4) {
	paramList.str("");
	paramList << l[1] << " " << l[2];
	age(paramList.str(), p);

	paramList.str("");
	paramList << l[3] << " " << l[4];
	logradius(paramList.str(), p);
      } else {
	paramList.str("");
	paramList << "0 " << l[1];
	age(paramList.str(), p);

	paramList.str("");
	paramList << l[2] << " " << l[3];
	logradius(paramList.str(), p);
      }	
      continue;
    }

    /* RTV log(radiusMin) log(radiusMax) minAge maxAge min_v0 max_v0
       = LOGRADIUS log(radiusMin) log(radiusMax); AGE minAge maxAge;
	 VELOCITY iso min_v0 max_v0; Q_D iso; Q_D iso; SUNCONE 0 180 */
    if (function == "rtv") {
      stringstream paramList;
      vector<string> l = StringConv(list).toVector<string>();

      paramList << l[0] << " " << l[1];
      logradius(paramList.str(), p);

      paramList.str("");
      paramList << l[2] << " " << l[3];
      age(paramList.str(), p);

      paramList.str("");
      paramList << "iso limit " << l[2] << " " << l[3];
      velocity(paramList.str(), p);
      q_d("iso", p);
      suncone("0 180", p);
      continue;
    }

    /* RTV_GRID steps log(radiusMin) log(radiusMax) minAge maxAge min_v0 max_v0
       = LOGRADIUS log(radiusMin) log(radiusMax) steps; AGE minAge maxAge steps;
	 VELOCITY iso range min_v0 max_v0; Q_D iso; SUNCONE 0 180 */
    if (function == "rtv_grid") {
      stringstream paramList;
      vector<string> l = StringConv(list).toVector<string>();

      paramList << l[1] << " " << l[2] << " " << l[0];
      logradius(paramList.str(), p);

      paramList.str("");
      paramList << l[3] << " " << l[4] << " " << l[0];
      age(paramList.str(), p);

      paramList.str("");
      paramList << "iso limit " << l[5] << " " << l[6] << " " << l[0];
      velocity(paramList.str(), p);
      q_d("iso", p);
      suncone("0 180", p);

      continue;
    }

    /* RTVL log(radiusMin) log(radiusMax) minAge maxAge min_v0 max_v0 axisLam axisBet minLat maxLat
       = cannot be reproduced with the modules (yet)
    if (function == "rtvl") {
      vector<double> input = toVector<double>(list);
      p.radiusDist.log10Min(input[0]);
      p.radiusDist.log10Max(input[1]);
      p.radiusDist.distribution(Distribution::LOG);

      p.ageDist.min(input[2]);
      p.ageDist.max(input[3]);
      p.ageDist.distribution(Distribution::LINEAR);

      p.vejGen.range();
      p.vejGen.vDist.min(input[4]);
      p.vejGen.vDist.max(input[5]);
      p.vejGen.vDist.distribution(Distribution::LINEAR);

      p.vejGen.axis(CoordTrans().longlatToVector(input[6] * M_PI / 180,
						 input[7] * M_PI / 180));
      
      p.vejGen.thetaDist.min(input[8] * M_PI / 180);
      p.vejGen.thetaDist.max(input[9] * M_PI / 180);
      p.vejGen.thetaDist.distribution(Distribution::ISOSPHERE);

      p.vejGen.phiDist.min(0);
      p.vejGen.phiDist.max(2 * M_PI);
      p.vejGen.phiDist.distribution(Distribution::LINEAR);

      continue;
    }
    */

    /* DEEPIMPACT maxAge dt vMin vMax thetaMin thetaMax
       = POLE 222.6 9.1;
         AGE (maxAge - dt / 3600) maxAge;
	 VELOCITY iso range vmin vmax;
	 LATITUDE thetaMin thetaMax;
	 Q_D iso;

	 maxAge in days
	 dt in seconds
	 vMin, vMax in km/s
	 thetaMin, thetaMax in degrees
     */
    if (function == "deepimpact") {
      stringstream paramList;
      vector<string> l = StringConv(list).toVector<string>();

      paramList << "222.6 9.1";
      pole(paramList.str(), p);

      paramList.str("");
      paramList << atof(l[0].c_str()) - atof(l[1].c_str()) / 3600 << 
	" " << l[0];
      age(paramList.str(), p);
      
      paramList.str("");
      paramList << "iso range " << l[2] << " " << l[3];
      velocity(paramList.str(), p);

      paramList.str("");
      paramList << l[4] << " " << l[5];
      latitude(paramList.str(), p);

      q_d("iso", p);
      continue;
    }

    cerr << "\nI don't know the particle function " << function << ".";
    _setup = false;
    return true;
  }

  _setup = true;
  if (DEBUG) cerr << "completed pfunctions::setup()\n";
  return false;
}

/** RADIUS min max [steps].  dn/d[log(radius)] proportional to 1. */
void pfunctions::radius(string list, particle &p) {
  vector<double> input = StringConv(list).toVector<double>();

  // the mandatory parameters: min max
  p.radiusDist.min(input[0]);
  p.radiusDist.max(input[1]);

  // the default distribution
  p.radiusDist.distribution(Distribution::LOG);

  if (input.size() == 3) {
    // we want gridded spacings
    p.radiusDist.nGridSteps(static_cast<int>(input[2]));
    p.radiusDist.distribution(Distribution::LOGGRID);
  }
}

/** LOGRADIUS min max [steps].  dn/d[log(radius)] proportional to 1. */
void pfunctions::logradius(string list, particle &p) {
  vector<double> input = StringConv(list).toVector<double>();

  // the mandatory parameters: min max
  p.radiusDist.log10Min(input[0]);
  p.radiusDist.log10Max(input[1]);

  // the default distribution
  p.radiusDist.distribution(Distribution::LOG);

  if (input.size() == 3) {
    // we want gridded spacings
    p.radiusDist.nGridSteps(static_cast<int>(input[2]));
    p.radiusDist.distribution(Distribution::LOGGRID);
  }
}

/** AGE min max [steps].  dn/d[age] proportional to 1.  min, max in
    days. */
void pfunctions::age(string list, particle &p) {
  vector<double> input = StringConv(list).toVector<double>();
  // the mandatory parameters: min max
  p.ageDist.min(input[0] * 86400);
  p.ageDist.max(input[1] * 86400);

  // the default distribution
  p.ageDist.distribution(Distribution::LINEAR);
	
  if (input.size() == 3) {
    // we want gridded spacings
    p.ageDist.nGridSteps(static_cast<int>(input[2]));
    p.ageDist.distribution(Distribution::GRID);
  }
}

/** RADIUSLAW u1.  v_ej proportional to radius^-u1. */
void pfunctions::radiuslaw(string list, particle &p) {
  p.vejGen.u1(atof(list.c_str()));
}

/** RHLAW u2.  v_ej proportional to rh^-u2. */
void pfunctions::rhlaw(string list, particle &p) {
  p.vejGen.u2(atof(list.c_str()));
}

/** Ejection velocity (well, actually speed):
    - VELOCITY ("iso"|"cos"|"temp") v0
    - VELOCITY ("iso"|"cos"|"temp") "normal" v0 sigma [absolute]
    - VELOCITY ("iso"|"cos"|"temp") "limit" v0min v0max [steps]
    - VELOCITY ("iso"|"cos"|"temp") "range" min max [steps] */
void pfunctions::velocity(string list, particle &p) {
  // get the first parameter as a string
  size_t found = list.find_first_of(" ");
  string type = list.substr(0, found);
  list.erase(0, found);

  // default is to assume "iso|cos|temp v0"
  if (type.find("iso") != string::npos) {
    p.vejGen.simpleActivity();
  } else if (type.find("cos") != string::npos) {
    p.vejGen.cosSimpleActivity();
  } else if (type.find("temp") != string::npos) {
    p.vejGen.tempSimpleActivity();
  } else {
    cerr << "VELOCITY parse error: Check your input type (" << type << ")\n";
    cerr << "                      Check your input list (" << list << ")\n";
  }

  // Is the word "range" next?
  bool range = false;
  found = list.find("range");
  if (found != string::npos) {
    range = true;
    list.erase(0, found+5);

    // use the range() functions
    if (type.find("iso") != string::npos) {
      p.vejGen.range();
    } else if (type.find("cos") != string::npos) {
      p.vejGen.cosRange();
    } else if (type.find("temp") != string::npos) {
      p.vejGen.tempRange();
    }
  }

  // Is the word "normal" next?
  bool normal = false;
  found = list.find("normal");
  if (found != string::npos) {
    normal = true;
    list.erase(0, found+6);

    // use the normal() functions
    if (type.find("iso") != string::npos) {
      p.vejGen.simpleActivityNormal();
    } else if (type.find("cos") != string::npos) {
      p.vejGen.cosSimpleActivityNormal();
    } else if (type.find("temp") != string::npos) {
      p.vejGen.tempSimpleActivityNormal();
    }
  }

  // Is the word "limit" next?
  bool limit = false;
  found = list.find("limit");
  if (found != string::npos) {
    limit = true;
    list.erase(0, found+5);

    // use the simpleActivityRange() functions
    if (type.find("iso") != string::npos) {
      p.vejGen.simpleActivityRange();
    } else if (type.find("cos") != string::npos) {
      p.vejGen.cosSimpleActivityRange();
    } else if (type.find("temp") != string::npos) {
      p.vejGen.tempSimpleActivityRange();
    }
  }

  // The remaining parameters are values
  vector<double> input = StringConv(list).toVector<double>();

  if (range) {
    p.vejGen.vDist.min(input[0]);
    p.vejGen.vDist.max(input[1]);
    if (input.size() == 3) {
      p.vejGen.vDist.distribution(Distribution::GRID);
      p.vejGen.vDist.nGridSteps(static_cast<int>(input[2]));
    } else {
      p.vejGen.vDist.distribution(Distribution::LINEAR);
    }
    return;
  }

  if (normal) {
    p.vejGen.v0(input[0]);
    p.vejGen.sigma(input[1]);
    p.vejGen.vDist.distribution(Distribution::NORMAL);
    return;
  }

  if (limit) {
    p.vejGen.v0min(input[0]);
    p.vejGen.v0(input[1]);
    if (input.size() == 3) {
      p.vejGen.vDist.distribution(Distribution::GRID);
      p.vejGen.vDist.nGridSteps(static_cast<int>(input[2]));
    } else {
      p.vejGen.vDist.distribution(Distribution::LINEAR);
    }
    return;
  }

  p.vejGen.vDist.distribution(Distribution::LINEAR);
  p.vejGen.v0(input[0]);
}

/** Q_D ("iso"|"cos"|"temp"). */
void pfunctions::q_d(string list, particle &p) {
  if (list.find("iso") != string::npos) {
    p.vejGen.thetaDist.distribution(Distribution::ISOSPHERE);
  } else if (list.find("cos") != string::npos) {
    p.vejGen.thetaDist.distribution(Distribution::COSSPHERE);
  } else if (list.find("temp") != string::npos) {
    p.vejGen.thetaDist.distribution(Distribution::TEMPSPHERE);
  } else {
    cerr << "Q_D parse error: Check your input string\n";
  }
  p.vejGen.phiDist.distribution(Distribution::LINEAR);
  return;
}

/** RHLIMIT max.  Dust will not be produced at rh > max. */
void pfunctions::rhlimit(string list, particle &p) {
  p.rhlimit(atof(list.c_str()) * _AU);
}

/** SUNCONE min max.  Dust will be ejected in a cone min and max
   degrees from the sun direction. */
void pfunctions::suncone(string list, particle &p) {
  vector<double> input = StringConv(list).toVector<double>();
  p.vejGen.thetaDist.min(input[0] * M_PI / 180);
  p.vejGen.thetaDist.max(input[1] * M_PI / 180);
  p.vejGen.phiDist.min(0);
  p.vejGen.phiDist.max(2 * M_PI);
}

/** POLE nPoleLambda nPoleBeta.  Set the nucleus north pole in
    ecliptic coordinates (degrees). */
void pfunctions::pole(string list, particle &p) {
  vector<double> input = StringConv(list).toVector<double>();
  longlat ll = {input[0], input[1]};
  p.pole(ll);
}

/** LATITUDE thetaMin thetaMax.  Set the ejection range distributed
    over a latitude range (degrees). */
void pfunctions::latitude(string list, particle &p) {
  vector<double> input = StringConv(list).toVector<double>();
  p.vejGen.thetaDist.min(input[0] * M_PI / 180);
  p.vejGen.thetaDist.max(input[1] * M_PI / 180);
}

/** JET longitude latitude angle period.  Eject grains out a
    precessing jet.  Requires POLE to be already set.  longitude and
    latitude are the direction of the jet in planetocentric
    coordinates (degrees) at t=0 (the time of observation); angle is
    the opening angle of the jet (degrees); and period is the period
    of the rotation (hours).  longitude = 0 is 90 +
    Pole_RightAscension from Right Ascension = 0, as defined by the
    IAU for small bodies.  For retrograde motion, set period < 0. */
void pfunctions::jet(string list, particle &p) {
  vector<double> input = StringConv(list).toVector<double>();
  double longitude = input[0];
  double latitude = input[1];
  double halfangle = input[2] / 2.0;
  p.vejGen.axis(p.pole());
  p.vejGen.thetaDist.min((90 - (latitude - halfangle)) * M_PI / 180);
  p.vejGen.thetaDist.max((90 - (latitude + halfangle)) * M_PI / 180);
  p.vejGen.phiDist.min((90 + p.pole().lambda + longitude - halfangle) * 
		       M_PI / 180 );
  p.vejGen.phiDist.max((90 + p.pole().lambda + longitude + halfangle) * 
		       M_PI / 180);
  p.vejGen.period(input[3] * 3600.0);  // Convert to seconds
  p.vejGen.jet();
}

/** COMPOSITION name.  Set the grain composition.  Grains composed of
    the default material, "geometric", respond to radiation pressure
    according to their physical cross section.  Available matarials
    (and their shorthand names):

      - geometric (g), 1 g/cm3
      - amorphous carbon (ac), 2.5 g/cm3
      - amorphous olivine 50 (ol50), 3.3 g/cm3, 50% Mg content

    If this list is updated, also update paramSet.cc, Physical.cc,
    Physical.h.
*/
void pfunctions::composition(string list, particle &p) {
  if ((list.find("geometric") != string::npos) || (list.compare("g") == 0)) {
    p.composition(GEOMETRIC);
    return;
  }
  if ((list.find("amorphous carbon") != string::npos) || (list.compare("ac") == 0)) {
    p.composition(AM_CARBON);
    return;
  }
  if ((list.find("amorphous olivine 50") != string::npos) || (list.compare("ol50") == 0)) {
    p.composition(AM_OLIVINE50);
    return;
  }
  /* No known material specified */
  cerr << "Unknown material (" << list << ") using default \"geometric\".\n";
  p.composition(GEOMETRIC);
}

/** BULKDENSITY rho.  Set the bulk material density [g/cm3].
    Currently, bulkdensity() only overrides the density of a geometric
    grain.  Overriding a real material's bulk density is not allowed.
*/
void pfunctions::bulkdensity(string list, particle &p) {
  if (p.composition().name == GEOMETRIC)
    p.composition().bulkdensity = atof(list.c_str());
}

/** FRACTALDIM D [a0].  The grain will be mixed with vacuum and have a
   fractally porous structure with a fractal dimension of D.  The
   vacuum inclusions have radii of 0.1 micron. */
void pfunctions::fractalDim(string list, particle &p) {
  p.fractalDim(atof(list.c_str()));
  p.minRadius(0.1);
}

/** Gets the absolute ranges in age and radius from a pfunction.  It
    should return UNLIMITED for parameters without limits.  Returns
    true on error. */
bool pfunctions::get_range(particle p, const int extrema_key,
			   double& extreme) {

  switch (extrema_key) {
  case RADIUS_MAX: {
    extreme = p.radiusDist.max();
    return false;
  }
  case RADIUS_MIN: {
    extreme = p.radiusDist.min();
    return false;
  }
  case AGE_MAX: {
    extreme = p.ageDist.max();
    return false;
  }
  case AGE_MIN: {
    extreme = p.ageDist.min();
    return false;
  }
  }

  return true;
}

/** Returns the pfunction name and parameters parsed from the
    parameter set.  Returns true if a new function is found. */
bool pfunctions::parse_parameters(string& pfunc, string& function,
				  string& list) {
  stringstream str;
  int i;
  char c[255];

  while (1) {
    // initialize str with all characters up to the next ";"
    i = pfunc.find_first_of(';');

    if (i == string::npos)
      return false;  // no function to extract

    str << pfunc.substr(0, i);
    pfunc.erase(0, i+1);  // remove everything up to and including ";"

    // get the function name and remove spaces
    function.clear();
    str >> function;
    while ((i = function.find(" ")) != string::npos)
      function.erase(i, 1);

    // we could have an empty command string (perhaps ";;")
    if (!function.empty())
      break;
  }

  // get the parameter list, if it exists
  str.getline(c, 255);
  if (str.fail())
    list = "";
  else
    list = c;

  // remove extra whitespace
  while ((i = list.find("  ")) != string::npos)
    list.erase(i, 1);
  while ((i = list.find_first_of(" ")) == 0)
    list.erase(i, 1);
  while ((i = list.find_last_of(" ")) == 0)
    list.erase(i, 1);

  return true;
}
