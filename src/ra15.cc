/***************************************************************************

  Interface between rundynamics and the RA15 integrator.

  Copyright (C) 2004-2007,2010 by Michael S. Kelley <msk@astro.umd.edu>

 ***************************************************************************/

#if HAVE_CONFIG_H
#  include <config.h>
#endif /* HAVE_CONFIG_H */

#include <sstream>
#include <cmath>
#include "ra15.h"
#include "Integrator.h"
#include "state.h"
#include "getxyz.h"

/** Advance a particle to t_final.  Throws outsideBox if the
    integration box is enabled and the particle has exceeded the
    integration boundary.

    \todo Allow backwards integrations for accuracy checking.
*/
void ra15::integrate() {
  double minStep, nFC, nTS, d2, et0, et;
  double box2;
  state final;
  stringstream str;

  box2 = parameters().box() * parameters().box();

  et0 = jd2et(parameters().obsDate());

  if (!_particle.parameters().isIntegrateXYZ()) {
    // Check the rh limit
    if ((_particle.istate().r.length() > _particle.rhlimit()) &&
	(_particle.rhlimit() > 0)) {
      // set the status string for run-time output
      str << _particle.beta() << " " << _particle.age() / 86400 << " outsideRhLimit\n";
      status = str.str();
      throw(outsideRhLimit);
    }

    final.r = _particle.istate().r;
    final.v = _particle.istate().v;
    final.t = 0;  // t=0 is the time of observation
  }

  // integrate to the micro-second
  if (abs(_particle.istate().t) > 1e-6) {
    double *r = final.r.dblarr();
    double *v = final.v.dblarr();

    if (parameters().box() > 0) {
      /* Check the comet-grain distance every age/10 seconds. */
      /* \todo Update to age/100 seconds? */
      int step = 0;
      double last_age = _particle.age();  // seconds
      double dAge = _particle.age() / 10;  // seconds
      double rc[3], vc[3];
      double nFCstep, nTSstep;

      nFC = nTS = 0;
      do {

	et = et0 - last_age;

	// take a step of dAge seconds
	everhart(r, v, dAge, et, _particle.beta(),
		 parameters().planets(), parameters().tolerance(),
		 &minStep, &nFCstep, &nTSstep,
		 parameters().planetLookUp(),
		 parameters().closeApproaches());
	final.r = r;
	final.v = v;
	last_age -= dAge;
	final.t = -last_age;

	nFC += nFCstep;
	nTS += nTSstep;
	step++;

	// Get the comet's position
	et = et0 - last_age;
	get_comet_xyz_et(parameters().comet().c_str(),
			 parameters().spkKernel().c_str(),
			 1, &et, rc, vc);

	// comet-grain distance^2
	d2 = (rc[0] - r[0]) * (rc[0] - r[0]) +
	  (rc[1] - r[1]) * (rc[1] - r[1]) +
	  (rc[2] - r[2]) * (rc[2] - r[2]);
      } while ((d2 <= box2) && (step < 10));
    } else {
      // take a step of -t seconds, where t is measured relative to
      // obsDate()
      et = et0 + _particle.istate().t;

      everhart(r, v, -_particle.istate().t, et,
      //      everhart(r, v, et0, et,
	       _particle.beta(), parameters().planets(),
	       parameters().tolerance(), &minStep, &nFC, &nTS,
	       parameters().planetLookUp(), parameters().closeApproaches());
    }
    final.r = r;
    final.v = v;
    final.t = 0;
  }

  _particle.fstate(final);

  // set the status string for run-time output
  str << _particle.radius() << " " << _particle.beta() << " " <<
    _particle.age() / 86400 << " " << minStep << " " << nFC;
  status = str.str();

  if ((parameters().box() > 0) && (d2 > box2)) {
    status += " outsideBox\n";
    throw(outsideBox);
  }

  status += "\n";
}
