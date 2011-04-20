/***************************************************************************

  A class that wraps the dynamical integrators.

  Copyright (C) 2008-2010 by Michael S. Kelley <msk@astro.umd.edu>

 ***************************************************************************/

#if HAVE_CONFIG_H
#  include <config.h>
#endif /* HAVE_CONFIG_H */

#include <queue>
#include "Integrator.h"
#include "pfunctions.h"
#include "particle.h"
#include "paramSet.h"
#include "rundynamics.h"
#include "xyzstream.h"

using namespace std;

/** Initializes the integrator. */
Integrator::Integrator() {
  _setup = false;
}

/** Initializes and sets up the integrator given the parameter set. */
Integrator::Integrator(const paramSet p) {
  _setup = false;
  setup(p);
}

/** Initializes and sets up the integrator given the parameter set
    stream. */
Integrator::Integrator(istream& inf) {
  _setup = false;
  setup(inf);
}

/** Return the current parameter set. */
paramSet Integrator::parameters() { return _parameters; }

/** Set the current parameters to integrate by. */
void Integrator::parameters(paramSet p) {
  if (DEBUG) cerr << "Setting parameters (Integrator)\n";
  _parameters = p;
  if (DEBUG) cerr << "Finished setting parameters (Integrator)\n";
}

/** Set the current parameters from a stream. */
void Integrator::parameters(istream& inf) {
  if (DEBUG) cerr << "Setting parameters (Integrator)\n";
  _parameters.loadParameters(inf);
  if (DEBUG) cerr << "Finished setting parameters (Integrator)\n";
}

/** Return the current particle. */
particle Integrator::p() { return _particle; }

/** Set the current particle. */
void Integrator::p(particle p) { _particle = p; }

/** Set up a particle with this parameter set.  Returns false for no
    errors. */
bool Integrator::setup(const paramSet p) {
  parameters(p);
  return setup();
}

/** Set up a particle with this parameter set stream.  Returns false
    for no errors. */
bool Integrator::setup(istream& inf) {
  parameters(inf);
  return setup();
}

/** Set up a particle.  Returns false for no errors. */
bool Integrator::setup() {
  if (_setup) return false;

  double age_max; // for planet lookup

  // open the output file and write the header
  if (DEBUG) cerr << "Integrator::setup()\n";
  xyzfile.xyzopen(parameters().outFile(), xyzstream::WRITE);
  xyzfile.writeHeader(parameters());
  xyzfile.initData(parameters());

  // particle needs the parameter set in order to generate inital
  // state vectors
  _particle.parameters(parameters());

  if (parameters().isSyndynes() || parameters().isMakeComet()) {
    // if there is a particle generator function, use it to set up our
    // particle, but not if we are in Integrate XYZ mode.
    if (parameters().pFunc().length()) {
      pfunctions pfunc;
      bool err = pfunc.setup(parameters(), _particle);
      if (err) return err;
      age_max = _particle.ageDist.max();
    }
  }

  if (parameters().isSyndynes()) {
    if (DEBUG) cerr << "  - Syndyne age, beta generation\n";	

    // if we are generating syndynes, setup age and beta, overriding
    // any pfunc parameters
    queue<double> radius, age;

    for (int i=0; i<parameters().beta().size(); i++) {
      for (int j=0; j<parameters().steps(); j++) {
	radius.push(0.57 / parameters().beta()[i] / _particle.composition().bulkdensity);
	age.push(parameters().nDays() * static_cast<double>(j) / 
		 (parameters().steps() - 1) * 86400);
	if (DEBUG) cerr << "    + age, radius = " << age.back() / 86400 << " days, " <<
		     radius.back() << " microns\n";
      }
    }

    _particle.ageDist.distribution(Distribution::SEQUENCE);
    _particle.ageDist.setSequence(age);
    _particle.radiusDist.distribution(Distribution::SEQUENCE);
    _particle.radiusDist.setSequence(radius);
    _nParticles = parameters().beta().size() * parameters().steps();
    age_max = parameters().nDays() * 86400;
  } else if (parameters().isMakeComet()) {
    // we are generating a Monte Carlo simulation
    _nParticles = parameters().nParticles();
  }

  _currentParticle = 0;  // we have not generated any particles

  // Initalize the planet lookup table, if needed.
  if (parameters().planetLookUp() && parameters().planets()) {
    cout << "Creating the planet lookup table...";
    cout.flush();
    planet_lookup_init(jd2et(parameters().obsDate() - age_max / 86400),
		       age_max, parameters().planets());

    cout << " done." << endl;
  }

  _setup = true;
  if (DEBUG) cerr << "Completed Integrator::setup()\n";
  return false;  // no error
}

/** Get one particle and integrate it to t_final.  Returns 0 if it is
    the last particle, 1 if the particle is good, 2 if the particle
    has left the integration box, 3 if the particle was ejected
    outside the rhLimit.

    \todo Use IntegrationFlags for return values.

*/
int Integrator::calculateOne() {
  // generate a new particle
  if (DEBUG) cerr << "Integrator::calculateOne()\n";
  if (!parameters().isIntegrateXYZ()) {
    _particle.next(_currentParticle);
  }

  // integrate its position to t_final
  try {
    integrate();
  } catch (IntegrationFlags flag) {
    // there is something wrong with this particle, do not save it
    _currentParticle++;
    if (_currentParticle == _nParticles) {
      if (parameters().isSyndynes() && (parameters().orbit() > 0))
	writeOrbit(_currentParticle);
      return 0;
    }
    if (flag == outsideBox) return 2;
    if (flag == outsideRhLimit) return 3;
    return 0;  // should not reach this point
  }

  // we are done with this particle
  xyzfile.writeParticle(_particle);

  // are we finished with integrating?
  _currentParticle++;
  if (_currentParticle == _nParticles) {
    if (parameters().isSyndynes() && (parameters().orbit() > 0))
      writeOrbit(_currentParticle);
    return 0;
  }
  return 1;
}

/** Write the comet's orbit to the file.  Returns true on error. */
bool Integrator::writeOrbit(long int n) {
  double r[3], v[3];
  state st;

  _particle.beta(-99);
  for (int i=0; i<parameters().steps(); i++) {
    double ageDays = parameters().orbit() * (static_cast<double>(i) /
					 (parameters().steps() - 1) - 0.5);
    double jd = parameters().obsDate() - ageDays;

    if (get_comet_xyz(parameters().comet().c_str(),
		      parameters().spkKernel().c_str(), 1, &jd, r, v)) {
      _particle.error = true;  // failed get_comet_xyz
      return true;
    }

    st.r = r;
    st.v = v;
    st.t = -ageDays * 86400;

    _particle.istate(st);
    _particle.fstate(st);
    _particle.generateLabel(n + i);
    xyzfile.writeParticle(_particle);
  }
  return false;
}
