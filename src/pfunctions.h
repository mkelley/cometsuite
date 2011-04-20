/***************************************************************************
  Copyright (C) 2005-2010 by Michael S. Kelley <msk@astro.umd.edu>

 ***************************************************************************/

#if !defined(__PFUNCTIONS)
#define __PFUNCTIONS 1

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <vector>
#include "paramSet.h"
#include "particle.h"
#include "rundynamics.h"

using namespace std;

/** \todo Change from a class to a set of functions? */
class pfunctions {
 public:
  // Constructors
  pfunctions();

  // Methods
  bool setup(paramSet, particle&);
  bool get_range(particle, const int, double&);

  void radius(string, particle&);
  void logradius(string, particle&);
  void age(string, particle&);
  void radiuslaw(string, particle&);
  void rhlaw(string, particle&);
  void velocity(string, particle&);
  void q_d(string, particle&);
  void rhlimit(string, particle&);
  void suncone(string, particle&);
  void pole(string, particle&);
  void latitude(string, particle&);
  void jet(string, particle&);
  void composition(string, particle&);
  void bulkdensity(string, particle&);
  void fractalDim(string, particle&);

 private:
  bool parse_parameters(string&, string&, string&);
  bool _setup;
};

#endif
