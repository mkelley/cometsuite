/***************************************************************************
  Copyright (C) 2008 by Michael S. Kelley <msk@astro.umd.edu>

 ***************************************************************************/

#if !defined(__STATE)
#define __STATE 1

#if HAVE_CONFIG_H
#  include <config.h>
#endif /* HAVE_CONFIG_H */

#include "Vector.h"

/** Units: km, km/s, seconds.  Time is usually relative to the
    observation date. */
struct state {
  Vector r;
  Vector v;
  double t;
};

#endif
