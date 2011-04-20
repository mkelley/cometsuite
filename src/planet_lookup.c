/***************************************************************************

  Creates and returns planet positions in a planet look up table.

  Copyright (C) 2004-2006,2008-2010 by Michael S. Kelley <msk@astro.umd.edu>

 ***************************************************************************/

#if HAVE_CONFIG_H
#  include <config.h>
#endif /* HAVE_CONFIG_H */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include "mskerr.h"
#include "rundynamics.h"

#define _DEBUG 0

/* variables to hold planet coordinates and dates */
static int plSetup = 0;
static float *planetxyz;
static double etStart = 0, etStep = 0;

int planet_lookup_init(double et, double nSeconds, int planets);
int planet_lookup(double et, int planets, double *r);

extern int get_planet_xyz(double et, int planets, double *r);

/******************************************************************************/
/** Initialize the planet lookup table.  Returns 1 on failure, 0 on
    success. */
int planet_lookup_init(double et, double nSeconds, int planets) {
  int j;
  long i;
  double now;
  double r[27];

  /* Init already completed, exit normally */
  if (plSetup != 0) return 0;

  planetxyz = (float *)calloc(_N_LOOKUP_DATES*9*3, sizeof(float));

  /* add 100 seconds of padding on either side of the look up table */
  etStart = et - 100;
  etStep = (fabs(nSeconds) + 200.0) / (_N_LOOKUP_DATES - 1);

  for (i=0; i<_N_LOOKUP_DATES; i++) {
    now = etStart + etStep * (double)i;

    if (get_planet_xyz_et(now, planets, r)) return 1;

    for (j=0; j<9; j++) {
	planetxyz[(i*9+j)*3  ] = r[j*3  ];
	planetxyz[(i*9+j)*3+1] = r[j*3+1];
	planetxyz[(i*9+j)*3+2] = r[j*3+2];
    }
  }

  plSetup = 1;

  return 0;
}

/******************************************************************************/
/** Get the planet positions via the lookup table.  Returns 1 on
   failure, 0 on success. */
int planet_lookup(double et, int planets, double *r) {
  int k;
  long i, j;
  double index;
  double wi, wj;
  double check[27];

  if (plSetup == 0) return 1;

  index = ((et - etStart) / etStep);
  /* index values */
  i = (long)floor(index);
  j = (long)ceil(index);

  /* we weight each index position by the fractional distance to the
     requested ephemeris time (linear interpolation) */
  wi = index - floor(index);
  wj = ceil(index) - index;

  if (i > _N_LOOKUP_DATES - 1 || j > _N_LOOKUP_DATES - 1) {
    fprintf(stderr, "Index too large, planet_lookup failed: et=%f etStart=%f etStep=%f i=%ld j=%ld\n", et, etStart, etStep, i, j);
    return 1;
  }
  if (i < 0 || j < 0) {
    fprintf(stderr, "Index too small, planet_lookup failed: et=%f etStart=%f etStep=%f i=%ld j=%ld\n", et, etStart, etStep, i, j);
    return 1;
  }

  if (i == j) {
    /* no interpolation necessary */
    for (k=0; k<9; k++) {
      r[k*3  ] = planetxyz[(i*9+k)*3  ];
      r[k*3+1] = planetxyz[(i*9+k)*3+1];
      r[k*3+2] = planetxyz[(i*9+k)*3+2];
    }
  } else {
    for (k=0; k<9; k++) {
      r[k*3  ] = wi * planetxyz[(i*9+k)*3  ] + wj * planetxyz[(j*9+k)*3  ];
      r[k*3+1] = wi * planetxyz[(i*9+k)*3+1] + wj * planetxyz[(j*9+k)*3+1];
      r[k*3+2] = wi * planetxyz[(i*9+k)*3+2] + wj * planetxyz[(j*9+k)*3+2];
    }
  }

  return 0;
}

