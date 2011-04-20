/***************************************************************************
  Copyright (C) 2006 Michael S. Kelley <msk@astro.umd.edu>
  
  ***************************************************************************/

#include <math.h>

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

double round(double);

double round(double x) {
  double y;
  y = floor(x);
  if ((x - y) < 0.5) return y;
  return y+1;
}
