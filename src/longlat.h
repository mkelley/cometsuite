/***************************************************************************
  Copyright (C) 2008 Michael S. Kelley <msk@astro.umd.edu>
  
  ***************************************************************************/

#if !defined(__LONGLAT)
#define __LONGLAT 1

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

struct longlat {
  double lambda, beta; // longitude and latitude
};

#endif
