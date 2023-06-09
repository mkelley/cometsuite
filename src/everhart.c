/***************************************************************************

  Transcoded from the RADU15.F integrator (Everhart 1985. in Dynamics
  of Comets: Their Origin and Evolution. A. Carusi and G. Valsecchi,
  eds. Astrophysics and Space Science Library 115 185).

  C version and enhancements: Copyright (C) 2005,2006,2007,2008,2009
  by Michael S Kelley <msk@astro.umd.edu>

 ***************************************************************************/

#if HAVE_CONFIG_H
#  include <config.h>
#endif /* HAVE_CONFIG_H */

#include <stdio.h>
#include <math.h>
#include <string.h>
#include "mskerr.h"
#include "rundynamics.h"

#define true 1
#define false 0

#define _DEBUG 0

/* number of differential equations */
#define _NDE 3

/* Default initial step size, 1 day is good, more than 10 days may
   affect the accuracy. */
#define _ISS 86400.0

/* Sun / planetsary system mass ratios (ssd.jpl.nasa.gov) DE405 */
/*  static const double massRatio[] = {   6023600.0,      /* (+/- 250.) merc */
/*  				       408523.71,     /* (+/- 0.06) ven */
/*  				       328900.5614,   /* (+/- 0.02) ear+moon */
/*  				      3098708.0,      /* (+/- 9.)   mars sys */
/*  				         1047.3486,   /* (+/- 0.0008) jup sys */
/*  				         3497.898,    /* (+/- 0.018)  sat sys */
/*  				        22902.98,     /* (+/- 0.03)   ura sys */
/*  				        19412.24,     /* (+/- 0.04)   nep sys */
/*  				    135200000.0 };    /* (+/- 0.07d8) plu sys */

/* GM for the planets via DE405 [km^3/s^2] (converted from au^3/day^2) */
static const double MUplanets[] = {     22032.080486417923,  /* merc */
				       324858.598826459784,  /* ven */
				       403503.233479087008,  /* ear+moon */
				        42828.314258067236,  /* mars sys */
				    126712767.857795968652,  /* jup sys */
				     37940626.061137281358,  /* sat sys */
				      5794549.007071875036,  /* ura sys */
				      6836534.063971338794,  /* nep sys */
					  981.600887707004 };/* plu sys */

/* Planet-comet distance required to be considered a close approach */
static const double CAdistance[] = { 0.1 * _AU, 0.1 * _AU, 0.1 * _AU,
				     0.1 * _AU, 1.0 * _AU, 1.0 * _AU,
				     1.0 * _AU, 1.0 * _AU, 0.1 * _AU };
/* The semi-major axis of each planet */
static const double CAplanet_a[] = { 0.387 * _AU, 0.723 * _AU,
				     1.000 * _AU, 1.523 * _AU,
				     5.203 * _AU, 9.537 * _AU,
				     19.191 * _AU, 30.068 * _AU,
				     39.481 * _AU };

/* h = t/tStep; fractional step sizes in terms of the total integration
   step size (Gaussian-Radau spacings sacled to the range [0,1] for
   integrating to order 15); the sum should be 3.7333333333333333 */
static const double h[] = {0.0, 0.05626256053692215, 0.18024069173689236,
			   0.35262471711316964, 0.54715362633055538,
			   0.73421017721541053, 0.88532094683909577,
			   0.97752061356128750};
static const double sr = 1.4;
static const double pw = 1.0 / 9.0;

int everhart(double *x, double *v, double tFinal, double et, double beta,
	     const int planets, double tol, double *minStep, double *DnFC,
	     double *DnSeq, const int planetLookUp, const int closeApproaches);
void initConstants(double *xc, double *vc, double *c, double *d, double *r);
int calcAccel(double *x, double *v, double t, double et, double beta,
	      const int planets, double *a, const int planetLookUp,
	      const int closeApproaches);

/******************************************************************************/
/** Integrate the position of a particle under solar radiation and
    gravity forces, and planetary perturbations.

    \todo Add relativistic approximations.
*/
int everhart(double *x, double *v, double tFinal, double et, double beta,
	     const int planets, double tol, double *minStep, double *DnFC,
	     double *DnSeq, const int planetLookUp, const int closeApproaches) {
  static double c[21], d[21], r[21];
  static double xc[8], vc[7];
  double a1[_NDE], aj[_NDE], y[_NDE], z[_NDE];
  double b[7][_NDE], g[7][_NDE], e[7][_NDE], bd[7][_NDE];
  double s[9];
  double tPrime, tStep, tTotal;
  double xstep, dir, temp, gk, hv;
  long nFC, nSeq;
  static int initDone = false;
  int first_sequence, final_sequence, constant_step, redo_first_sequence;
  int i, j, k, l, m, nIter, nCount;

  if (_DEBUG) {
    _ERR("Input Variables:\n");
    _V_ERR("   x", x); _NL;
    _V_ERR("   v", v); _NL;
    _E_ERR("   tFinal", tFinal); _E_ERR(" et", et); _E_ERR(" beta", beta); _NL;
    _I_ERR("   planets", planets); _E_ERR(" tol", tol); _NL; _NL;
  }

  xstep = 1e3;

  /* zero some arrays */
  memset(b, 0, sizeof(double)*7*_NDE);
  memset(bd, 0, sizeof(double)*7*_NDE);
  memset(g, 0, sizeof(double)*7*_NDE);
  memset(e, 0, sizeof(double)*7*_NDE);

  if (!initDone) {
    initConstants(xc, vc, c, d, r);
    initDone = true;
  }

  /* Now that the constants are initialized, make an estimate of tPrime */
  dir = (tFinal<0.0)?-1.0:1.0;
  if (tol < 0.0) {
    constant_step = true;
    tPrime = xstep;
  } else {
    if (tol == 0.0) tol = 1e-8;
    constant_step = false;
    /*    tPrime = 0.1 * dir; */
    tPrime = _ISS * dir;
  }

  if ((tPrime / tFinal) > 0.5) tPrime = 0.5 * tFinal;

  nCount = 0;
  first_sequence = true;
  final_sequence = false;
  if (_DEBUG) printf("%15s %9s %24s %24s %12s %12s %24s\n",
		     "Function calls", "Sequences", "x[0]",
		     "v[0]", "tStep", "tTotal", "tFinal");

  do { /* while(1) */
    do { /* while(redo_first_sequence) */
      if (first_sequence) {
	nSeq = 0;
	nIter = 6;
	tTotal = 0.0;
	*minStep = 0.0;
	if (calcAccel(x, v, 0.0, et, beta, planets, a1, planetLookUp,
		      closeApproaches))
	  return 1;
	nFC = 1;
      }

      tStep = tPrime;
      if (nSeq % 1000 == 0) {
	if (_DEBUG) printf("%15ld %9ld %24.16e %24.16e %12e %12e %24.16e\n",
			   nFC, nSeq, x[0], v[0], tStep, tTotal, tFinal);
      }

      /* Find new g values from predicted b values */
      for (k=0; k<_NDE; k++) {
	g[0][k] = d[15]*b[6][k] + d[10]*b[5][k] + d[6]*b[4][k] + d[3]*b[3][k] + d[1]*b[2][k] 
	  + d[0]*b[1][k] + b[0][k];
	g[1][k] = d[16]*b[6][k] + d[11]*b[5][k] + d[7]*b[4][k] + d[4]*b[3][k] + d[2]*b[2][k] +
	  b[1][k];
	g[2][k] = d[17]*b[6][k] + d[12]*b[5][k] + d[8]*b[4][k] + d[5]*b[3][k] +      b[2][k];
	g[3][k] = d[18]*b[6][k] + d[13]*b[5][k] + d[9]*b[4][k] +      b[3][k];
	g[4][k] = d[19]*b[6][k] + d[14]*b[5][k] +      b[4][k];
	g[5][k] = d[20]*b[6][k] +       b[5][k];
	g[6][k] =       b[6][k];
      }

      for (m=1; m<=nIter; m++) {
	for (j=1; j<8; j++) {
	  s[0] = tStep * h[j];
	  s[1] = s[0] * s[0];
	  s[2] = s[1] * h[j];         s[1] = s[1] * xc[0];
	  s[3] = s[2] * h[j];         s[2] = s[2] * xc[1];
	  s[4] = s[3] * h[j];         s[3] = s[3] * xc[2];
	  s[5] = s[4] * h[j];         s[4] = s[4] * xc[3];
	  s[6] = s[5] * h[j];         s[5] = s[5] * xc[4];
	  s[7] = s[6] * h[j];         s[6] = s[6] * xc[5];
	  s[8] = s[7] * h[j] * xc[7]; s[7] = s[7] * xc[6];

	  for (k=0; k<_NDE; k++) {
	    y[k] = x[k] + 
	         v[k] * s[0] +
	        a1[k] * s[1] +
	      b[0][k] * s[2] +
	      b[1][k] * s[3] +
	      b[2][k] * s[4] +
	      b[3][k] * s[5] +
	      b[4][k] * s[6] +
	      b[5][k] * s[7] +
	      b[6][k] * s[8];
	  }

	  s[0] = tStep * h[j];
	  s[1] = s[0] * h[j];
	  s[2] = s[1] * h[j];         s[1] = s[1] * vc[0];
	  s[3] = s[2] * h[j];         s[2] = s[2] * vc[1];
	  s[4] = s[3] * h[j];         s[3] = s[3] * vc[2];
	  s[5] = s[4] * h[j];         s[4] = s[4] * vc[3];
	  s[6] = s[5] * h[j];         s[5] = s[5] * vc[4];
	  s[7] = s[6] * h[j] * vc[6]; s[6] = s[6] * vc[5];
	  
	  for (k=0; k<_NDE; k++) {
	    z[k] = v[k] + 
	        a1[k] * s[0] +
	      b[0][k] * s[1] +
	      b[1][k] * s[2] +
	      b[2][k] * s[3] +
	      b[3][k] * s[4] +
	      b[4][k] * s[5] +
	      b[5][k] * s[6] +
	      b[6][k] * s[7];
	  }	  

	  if (calcAccel(y, z, tTotal + h[j] * tStep, et, beta, planets, aj,
			planetLookUp, closeApproaches)) return 1;
	  nFC++;

	  for (k=0; k<_NDE; k++) {
	    temp = g[j-1][k];
	    gk = (aj[k] - a1[k]) / h[j];

	    switch(j) {
	    case 1 : g[0][k] =       gk; break;
	    case 2 : g[1][k] =      (gk-g[0][k])*r[ 0]; break;
	    case 3 : g[2][k] =     ((gk-g[0][k])*r[ 1]-g[1][k])*r[ 2]; break;
	    case 4 : g[3][k] =    (((gk-g[0][k])*r[ 3]-g[1][k])*r[ 4]-g[2][k])*r[ 5]; break;
	    case 5 : g[4][k] =   ((((gk-g[0][k])*r[ 6]-g[1][k])*r[ 7]-g[2][k])*r[ 8]-g[3][k])*
		                 r[ 9]; break;
	    case 6 : g[5][k] =  (((((gk-g[0][k])*r[10]-g[1][k])*r[11]-g[2][k])*r[12]-g[3][k])*
				 r[13]-g[4][k])*r[14]; break;
	    case 7 : g[6][k] = ((((((gk-g[0][k])*r[15]-g[1][k])*r[16]-g[2][k])*r[17]-g[3][k])*
				 r[18]-g[4][k])*r[19]-g[5][k])*r[20]; break;
	    }

	    temp = g[j-1][k] - temp;
	    b[j-1][k] = b[j-1][k] + temp;

	    switch(j) {
	    case 2 : b[0][k] = b[0][k] + c[0] * temp;
	             break;
	    case 3 : b[0][k] = b[0][k] + c[1] * temp;
	             b[1][k] = b[1][k] + c[2] * temp;
		     break;
	    case 4 : b[0][k] = b[0][k] + c[3] * temp;
	             b[1][k] = b[1][k] + c[4] * temp;
	             b[2][k] = b[2][k] + c[5] * temp;
		     break;
	    case 5 : b[0][k] = b[0][k] + c[6] * temp;
	             b[1][k] = b[1][k] + c[7] * temp;
	             b[2][k] = b[2][k] + c[8] * temp;
	             b[3][k] = b[3][k] + c[9] * temp;
		     break;
	    case 6 : b[0][k] = b[0][k] + c[10] * temp;
	             b[1][k] = b[1][k] + c[11] * temp;
	             b[2][k] = b[2][k] + c[12] * temp;
	             b[3][k] = b[3][k] + c[13] * temp;
	             b[4][k] = b[4][k] + c[14] * temp;
		     break;
	    case 7 : b[0][k] = b[0][k] + c[15] * temp;
	             b[1][k] = b[1][k] + c[16] * temp;
	             b[2][k] = b[2][k] + c[17] * temp;
	             b[3][k] = b[3][k] + c[18] * temp;
	             b[4][k] = b[4][k] + c[19] * temp;
	             b[5][k] = b[5][k] + c[20] * temp;
		     break;
	    }
	  } /* End k loop */
	} /* End j loop */
      } /* End nIter loop */

      if (!constant_step) {
	/* Sequence size control */
	hv = 0.0;
	for (k=0; k<_NDE; k++) {
	  hv = (hv>fabs(b[6][k]))?hv:fabs(b[6][k]);
	}
	hv *= xc[7] / pow(fabs(tStep), 7);
      }

      redo_first_sequence = false;
      if (first_sequence) {
	if (constant_step) {
	  tPrime = xstep;
	} else {
	  tPrime = dir * pow(tol / hv, pw);
	  if ((tPrime / tStep) <= 1.0) {
	    /* If the new tPrime is smaller than the last tStep */
	    /* then restart with: */
	    tPrime = 0.8 * tPrime;
	    nCount++;
	    if (nCount > 1) if (_DEBUG) printf("%15d %24.16e %24.16e\n", nCount,tStep,tPrime);

	    if (nCount > 10) {
	      fprintf(stderr, "Too many time refinements!\n");
	      return 1;
	    }
	  }
	}
	if (!redo_first_sequence) first_sequence = false;
      }
    } while(redo_first_sequence);

    /* Find the new x and v values */
    s[0] = tStep;
    s[1] = s[0] * s[0];
    s[2] = s[1] * xc[1];
    s[3] = s[1] * xc[2];
    s[4] = s[1] * xc[3];
    s[5] = s[1] * xc[4];
    s[6] = s[1] * xc[5];
    s[7] = s[1] * xc[6];
    s[8] = s[1] * xc[7];
    s[1] = s[1] * xc[0];

    for (k=0; k<_NDE; k++) {
      x[k] = x[k] + 
	   v[k] * s[0] +
	  a1[k] * s[1] +
	b[0][k] * s[2] +
	b[1][k] * s[3] +
	b[2][k] * s[4] +
	b[3][k] * s[5] +
	b[4][k] * s[6] +
	b[5][k] * s[7] +
	b[6][k] * s[8];
    }
      
    s[0] = tStep;
    s[1] = s[0] * vc[0];
    s[2] = s[0] * vc[1];
    s[3] = s[0] * vc[2];
    s[4] = s[0] * vc[3];
    s[5] = s[0] * vc[4];
    s[6] = s[0] * vc[5];
    s[7] = s[0] * vc[6];
  
    for (k=0; k<_NDE; k++) {
      v[k] = v[k] + 
	a1[k] * s[0] +
	b[0][k] * s[1] +
	b[1][k] * s[2] +
	b[2][k] * s[3] +
	b[3][k] * s[4] +
	b[4][k] * s[5] +
	b[5][k] * s[6] +
	b[6][k] * s[7];
    }

    tTotal += tStep;
    if ((*minStep > tStep) || (*minStep == 0.0)) *minStep = tStep;
    nSeq++;

    /* If we are done, then return */
    if (final_sequence) {
      if (_DEBUG) printf("%15ld %9ld %24.16e %24.16e %12e %12e %24.16e\n",
			 nFC, nSeq, x[0], v[0], tStep, tTotal, tFinal);
      *DnFC = (double)nFC;
      *DnSeq = (double)nSeq;
      return 0;
    }

    /* Control the size of the next sequence and adjust the last sequence */
    /* to exactly cover the integration span. */
    if (constant_step) {
      tPrime = xstep;
    } else {
      tPrime = dir * pow(tol / hv, pw);
      if ((tPrime / tStep) > sr) tPrime = tStep * sr;
    }

    if ((dir * (tTotal + tPrime)) >= (dir * tFinal - 1e-8)) {
      tPrime = tFinal - tTotal;
      final_sequence = true;
    }

    /* Get the acceleration at the begining of the next sequence */
    if (calcAccel(x, v, tTotal, et, beta, planets, a1, planetLookUp,
		  closeApproaches))
      return 1;
    nFC++;

    /* Predict b values for the next step.  Values from the preceeding */
    /* sequence were saved in the e matrix.  The correction bd is
       applied below */
    s[0] = tPrime / tStep;
    s[1] = s[0] * s[0];
    s[2] = s[1] * s[0];
    s[3] = s[2] * s[0];
    s[4] = s[3] * s[0];
    s[5] = s[4] * s[0];
    s[6] = s[5] * s[0];

    for (k=0; k<_NDE; k++) {
      if (nSeq != 1) {
	for (j=0; j<7; j++) {
	  bd[j][k] = b[j][k] - e[j][k];
	}
      }

      e[0][k] = s[0]*( 7.0*b[6][k] +  6.0*b[5][k] +  5.0*b[4][k] + 4.0*b[3][k] + 3.0*b[2][k] + 
		       2.0*b[1][k] + b[0][k]);
      e[1][k] = s[1]*(21.0*b[6][k] + 15.0*b[5][k] + 10.0*b[4][k] + 6.0*b[3][k] + 3.0*b[2][k] +
		      b[1][k]);
      e[2][k] = s[2]*(35.0*b[6][k] + 20.0*b[5][k] + 10.0*b[4][k] + 4.0*b[3][k] +     b[2][k]);
      e[3][k] = s[3]*(35.0*b[6][k] + 15.0*b[5][k] +  5.0*b[4][k] +     b[3][k]);
      e[4][k] = s[4]*(21.0*b[6][k] +  6.0*b[5][k] +      b[4][k]);
      e[5][k] = s[5]*( 7.0*b[6][k] +      b[5][k]);
      e[6][k] = s[6]*(     b[6][k]);

      for (l=0; l<7; l++) b[l][k] = e[l][k] + bd[l][k];
    }

    /* Two iterations for every sequence */
    nIter = 2;
  } while(1);

  return 0;
}

/******************************************************************************/
/* initConstants()
   One time initialization of constants.
*/
void initConstants(double *xc, double *vc, double *c, double *d, double *r) {
  int k, l, la, lb, lc, ld, le;
  int nw[] = {-1, -1, 0, 2, 5, 9, 14, 20};

  memset(c, 0, sizeof(double)*21);
  memset(d, 0, sizeof(double)*21);
  memset(r, 0, sizeof(double)*21);

  /* Prepare the constants */
  xc[0] = 1.0 / 2.0;
  xc[1] = 1.0 / 6.0;
  xc[2] = 1.0 / 12.0;
  xc[3] = 1.0 / 20.0;
  xc[4] = 1.0 / 30.0;
  xc[5] = 1.0 / 42.0;
  xc[6] = 1.0 / 56.0;
  xc[7] = 1.0 / 72.0;

  vc[0] = 1.0 / 2.0;
  vc[1] = 1.0 / 3.0;
  vc[2] = 1.0 / 4.0;
  vc[3] = 1.0 / 5.0;
  vc[4] = 1.0 / 6.0;
  vc[5] = 1.0 / 7.0;
  vc[6] = 1.0 / 8.0;

  c[0] = -h[1];
  d[0] = h[1];
  r[0] = 1.0 / (h[2] - h[1]); 
 la = 0;
  lc = 0;

  for (k=2; k<7; k++) {
    lb = la;
    la = lc + 1;
    lc = nw[k+1];

    c[la] = -h[k] * c[lb];
    c[lc] = c[la-1] - h[k];
    d[la] = h[1] * d[lb];
    d[lc] = -c[lc];
    r[la] = 1.0 / (h[k+1] - h[1]);
    r[lc] = 1.0 / (h[k+1] - h[k]);
    if (k > 2) {
      for (l=2; l<k; l++) {
	ld = la + l - 1;
	le = lb + l - 2;
	c[ld] = c[le] - h[k] * c[le+1];
	d[ld] = d[le] + h[l] * d[le+1];
	r[ld] = 1.0 / (h[k+1] - h[l]);
      }
    }
  }
  return;
}

/******************************************************************************/
/* calcAccel()
   Syndyne force + planetary perturbations.
*/
int calcAccel(double *RHc, double *v, double t, double et, double beta,
	      const int planets, double *a, const int planetLookUp,
	      const int closeApproaches) {
  double RHp[27], rhat[3], RHtemp[27];
  double RHc1, RHc2, RHc3, RCp1, RHp3, RCp3, vr, min;
  double now;
  int i, j, pBit, CAplanet;

  /* Comet heliocentric distance and (the same)^2 and ^3 */
  RHc1 = sqrt(RHc[0]*RHc[0] + RHc[1]*RHc[1] + RHc[2]*RHc[2]);
  RHc2 = pow(RHc1, 2);
  RHc3 = pow(RHc1, 3);

  /* Acceleration due to the Sun */
  a[0] = -_MU * RHc[0] / RHc3;
  a[1] = -_MU * RHc[1] / RHc3;
  a[2] = -_MU * RHc[2] / RHc3;

  /* Add radiation forces (including Poynting-Robertson drag) */
  if (beta != 0.0) {
    rhat[0] = RHc[0] / RHc1;
    rhat[1] = RHc[1] / RHc1;
    rhat[2] = RHc[2] / RHc1;

    /* v dot rhat = radial velocity */
    vr = v[0]*rhat[0] + v[1]*rhat[1] + v[2]*rhat[2];

    a[0] += beta * _MU / RHc2 * ((1.0 - vr / _C) * rhat[0] - v[0] / _C);
    a[1] += beta * _MU / RHc2 * ((1.0 - vr / _C) * rhat[1] - v[1] / _C);
    a[2] += beta * _MU / RHc2 * ((1.0 - vr / _C) * rhat[2] - v[2] / _C);
  }

  /* Planetary perturbations */
  if (planets) {
    now = et + t;

    if (planetLookUp) {
      if (planet_lookup(now, planets, RHp)) return 1;

      if (closeApproaches) {
        /* Determine which planet we could be close to */
        min = 60 * _AU;
        for (i=0; i<9; i++) {
          if (min > RHc1 - CAplanet_a[i]) {
            j = i;
            min = RHc1 - CAplanet_a[i];
          }
        }

        /* How close are we? */
        RCp1 = sqrt((RHc[0] - RHp[j*3  ])*(RHc[0] - RHp[j*3  ]) +
                    (RHc[1] - RHp[j*3+1])*(RHc[1] - RHp[j*3+1]) +
                    (RHc[2] - RHp[j*3+2])*(RHc[2] - RHp[j*3+2]));

        if (RCp1 < CAdistance[j]) {
          if (get_planet_xyz_et(now, (int)pow(2, j), RHtemp)) return 1;
          RHp[j*3  ] = RHtemp[j*3  ];
          RHp[j*3+1] = RHtemp[j*3+1];
          RHp[j*3+2] = RHtemp[j*3+2];
          i = 9;
        }
      }
    } else {
      if (get_planet_xyz_et(now, planets, RHp)) return 1;
    }

    for (i=0; i<9; i++) {
      if (i == 0) pBit = 1; else pBit *= 2;
      if (planets & pBit) {	/* There must be a planet to include */
	/* Planet heliocentric distance cubed */
	RHp3 = pow(sqrt(RHp[i*3  ]*RHp[i*3  ] +
			RHp[i*3+1]*RHp[i*3+1] + 
			RHp[i*3+2]*RHp[i*3+2]), 3);
	/* Planet cometocentric distance cubed */
	RCp3 = pow(sqrt((RHc[0] - RHp[i*3  ])*(RHc[0] - RHp[i*3  ]) +
			(RHc[1] - RHp[i*3+1])*(RHc[1] - RHp[i*3+1]) +
			(RHc[2] - RHp[i*3+2])*(RHc[2] - RHp[i*3+2])), 3);

	/* Add in planetary attraction, the first term is the direct accel, */
	/* the second term is the indirect accel because the Sun at the */
	/* origin is not the center of mass of the system */
	a[0] += MUplanets[i] * ((RHp[i*3  ] - RHc[0]) / RCp3 - RHp[i*3  ] / RHp3);
	a[1] += MUplanets[i] * ((RHp[i*3+1] - RHc[1]) / RCp3 - RHp[i*3+1] / RHp3);
	a[2] += MUplanets[i] * ((RHp[i*3+2] - RHc[2]) / RCp3 - RHp[i*3+2] / RHp3);
      }
    } /* end for */
  } /* end if (planets) */

  return 0;
}
