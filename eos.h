/***********************************************************************
 *
 * eos.h
 *
 * Copyright (c) 2015 Ulf D. Schiller <ulf@lattice-boltzmann.de>
 * All rights reserved.
 *
 ***********************************************************************/

#define EOS_GAUSS

#ifndef EOS_H
#define EOS_H

#include <stdlib.h>

#define SCALE (1.0)

/* phase separation densities */
#define RHO_LOW  (1.049171*pow(SCALE,3))  // 1.049171 1.113113 1.481619
#define RHO_HIGH (1.248991*pow(SCALE,3))  // 1.248991 1.541998 1.886775
#define RHO_MEAN ((RHO_LOW+RHO_HIGH)/2.0) // 1.149081 1.327556 1.684197

#define R1 (0.5*pow(SCALE,3)) /* densities */
#define R2 (1.0*pow(SCALE,3))
#define R3 (2.0*R2-R1)

/* constants that define equation of state: p/rho */
#if defined EOS_PSI
#define A0 (1.2/pow(SCALE,3)) /* specific volume */
#define S1 0.6                /* reference speed of sound squared */
#endif

#ifdef EOS_GAUSS
#define S1 0.6                /* baseline speed of sound squared */
#define S2 0.75
#define W  0.25*(R2-R1)
#endif

#ifdef EOS_VDW
#define A0 (9./49.)
#define B0 (2./21.)
#define T0 0.57
#endif

/***********************************************************************/

#if defined EOS_PSI

static double psi(double rho) {
  return A0*sin(M_PI*(rho-R1)/(R2-R1));
}

static double dpsi(double rho) {
  return A0*M_PI/(R2-R1)*cos(M_PI*(rho-R1)/(R2-R1));
}

static double intpsi(double rho) {
  return A0*(R2-R1)/M_PI*(1.0 - cos(M_PI*(rho-R1)/(R2-R1)));
}

/***********************************************************************/

static double eq_state(double rho) {
  double cs2;
  
  if (rho <= R1 || rho >= R3) {
    fprintf(stderr, "Don't like to be here with rho = %f\n", rho);
    cs2 = S1;
  } else if (rho < R3) {
    cs2 = S1 * exp(intpsi(rho));
  } else {
    printf ("I've lost myself in the EOS (rho=%f)!\n", rho);
    exit (0);
  }

  return cs2;
    
}

/***********************************************************************/

static double derP(double rho) {
  double dp;

  if (rho <= R1 || rho >= R3) {
    dp = S1;
  } else {
    /* dp = eq_state(rho) + rho*derS(rho); */
    dp = eq_state(rho)*(1. + rho*psi(rho));
  }

  return dp;

}
  
/***********************************************************************/

static double der2P(double rho) {
  double d2p;
  
  if (rho <= R1 || rho >= R3) {
    d2p = 0.0;
  } else {
    /* d2p = 2.0*derS(rho) + rho*der2S(rho); */
    d2p = eq_state(rho)*(psi(rho)*(2. + rho*psi(rho)) + rho*dpsi(rho));
  }

  return d2p;

}

#elif defined EOS_GAUSS

/***********************************************************************/

static double eq_state(double rho) {
  double cs2;

  cs2 = (S2 - S1) * exp( -(rho-R2)*(rho-R2)/(2.*W*W) ) + S1;

  return cs2;
}

static double derS(double rho) {
  double ds;

  ds = -(rho-R2)/(W*W)*(eq_state(rho) - S1);

  return ds;
}

static double der2S(double rho) {
  double d2s;

  d2s = -(eq_state(rho) - S1 + (rho-R2)*derS(rho))/(W*W);

  return d2s;
}

static double derP(double rho) {
  double dp;

  dp = eq_state(rho) + rho*derS(rho);

  return dp;
}

static double der2P(double rho) {
  double d2p;

  d2p = 2.*derS(rho) + rho*der2S(rho);

  return d2p;
}

#else

/***********************************************************************/

/* van der Waals equation of state */

static double eq_state(double rho) {
  double cs2;

  cs2= T0/(1.-B0*rho)-A0*rho;

  return cs2;

}

static double derS(double rho) {
  double ds;

  ds = T0*B0/((1.-B0*rho)*(1.-B0*rho)) - A0;

  return ds;
}

static double der2S(double rho) {
  double d2s;

  d2s = 2.*T0*B0*B0/((1.-B0*rho)*(1.-B0*rho)*(1.-B0*rho));

  return d2s;
}

static double derP(double rho) {
  double dp;

  dp = eq_state(rho) + rho*derS(rho);

  return dp;
}

static double der2P(double rho) {
  double d2p;

  d2p = 2.*derS(rho) + rho*der2S(rho);

  return d2p;
}

#endif

/***********************************************************************/

inline static void write_eos() {
  FILE* file;
  double rho, cs2, p, dp, d2p, psi;

  file = fopen("eos.dat", "w");

  for (rho=R1; rho<=R3; rho+=(R3-R1)/100) {

    cs2 = eq_state(rho);
    p = rho*cs2;
    dp = derP(rho);
    d2p = der2P(rho);
    psi = (derP(rho)/cs2 - 1.)/rho;

    fprintf(file, "%f %f %f %f %f %f\n", rho, cs2, p, dp, d2p, psi);

  }

  fclose(file);

}

/***********************************************************************/

#endif
