/***********************************************************************
 *
 * d2q21.c
 *
 * Copyright (c) 2009-2015 Ulf D. Schiller <ulf@lattice-boltzmann.de>
 * All rights reserved.
 *
 * Please cite the following publications when using this code:
 *
 * U. D. Schiller 
 * "A unified operator splitting approach for multi-scale
 * fluid–particle coupling in the lattice Boltzmann method"
 * Comp. Phys. Comm. 185, 2586-2597 (2014)
 * http://dx.doi.org/10.1016/j.cpc.2014.06.005
 * 
 * B. Duenweg, U. D. Schiller, A. J. C. Ladd
 * "Statistical mechanics of the fluctuating lattice Boltzmann equation"
 * Phys. Rev. E 76, 036704 (2007)
 * http://dx.doi.org/10.1103/PhysRevE.76.036704
 *
 ***********************************************************************/

#include <time.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "lis.h"
#include "derivFD.h"
#include "d2q21.h"
#include "mlb.h"
#include "eos.h"

/***********************************************************************/

/* temporary/secondary grid */
#define PFI *FI(NVEL,WGRID)
#define FI(i,x) fi[i][x]

/***********************************************************************/

static double *lbf = NULL;
static double PFI;

/***********************************************************************/

#define TOLERANCE 1.e-6

void lb_check_halo() {
  int i, x, y, xt, yt;
  double *f1, *f2, *m1, *m2;

  /***************
   * X direction *
   ***************/

  for (y=0; y<lblattice.halo_grid[1]; ++y) {
    for (x=0; x<lblattice.halo_size[0]; ++x) {
      xt = lblattice.grid[0] + x;
      f1 = lbf + (x*lblattice.stride[0]+y)*lbmodel.n_vel;
      f2 = lbf + (xt*lblattice.stride[0]+y)*lbmodel.n_vel;
      m1 = f1 + lblattice.halo_grid_volume*lbmodel.n_vel;
      m2 = f2 + lblattice.halo_grid_volume*lbmodel.n_vel;
      for (i=0; i<lbmodel.n_vel; ++i) {
	if (fabs(f2[i]-f1[i]) > TOLERANCE) {
	  fprintf(stderr, "f-Halo error at (%d,%d) to (%d,%d) [%d]\n", x, y, xt, y, i);
	  return;
	}
	if (fabs(m2[i]-m1[i]) > TOLERANCE) {
	  fprintf(stderr, "m-Halo error at (%d,%d) to (%d,%d) [%d]\n", x, y, xt, y, i);
	  return;
	}
      }
    }
    for (x=lblattice.grid[0]+lblattice.halo_size[0]; x<lblattice.halo_grid[0]; ++x) {
      xt = x - lblattice.grid[0];
      f1 = lbf + (x*lblattice.stride[0]+y)*lbmodel.n_vel;
      f2 = lbf + (xt*lblattice.stride[0]+y)*lbmodel.n_vel;
      m1 = f1 + lblattice.halo_grid_volume*lbmodel.n_vel;
      m2 = f2 + lblattice.halo_grid_volume*lbmodel.n_vel;
      for (i=0; i<lbmodel.n_vel; ++i) {
	if (fabs(f2[i]-f1[i]) > TOLERANCE) {
	  fprintf(stderr, "f-Halo error at (%d,%d) to (%d,%d) [%d]\n", x, y, xt, y, i);
	  return;
	}
	if (fabs(m2[i]-m1[i]) > TOLERANCE) {
	  fprintf(stderr, "m-Halo error at (%d,%d) to (%d,%d) [%d]\n", x, y, xt, y, i);
	  return;
	}
      }
    }
  }

  /***************
   * Y direction *
   ***************/

  for (x=0; x<lblattice.halo_grid[0]; ++x) {
    for (y=0; y<lblattice.halo_size[1]; ++y) {
      yt = lblattice.grid[1] + y;
      f1 = lbf + (x*lblattice.stride[0]+y)*lbmodel.n_vel;
      f2 = lbf + (x*lblattice.stride[0]+yt)*lbmodel.n_vel;
      m1 = f1 + lblattice.halo_grid_volume*lbmodel.n_vel;
      m2 = f2 + lblattice.halo_grid_volume*lbmodel.n_vel;
      for (i=0; i<lbmodel.n_vel; ++i) {
	if (fabs(f2[i]-f1[i]) > TOLERANCE) {
	  fprintf(stderr, "f-Halo error at (%d,%d) to (%d,%d) [%d]\n", x, y, x, yt, i);
	  return;
	}
	if (fabs(m2[i]-m1[i]) > TOLERANCE) {
	  fprintf(stderr, "m-Halo error at (%d,%d) to (%d,%d) [%d]\n", x, y, x, yt, i);
	  return;
	}
      }
    }
    for (y=lblattice.grid[1]+lblattice.halo_size[1]; y<lblattice.halo_grid[1]; ++y) {
      yt = y - lblattice.grid[1];
      f1 = lbf + (x*lblattice.stride[0]+y)*lbmodel.n_vel;
      f2 = lbf + (x*lblattice.stride[0]+yt)*lbmodel.n_vel;
      m1 = f1 + lblattice.halo_grid_volume*lbmodel.n_vel;
      m2 = f2 + lblattice.halo_grid_volume*lbmodel.n_vel;
      for (i=0; i<lbmodel.n_vel; ++i) {
        if (fabs(f2[i]-f1[i]) > TOLERANCE) {
          fprintf(stderr, "f-Halo error at (%d,%d) to (%d,%d) [%d]\n", x, y, x, yt, i);
	  return;
        }
	if (fabs(m2[i]-m1[i]) > TOLERANCE) {
	  fprintf(stderr, "m-Halo error at (%d,%d) to (%d,%d) [%d]\n", x, y, x, yt, i);
	  return;
	}
      }
    }
  }

}

/***********************************************************************/

void lb_halo_copy() {

  int x;
  int lsrc, rsrc, ldst, rdst, size;

  int totalsize = lbmodel.n_vel*lblattice.halo_grid_volume;

  int xstride = lbmodel.n_vel*lblattice.stride[0];
  int ystride = lbmodel.n_vel*lblattice.stride[1];

  /***************
   * X direction *
   ***************/

  size = lblattice.halo_size[0]*xstride;

  ldst = 0;
  lsrc = size;
  rsrc = lblattice.grid[0]*xstride;
  rdst = rsrc + size;

  memcpy(lbf+ldst, lbf+rsrc, size*sizeof(*lbf));
  memcpy(lbf+rdst, lbf+lsrc, size*sizeof(*lbf));

  /* also copy the modes */
  memcpy(lbf+ldst+totalsize, lbf+rsrc+totalsize, size*sizeof(*lbf));
  memcpy(lbf+rdst+totalsize, lbf+lsrc+totalsize, size*sizeof(*lbf));

  /***************
   * Y direction *
   ***************/

  size = lblattice.halo_size[1]*ystride;

  ldst = 0;
  lsrc = size;
  rsrc = lblattice.grid[1]*ystride;
  rdst = rsrc + size;

  for (x=0; x<lblattice.halo_grid[0]; x++) {

    /* populations are handled directly in the update loop */
    /* but only certain velocity directions??? */
    memcpy(lbf+ldst, lbf+rsrc, size*sizeof(*lbf));
    memcpy(lbf+rdst, lbf+lsrc, size*sizeof(*lbf));

    /* only copy the modes */
    memcpy(lbf+ldst+totalsize, lbf+rsrc+totalsize, size*sizeof(*lbf));
    memcpy(lbf+rdst+totalsize, lbf+lsrc+totalsize, size*sizeof(*lbf));

    ldst += xstride;
    lsrc += xstride;
    rsrc += xstride;
    rdst += xstride;

  }

}

/***********************************************************************/

static void lb_calc_local_moments(double *f) {

  int i;
  double *m = f + lblattice.halo_grid_volume*lbmodel.n_vel;
  double rho, p, dp, d2p;

  for (i=0; i<lbmodel.n_vel; ++i) {
    m[i] = 0.0;
  }

  for (i=0; i<lbmodel.n_vel; ++i) {
    m[0] += f[i];
    m[1] += f[i]*lbmodel.c[i][0];
    m[2] += f[i]*lbmodel.c[i][1];
    m[3] += f[i]*lbmodel.c[i][0]*lbmodel.c[i][0];
    m[4] += f[i]*lbmodel.c[i][1]*lbmodel.c[i][1];
    m[5] += f[i]*lbmodel.c[i][0]*lbmodel.c[i][1];
  }

  /* evaluate equation of state */
  rho = m[0];
  p   = rho*eq_state(rho);
  dp  = derP(rho);
  d2p = der2P(rho);

  LB_Moments *lbmom = (LB_Moments *)m;

  /* pressure and derivatives */
  lbmom->p     = p;
  lbmom->dp    = dp;
  lbmom->pmrdp = p - rho*dp;
  lbmom->rd2p  = rho*d2p;

  /* force */
  lbmom->force[0] = 0.0;
  lbmom->force[1] = 0.0;

  /* velocity */
  lbmom->u[0] = m[1]/m[0];
  lbmom->u[1] = m[2]/m[0];

  /* correction current */
  //lbmom->jcorr[0] = 0.0;
  //lbmom->jcorr[1] = 0.0;

}

/***********************************************************************/

void lb_calc_moments() {
  int x, y, xl, xh, yl, yh, xoff;
  double *f = lbf;

  xl = lblattice.halo_size[0];
  xh = lblattice.halo_size[0] + lblattice.grid[0];
  yl = lblattice.halo_size[1];
  yh = lblattice.halo_size[1] + lblattice.grid[1];

  xoff = lblattice.stride[0] - (yh - yl);

  f += lbmodel.n_vel*(xl*lblattice.stride[0]+yl);

  for (x=xl; x<xh; ++x, f+=lbmodel.n_vel*xoff) {
    for (y=yl; y<yh; ++y, f+=lbmodel.n_vel) {
      lb_calc_local_moments(f);
    }
  }

}

/***********************************************************************/

static void lb_calc_equilibrium(double *f_eq, double *f) {

  int i;
  double w[lbmodel.n_vel];
  double rho, *u, uc, u2, cs2;
  LB_Moments *m = (LB_Moments *)(f + lblattice.halo_grid_volume*lbmodel.n_vel);

  rho  = m->rho;
  //u[0] = (m[1] + 0.5*force[0])/rho;
  //u[1] = (m[2] + 0.5*force[1])/rho;
  u    = m->u; // m + 12

  cs2 = eq_state(rho);
  lb_weights(w, cs2);

  u2   = (u[0]*u[0] + u[1]*u[1])/cs2;

  for (i=0; i<lbmodel.n_vel; ++i) {
    uc = (u[0]*lbmodel.c[i][0] + u[1]*lbmodel.c[i][1])/cs2;
    f_eq[i] = w[i]*rho*(1.0 + uc + 0.5*uc*uc - 0.5*u2);
  }

}

/***********************************************************************/

static void lb_bulk_collisions(double *f) {

  int i;
  double f_eq[lbmodel.n_vel];
  
  lb_calc_equilibrium(f_eq, f);

  for (i=0; i<lbmodel.n_vel; ++i) {
    f[i] += (lbpar.gamma - 1.) * ( f[i] - f_eq[i] );
  }
  
}

/***********************************************************************/

static void lb_collisions(double *f, int x, int y) {

  lb_bulk_collisions(f);

  mlb_interface_collisions(f);

  mlb_correction_collisions(f);

  int i;
  double rho, w[lbmodel.n_vel];
  for (i=0; i<lbmodel.n_vel; ++i) {
    if (f[i] < 0.) {
      rho = (f+lbmodel.n_vel*lblattice.halo_grid_volume)[0];
      lb_weights(w, eq_state(rho));
      fprintf(stderr, "Negative populations at (%d,%d) f[%d] = %f rho = %f w[%d] = %f\n", x, y, i, f[i], rho, i, w[i]);
      exit(-1);
    }
  }

}

/***********************************************************************/

static void lb_stream(double *f, double PFI, int x, int y) {

  int xc, xp, xm, xp2, xm2, xp4, xm4;
  xc = x%WGRID;
  xp  = (x+1)%WGRID; xm  = (x-1+WGRID)%WGRID;
  xp2 = (x+2)%WGRID; xm2 = (x-2+WGRID)%WGRID;
  xp4 = (x+4)%WGRID; xm4 = (x-4+WGRID)%WGRID;

  FI( 0, xc )[y]   = f[0];
  FI( 1, xp )[y]   = f[1];
  FI( 2, xm )[y]   = f[2];
  FI( 3, xc )[y+1] = f[3];
  FI( 4, xc )[y-1] = f[4];
  FI( 5, xp )[y+1] = f[5];
  FI( 6, xm )[y-1] = f[6];
  FI( 7, xm )[y+1] = f[7];
  FI( 8, xp )[y-1] = f[8];
  FI( 9, xp2)[y]   = f[9];
  FI(10, xm2)[y]   = f[10];
  FI(11, xc )[y+2] = f[11];
  FI(12, xc )[y-2] = f[12];
  FI(13, xp2)[y+2] = f[13];
  FI(14, xm2)[y-2] = f[14];
  FI(15, xm2)[y+2] = f[15];
  FI(16, xp2)[y-2] = f[16];
  FI(17, xp4)[y]   = f[17];
  FI(18, xm4)[y]   = f[18];
  FI(19, xc )[y+4] = f[19];
  FI(20, xc )[y-4] = f[20];
  
}

/***********************************************************************/

static void lb_pbc_y(double PFI, int x) {

  int yl, yh;
  yl = lblattice.halo_size[1];
  yh = lblattice.halo_size[1] + lblattice.grid[1] - 1;

  int xc, xp, xm, xp2, xm2;
  xc  = x%WGRID;
  xp  = (x+1)%WGRID; xm  = (x-1+WGRID)%WGRID;
  xp2 = (x+2)%WGRID; xm2 = (x-2+WGRID)%WGRID;

  /* periodic boundary conditions in y */
  FI( 3, xc )[yl]   = FI( 3, xc )[yh+1];
  FI( 5, xp )[yl]   = FI( 5, xp )[yh+1];
  FI( 7, xm )[yl]   = FI( 7, xm )[yh+1];
  FI( 4, xc )[yh]   = FI( 4, xc )[yl-1];
  FI( 6, xm )[yh]   = FI( 6, xm )[yl-1];
  FI( 8, xp )[yh]   = FI( 8, xp )[yl-1];

  FI(11, xc )[yl]   = FI(11, xc )[yh+1];
  FI(13, xp2)[yl]   = FI(13, xp2)[yh+1];
  FI(15, xm2)[yl]   = FI(15, xm2)[yh+1];
  FI(12, xc )[yh]   = FI(12, xc )[yl-1];
  FI(14, xm2)[yh]   = FI(14, xm2)[yl-1];
  FI(16, xp2)[yh]   = FI(16, xp2)[yl-1];
  FI(11, xc )[yl+1] = FI(11, xc )[yh+2];
  FI(13, xp2)[yl+1] = FI(13, xp2)[yh+2];
  FI(15, xm2)[yl+1] = FI(15, xm2)[yh+2];
  FI(12, xc )[yh-1] = FI(12, xc )[yl-2];
  FI(14, xm2)[yh-1] = FI(14, xm2)[yl-2];
  FI(16, xp2)[yh-1] = FI(16, xp2)[yl-2];

  FI(19, xc )[yl]   = FI(19, xc )[yh+1];
  FI(19, xc )[yl+1] = FI(19, xc )[yh+2];
  FI(19, xc )[yl+2] = FI(19, xc )[yh+3];
  FI(19, xc )[yl+3] = FI(19, xc )[yh+4];
  FI(20, xc )[yh]   = FI(20, xc )[yl-1];
  FI(20, xc )[yh-1] = FI(20, xc )[yl-2];
  FI(20, xc )[yh-2] = FI(20, xc )[yl-3];
  FI(20, xc )[yh-3] = FI(20, xc )[yl-4];

}

/***********************************************************************/

static void lb_write_back(double *f, double PFI, int x, int y) {

  int xc = x%WGRID;

  f[ 0] = FI( 0, xc)[y];
  f[ 1] = FI( 1, xc)[y];
  f[ 2] = FI( 2, xc)[y];
  f[ 3] = FI( 3, xc)[y];
  f[ 4] = FI( 4, xc)[y];
  f[ 5] = FI( 5, xc)[y];
  f[ 6] = FI( 6, xc)[y];
  f[ 7] = FI( 7, xc)[y];
  f[ 8] = FI( 8, xc)[y];
  f[ 9] = FI( 9, xc)[y];
  f[10] = FI(10, xc)[y];
  f[11] = FI(11, xc)[y];
  f[12] = FI(12, xc)[y];
  f[13] = FI(13, xc)[y];
  f[14] = FI(14, xc)[y];
  f[15] = FI(15, xc)[y];
  f[16] = FI(16, xc)[y];
  f[17] = FI(17, xc)[y];
  f[18] = FI(18, xc)[y];
  f[19] = FI(19, xc)[y];
  f[20] = FI(20, xc)[y];

}

/***********************************************************************/

static void lb_read(double *f, double PFI, int x, int y) {

  lb_collisions(f, x, y);
  lb_stream(f, fi, x, y);

}

/***********************************************************************/

static void lb_write(double *f, double PFI, int x, int y) {

  lb_write_back(f, fi, x, y);
  lb_calc_local_moments(f);

}

/***********************************************************************/

static void lb_read_column(double *f, double PFI, int x) {
  int y, yl, yh;

  yl = lblattice.halo_size[1];
  yh = lblattice.halo_size[1] + lblattice.grid[1] - 1;

  for (y=yl, f+=yl*lbmodel.n_vel; y<=yh; ++y, f+=lbmodel.n_vel) {
    lb_read(f, fi, x, y);
  }

  lb_pbc_y(fi, x);

}

/***********************************************************************/

static void lb_write_column(double *f, double PFI, int x) {
  int y, yl, yh;

  yl = lblattice.halo_size[1];
  yh = lblattice.halo_size[1] + lblattice.grid[1] - 1;

  for (y=yl, f+=yl*lbmodel.n_vel; y<=yh; ++y, f+=lbmodel.n_vel) {
    lb_write(f, fi, x, y);
  }
  
}

/***********************************************************************/

static void lb_collide_stream(double *f) {
  int x, xl, xh;
  int xstride = lblattice.stride[0]*lbmodel.n_vel;

  lb_halo_copy(); /* need up to date moments in halo */

  xl = lblattice.halo_size[0]-VMAX;
  xh = lblattice.halo_size[0]+lblattice.grid[0]+VMAX;

  /* Columns in the lower range will be read only */
  for (x=xl, f+=xl*xstride; x<xl+2*VMAX; ++x, f+=xstride) {
    lb_read_column(f, fi, x);
  }

  /* Collide and stream column x, read back column x-MAXV
   * x-MAXV can be overwritten and all info is available now */
  for (x=xl+2*VMAX; x<xh; ++x, f+=xstride) {
    lb_read_column(f, fi, x);
    lb_write_column(f-VMAX*xstride, fi, x-VMAX);
  }

}

/***********************************************************************/

static void lb_update(double *f) {
  double *m = f + lblattice.halo_grid_volume*lbmodel.n_vel;

  mlb_correction_current(m);

  lb_collide_stream(f);

  lb_calc_moments();

}

/***********************************************************************/

static void lb_init_interface() {
  int x, y, i;
  double rho, cs2, w[lbmodel.n_vel];
  double dx;
  double *f = lbf;

  for (x=0; x<lblattice.halo_grid[0]; ++x) {
    for (y=0; y<lblattice.halo_grid[1]; ++y, f+=lbmodel.n_vel) {
      dx = (double)(x - lblattice.halo_grid[0]/2);
      rho = 0.5*(RHO_HIGH-RHO_LOW)*(cos(2.*M_PI/lblattice.grid[0]*dx)) + RHO_MEAN;
      for (i=0; i<lbmodel.n_vel; ++i) {
	cs2 = eq_state(rho);
	lb_weights(w, cs2);
	f[i] = w[i]*rho;
      }
      lb_calc_local_moments(f);
    }
  }

}

/***********************************************************************/

inline static void lb_init_droplet() {
  int x, y, i;
  double rho, cs2, w[lbmodel.n_vel];
  double dx, dy, rc, re, a, b, iw, cos;
  double *f = lbf;

  iw = 0.1;//*lblattice.grid[0]/20.0;
  a  = 5.0;//*lblattice.grid[0]/20.0;
  b  = 5.0;//*lblattice.grid[0]/20.0;

  for (x=0; x<lblattice.halo_grid[0]; ++x) {
    for (y=0; y<lblattice.halo_grid[1]; ++y, f+=lbmodel.n_vel) {
      dx = (double)(x - lblattice.halo_grid[0]/2);
      dy = (double)(y - lblattice.halo_grid[1]/2);
      rc = sqrt(dx*dx + dy*dy);
      if (rc==0.) {
	cos = 1.;
      } else {
	cos = dx/rc;
      }
      re = a*b/sqrt(a*a + (b*b-a*a)*cos*cos);
      rho = 0.5*(RHO_HIGH-RHO_LOW)*tanh((re-rc)/iw) + (RHO_HIGH+RHO_LOW)/2;
      for (i=0; i<lbmodel.n_vel; ++i) {
	cs2 = eq_state(rho);
	lb_weights(w, cs2);
	f[i] = w[i]*rho;
      }
      lb_calc_local_moments(f);
    }
  }

}

/***********************************************************************/

inline static void lb_init_random() {
}

/***********************************************************************/

static void lb_init_fluid() {

  lb_init_interface();

}

/***********************************************************************/

static void lb_init_lattice(int *grid) {
  int i, x, hgrid[lbmodel.n_dim], hsize[lbmodel.n_dim];

  lblattice.grid[0] = grid[0];
  lblattice.grid[1] = grid[1];

  lblattice.halo_size[0] = hsize[0] = HALO;
  lblattice.halo_size[1] = hsize[1] = HALO;

  lblattice.halo_grid[0] = hgrid[0] = lblattice.grid[0] + 2*hsize[0];
  lblattice.halo_grid[1] = hgrid[1] = lblattice.grid[1] + 2*hsize[1];

  lblattice.stride[1] = 1;
  lblattice.stride[0] = hgrid[1];

  lblattice.halo_grid_volume = hgrid[0]*hgrid[1];

  lbf = calloc(2*lblattice.halo_grid_volume*lbmodel.n_vel, sizeof(*lbf));

  for (i=0; i<lbmodel.n_vel; ++i) {
    for (x=0; x<WGRID; ++x) {
      FI(i,x) = calloc(lblattice.halo_grid[1],sizeof(*FI(0,0)));
    }
    lblattice.nb_offset[i] = (int)(lbmodel.c[i][0])*hgrid[1]+(int)(lbmodel.c[i][1]);
  }

}

/***********************************************************************/

static void lb_init_parameters(double rho, double gamma, double kappa) {

  lbpar.gamma = gamma;
  lbpar.kappa = kappa;

}

/***********************************************************************/

void lb_init(int *grid, double rho, double gamma, double kappa) {

  lb_init_parameters(rho, gamma, kappa);

  lb_init_lattice(grid);

  lb_init_fluid();

}

/***********************************************************************/

void lb_finalize() {
  int i, x;

  for (i=0; i<lbmodel.n_vel; ++i) {
    for (x=0; x<WGRID; ++x) {
      free(FI(i,x));
    }
  }    

  free(lbf);

}

/***********************************************************************/

void lb_mass_mom(int i) {
  int x, y, xl, xh, yl, yh, xoff;
  double rho, j[lbmodel.n_dim], j_eff[lbmodel.n_dim];
  double *u, *m = lbf + lblattice.halo_grid_volume*lbmodel.n_vel;

  xl = lblattice.halo_size[0];
  xh = lblattice.halo_size[0] + lblattice.grid[0];
  yl = lblattice.halo_size[1];
  yh = lblattice.halo_size[1] + lblattice.grid[1];

  xoff = lblattice.stride[0] - (yh - yl);

  m += lbmodel.n_vel*(lblattice.stride[0]*xl+yl);

  rho = j[0] = j[1] = j_eff[0] = j_eff[1] = 0.0;
  for (x=xl; x<xh; ++x, m+=lbmodel.n_vel*xoff) {
    for (y=yl; y<yh; ++y, m+=lbmodel.n_vel) {
      rho  += m[0];
      j[0] += m[1];
      j[1] += m[2];
      u = ((LB_Moments *)m)->u;
      j_eff[0] += m[0]*u[0];
      j_eff[1] += m[0]*u[1];
    }
  }

  fprintf(stderr, "#%d: rho = %f j = (%f,%f) j_eff=(%f,%f)\n", i, rho, j[0], j[1], j_eff[0], j_eff[1]);

}

/***********************************************************************/

void write_profile(char filename[], int write_halo) {
  int x, y, xl, xh, yl, yh, xoff;
  double rho, j[lbmodel.n_dim], *u;
  double *m = lbf + lblattice.halo_grid_volume*lbmodel.n_vel;
  FILE *file;

  if (write_halo) {
    lb_halo_copy();
    xl = 0;
    xh = lblattice.halo_grid[0];
    yl = 0;
    yh = lblattice.halo_grid[1];
  } else {
    xl = lblattice.halo_size[0];
    xh = lblattice.halo_size[0] + lblattice.grid[0];
    yl = lblattice.halo_size[1];
    yh = lblattice.halo_size[1] + lblattice.grid[1];
  }

  xoff = lblattice.stride[0] - (yh - yl);

  m += lbmodel.n_vel*(lblattice.stride[0]*xl+yl);

  file = fopen(filename,"w");
  for (x=xl; x<xh; ++x, m+=lbmodel.n_vel*xoff) {
    for (y=yl; y<yh; ++y, m+=lbmodel.n_vel) {
      rho  = m[0];
      j[0] = m[1];
      j[1] = m[2];
      u = ((LB_Moments *)m)->u;
      fprintf(file,"%f %f %f %f %f %f %f\n",(double)x-xl,(double)y-yl,rho,j[0]/rho,j[1]/rho,u[0],u[1]);
    }
  }
  fclose(file);

}

/***********************************************************************/

#if 1
int main(int argc, char *argv[]) {
  int i, n_steps, grid[lbmodel.n_dim], vol;
  double rho, gamma, kappa;
  double start, finish, elapsed, mups;

  if (argc!=3) {
    fprintf(stderr, "Usage: ./run <kappa> <nsteps>\n");
    return -1;
  }	   

  lis_initialize(&argc, &argv);

  n_steps = atoi(argv[2]);

  grid[0] = 100;
  grid[1] = 20;

  vol = grid[0]*grid[1];

  rho   = 1.0;
  gamma = 0.0;
  kappa = atof(argv[1]);

  write_eos();

  char filename[1024];
  sprintf(filename, "profile_k%.03f.dat", kappa);

  lb_init(grid,rho,gamma,kappa); lb_mass_mom(0);

  fprintf(stdout, "Running  %d iterations\n", n_steps); fflush(stdout);

  start = (double) clock();
  for (i=0; i<n_steps; ++i) {
    lb_update(lbf);
    lb_mass_mom(i+1);
    write_profile(filename, 0);
  }
  finish = (double) clock();

  elapsed = (finish-start)/CLOCKS_PER_SEC;
  mups = vol*n_steps/elapsed/1e6;

  fprintf(stdout, "Elapsed time: %.3f s (%.3e MUPS)\n", elapsed, mups); fflush(stdout); 

  write_profile(filename, 0);

  lb_finalize();

  lis_finalize();

  return EXIT_SUCCESS;

}
#endif

/***********************************************************************/
