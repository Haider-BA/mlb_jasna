/***********************************************************************
 *
 * mlb.c
 *
 * Copyright (c) 2015 Ulf D. Schiller <ulf@lattice-boltzmann.de>
 * All rights reserved.
 *
 ***********************************************************************/

#include <math.h>
#include <stdio.h>
#include "mkl_service.h"
#include "mkl_trans.h"
#include "mkl_blas.h"
#include "mkl_lapack.h"
#include "derivFD.h"
#include "eos.h"

#define TOLERANCE 1.e-3

#define GAUSS_SEIDEL
#define MINRES

/***********************************************************************/

void mlb_calc_force(double *force, LB_Moments *m, int x, int y) {
  const double *w = lbmodel.fd_weights[3];
  const double (*c)[lbmodel.n_dim] = lbmodel.c;
  int i;
  double rho, nb_rho;

  force[0] = force[1] = 0.0;

  rho = m->rho;

  for (i=0; i<lbmodel.n_fd; ++i) {
    nb_rho = ((double *)m)[lblattice.nb_offset[i]*lbmodel.n_vel];
    force[0] += w[i]*c[i][0]*nb_rho;
    force[1] += w[i]*c[i][1]*nb_rho;
  }

  //thirdDer(force, &(m->rho));

  force[0] *= lbpar.kappa*rho;
  force[1] *= lbpar.kappa*rho;

}

/***********************************************************************/

static void mlb_calc_current(double *jc, LB_Moments *m, int x, int y) {
  int i, j;
  double Dp[lbmodel.n_dim],
    Dpmrdp[lbmodel.n_dim],
    Du[lbmodel.n_dim][lbmodel.n_dim],
    D2u[lbmodel.n_dim][lbmodel.n_dim][lbmodel.n_dim];

  double *p     = &(m->p);     // m + 6
  double *pmrdp = &(m->pmrdp); // m + 8
  double *u     = m->u;        // m + 12

  firstDer(Dp, p);
  firstDer(Dpmrdp, pmrdp);
  for (i=0; i<lbmodel.n_dim; ++i) {
    firstDer(Du[i], &u[i]);
    secDerAB(D2u[i], &u[i]);
  }

  for (i=0; i<lbmodel.n_dim; ++i) {
    jc[i] = 0.0;
    for (j=0; j<lbmodel.n_dim; ++j) {
      jc[i] += Dpmrdp[i] * Du[j][j];
      jc[i] += Dp[j] * (Du[i][j] + Du[j][i]);
      jc[i] += (*pmrdp + *p) * D2u[j][i][j];
      jc[i] += *p * D2u[i][j][j];
    }
    jc[i] /= -12.;
  }

}

/***********************************************************************/

inline static void mlb_matrix_current(double *jc, double **M, LB_Moments *m, int x, int y) {
  int i, j, k, ind, nbi;
  double rho, *u, *u_nb;

  ind = 2*(x*lblattice.stride[0] + y);
  rho = m->rho;
  u   = m->u;

  for (j=0; j<lbmodel.n_dim; ++j) {

    jc[j] = 0.;

    for (i=0; i<lbmodel.n_fd; ++i) {

      nbi = ind + lblattice.nb_offset[i]*2;

      u_nb = u + lblattice.nb_offset[i]*lbmodel.n_vel;

      for (k=0; k<lbmodel.n_dim; ++k) {

	jc[j] += M[ind+j][nbi+k] * u_nb[k];

      }

    }

    jc[j] += rho*u[j];

  }

}

/***********************************************************************/

inline static void mlb_construct_matrix(double **M, double *b, double *m0) {
  int i, j, k, x, y, xl, xh, yl, yh, xoff;
  int ind, nbi;
  double *m, *rho, *q, *f, *p, *pmrdp;
  double Dp[lbmodel.n_dim], Dpmrdp[lbmodel.n_dim];

  double first[lbmodel.n_fd][lbmodel.n_dim],
    second[lbmodel.n_fd][lbmodel.n_dim][lbmodel.n_dim];

  firstDer_coeff(first);
  secDer_coeff(second);

  xl = lblattice.halo_size[0];
  xh = lblattice.halo_size[0] + lblattice.grid[0];
  yl = lblattice.halo_size[1];
  yh = lblattice.halo_size[1] + lblattice.grid[1];

  xoff = lblattice.stride[0] - (yh - yl);

  m = m0 + lbmodel.n_vel*(xl*lblattice.stride[0]+yl);

  for (x=xl; x<xh; ++x, m+=lbmodel.n_vel*xoff) {
    for (y=yl; y<yh; ++y, m+=lbmodel.n_vel) {

      ind = 2*(x*lblattice.stride[0] + y);

      rho   = &(((LB_Moments *)m)->rho);   // m + 0
      q     =  (((LB_Moments *)m)->j);     // m + 1
      p     = &(((LB_Moments *)m)->p);     // m + 6
      pmrdp = &(((LB_Moments *)m)->pmrdp); // m + 8
      f     =  (((LB_Moments *)m)->force); // m + 10

      firstDer(Dp, p);
      firstDer(Dpmrdp, pmrdp);

      for (j=0; j<lbmodel.n_dim; ++j) {

	for (i=0; i<lbmodel.n_fd; ++i) {

	  nbi  = ind + lblattice.nb_offset[i]*2;

	  for (k=0; k<lbmodel.n_dim; ++k) {
	    M[ind+j][nbi+k] = 0.;
	  }

	  for (k=0; k<lbmodel.n_dim; ++k) {
	    M[ind+j][nbi+k] += Dpmrdp[j] * first[i][k];
	    M[ind+j][nbi+j] += Dp[k] * first[i][k];
	    M[ind+j][nbi+k] += Dp[k] * first[i][j];
	    M[ind+j][nbi+k] += (*pmrdp + *p) * second[i][j][k];
	    M[ind+j][nbi+j] += *p * second[i][k][k];
	  }

	  for (k=0; k<lbmodel.n_dim; ++k) {
	    M[ind+j][nbi+k] /= -12.;
	  }

	}

	M[ind+j][ind+j] -= *rho;

	b[ind+j] = -(q[j] + 0.5*f[j]);

      }

    }
  }

}

/***********************************************************************/

#ifdef GAUSS_SEIDEL
static void gauss_seidel(double **M, double *b, double *phi) {
  int i, j, k, x, y, xl, xh, yl, yh, ind, nbi, niter=0;
  double sigma, d, dmax;

  xl = lblattice.halo_size[0];
  xh = lblattice.halo_size[0] + lblattice.grid[0];
  yl = lblattice.halo_size[1];
  yh = lblattice.halo_size[1] + lblattice.grid[1];

  /* Gauss-Seidel iteration */

  do {

    dmax = 0.0; ++niter;

    for (x=xl; x<xh; ++x) {
      for (y=yl; y<yh; ++y) {
	ind = 2*(x*lblattice.stride[0] + y);
	for (j=0; j<lbmodel.n_dim; ++j) {
	  sigma = 0.;
	  //for (x2=xl; x2<xh; ++x2) {
	  //  for (y2=yl; y2<yh; ++y2) {
	  //    nbi = 2*(x2*lblattice.stride[0] + y2);
	  for (i=0; i<lbmodel.n_fd; ++i) {
	      nbi = ind + lblattice.nb_offset[i]*2;
	      for (k=0; k<lbmodel.n_dim; ++k) {
		sigma += M[ind+j][nbi+k]*phi[nbi+k];
	      }
	      //}
	  }
	  //if (fabs(M[ind+i][ind+i]) < 1.e-9) {
	  //  fprintf(stderr, "Singular matrix (%d,%d)[%d]=%f!\n",x,y,i,M[ind+i][ind+i]);
	  //}
	  sigma = (b[ind+j] - sigma)/M[ind+j][ind+j];
	  d = fabs(sigma);
	  if (d > dmax) dmax = d;
	  phi[ind+j] += sigma;
	}
      }
    }

    //for (x=xl; x<xh; ++x) {
    //  for (y=yl; y<yh; ++y) {
    //	ind = 2*(x*lblattice.stride[0] + y);
    //	for (j=0; j<lbmodel.n_dim; ++j) {
    //	  phi[ind+j] = phi_new[ind+j];
    //	}
    //  }
    //}

#if 0
    for (i=il; i<ih; ++i) {
      sigma = 0.;
      for (j=il; j<ih; ++j) {
	sigma += M[i][j]*phi[j];
      }
      sigma -= M[i][i]*phi[i]; // subtract i=j contribution from loop above
      if (fabs(M[i][i]) < 1.e-9) {
	fprintf(stderr, "Singular matrix (%d,%d)=%f!\n",i,i,M[i][i]);
      }
      sigma = (b[i] - sigma)/M[i][i];
      d = fabs(phi[i] - sigma);
      if (d > dmax) dmax = d;
      phi[i] = sigma;
    }
#endif

  } while (dmax > TOLERANCE);

  fprintf(stderr, "Gauss-Seidel converged after %d iteration(s).\n", niter);

}
#endif

/***********************************************************************/

#ifdef MINRES
static void minres(double **M, double *b, double *phi) {
  int j, k, x, y, x2, y2, xl, xh, yl, yh, ind, nbi, niter=0;
  double *ptmp=NULL;

  xl = lblattice.halo_size[0];
  xh = lblattice.halo_size[0] + lblattice.grid[0];
  yl = lblattice.halo_size[1];
  yh = lblattice.halo_size[1] + lblattice.grid[1];

  /* MINRES */

  double *r  = malloc(lblattice.halo_grid_volume*2*sizeof(*r));
  double *p0 = malloc(lblattice.halo_grid_volume*2*sizeof(*r));
  double *p1 = malloc(lblattice.halo_grid_volume*2*sizeof(*r));
  double *p2 = malloc(lblattice.halo_grid_volume*2*sizeof(*r));
  double *s0 = malloc(lblattice.halo_grid_volume*2*sizeof(*r));
  double *s1 = malloc(lblattice.halo_grid_volume*2*sizeof(*r));
  double *s2 = malloc(lblattice.halo_grid_volume*2*sizeof(*r));

  double alpha, beta, r2, rdots, sdots, s0dots1, s1dots1, s0dots2, s2dots2;

  r2 = 0.;
  for (x=xl; x<xh; ++x) {
    for (y=yl; y<yh; ++y) {
      ind = 2*(x*lblattice.stride[0]+y);
      for (j=0; j<lbmodel.n_dim; ++j) {
	r[ind+j] = b[ind+j];
	for (x2=xl; x2<xh; ++x2) {
	  for (y2=yl; y2<yh; ++y2) {
	    nbi = 2*(x2*lblattice.stride[0]+y2);
	    for (k=0; k<lbmodel.n_dim; ++k) {
	      r[ind+j] -= M[ind+j][nbi+k]*phi[nbi+k];
	    }
	  }
	}
	p0[ind+j] = r[ind+j];
	r2 += r[ind+j]*r[ind+j];
      }
    }
  }
  sdots = 0.;
  for (x=xl; x<xh; ++x) {
    for (y=yl; y<yh; ++y) {
      ind = 2*(x*lblattice.stride[0]+y);
      for (j=0; j<lbmodel.n_dim; ++j) {
	s0[ind+j] = 0.;
	for (x2=xl; x2<xh; ++x2) {
	  for (y2=yl; y2<yh; ++y2) {
	    nbi = 2*(x2*lblattice.stride[0]+y2);
	    for (k=0; k<lbmodel.n_dim; ++k) {
	      s0[ind+j] += M[ind+j][nbi+k]*p0[nbi+k];
	    }
	  }
	}
	sdots += s0[ind+j]*s0[ind+j];
      }
    }
  }

  while (r2 > TOLERANCE*TOLERANCE) {

    ++niter;

    ptmp = p2; p2 = p1; p1 = p0; p0 = ptmp;
    ptmp = s2; s2 = s1; s1 = s0; s0 = ptmp;

    rdots = sdots = 0.;
    for (x=xl; x<xh; ++x) {
      for (y=yl; y<yh; ++y) {
	ind = 2*(x*lblattice.stride[0]+y);
	for (j=0; j<lbmodel.n_dim; ++j) {
	  rdots +=  r[ind+j]*s1[ind+j];
	  sdots += s1[ind+j]*s1[ind+j];
	}
      }
    }
    alpha = rdots/sdots;

    r2 = 0.;
    for (x=xl; x<xh; ++x) {
      for (y=yl; y<yh; ++y) {
	ind = 2*(x*lblattice.stride[0]+y);
	for (j=0; j<lbmodel.n_dim; ++j) {
	  phi[ind+j] += alpha*p1[ind+j];
	  r[ind+j]   -= alpha*s1[ind+j];
	  r2 += r[ind+j]*r[ind+j];
	}
      }
    }

    if (r2 < TOLERANCE*TOLERANCE) break;

    s0dots1 = s1dots1 = s0dots2 = s2dots2 = 0.;
    for (x=xl; x<xh; ++x) {
      for (y=yl; y<yh; ++y) {
	ind = 2*(x*lblattice.stride[0]+y);
	for (j=0; j<lbmodel.n_dim; ++j) {
	  p0[ind+j] = s1[ind+j];
	  s0[ind+j] = 0.;
	  for (x2=xl; x2<xh; ++x2) {
	    for (y2=yl; y2<yh; ++y2) {
	      nbi = 2*(x2*lblattice.stride[0]+y2);
	      for (k=0; k<lbmodel.n_dim; ++k) {
		s0[ind+j] += M[ind+j][nbi+k]*s1[nbi+k];
	      }
	    }
	  }
	  s0dots1 += s0[ind+j]*s1[ind+j];
	  s0dots2 += s0[ind+j]*s2[ind+j];
	  s1dots1 += s1[ind+j]*s1[ind+j];
	  s2dots2 += s2[ind+j]*s2[ind+j];
	}
      }
    }
    //s0dots1 = s1dots1 = 0.;
    //for (x=xl; x<xh; ++x) {
    //  for (y=yl; y<yh; ++y) {
    //	ind = 2*(x*lblattice.stride[0]+y);
    //	for (j=0; j<lbmodel.n_dim; ++j) {
    //	  s0dots1 += s0[ind+j]*s1[ind+j];
    //	  s1dots1 += s1[ind+j]*s1[ind+j];
    //	}
    //  }
    //}
    beta = s0dots1/s1dots1;
    for (x=xl; x<xh; ++x) {
      for (y=yl; y<yh; ++y) {
	ind = 2*(x*lblattice.stride[0]+y);
	for (j=0; j<lbmodel.n_dim; ++j) {
	  p0[ind+j] -= beta*p1[ind+j];
	  s0[ind+j] -= beta*s1[ind+j];
	}
      }
    }
    if (niter > 1) {
      //s0dots1 = s1dots1 = 0.;
      //for (x=xl; x<xh; ++x) {
      //	for (y=yl; y<yh; ++y) {
      //	  ind = 2*(x*lblattice.stride[0]+y);
      //	  for (j=0; j<lbmodel.n_dim; ++j) {
      //	    s0dots1 += s0[ind+j]*s2[ind+j];
      //	    s1dots1 += s2[ind+j]*s2[ind+j];
      //	  }
      //	}
      //}
      beta = s0dots2/s2dots2;
      for (x=xl; x<xh; ++x) {
	for (y=yl; y<yh; ++y) {
	  ind = 2*(x*lblattice.stride[0]+y);
	  for (j=0; j<lbmodel.n_dim; ++j) {
	    p0[ind+j] -= beta*p2[ind+j];
	    s0[ind+j] -= beta*s2[ind+j];
	  }
	}
      }
    }

  }

  fprintf(stderr, "MINRES converged after %d iteration(s).\n", niter);

  free(r);
  free(p0);
  free(p1);
  free(p2);
  free(s0);
  free(s1);
  free(s2);

}
#endif

/***********************************************************************/

static void minres2(double **M, double *b, double *sol) {
  int j, k, x, y, x2, y2, xl, xh, yl, yh, ind, nbi, niter=0;
  double r2, alfa, beta, beta1, g1, g2, d1, d2, e1, e2, c, s, phi, tau;
  double *ptmp=NULL;

  xl = lblattice.halo_size[0];
  xh = lblattice.halo_size[0] + lblattice.grid[0];
  yl = lblattice.halo_size[1];
  yh = lblattice.halo_size[1] + lblattice.grid[1];

  double *v, *v1, *v2, *w, *w1, *w2;
  v  = malloc(lblattice.halo_grid_volume*2*sizeof(*v));
  v1 = malloc(lblattice.halo_grid_volume*2*sizeof(*v));
  v2 = malloc(lblattice.halo_grid_volume*2*sizeof(*v));
  w  = malloc(lblattice.halo_grid_volume*2*sizeof(*w));
  w1 = malloc(lblattice.halo_grid_volume*2*sizeof(*w));
  w2 = malloc(lblattice.halo_grid_volume*2*sizeof(*w));

  beta = 0.;
  for (x=xl; x<xh; ++x) {
    for (y=yl; y<yh; ++y) {
      ind = 2*(x*lblattice.stride[0]+y);
      for (j=0; j<lbmodel.n_dim; ++j) {
	v[ind+j] = b[ind+j];
	for (x2=xl; x2<xh; ++x2) {
	  for (y2=yl; y2<yh; ++y2) {
	    nbi = 2*(x2*lblattice.stride[0]+y2);
	    for (k=0; k<lbmodel.n_dim; ++k) {
	      v[ind+j] -= M[ind+j][nbi+k]*sol[ind+k];
	    }
	  }
	}
	beta += v[ind+j]*v[ind+j];
	v1[ind+j] = 0.;
	w[ind+j]  = 0.;
	w1[ind+j] = 0.;
      }
    }
  }
  beta  = sqrt(beta);

  for (x=xl; x<xh; ++x) {
    for (y=yl; y<yh; ++y) {
      ind = 2*(x*lblattice.stride[0]+y);
      for (j=0; j<lbmodel.n_dim; ++j) {
	v[ind+j] /= beta;
      }
    }
  }

  r2  = beta;
  tau = phi = beta;
  d1  = 0.;
  c   = -1.;
  s   = 0.;

  while (r2 > TOLERANCE) {

    ++niter;

    ptmp = v2; v2 = v1; v1 = v; v = ptmp;
    ptmp = w2; w2 = w1; w1 = w; w = ptmp;

    beta1 = beta;
    e2    = e1;

    /* Lanczos recurrence */
    alfa = 0.0;
    for (x=xl; x<xh; ++x) {
      for (y=yl; y<yh; ++y) {
	ind = 2*(x*lblattice.stride[0]+y);
	for (j=0; j<lbmodel.n_dim; ++j) {
	  v1[ind+j] /= beta1;
	  v[ind+j] = 0.;
	  for (x2=xl; x2<xh; ++x2) {
	    for (y2=yl; y2<yh; ++y2) {
	      nbi = 2*(x2*lblattice.stride[0]+y2);
	      for (k=0; k<lbmodel.n_dim; ++k) {
		v[ind+j] += M[ind+j][nbi+k]*v1[nbi+k];  /* r = A*v */
	      }
	    }
	  }
	  alfa += v1[ind+j]*v[ind+j];                   /* alfa = v'*r */
	}
      }
    }

    for (x=xl; x<xh; ++x) {
      for (y=yl; y<yh; ++y) {
	ind = 2*(x*lblattice.stride[0]+y);
	for (j=0; j<lbmodel.n_dim; ++j) {
	  v[ind+j] -= alfa*v1[ind+j] + beta1*v2[ind+j]; /* Lanczos update */
	  beta += v[ind+j]*v[ind+j];                    /* beta = v'*v */
	}
      }
    }
    beta  = sqrt(beta);

    /* QR factorisation */

    /* left orthogonalization by previous Givens rotation */
    d2 =   c * d1 + s * alfa;
    g1 =   s * d1 - c * alfa;
    e1 =            s * beta;
    d1 =          - c * beta;

    /* new Givens rotation for subdiagonal */
    g2  = sqrt( g1*g1 + beta*beta );
    c   = g1 / g2;
    s   = beta / g2;
    tau = c * phi;
    phi = s * phi;

    /* update of solution */
    for (x=xl; x<xh; ++x) {
      for (y=yl; y<yh; ++y) {
	ind = 2*(x*lblattice.stride[0]+y);
	for (j=0; j<lbmodel.n_dim; ++j) {
	  w[ind+j] = v1[ind+j] - d2*w1[ind+j] - e2*w2[ind+j];
	  w[ind+j] /= g2;
	  sol[ind+j] += tau*w[ind+j];
	}
      }
    }

    r2 *= s;

  }

  fprintf(stderr, "MINRES2 converged after %d iteration(s).\n", niter);

  free(v);
  free(v1);
  free(v2);
  free(w);
  free(w1);
  free(w2);

}

/***********************************************************************/

static void mkl_factorise(double **M, double *b, double *phi) {
  int j, x, y, xl, xh, yl, yh, ind;

  xl = lblattice.halo_size[0];
  xh = lblattice.halo_size[0] + lblattice.grid[0];
  yl = lblattice.halo_size[1];
  yh = lblattice.halo_size[1] + lblattice.grid[1];

  int i, *ipiv, info;
  int n     = 2*(yh - yl)*(xh - xl);
  int lda   = n;
  int ldb   = n;
  int nrhs  = 1;
  int lwork = 64*n*sizeof(double);

  double **a, *bb, *work;
  char uplo = 'L';

  a     = malloc(n*sizeof(*a));
  ipiv  = malloc(n*sizeof(*ipiv));
  *a    = mkl_malloc(n*n*sizeof(**a),16);
  bb    = mkl_malloc(n*sizeof(*bb),16);
  work  = mkl_malloc(64*n*sizeof(*work),16);

  for (i=0; i<n; ++i) a[i] = a[0] + i*n;

  int k, x0, y0, x2, y2, nbi;
  for (x0=0, x=xl; x<xh; ++x) {
    for (y=yl; y<yh; ++y) {
      ind = 2*(x*lblattice.stride[0]+y);
      for (j=0; j<lbmodel.n_dim; ++j, ++x0) {
	for (y0=0, x2=xl; x2<xh; ++x2) {
	  for (y2=yl; y2<yh; ++y2) {
	    nbi = 2*(x2*lblattice.stride[0]+y2);
	    for (k=0; k<lbmodel.n_dim; ++k, ++y0) {
	      a[x0][y0] = M[ind+j][nbi+k];
	    }
	  }
	}
	bb[x0] = b[ind+j];
      }
    }
  }

  dsytrf(&uplo, &n, *a, &lda, ipiv, work, &lwork, &info);
  dsytrs(&uplo, &n, &nrhs, *a, &lda, ipiv, bb, &ldb, &info);

   for (x0=0, x=xl; x<xh; ++x) {
    for (y=yl; y<yh; ++y) {
      ind = 2*(x*lblattice.stride[0]+y);
      for (j=0; j<lbmodel.n_dim; ++j, ++x0) {
	phi[ind+j] = bb[x0];
      }
    }
  }

  mkl_free(work);
  mkl_free(bb);
  mkl_free(*a);
  free(ipiv);
  free(a);

}

/***********************************************************************/

static void mlb_solve_matrix(double **M, double *b, double *m0) {
  int j, x, y, xl, xh, yl, yh, xoff, ind;
  double *m, *phi0, *phi, *phi_new;

  phi0 = malloc(lblattice.halo_grid_volume*2*sizeof(*phi));
  phi = malloc(lblattice.halo_grid_volume*2*sizeof(*phi));
  phi_new = malloc(lblattice.halo_grid_volume*2*sizeof(*phi));

  xl = lblattice.halo_size[0];
  xh = lblattice.halo_size[0] + lblattice.grid[0];
  yl = lblattice.halo_size[1];
  yh = lblattice.halo_size[1] + lblattice.grid[1];

  xoff = lblattice.stride[0] - (yh - yl);

  /* load initial guess from LB moments */
  m = m0 + lbmodel.n_vel*(xl*lblattice.stride[0]+yl);
  for (x=xl; x<xh; ++x, m+=lbmodel.n_vel*xoff) {
    for (y=yl; y<yh; ++y, m+=lbmodel.n_vel) {
      ind = 2*(x*lblattice.stride[0] + y);
      for (j=0; j<lbmodel.n_dim; ++j) {
	phi0[ind+j] = 0.0;//((LB_Moments *)m)->u[j];
      }
    }
  }

  memcpy(phi, phi0, lblattice.halo_grid_volume*2*sizeof(*phi));
  gauss_seidel(M, b, phi);

  memcpy(phi, phi0, lblattice.halo_grid_volume*2*sizeof(*phi));
  minres(M, b, phi);

  memcpy(phi_new, phi0, lblattice.halo_grid_volume*2*sizeof(*phi));
  minres2(M, b, phi_new);

  memcpy(phi, phi0, lblattice.halo_grid_volume*2*sizeof(*phi));
  mkl_factorise(M, b, phi);

  /* copy solution to LB moments */
  m = m0 + lbmodel.n_vel*(xl*lblattice.stride[0]+yl);
  for (x=xl; x<xh; ++x, m+=lbmodel.n_vel*xoff) {
    for (y=yl; y<yh; ++y, m+=lbmodel.n_vel) {
      ind = 2*(x*lblattice.stride[0] + y);
      for (j=0; j<lbmodel.n_dim; ++j) {
	fprintf(stderr, "phi[%d] = %f phi_new[%d] = %f u(%d,%d)[%d] = %f\n",ind+j,phi[ind+j],ind+j,phi_new[ind+j],x,y,j,((LB_Moments *)m)->u[j]);
	double rho = ((LB_Moments *)m)->rho;
	double *q  = ((LB_Moments *)m)->j;
	double *u  = ((LB_Moments *)m)->u;
	double *g  = ((LB_Moments *)m)->force;
	((LB_Moments *)m)->u[j]     = phi[ind+j];
	((LB_Moments *)m)->jcorr[j] = rho*u[j] - q[j] - 0.5*g[j];
      }
    }
  }

  free(phi);
  free(phi_new);

}

/***********************************************************************/

static void ic_read(double *dmax, double **M, LB_Moments *m, int x, int y) {
  double jnew[lbmodel.n_dim], *jc = m->jcorr, d;
  double jc2[lbmodel.n_dim];

  jnew[0] = jnew[1] = 0.0;

  mlb_calc_current(jnew, m, x, y);

  mlb_matrix_current(jc2, M, m, x, y);

  //fprintf(stderr, "jnew = (%.9f,%.9f)\tjmatrix = (%.9f,%.9f)\n",jnew[0],jnew[1],jc2[0],jc2[1]);

  d = fabs(jnew[0] - jc[0]);
  if (d > *dmax) *dmax = d;
  d = fabs(jnew[1] - jc[1]);
  if (d > *dmax) *dmax = d;

  jc[0] = jnew[0];
  jc[1] = jnew[1];

  //double norm = jc[0]*jc[0]+jc[1]*jc[1];
  //if (norm > 1e3) {
  //  fprintf(stderr, "Warning! Large correction current jc = (%f,%f)\n", jc[0],jc[1]);
  //}

}

/***********************************************************************/

static void ic_write(LB_Moments *m, int x, int y) {

  double rho = m->rho;
  double *j  = m->j;
  double *g  = m->force;
  double *u  = m->u;
  double *jc = m->jcorr;

  u[0] = (j[0] + 0.5*g[0] + jc[0])/rho;
  u[1] = (j[1] + 0.5*g[1] + jc[1])/rho;

}

/***********************************************************************/

static void ic_read_column(double *dmax, double **M, double *m, int x) {
  int y, yl, yh;

  yl = lblattice.halo_size[1];
  yh = lblattice.halo_size[1] + lblattice.grid[1] - 1;

  for (y=yl, m+=yl*lbmodel.n_vel; y<=yh; ++y, m+=lbmodel.n_vel) {
    ic_read(dmax, M, (LB_Moments *)m, x, y);
  }

}

/***********************************************************************/

static void ic_write_column(double *m, int x) {
  int y, yl, yh;

  yl = lblattice.halo_size[1];
  yh = lblattice.halo_size[1] + lblattice.grid[1] - 1;

  for (y=yl, m+=yl*lbmodel.n_vel; y<=yh; ++y, m+=lbmodel.n_vel) {
    ic_write((LB_Moments *)m, x, y);
  }

}

/***********************************************************************/

static void mlb_init_current(double *m) {
  int x, y, xl, xh, yl, yh, xoff;
  double rho, p, dp, d2p, force[lbmodel.n_dim];

  xl = lblattice.halo_size[0];
  xh = lblattice.halo_size[0] + lblattice.grid[0];
  yl = lblattice.halo_size[1];
  yh = lblattice.halo_size[1] + lblattice.grid[1];

  xoff = lblattice.stride[0] - (yh - yl);

  m += lbmodel.n_vel*(xl*lblattice.stride[0]+yl);

  for (x=xl; x<xh; ++x, m+=lbmodel.n_vel*xoff) {
    for (y=yl; y<yh; ++y, m+=lbmodel.n_vel) {

      LB_Moments *lbmom = (LB_Moments *)m;

      rho = lbmom->rho;
      p   = rho*eq_state(rho);
      dp  = derP(rho);
      d2p = der2P(rho);

      /* Step 1 */
      lbmom->p     = p;
      lbmom->dp    = dp;
      lbmom->pmrdp = p - rho*dp;
      lbmom->rd2p  = rho*d2p;

      /* Step 2 */
      mlb_calc_force(force, lbmom, x, y);
      lbmom->force[0] = force[0];
      lbmom->force[1] = force[1];

      /* correction current */
      //lbmom->jcorr[0] = 0.0; /* correction current */
      //lbmom->jcorr[1] = 0.0;

      /* Step 3 */
      lbmom->u[0] = (lbmom->j[0] + 0.5*force[0] + lbmom->jcorr[0])/rho;
      lbmom->u[1] = (lbmom->j[1] + 0.5*force[1] + lbmom->jcorr[1])/rho;

    }
  }

}

/***********************************************************************/

void mlb_correction_current(double *m0) {
  double *m = m0;
  int i, x, xl, xh;
  int xstride = lblattice.stride[0]*lbmodel.n_vel;

  int niter = 0;
  double dmax = 0.0;

  double ** M = malloc(lblattice.halo_grid_volume*2*sizeof(*M));
  double *b   = calloc(lblattice.halo_grid_volume*2,sizeof(*b));
  *M = calloc(lblattice.halo_grid_volume*2*lblattice.halo_grid_volume*2,sizeof(**M));
  for (i=0; i<lblattice.halo_grid_volume*2; ++i) {
    M[i] = *M + i*lblattice.halo_grid_volume*2;
  }

  lb_halo_copy(); /* need up to date densities in halo */

  mlb_init_current(m0);

  mlb_construct_matrix(M, b, m0);

  do {

    m = m0; dmax = 0.0; ++niter;

    lb_halo_copy(); /* need up to date currents in halo */

    //fprintf(stderr, "Starting iteration #%d of implicit algorithm...\n", niter);

    /* no information is `streaming' so we loop the internal region */
    xl = lblattice.halo_size[0];
    xh = lblattice.halo_size[0]+lblattice.grid[0]+VMAX;

    /* Columns in the lower range will read only */
    for (x=xl, m+=xl*xstride; x<xl+VMAX; ++x, m+=xstride) {
      ic_read_column(&dmax, M, m, x);
    }

    /* Calculate the correction current but do not overwrite the column yet
     * x-MAXV can be overwritten with the updated current */
    for (x=xl+VMAX; x<xh-VMAX; ++x, m+=xstride) {
      ic_read_column(&dmax, M, m, x);
      ic_write_column(m-VMAX*xstride, x-VMAX);
    }

    /* Columns in the higher range will write only */
    for (x=xh-VMAX; x<xh; ++x, m+=xstride) {
      ic_write_column(m-VMAX*xstride, x-VMAX);
    }

    //fprintf(stderr, "Iteration #%d: dmax = %f\n", niter, dmax);

  } while (dmax > TOLERANCE);

  fprintf(stderr, "Implicit algorithm converged after %d iteration(s).\n", niter);

  mlb_solve_matrix(M, b, m0);

  free(*M);
  free(M);
  free(b);

}
  
/***********************************************************************/

void mlb_interface_collisions(double *f) {
  int i;
  LB_Moments *m = (LB_Moments *)(f + lblattice.halo_grid_volume*lbmodel.n_vel);
  double cs2, fc, *force = m->force;
  double w[lbmodel.n_vel];

  cs2 = eq_state(RHO_MEAN);
  lb_weights(w, cs2);

  for (i=0; i<lbmodel.n_vel; ++i) {
    fc = lbmodel.c[i][0]*force[0] + lbmodel.c[i][1]*force[1];
    f[i] += 0.5*(1. + lbpar.gamma)*w[i]/cs2*fc;
  }

}

/***********************************************************************/

inline static double delta(int a, int b){
	if (a == b) return 1.0;
	else return 0.0;
}

/***********************************************************************/

inline static void mlb_calc_sigma(double Sigma[][lbmodel.n_dim], LB_Moments *m) {
  const double gamma = lbpar.gamma;
  int i, j;

  double Dr[lbmodel.n_dim],
    Dp[lbmodel.n_dim],
    Dpdu[lbmodel.n_dim],
    Dpmrdp[lbmodel.n_dim],
    Dpmrdpdivu[lbmodel.n_dim],
    Drdpudu[lbmodel.n_dim][lbmodel.n_dim],
    Du[lbmodel.n_dim][lbmodel.n_dim],
    D2p[lbmodel.n_dim][lbmodel.n_dim],
    D2u[lbmodel.n_dim][lbmodel.n_dim][lbmodel.n_dim],
    divu, divj;

  double *rho   = &(m->rho);   // m
  double *p     = &(m->p);     // m + 6
  double *dp    = &(m->dp);    // m + 7
  double *pmrdp = &(m->pmrdp); // m + 8
  double *rd2p  = &(m->rd2p);  // m + 9
  double *u     = m->u;        // m + 12;

  firstDer(Dr, rho);
  firstDer(Dp, p);
  firstDer(Dpmrdp, pmrdp);
  firstDer(Du[0], &u[0]);
  firstDer(Du[1], &u[1]);
  secDerAB(D2p, p);
  secDerAB(D2u[0], &u[0]);
  secDerAB(D2u[1], &u[1]);

  divu = divj = 0.;
  for (i=0; i<lbmodel.n_dim; ++i) {
    divu += Du[i][i];
    divj += Dr[i]*u[i] + *rho*Du[i][i];
  }

  /* Step 10 and 11 */
  for (i=0; i<lbmodel.n_dim; ++i) {
    Dpdu[i] = ( Dp[0] * ( Du[0][i] + Du[i][1] )
		+ Dp[1] * ( Du[1][i] + Du[i][1] )
		+ *p * ( D2u[0][0][i] + D2u[1][1][i]
			 + D2u[i][0][0] + D2u[i][1][1] ) );
    Dpmrdpdivu[i] = Dpmrdp[i]*divu - *pmrdp*(D2u[0][i][0] + D2u[1][i][1]);
    for (j=0; j<lbmodel.n_dim; ++j) {
      Drdpudu[i][j] = ( D2p[i][j]/(*rho) - Dr[i]*Dp[j]/(*rho**rho)
		       + Du[0][i]*Du[j][0] + Du[1][i]*Du[j][1]
		       + u[0]*D2u[j][i][0] + u[1]*D2u[j][i][1] );
    }
  }

  /* calculate Sigma */
  for (i=0; i<lbmodel.n_dim; ++i) {
    for (j=0; j<lbmodel.n_dim; ++j) {
      Sigma[i][j] = ( -0.25*(gamma+1.)*(gamma+1.)/(gamma-1.)
		      * ( u[i]*Dpmrdpdivu[j] + u[j]*Dpmrdpdivu[i]
			  + u[i]*Dpdu[j] + u[j]*Dpdu[i] )
		      - (gamma*gamma + 4.*gamma + 1.)/(gamma-1.)/6.
		      * ( *dp * divj * ( Du[i][j] + Du[j][i] )
			  + *p * ( Drdpudu[i][j] + Drdpudu[j][i] ) ) );
    }
    Sigma[i][i] -= ( (gamma*gamma + 4.*gamma + 1.)/(gamma-1.)/6.
		     * ( *pmrdp * ( Drdpudu[0][0] + Drdpudu[1][1] )
			 - *rd2p * divu * divj ) );
  }

}

/***********************************************************************/

inline static void mlb_calc_xi(double Xi[][lbmodel.n_dim][lbmodel.n_dim],
			       LB_Moments *m) {
  const double gamma = lbpar.gamma;
  int i, j, k;

  double Dr[lbmodel.n_dim],
    Dp[lbmodel.n_dim],
    Dpr[lbmodel.n_dim],
    Du[lbmodel.n_dim][lbmodel.n_dim],
    Dpuu[lbmodel.n_dim][lbmodel.n_dim][lbmodel.n_dim],
    divu;

  double *rho   = &(m->rho);   // m
  double *p     = &(m->p);     // m + 6
  double *pmrdp = &(m->pmrdp); // m + 8
  double *u     = m->u;        // m + 12
  double *jcorr = m->jcorr;    // m + 14

  double cs2 = eq_state(RHO_MEAN);

  //fprintf(stderr, "cs2= %f p/rho = %f\t%f\n",cs2,*p/(*rho),cs2-*p/(*rho));

  firstDer(Dr, rho);
  firstDer(Dp, p);
  firstDer(Du[0], &u[0]);
  firstDer(Du[1], &u[1]);

  /* Step 10 and 11 */
  divu = 0.;
  for (i=0; i<lbmodel.n_dim; ++i) {
    divu += Du[i][i];
    Dpr[i] = Dp[i]/(*rho) - *p/(*rho**rho)*Dr[i];
    for (j=0; j<lbmodel.n_dim; ++j) {
      for (k=0; k<lbmodel.n_dim; ++k) {
	Dpuu[i][j][k] = Dp[i]*u[j]*u[k] + *p*Du[j][i]*u[k] + *p*u[j]*Du[k][i];
      }
    }
  }

  /* calculate Xi */
  for (i=0; i<lbmodel.n_dim; ++i) {
    for (j=0; j<lbmodel.n_dim; ++j) {
      for (k=0; k<lbmodel.n_dim; ++k) {
	Xi[i][j][k] = ( (gamma*gamma + 4.*gamma + 1.)/(gamma + 1.)/3.
			* ( Dpuu[k][i][j] + Dpuu[i][j][k] + Dpuu[j][k][i]
			    + delta(i,j) * ( *p*Dpr[k] + *pmrdp*u[k]*divu )
			    + delta(j,k) * ( *p*Dpr[i] + *pmrdp*u[i]*divu )
			    + delta(k,i) * ( *p*Dpr[j] + *pmrdp*u[j]*divu ) ) );
	Xi[i][j][k] -= ( (gamma - 1.) * cs2
			* ( delta(i,j)*jcorr[k]
			    + delta(j,k)*jcorr[i]
			    + delta(k,i)*jcorr[j] ) );
      }
    }
  }

}

/***********************************************************************/

void mlb_correction_collisions(double *f) {
  const double (*c)[lbmodel.n_dim] = lbmodel.c;
  LB_Moments *m = (LB_Moments *)(f + lblattice.halo_grid_volume*lbmodel.n_vel);
  int i, j, k, l;
  double cs2, jc, sc, xc;
  double w[lbmodel.n_vel];
  double Sigma[lbmodel.n_dim][lbmodel.n_dim];
  double Xi[lbmodel.n_dim][lbmodel.n_dim][lbmodel.n_dim];

  double *jcorr = m->jcorr; // m + 14

  cs2 = eq_state(RHO_MEAN);
  lb_weights(w, cs2);

  mlb_calc_sigma(Sigma, m);
  mlb_calc_xi(Xi, m);

  for (i=0; i<lbmodel.n_vel; ++i) {
    jc = sc = xc = 0.;
    for (j=0; j<lbmodel.n_dim; ++j) {
      jc += jcorr[j]*c[i][j];
      for (k=0; k<lbmodel.n_dim; ++k) {
	sc += Sigma[j][k]*(c[i][j]*c[i][k] - cs2 * delta(j,k));
	for (l=0; l<lbmodel.n_dim; ++l) {
	  xc += Xi[j][k][l] * ( c[i][j]*c[i][k]*c[i][l] - cs2 * delta(j,k) * c[i][l]- cs2 * delta(k,l) * c[i][j] - cs2 * delta(l,j) * c[i][k] );
	}
      }
    }
    f[i] += (lbpar.gamma - 1.)*w[i]/cs2*jc;
    f[i] += w[i]/(2.*cs2*cs2)*sc;
    f[i] += w[i]/(6.*cs2*cs2*cs2)*xc;
  }

}
   
   /* see eos.h
double eq_state(double rho){
	double sigma=0.0;
	double integral=0.0;
	
	if(rho <= R1 || rho >= R3 ){
		integral = 0.0;
	}
	else if(R1 < rho && rho < R3){
		integral = 2*A*(R2-R1)/M_PI*sin(0.5*M_PI*(rho-2*R1+R0)/(R2-R1))*sin(0.5*M_PI*(rho-R0)/(R2-R1));
	} else {
		printf ("I've lost myself in the EOS!\n");
		exit (0);
	}
	sigma= S1 * exp(integral);
	return sigma;
}
   */
   
/***********************************************************************/
