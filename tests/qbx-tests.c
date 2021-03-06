/* This file is part of QBX, a library for Quadrature By Expansion
 *
 * Copyright (C) 2020 Michael Carley
 *
 * QBX is free software: you can redistribute it and/or modify it
 * under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.  QBX is distributed in the
 * hope that it will be useful, but WITHOUT ANY WARRANTY; without even
 * the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
 * PURPOSE.  See the GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with QBX.  If not, see <https://www.gnu.org/licenses/>.
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#include <math.h>

#include <glib.h>

#include <qbx.h>

#include <blaswrap.h>

#include "qbx-private.h"

gchar *tests[] = {"planar_test",
		  "closest_point",
		  "self_test",
		  "off_test",
		  "koornwinder",
		  "koornwinder_orthogonality",
		  "blas",
		  "koornwinder_interpolation",
		  "foreach",
		  "target_specific",
		  "weights",
		  "matrix_weights",
		  "off_weights",
		  "adaptive",
		  ""} ;

GTimer *timer ;

static gint parse_test(gchar *arg)

{
  gint i = 0 ;
  
  while ( strlen(tests[i]) != 0) {
    if ( !strcmp(tests[i], arg) ) return i ;
    i ++ ;
  }

  return -1 ;
}

static gint read_element(FILE *f, gdouble *xe, gint *ne, gint *xstr)

{
  gint i ;

  *xstr = 4 ;
  fscanf(f, "%d", ne) ;

  for ( i = 0 ; i < (*ne) ; i ++ ) {
    fscanf(f, "%lg", &(xe[i*(*xstr)+0])) ;
    fscanf(f, "%lg", &(xe[i*(*xstr)+1])) ;
    fscanf(f, "%lg", &(xe[i*(*xstr)+2])) ;
  }
  
  return 0 ;
}

static gint planar_triangle_shape_test(gdouble *xe, gint xstr, gint ne,
				       gint nq, gint N,
				       gdouble *x0, gdouble s0, gdouble t0,
				       gdouble rc, gint depth, gdouble tol,
				       gint nx)

{
  gdouble *q, C[8192], xc[3], c[3], n[3], J, fe[128], x[3], f[64], g[64] ;
  gdouble ee[16], w, rcopt, r0 ;
  gint order, str, i, Nopt, smax, sopt, fstr, nf ;
  
  g_assert(ne == 3) ;
  nf = ne ;
  
  smax = 7 ;
  
  nf = 3 ; fstr = 5 ;

  str = 2*nf ;
  qbx_quadrature_select(nq, &q, &order) ;

  if ( x0 == NULL ) {
    qbx_element_point_3d(xe, xstr, ne, s0, t0, xc, n, &J, NULL) ;
    x[0] = xc[0] ; x[1] = xc[1] ; x[2] = xc[2] ;    
  } else {
    qbx_element_point_3d(xe, xstr, ne, 0.25, 0.25, xc, n, &J, NULL) ;
    x[0] = x0[0] ; x[1] = x0[1] ; x[2] = x0[2] ;
  }

  w = sqrt(4*qbx_element_area(xe, xstr, ne, q, nq)/M_PI) ;
  /* qbx_quadrature_optimal(w, 0.5*w/(1 << smax), 4.0*w/(1 << smax), */
  /* 			 16, order, nq, 20, smax, tol, */
  /* 			 &rcopt, &Nopt, &sopt) ; */

  r0 = qbx_vector_distance(x, xc) ;
  if ( r0 < 1e-12 ) {
    qbx_quadrature_optimal_points(w, 0.5*w/(1 << smax), 4.0*w/(1 << smax),
				  16, order, nq, 20, smax, tol,
				  xc, n, x,
				  &rcopt, &Nopt, &sopt) ;
  } else {
    qbx_quadrature_optimal_points(w, 0.125*r0, 2*r0,
				  16, order, nq, 20, smax, tol,
				  xc, n, x,
				  &rcopt, &Nopt, &sopt) ;
  }

  c[0] = xc[0] + rc*n[0] ; 
  c[1] = xc[1] + rc*n[1] ; 
  c[2] = xc[2] + rc*n[2] ; 
  
  fprintf(stderr, "planar triangle test\n") ;
  fprintf(stderr, "====================\n") ;

  fprintf(stderr, "quadrature: %d points, %dth order\n", nq, order) ;
  fprintf(stderr, "subdivision depth: %d %d\n", depth, sopt) ;
  fprintf(stderr, "N: %d %d\n", N, Nopt) ;
  fprintf(stderr, "rc: %lg %lg\n", rc, rcopt) ;
  fprintf(stderr, "xc: %lg %lg %lg\n", xc[0], xc[1], xc[2]) ;
  fprintf(stderr, "n:  %lg %lg %lg\n", n[0], n[1], n[2]) ;
  fprintf(stderr, "c:  %lg %lg %lg\n", c[0], c[1], c[2]) ;
  
  /*set nodal quantities to shape function*/
  memset(fe, 0, 3*fstr*sizeof(gdouble)) ;
  fe[0*fstr+0] = 1.0 ; fe[1*fstr+1] = 1.0 ; fe[2*fstr+2] = 1.0 ;

  fprintf(stderr, "generating expansion, t=%lg\n",
	  g_timer_elapsed(timer, NULL)) ;
  qbx_expansion_make_laplace_adaptive(xe, xstr, ne, fe, fstr, nf,
				      q, nq, order, c, rc, N,
				      C, str, depth, tol, w,
				      TRUE, TRUE) ;
  qbx_expansion_make_laplace_adaptive(xe, xstr, ne, fe, fstr, nf,
  				      q, nq, order, c, rc, N,
  				      &(C[nf]), str, depth, tol, w,
  				      FALSE, FALSE) ;
  fprintf(stderr, "expansion generated, t=%lg\n",
	  g_timer_elapsed(timer, NULL)) ;

  /*evaluate on surface and find error*/
  memset(f, 0, 2*ne*sizeof(gdouble)) ;
  qbx_expansion_eval_laplace(c, N, C, str, x, f, ee) ;
  newman_tri_shape(x, &(xe[xstr*0]), &(xe[xstr*1]), &(xe[xstr*2]), NULL, 0,
		   &(g[0]), &(g[ne])) ;
  for ( i = 0 ; i < 2*ne ; i ++ ) g[i] *= -1.0/4.0/M_PI ;
  if ( s0 != G_MAXDOUBLE ) {
    g[ne  ] = (1 - s0 - t0)/2 ; g[ne+1] = s0/2 ; g[ne+2] = t0/2 ;
  }
  
  fprintf(stderr, "expansion:") ;
  for ( i = 0 ; i < nf ; i ++ ) fprintf(stderr, " %+lg", f[i]) ;
  for ( i = 0 ; i < nf ; i ++ ) fprintf(stderr, " %+lg", f[nf+i]) ;
  fprintf(stderr, "\n") ;
  fprintf(stderr, "planar:   ") ;
  for ( i = 0 ; i < 2*ne ; i ++ ) fprintf(stderr, " %+lg", g[i]) ;
  fprintf(stderr, "\n") ;
  fprintf(stderr, "error:    ") ;
  for ( i = 0 ; i < nf ; i ++ )
    fprintf(stderr, " %+lg", fabs(f[i]-g[i])) ;
  for ( i = 0 ; i < nf ; i ++ )
    fprintf(stderr, " %+lg", fabs(f[nf+i]-g[ne+i])) ;
  fprintf(stderr, "\n") ;

  return 0 ;
}

static gint self_test(gdouble *xe, gint xstr, gint ne,
		      gint nq, gint N,
		      gdouble *x0, gdouble s0, gdouble t0,
		      gdouble rc, gint depth, gdouble tol,
		      gint nx)

{
  gdouble work[8192], n[3], J, x[3], f[64], g[64], *q, fs[8], gs[8], t ;
  gint order, i, nf, fstr ;
  
  nf = ne ; fstr = ne ;

  qbx_quadrature_select(nq, &q, &order) ;

  fprintf(stderr, "triangle self-point test\n") ;
  fprintf(stderr, "========================\n") ;

  fprintf(stderr, "quadrature: %d points, %dth order\n", nq, order) ;
  fprintf(stderr, "subdivision depth: %d\n", depth) ;
  fprintf(stderr, "N: %d\n", N) ;

  fprintf(stderr, "evaluating integrals, t=%lg\n",
	  (t = g_timer_elapsed(timer, NULL))) ;
  qbx_triangle_laplace_self_quad(xe, xstr, ne, s0, t0, FALSE, q, nq, order,
				 N, depth, tol,
				 &(f[0]), fstr, &(f[ne]), fstr, work) ;
  qbx_triangle_laplace_self_quad(xe, xstr, ne, s0, t0, TRUE, q, nq, order,
				 N, depth, tol,
				 &(f[2*ne]), fstr, &(f[2*ne+ne]), fstr, work) ;
  fprintf(stderr, "integrals evaluated, t=%lg (%lg)\n",
	  g_timer_elapsed(timer, NULL), g_timer_elapsed(timer, NULL)-t) ;

  qbx_element_point_3d(xe, xstr, ne, s0, t0, x, n, &J, NULL) ;
  newman_tri_shape(x, &(xe[xstr*0]), &(xe[xstr*1]), &(xe[xstr*2]), NULL, 0,
		   &(g[0]), &(g[ne])) ;
  for ( i = 0 ; i < 2*ne ; i ++ ) g[i] *= -1.0/4.0/M_PI ;
  if ( s0 != G_MAXDOUBLE ) {
    qbx_element_shape_3d(ne, s0, t0, &(g[nf]), NULL, NULL, NULL, NULL, NULL) ;
    for ( i = 0 ; i < ne ; i ++ ) g[nf+i] *= 0.5 ;
  }

  memset(fs, 0, ne*sizeof(gdouble)) ;
  memset(gs, 0, ne*sizeof(gdouble)) ;

  for ( i = 0 ; i < ne ; i ++ ) {
    fs[0] += f[i]*1.0 ;
    fs[1] += f[i]*xe[i*xstr+0] ;
    fs[2] += f[i]*xe[i*xstr+1] ;
    fs[3] += f[2*ne+i]*1.0 ;
    fs[4] += f[2*ne+i]*xe[i*xstr+0] ;
    fs[5] += f[2*ne+i]*xe[i*xstr+1] ;
    gs[0] += g[i]*1.0 ;
    gs[1] += g[i]*xe[i*xstr+0] ;
    gs[2] += g[i]*xe[i*xstr+1] ;
  }
  
  fprintf(stderr, "single layer\n") ;
  fprintf(stderr, "expansion:") ;
  for ( i = 0 ; i < nf ; i ++ ) fprintf(stderr, " %+lg", f[i]) ;
  fprintf(stderr, "\n") ;
  fprintf(stderr, "internal: ") ;
  for ( i = 0 ; i < nf ; i ++ ) fprintf(stderr, " %+lg", f[2*ne+i]) ;
  fprintf(stderr, "\n") ;
  fprintf(stderr, "planar:   ") ;
  for ( i = 0 ; i < nf ; i ++ ) fprintf(stderr, " %+lg", g[i]) ;
  fprintf(stderr, "\n") ;
  fprintf(stderr, "error:    ") ;
  for ( i = 0 ; i < nf ; i ++ ) fprintf(stderr, " %+lg", fabs(f[i]-g[i])) ;
  fprintf(stderr, "\n") ;
  fprintf(stderr, "source:   ") ;
  for ( i = 0 ; i < 3 ; i ++ ) fprintf(stderr, " %+lg", fs[i]) ;
  fprintf(stderr, "\n") ;
  fprintf(stderr, "          ") ;
  for ( i = 0 ; i < 3 ; i ++ ) fprintf(stderr, " %+lg", gs[i]) ;
  fprintf(stderr, "\n") ;
  fprintf(stderr, "          ") ;
  for ( i = 0 ; i < 3 ; i ++ ) fprintf(stderr, " %+lg", fabs(gs[i]-fs[i])) ;
  fprintf(stderr, "\n") ;
  
  memset(fs, 0, ne*sizeof(gdouble)) ;
  memset(gs, 0, ne*sizeof(gdouble)) ;

  for ( i = 0 ; i < ne ; i ++ ) {
    fs[0] += f[nf+i]*1.0 ;
    fs[1] += f[nf+i]*xe[i*xstr+0] ;
    fs[2] += f[nf+i]*xe[i*xstr+1] ;
    gs[0] += g[nf+i]*1.0 ;
    gs[1] += g[nf+i]*xe[i*xstr+0] ;
    gs[2] += g[nf+i]*xe[i*xstr+1] ;
  }

  fprintf(stderr, "double layer\n") ;
  fprintf(stderr, "expansion:") ;
  for ( i = 0 ; i < nf ; i ++ ) fprintf(stderr, " %+lg", f[nf+i]) ;
  fprintf(stderr, "\n") ;
  fprintf(stderr, "internal: ") ;
  for ( i = 0 ; i < nf ; i ++ ) fprintf(stderr, " %+lg", f[2*ne+nf+i]) ;
  fprintf(stderr, "\n") ;
  fprintf(stderr, "planar:   ") ;
  for ( i = 0 ; i < nf ; i ++ ) fprintf(stderr, " %+lg", g[nf+i]) ;
  fprintf(stderr, "\n") ;
  fprintf(stderr, "error:    ") ;
  for ( i = 0 ; i < nf ; i ++ )
    fprintf(stderr, " %+lg", fabs(f[nf+i]-g[nf+i])) ;
  fprintf(stderr, "\n") ;
  fprintf(stderr, "source:   ") ;
  for ( i = 0 ; i < 3 ; i ++ ) fprintf(stderr, " %+lg", fs[i]) ;
  fprintf(stderr, "\n") ;
  fprintf(stderr, "          ") ;
  for ( i = 0 ; i < 3 ; i ++ ) fprintf(stderr, " %+lg", gs[i]) ;
  fprintf(stderr, "\n") ;
  fprintf(stderr, "          ") ;
  for ( i = 0 ; i < 3 ; i ++ ) fprintf(stderr, " %+lg", fabs(gs[i]-fs[i])) ;
  fprintf(stderr, "\n") ;

  return 0 ;
}

static gint target_specific_test(gdouble *xe, gint xstr, gint ne,
				 gint nq, gint N,
				 gdouble *x0, gdouble s0, gdouble t0,
				 gdouble rc, gint depth, gdouble tol,
				 gint nx)

{
  gdouble n[3], J, x[3], f[64], g[64], *q, fs[8], gs[8], t, w ;
  gdouble rcopt ;
  gint order, i, nf, fstr, Nopt, sopt, smax ;
  
  nf = ne ; fstr = ne ; smax = 10 ;

  qbx_quadrature_select(nq, &q, &order) ;

  w = sqrt(4*qbx_element_area(xe, xstr, ne, q, nq)/M_PI) ;
  qbx_element_point_3d(xe, xstr, ne, s0, t0, x, n, &J, NULL) ;
  qbx_quadrature_optimal_points(w, 0.5*w/(1 << smax), 4.0*w/(1 << smax),
				16, order, nq, 20, smax, tol,
				x, n, x,
				&rcopt, &Nopt, &sopt) ;
  
  fprintf(stderr, "target specific self-point test\n") ;
  fprintf(stderr, "===============================\n") ;
  fprintf(stderr, "quadrature: %d points, %dth order\n", nq, order) ;
  fprintf(stderr, "subdivision depth: %d %d\n", depth, sopt) ;
  fprintf(stderr, "N: %d %d\n", N, Nopt) ;
  fprintf(stderr, "rc: %lg %lg\n", rc, rcopt) ;
  /* fprintf(stderr, "xc: %lg %lg %lg\n", xc[0], xc[1], xc[2]) ; */
  /* fprintf(stderr, "n:  %lg %lg %lg\n", n[0], n[1], n[2]) ; */
  /* fprintf(stderr, "c:  %lg %lg %lg\n", c[0], c[1], c[2]) ; */

  fprintf(stderr, "evaluating integrals, t=%lg\n",
	  (t = g_timer_elapsed(timer, NULL))) ;
  memset(f, 0, 64*sizeof(gdouble)) ;
  
  qbx_laplace_ts_integrate(xe, xstr, ne, q, nq, order, rc, N,
			   s0, t0, TRUE, &(f[0]), fstr, depth, tol, w) ;
  for ( i = 0 ; i < ne ; i ++ ) f[ne+i] *= -1 ;
  qbx_laplace_ts_integrate(xe, xstr, ne, q, nq, order, rc, N,
			   s0, t0, FALSE, &(f[0]), fstr, depth, tol, w) ;
  for ( i = 0 ; i < 2*ne ; i ++ ) f[i] *= 0.5 ;

  fprintf(stderr, "integrals evaluated, t=%lg (%lg)\n",
	  g_timer_elapsed(timer, NULL), g_timer_elapsed(timer, NULL)-t) ;

  qbx_element_point_3d(xe, xstr, ne, s0, t0, x, n, &J, NULL) ;
  newman_tri_shape(x, &(xe[xstr*0]), &(xe[xstr*1]), &(xe[xstr*2]), NULL, 0,
		   &(g[0]), &(g[ne])) ;
  for ( i = 0 ; i < 2*ne ; i ++ ) g[i] *= -1.0/4.0/M_PI ;
  if ( s0 != G_MAXDOUBLE ) {
    qbx_element_shape_3d(ne, s0, t0, &(g[nf]), NULL, NULL, NULL, NULL, NULL) ;
    for ( i = 0 ; i < ne ; i ++ ) g[nf+i] *= 0.5 ;
  }

  memset(fs, 0, ne*sizeof(gdouble)) ;
  memset(gs, 0, ne*sizeof(gdouble)) ;

  for ( i = 0 ; i < ne ; i ++ ) {
    fs[0] += f[i]*1.0 ;
    fs[1] += f[i]*xe[i*xstr+0] ;
    fs[2] += f[i]*xe[i*xstr+1] ;
    fs[3] += f[2*ne+i]*1.0 ;
    fs[4] += f[2*ne+i]*xe[i*xstr+0] ;
    fs[5] += f[2*ne+i]*xe[i*xstr+1] ;
    gs[0] += g[i]*1.0 ;
    gs[1] += g[i]*xe[i*xstr+0] ;
    gs[2] += g[i]*xe[i*xstr+1] ;
  }
  
  fprintf(stderr, "single layer average\n") ;
  fprintf(stderr, "expansion:") ;
  for ( i = 0 ; i < nf ; i ++ ) fprintf(stderr, " %+lg", f[i]) ;
  fprintf(stderr, "\n") ;
  /* fprintf(stderr, "internal: ") ; */
  /* for ( i = 0 ; i < nf ; i ++ ) fprintf(stderr, " %+lg", f[2*ne+i]) ; */
  /* fprintf(stderr, "\n") ; */
  fprintf(stderr, "planar:   ") ;
  for ( i = 0 ; i < nf ; i ++ ) fprintf(stderr, " %+lg", g[i]) ;
  fprintf(stderr, "\n") ;
  fprintf(stderr, "error:    ") ;
  for ( i = 0 ; i < nf ; i ++ ) fprintf(stderr, " %+lg", fabs(f[i]-g[i])) ;
  fprintf(stderr, "\n") ;
  fprintf(stderr, "source:   ") ;
  for ( i = 0 ; i < 3 ; i ++ ) fprintf(stderr, " %+lg", fs[i]) ;
  fprintf(stderr, "\n") ;
  fprintf(stderr, "          ") ;
  for ( i = 0 ; i < 3 ; i ++ ) fprintf(stderr, " %+lg", gs[i]) ;
  fprintf(stderr, "\n") ;
  fprintf(stderr, "          ") ;
  for ( i = 0 ; i < 3 ; i ++ ) fprintf(stderr, " %+lg", fabs(gs[i]-fs[i])) ;
  fprintf(stderr, "\n") ;
  
  memset(fs, 0, ne*sizeof(gdouble)) ;
  memset(gs, 0, ne*sizeof(gdouble)) ;

  for ( i = 0 ; i < ne ; i ++ ) {
    fs[0] += f[nf+i]*1.0 ;
    fs[1] += f[nf+i]*xe[i*xstr+0] ;
    fs[2] += f[nf+i]*xe[i*xstr+1] ;
    gs[0] += g[nf+i]*1.0 ;
    gs[1] += g[nf+i]*xe[i*xstr+0] ;
    gs[2] += g[nf+i]*xe[i*xstr+1] ;
  }

  fprintf(stderr, "double layer jump\n") ;
  fprintf(stderr, "expansion:") ;
  for ( i = 0 ; i < nf ; i ++ ) fprintf(stderr, " %+lg", f[nf+i]) ;
  fprintf(stderr, "\n") ;
  /* fprintf(stderr, "internal: ") ; */
  /* for ( i = 0 ; i < nf ; i ++ ) fprintf(stderr, " %+lg", f[2*ne+nf+i]) ; */
  /* fprintf(stderr, "\n") ; */
  fprintf(stderr, "exact:    ") ;
  for ( i = 0 ; i < nf ; i ++ ) fprintf(stderr, " %+lg", g[nf+i]) ;
  fprintf(stderr, "\n") ;
  fprintf(stderr, "error:    ") ;
  for ( i = 0 ; i < nf ; i ++ )
    fprintf(stderr, " %+lg", fabs(f[nf+i]-g[nf+i])) ;
  fprintf(stderr, "\n") ;
  fprintf(stderr, "source:   ") ;
  for ( i = 0 ; i < 3 ; i ++ ) fprintf(stderr, " %+lg", fs[i]) ;
  fprintf(stderr, "\n") ;
  fprintf(stderr, "          ") ;
  for ( i = 0 ; i < 3 ; i ++ ) fprintf(stderr, " %+lg", gs[i]) ;
  fprintf(stderr, "\n") ;
  fprintf(stderr, "          ") ;
  for ( i = 0 ; i < 3 ; i ++ ) fprintf(stderr, " %+lg", fabs(gs[i]-fs[i])) ;
  fprintf(stderr, "\n") ;

  return 0 ;
}

static gint weight_test(gdouble *xe, gint xstr, gint ne,
			gint nq, gint N,
			gdouble *x0, gdouble s0, gdouble t0,
			gdouble rc, gint depth, gdouble tol,
			gint nx)

{
  gdouble n[3], J, x[3], f[64], g[64], *q, fs[8], gs[8], t, w, sc ;
  gdouble rcopt, K[453*453], wt[453], L[16] ;
  gint order, i, nf, j, fstr, Nopt, sopt, smax, Nk ;
  
  nf = ne ; fstr = ne ; smax = 10 ;

  qbx_quadrature_select(nq, &q, &order) ;

  w = sqrt(4*qbx_element_area(xe, xstr, ne, q, nq)/M_PI) ;
  qbx_element_point_3d(xe, xstr, ne, s0, t0, x, n, &J, NULL) ;
  qbx_quadrature_optimal_points(w, 0.5*w/(1 << smax), 4.0*w/(1 << smax),
				16, order, nq, 20, smax, tol,
				x, n, x,
				&rcopt, &Nopt, &sopt) ;
  
  fprintf(stderr, "Koornwinder weighted self-point test\n") ;
  fprintf(stderr, "====================================\n") ;
  fprintf(stderr, "quadrature: %d points, %dth order\n", nq, order) ;
  fprintf(stderr, "subdivision depth: %d %d\n", depth, sopt) ;
  fprintf(stderr, "N: %d %d\n", N, Nopt) ;
  fprintf(stderr, "rc: %lg %lg\n", rc, rcopt) ;

  qbx_quadrature_select(nq, &q, &order) ;

  fprintf(stderr, "generating Koornwinder matrix, t=%lg\n",
	  (t = g_timer_elapsed(timer, NULL))) ;
  Nk = qbx_koornwinder_interp_matrix(q, nq, K) ;
  
  qbx_laplace_ts_self_weights(xe, xstr, ne, q, nq, order, K, Nk, rc, N,
			      s0, t0, TRUE, wt, &(wt[nq]), depth, tol, w) ;
  sc = -1.0 ; i = 1 ;
  blaswrap_dscal(nq, sc, &(wt[nq]), i) ;

  qbx_laplace_ts_self_weights(xe, xstr, ne, q, nq, order, K, Nk, rc, N,
			      s0, t0, FALSE, wt, &(wt[nq]), depth, tol, w) ;
  sc = 0.5 ;
  blaswrap_dscal(nq, sc, wt, i) ;
  blaswrap_dscal(nq, sc, &(wt[nq]), i) ;

  fprintf(stderr, "evaluating integrals, t=%lg\n",
	  (t = g_timer_elapsed(timer, NULL))) ;
  memset(f, 0, 64*sizeof(gdouble)) ;
  memset(g, 0, 64*sizeof(gdouble)) ;
  
  qbx_laplace_ts_integrate(xe, xstr, ne, q, nq, order, rc, N,
			   s0, t0, TRUE, &(f[0]), fstr, depth, tol, w) ;
  memset(f, 0, ne*sizeof(QBX_REAL)) ;
  sc = -1.0 ; i = 1 ;
  blaswrap_dscal(nq, sc, &(f[ne]), i) ;
  qbx_laplace_ts_integrate(xe, xstr, ne, q, nq, order, rc, N,
  			   s0, t0, FALSE, &(f[0]), fstr, depth, tol, w) ;
  sc = 0.5 ;
  blaswrap_dscal(nq, sc, &(f[ne]), i) ;

  for ( i = 0 ; i < nq ; i ++ ) {
    qbx_element_shape_3d(ne, q[3*i+0], q[3*i+1], L,
			 NULL, NULL, NULL, NULL, NULL) ;
    for ( j = 0 ; j < ne ; j ++ ) {
      g[   j] += wt[0*nq+i]*L[j] ;
      g[ne+j] += wt[1*nq+i]*L[j] ;
    }
  }

  fprintf(stderr, "integrals evaluated, t=%lg (%lg)\n",
	  g_timer_elapsed(timer, NULL), g_timer_elapsed(timer, NULL)-t) ;
  
  qbx_element_shape_3d(ne, s0, t0, gs, NULL, NULL, NULL, NULL, NULL) ;
  sc = 0.5 ; i = 1 ;
  blaswrap_dscal(ne, sc, gs, i) ;

  fprintf(stderr, "single layer average\n") ;
  fprintf(stderr, "expansion:") ;
  for ( i = 0 ; i < nf ; i ++ ) fprintf(stderr, " %+lg", f[i]) ;
  fprintf(stderr, "\n") ;
  fprintf(stderr, "weights:  ") ;
  for ( i = 0 ; i < nf ; i ++ ) fprintf(stderr, " %+lg", g[i]) ;
  fprintf(stderr, "\n") ;
  fprintf(stderr, "error:    ") ;
  for ( i = 0 ; i < nf ; i ++ ) fprintf(stderr, " %+lg", fabs(f[i]-g[i])) ;
  fprintf(stderr, "\n") ;
  
  fprintf(stderr, "double layer jump\n") ;
  fprintf(stderr, "expansion:") ;
  for ( i = 0 ; i < nf ; i ++ ) fprintf(stderr, " %+lg", f[nf+i]) ;
  fprintf(stderr, "\n") ;
  fprintf(stderr, "jump:     ") ;
  for ( i = 0 ; i < nf ; i ++ ) fprintf(stderr, " %+lg", gs[  i]) ;
  fprintf(stderr, "\n") ;
  fprintf(stderr, "weights:  ") ;
  for ( i = 0 ; i < nf ; i ++ ) fprintf(stderr, " %+lg", g[nf+i]) ;
  fprintf(stderr, "\n") ;
  fprintf(stderr, "error:    ") ;
  for ( i = 0 ; i < nf ; i ++ )
    fprintf(stderr, " %+lg", fabs(g[nf+i]-gs[i])) ;
  fprintf(stderr, "\n") ;

  return 0 ;
}

static gint off_weight_test(gdouble *xe, gint xstr, gint ne,
			    gint nq, gint N,
			    gdouble *x0,
			    gint depth, gdouble tol,
			    gint nx)

{
  gdouble n[3], J, x[3], f[64], g[64], *q, fs[8], gs[8], t, w, sc ;
  gdouble rcopt, K[453*453], wd[453], ws[453], L[16], s[453], work[1024] ;
  gint order, i, nf, j, fstr, Nopt, sopt, smax, Nk, isrc = 1, i1 = 1 ;
  
  nf = ne ; fstr = ne ;

  qbx_quadrature_select(nq, &q, &order) ;

  /*set up a sample source vector*/
  for ( i = 0 ; i < nq ; i ++ ) {
    qbx_element_shape_3d(ne, q[3*i+0], q[3*i+1], L,
			 NULL, NULL, NULL, NULL, NULL) ;
    s[i] = L[isrc] ;
  }

  w = sqrt(4*qbx_element_area(xe, xstr, ne, q, nq)/M_PI) ;
  
  fprintf(stderr, "Koornwinder weighted off point test\n") ;
  fprintf(stderr, "===================================\n") ;
  fprintf(stderr, "quadrature: %d points, %dth order\n", nq, order) ;
  fprintf(stderr, "subdivision depth: %d\n", depth) ;
  fprintf(stderr, "target: %lg %lg %lg\n", x0[0], x0[1], x0[2]) ;
  
  qbx_quadrature_select(nq, &q, &order) ;

  fprintf(stderr, "generating Koornwinder matrix, t=%lg\n",
	  (t = g_timer_elapsed(timer, NULL))) ;
  Nk = qbx_koornwinder_interp_matrix(q, nq, K) ;

  qbx_laplace_weights(xe, xstr, ne, q, nq, order, K, Nk,
		      x0, ws, wd, depth, tol, w) ;
  sc = 1.0 ;
  f[0] = blaswrap_ddot(nq, ws, i1, s, i1) ;
  f[1] = blaswrap_ddot(nq, wd, i1, s, i1) ;

  memset(g, 0, 64*sizeof(gdouble)) ;
  qbx_triangle_laplace_quad(xe, xstr, ne, x0, q, nq, order,
			    0, depth, tol,
			    &(g[0]), 1, &(g[ne]), 1, work) ;

  fprintf(stderr, "single layer\n") ;
  fprintf(stderr, "%lg %lg (%e)\n",
	  f[0], g[isrc], fabs(f[0]-g[isrc])) ;
  fprintf(stderr, "double layer\n") ;
  fprintf(stderr, "%lg %lg (%e)\n",
	  f[1], g[ne+isrc], fabs(f[1]-g[ne+isrc])) ;

  return 0 ;
}

static gint matrix_weight_test(gdouble *xe, gint xstr, gint ne,
			       gint nq, gint N,
			       gdouble *x0, gdouble s0, gdouble t0,
			       gdouble rc, gint depth, gdouble tol,
			       gint nx)

{
  gdouble n[3], J, x[3], f[64], g[64], *q, t, w, al, bt ;
  gdouble rcopt, K[453*453], wt[453], L[16], As[453*453], Ad[453*453] ;
  gdouble s[453], fs[453], fd[453], ed ;
  gint order, i, nf, j, fstr, Nopt, sopt, smax, Nk, isrc, one = 1 ;
  
  smax = 10 ;

  isrc = 2 ;
  
  qbx_quadrature_select(nq, &q, &order) ;

  w = sqrt(4*qbx_element_area(xe, xstr, ne, q, nq)/M_PI) ;
  qbx_element_point_3d(xe, xstr, ne, s0, t0, x, n, &J, NULL) ;
  qbx_quadrature_optimal_points(w, 0.5*w/(1 << smax), 4.0*w/(1 << smax),
				16, order, nq, 20, smax, tol,
				x, n, x,
				&rcopt, &Nopt, &sopt) ;
  
  fprintf(stderr, "element layer potential matrices test\n") ;
  fprintf(stderr, "=====================================\n") ;
  fprintf(stderr, "quadrature: %d points, %dth order\n", nq, order) ;
  fprintf(stderr, "subdivision depth: %d %d\n", depth, sopt) ;
  fprintf(stderr, "N: %d %d\n", N, Nopt) ;
  fprintf(stderr, "rc: %lg %lg\n", rc, rcopt) ;

  /*set up a sample source vector*/
  for ( i = 0 ; i < nq ; i ++ ) {
    qbx_element_shape_3d(ne, q[3*i+0], q[3*i+1], L,
			 NULL, NULL, NULL, NULL, NULL) ;
    s[i] = L[isrc] ;
  }
  
  qbx_quadrature_select(nq, &q, &order) ;

  fprintf(stderr, "generating Koornwinder matrix, t=%lg\n",
	  (t = g_timer_elapsed(timer, NULL))) ;
  Nk = qbx_koornwinder_interp_matrix(q, nq, K) ;

  fprintf(stderr, "generating layer potential matrices, t=%lg\n",
	  (t = g_timer_elapsed(timer, NULL))) ;

  qbx_triangle_laplace_ts_self_matrix(xe, xstr, ne, q, nq, order, K, Nk, rc, N,
				      As, Ad, depth, tol, w) ;

  fprintf(stderr, "layer potential matrices generated, t=%lg\n",
	  (t = g_timer_elapsed(timer, NULL))) ;

  fprintf(stderr, "matrix evaluation of layer potentials, t=%lg\n",
	  (t = g_timer_elapsed(timer, NULL))) ;
  al = 1.0 ; bt = 0.0 ;
  blaswrap_dgemv(FALSE, nq, nq, al, As, nq, s, one, bt, fs, one) ;
  blaswrap_dgemv(FALSE, nq, nq, al, Ad, nq, s, one, bt, fd, one) ;

  ed = 0.0 ;
  for ( i = 0 ; i < nq ; i ++ ) {
    fprintf(stdout,
	    "%lg %lg %lg %lg %lg (%e)\n",
	    q[3*i+0], q[3*i+1], fs[i], fd[i], 0.5*s[i], fabs(fd[i]-0.5*s[i])) ;
    ed = MAX(ed, fabs(fd[i]-0.5*s[i])) ;
  }

  fprintf(stderr, "double layer maximum error: %lg\n", ed) ;
  
  return 0 ;
}

static gint off_test(gdouble *xe, gint xstr, gint ne,
		     gint nq, gint N,
		     gdouble *x0, gdouble s0, gdouble t0,
		     gdouble rc, gint depth, gdouble tol,
		     gint nx)

{
  gdouble work[8192], n[3], J, x[3], f[64], g[64], *q, fs[8], gs[8] ;
  gint order, i, nf ;
  
  nf = ne ;

  qbx_quadrature_select(nq, &q, &order) ;

  fprintf(stderr, "triangle off-element point test\n") ;
  fprintf(stderr, "===============================\n") ;

  fprintf(stderr, "quadrature: %d points, %dth order\n", nq, order) ;
  fprintf(stderr, "subdivision depth: %d\n", depth) ;
  fprintf(stderr, "N: %d\n", N) ;

  qbx_element_point_3d(xe, xstr, ne, s0, t0, x, n, &J, NULL) ;
  qbx_vector_shift(x,x,n,rc) ;

  fprintf(stderr, "x: %lg %lg %lg\n", x[0], x[1], x[2]) ;
  
  fprintf(stderr, "evaluating integrals, t=%lg\n",
	  g_timer_elapsed(timer, NULL)) ;
  qbx_triangle_laplace_quad(xe, xstr, ne, x, q, nq, order,
			    N, depth, tol,
			    &(f[0]), 1, &(f[ne]), 1, work) ;
  fprintf(stderr, "integrals evaluated, t=%lg\n",
	  g_timer_elapsed(timer, NULL)) ;

  newman_tri_shape(x, &(xe[xstr*0]), &(xe[xstr*1]), &(xe[xstr*2]), NULL, 0,
		   &(g[0]), &(g[ne])) ;
  for ( i = 0 ; i < 2*ne ; i ++ ) g[i] *= -1.0/4.0/M_PI ;

  memset(fs, 0, ne*sizeof(gdouble)) ;
  memset(gs, 0, ne*sizeof(gdouble)) ;

  for ( i = 0 ; i < ne ; i ++ ) {
    fs[0] += f[i]*1.0 ;
    fs[1] += f[i]*xe[i*xstr+0] ;
    fs[2] += f[i]*xe[i*xstr+1] ;
    gs[0] += g[i]*1.0 ;
    gs[1] += g[i]*xe[i*xstr+0] ;
    gs[2] += g[i]*xe[i*xstr+1] ;
  }
  
  fprintf(stderr, "single layer\n") ;
  fprintf(stderr, "expansion:") ;
  for ( i = 0 ; i < nf ; i ++ ) fprintf(stderr, " %+lg", f[i]) ;
  fprintf(stderr, "\n") ;
  fprintf(stderr, "planar:   ") ;
  for ( i = 0 ; i < nf ; i ++ ) fprintf(stderr, " %+lg", g[i]) ;
  fprintf(stderr, "\n") ;
  fprintf(stderr, "error:    ") ;
  for ( i = 0 ; i < nf ; i ++ ) fprintf(stderr, " %+lg", fabs(f[i]-g[i])) ;
  fprintf(stderr, "\n") ;
  fprintf(stderr, "source:   ") ;
  for ( i = 0 ; i < 3 ; i ++ ) fprintf(stderr, " %+lg", fs[i]) ;
  fprintf(stderr, "\n") ;
  fprintf(stderr, "          ") ;
  for ( i = 0 ; i < 3 ; i ++ ) fprintf(stderr, " %+lg", gs[i]) ;
  fprintf(stderr, "\n") ;
  fprintf(stderr, "          ") ;
  for ( i = 0 ; i < 3 ; i ++ ) fprintf(stderr, " %+lg", fabs(gs[i]-fs[i])) ;
  fprintf(stderr, "\n") ;
  
  memset(fs, 0, ne*sizeof(gdouble)) ;
  memset(gs, 0, ne*sizeof(gdouble)) ;

  for ( i = 0 ; i < ne ; i ++ ) {
    fs[0] += f[nf+i]*1.0 ;
    fs[1] += f[nf+i]*xe[i*xstr+0] ;
    fs[2] += f[nf+i]*xe[i*xstr+1] ;
    gs[0] += g[nf+i]*1.0 ;
    gs[1] += g[nf+i]*xe[i*xstr+0] ;
    gs[2] += g[nf+i]*xe[i*xstr+1] ;
  }

  fprintf(stderr, "double layer\n") ;
  fprintf(stderr, "expansion:") ;
  for ( i = 0 ; i < nf ; i ++ ) fprintf(stderr, " %+lg", f[nf+i]) ;
  fprintf(stderr, "\n") ;
  fprintf(stderr, "planar:   ") ;
  for ( i = 0 ; i < nf ; i ++ ) fprintf(stderr, " %+lg", g[nf+i]) ;
  fprintf(stderr, "\n") ;
  fprintf(stderr, "error:    ") ;
  for ( i = 0 ; i < nf ; i ++ ) fprintf(stderr, " %+lg", fabs(f[nf+i]-g[nf+i])) ;
  fprintf(stderr, "\n") ;
  fprintf(stderr, "source:   ") ;
  for ( i = 0 ; i < 3 ; i ++ ) fprintf(stderr, " %+lg", fs[i]) ;
  fprintf(stderr, "\n") ;
  fprintf(stderr, "          ") ;
  for ( i = 0 ; i < 3 ; i ++ ) fprintf(stderr, " %+lg", gs[i]) ;
  fprintf(stderr, "\n") ;
  fprintf(stderr, "          ") ;
  for ( i = 0 ; i < 3 ; i ++ ) fprintf(stderr, " %+lg", fabs(gs[i]-fs[i])) ;
  fprintf(stderr, "\n") ;

  return 0 ;
}

static gint element_closest_point_test(gdouble *xe, gint xstr, gint ne,
				       gdouble s0,  gdouble t0,
				       gdouble *x0, gdouble rc)

{
  gdouble xc[3], n[3], xf[3], J, s, t, xn[3], r[3], kn, km ;
  gint ni ;
  
  fprintf(stderr, "closest point test\n") ;
  fprintf(stderr, "===================\n") ;

  fprintf(stderr, "element: %d nodes\n", ne) ;
  fprintf(stderr, "rc: %lg\n", rc) ;
  fprintf(stderr, "s,t: %lg %lg\n", s0, t0) ;

  /*closest point on element for test purposes*/
  qbx_element_point_3d(xe, xstr, ne, s0, t0, xc, n, &J, NULL) ;  
  fprintf(stderr, "xc: %lg %lg %lg\n", xc[0], xc[1], xc[2]) ;
  fprintf(stderr, "n:  %lg %lg %lg\n", n[0], n[1], n[2]) ;
  qbx_triangle_curvature(xe, xstr, ne, s0, t0, &kn, &km) ;
  fprintf(stderr, "kn: %lg\n", kn) ;
  
  xf[0] = xc[0] + n[0]*rc ;
  xf[1] = xc[1] + n[1]*rc ;
  xf[2] = xc[2] + n[2]*rc ;

  fprintf(stderr, "xf: %lg %lg %lg\n", xf[0], xf[1], xf[2]) ;

  /* s = 1/3.0 ; t = 1/3.0 ; */
  /* s = t = 0.0 ; */
  s = 1.0 ; t = 0.0 ;
  /* s = s0 ; t = t0+0.000001 ; */
  ni = qbx_element_nearest_point(xe, xstr, ne, xf, &s, &t, xn, 1e-9, 256, &kn) ;
  
  qbx_element_point_3d(xe, xstr, ne, s, t, xc, n, &J, NULL) ;  

  r[0] = xf[0] - xc[0] ; r[1] = xf[1] - xc[1] ; r[2] = xf[2] - xc[2] ;
  fprintf(stderr, "s,t: %lg %lg\n", s, t) ;
  fprintf(stderr, "ni:  %d\n", ni) ;
  fprintf(stderr, "x0:  %lg %lg %lg\n", xc[0], xc[1], xc[2]) ;
  fprintf(stderr, "kn:  %lg\n", kn) ;
  fprintf(stderr, "n:   %lg %lg %lg\n", n[0], n[1], n[2]) ;
  fprintf(stderr, "r:   %lg %lg %lg\n", r[0], r[1], r[2]) ;
  fprintf(stderr, "r.n: %lg (%lg)\n",
	  qbx_vector_scalar(r,n), qbx_vector_length(r)) ;
  
  return 0 ;
}
				       
static gint koornwinder_test(gint N, gdouble u, gdouble v)

{
  gdouble Knm[1024] ;
  gint n, m, str ;

  fprintf(stderr, "koornwinder test\n") ;
  fprintf(stderr, "================\n") ;

  fprintf(stderr, "N = %d\n", N) ;
  fprintf(stderr, "(u,v) = (%lg,%lg)\n", u, v) ;

  str = 3 ;
  
  qbx_koornwinder_nm(N, u, v, str, 8192, Knm) ;

  for ( n = 0 ; n <= N ; n ++ ) {
    for ( m = 0 ; m <= n ; m ++ ) {
      fprintf(stdout,
	      "%d %d %1.16e %1.16e %1.16e\n",
	      n, m, u, v, Knm[str*(n*(n+1)/2+m)]) ;
    }
  }

  return 0 ;
}

static gint koornwinder_orthogonality_test(gint N)

{
  gint nq, order, i, n1, m1, n2, m2, idx1, idx2, str ;
  gdouble s, t, Knm[4096], *q, I, w, tol ;

  str = 3 ;
  
  tol = 1e-12 ;
  
  fprintf(stderr, "koornwinder orthogonality test\n") ;
  fprintf(stderr, "==============================\n") ;

  fprintf(stderr, "N = %d\n", N) ;

  nq = 453 ;
  
  qbx_quadrature_select(nq, &q, &order) ;

  for ( n1 = 0 ; n1 <= N ; n1 ++ ) {
    for ( m1 = 0 ; m1 <= n1 ; m1 ++ ) {
      for ( n2 = 0 ; n2 <= N ; n2 ++ ) {
	for ( m2 = 0 ; m2 <= n2 ; m2 ++ ) {
	  idx1 = n1*(n1+1)/2 + m1 ; 
	  idx2 = n2*(n2+1)/2 + m2 ; 
  
	  I = 0.0 ;
	  for ( i = 0 ; i < nq ; i ++ ) {
	    s = q[3*i+0] ; 
	    t = q[3*i+1] ;
	    w = q[3*i+2] ; 
	    
	    qbx_koornwinder_nm(N, s, t, str, 8192, Knm) ;
	    
	    I += Knm[str*idx1]*Knm[str*idx2]*w ;
	  }

	  fprintf(stderr, "%d %d %d %d %e ", n1, m1, n2, m2, I) ;
	  if ( n1 != n2 || m1 != m2 ) {
	    if ( fabs(I) > tol ) 
	      fprintf(stderr, "FAIL\n") ;
	    else
	      fprintf(stderr, "PASS\n") ;
	  } else {
	    if ( fabs(I-1.0) > tol ) 
	      fprintf(stderr, "FAIL\n") ;
	    else
	      fprintf(stderr, "PASS\n") ;
	  }	  
	}
      }
    }
  }
  return 0 ;
}

static gint koornwinder_interpolation_test(gint N)

{
  gint nq, order, i, i1 = 1 ;
  gdouble s, t, Knm[32768], A[4*65536], *q, f, fr, fi[512], al, bt, c[512] ;
  
  fprintf(stderr, "koornwinder interpolation test\n") ;
  fprintf(stderr, "==============================\n") ;

  nq = 85 ; 
  
  qbx_quadrature_select(nq, &q, &order) ;

  N = qbx_koornwinder_interp_matrix(q, nq, A) ;

  fprintf(stderr, "Knm N max: %d\n", N) ;
  
  for ( i = 0 ; i < nq ; i ++ ) {
    s = q[i*3+0] ; t = q[i*3+1] ;
    /* fi[i] = 3.0*s*t - t*t ; */
    fi[i] = sin(2.0*M_PI*s*t/8) ;
  }

  al = 1.0 ; bt = 0.0 ;
  blaswrap_dgemv(FALSE, nq, nq, al, A, nq, fi, i1, bt, c, i1) ;

  for ( s = 0 ; s <= 1.0 ; s += 0.1 ) {
    for ( t = 0 ; t <= 1.0-s ; t += 0.1 ) {
      qbx_koornwinder_nm(N, s, t, 1, nq, Knm) ;

      f = blaswrap_ddot(nq, c, i1, Knm, i1) ;
      /* fr = 3.0*s*t - t*t ; */
      fr = sin(2.0*M_PI*s*t/8) ;
  
      fprintf(stderr, "%lg %lg %lg %lg (%lg)\n", s, t, fr, f, fabs(fr-f)) ;
    }
  }

  return 0 ;
}

static gint blas_tests(gint N)

{
  gint i, j, stra, strx, stry, nr, nc, i1 = 1 ;
  gdouble A[8192], x[256], y[256], yref[256], al, bt, d, dref ;

  fprintf(stderr, "BLAS test\n") ;
  fprintf(stderr, "=========\n") ;

  nr = N ; nc = nr + 3 ;

  al = -1.23 ; bt = 0.7 ;
  
  stra = nc+5 ; strx = 2 ; stry = 5 ;
  /* stra = nc ; strx = 1 ; stry = 1 ; */
  
  fprintf(stderr, "A: %dx%d, stride %d\n", nr, nc, stra) ;
  fprintf(stderr, "x: %d elements, stride %d\n", nc, strx) ;
  fprintf(stderr, "y: %d elements, stride %d\n", nr, stry) ;
  fprintf(stderr, "al = %lg; bt = %lg\n", al, bt) ;
  
  for ( i = 0 ; i < nr ; i ++ ) {
    for ( j = 0 ; j < nc ; j ++ ) {
      A[i*stra+j] = (gdouble)(i+3)*(j+1) ;
    }
  }
  for ( j = 0 ; j < nc ; j ++ ) {
    y[stry*j] = -(gdouble)((j-0.3)*0.7) ;
    x[strx*j] = -(gdouble)((j+5.3)*0.7) ;
  }

  for ( i = 0 ; i < nr ; i ++ ) {
    yref[i] = bt*y[i*stry] ;
    for ( j = 0 ; j < nc ; j ++ ) {
      yref[i] += al*A[i*nc+j]*x[j*strx] ;
    }
  }

  blaswrap_dgemv(FALSE, nr, nc, al, A, stra, x, strx, bt, y, stry) ;

  for ( i = 0 ; i < nr ; i ++ ) {
    fprintf(stderr, "%lg ", yref[i]) ;
  }
  fprintf(stderr, "\n") ;
  for ( i = 0 ; i < nr ; i ++ ) {
    fprintf(stderr, "%lg ", y[stry*i]) ;
  }
  fprintf(stderr, "\n") ;
  dref = 0.0 ;
  for ( i = 0 ; i < nr ; i ++ ) {
    fprintf(stderr, "%lg ", fabs(y[stry*i]-yref[i])) ;
    dref += y[stry*i]*yref[i] ;
  }
  fprintf(stderr, "\n") ;

  d = blaswrap_ddot(nr, y, stry, yref, i1) ;

  fprintf(stderr, "dot: %lg %lg (%lg)\n", d, dref, fabs(d-dref)) ;
  
  return 0 ;
}

static gint foreach_test_func(gdouble s, gdouble t, gdouble w,
			       gdouble *x, gdouble *y, gdouble *n,
			       gint N,
			       gdouble *fq, gint fstr, gint nf,
			       gpointer data[])
{
  gdouble *xe = data[0] ;
  gint xstr = *((gint *)(data[1])) ;
  gint ne = *((gint *)(data[2])) ;
  gdouble L[32], Ls[32], Lt[32], J, yc[3], nc[3], err ;
  
  /* fprintf(stdout, "%e %e %e\n", y[0], y[1], y[2]) ; */
  qbx_element_shape_3d(ne, s, t, L, Ls, Lt, NULL, NULL, NULL) ;
  qbx_element_point_interp_3d(xe, xstr, ne, L, Ls, Lt, yc, nc, &J, NULL) ;

  err = qbx_vector_distance(y, yc) ;  
  fprintf(stdout, "%e %e %e %e ", y[0], y[1], y[2], err) ;
  err = qbx_vector_distance(n, nc) ;
  fprintf(stdout, "%e\n", err) ;

  fq[0] += w ;
  
  return 0 ;
}

static gint foreach_test(gdouble *xe, gint xstr, gint ne,
			  gint nq, gint N,
			  gdouble *x0, gdouble s0, gdouble t0,
			  gdouble rc, gint depth, gdouble tol,
			  gint nx)

{
  gdouble st[32], *q, f[512], w, xt[3], J, n[3], c[3], rcopt ;
  gint oq, smax, sopt, Nopt ;
  gpointer data[4] ;
  
  fprintf(stderr, "adaptive subdivision test\n") ;
  fprintf(stderr, "=========================\n") ;

  smax = depth ;
  
  st[0*2+0] = 0.0  ; st[0*2+1] = 0.0 ; 
  st[1*2+0] = 1.0  ; st[1*2+1] = 0.0 ; 
  st[2*2+0] = 0.0  ; st[2*2+1] = 1.0 ; 
  st[3*2+0] = 0.5  ; st[3*2+1] = 0.0 ; 
  st[4*2+0] = 0.5  ; st[4*2+1] = 0.5 ; 
  st[5*2+0] = 0.0  ; st[5*2+1] = 0.5 ; 

  /* nq = 25 ; */
  qbx_quadrature_select(nq, &q, &oq) ;
  w = sqrt(4*qbx_element_area(xe, xstr, ne, q, nq)/M_PI) ;

  /*surface target point*/
  qbx_element_point_3d(xe, xstr, ne, s0, t0, xt, n, &J, NULL) ;
  /* x[0] = xc[0] ; x[1] = xc[1] ; x[2] = xc[2] ;     */

  c[0] = xt[0] + rc*n[0] ; 
  c[1] = xt[1] + rc*n[1] ; 
  c[2] = xt[2] + rc*n[2] ; 

  qbx_quadrature_optimal_points(w, 0.5*w/(1 << smax), 4.0*w/(1 << smax),
				16, oq, nq, 20, smax, tol,
				xt, n, xt,
				&rcopt, &Nopt, &sopt) ;

  fprintf(stderr, "quadrature: %d points, %dth order\n", nq, oq) ;
  fprintf(stderr, "subdivision depth: %d %d\n", depth, sopt) ;
  fprintf(stderr, "tol: %e\n", tol) ;
  fprintf(stderr, "N:  %d %d\n", N, Nopt) ;
  fprintf(stderr, "rc: %lg %lg\n", rc, rcopt) ;
  fprintf(stderr, "xt: %lg %lg %lg\n", xt[0], xt[1], xt[2]) ;
  fprintf(stderr, "n:  %lg %lg %lg\n", n[0], n[1], n[2]) ;
  fprintf(stderr, "c:  %lg %lg %lg\n", c[0], c[1], c[2]) ;

  f[0] = 0.0 ;
  
  data[0] = xe ; data[1] = &xstr ; data[2] = &ne ;
  qbx_triangle_adaptive(foreach_test_func,
			xe, xstr, ne, xe, xstr, st, w, c, N, depth,
			q, nq, oq, tol,
			data) ;
			
  fprintf(stderr, "A:  %lg %lg (%e)\n",
	  f[0], qbx_element_area(xe, xstr, ne, q, nq),
	  fabs(f[0] - qbx_element_area(xe, xstr, ne, q, nq))) ;
  
  return 0 ;
}

gint adaptive_quad_func(gdouble s, gdouble t, gdouble w,
			gdouble *y, gdouble *n,
			gdouble *quad, gint nq, gpointer data[])
{
  gdouble *x = data[0] ;
  gdouble R, dR ;

  R = qbx_vector_distance(x, y) ;

  dR = qbx_vector_diff_scalar(x, y, n)/R/R ;  
  w *= 0.25*M_1_PI/R ;
  
  quad[0] += w*(1 - s - t)  ;
  quad[1] += w*(    s    )  ;
  quad[2] += w*(        t)  ;
  quad[3] += w*(1 - s - t)*dR ;
  quad[4] += w*(    s    )*dR ;
  quad[5] += w*(        t)*dR ;
  
  return 0 ;
}

gint adaptive_quad_test(gdouble *xe, gint xstr, gint ne,
			gint nq, gint N,
			gdouble *x,
			gint depth, gdouble tol,
			gint nx)

{
  gdouble st[32], *q, f[512], J, n[3], g[8], dg[8], work[1024], t ;
  gint oq, nc, i ;
  gpointer data[4] ;

  nc = 2*ne ;
  
  fprintf(stderr, "adaptive quadrature test\n") ;
  fprintf(stderr, "========================\n") ;
  fprintf(stderr, "x = %lg %lg %lg\n", x[0], x[1], x[2]) ;
  fprintf(stderr, "tol = %lg\n", tol) ;
  
  st[0*2+0] = 0.0  ; st[0*2+1] = 0.0 ; 
  st[1*2+0] = 1.0  ; st[1*2+1] = 0.0 ; 
  st[2*2+0] = 0.0  ; st[2*2+1] = 1.0 ; 
  st[3*2+0] = 0.5  ; st[3*2+1] = 0.0 ; 
  st[4*2+0] = 0.5  ; st[4*2+1] = 0.5 ; 
  st[5*2+0] = 0.0  ; st[5*2+1] = 0.5 ; 

  qbx_quadrature_select(nq, &q, &oq) ;

  data[0] = x ;
  fprintf(stderr, "starting integration, t=%lg\n",
	  t = g_timer_elapsed(timer, NULL)) ;  
  qbx_adaptive_quad_tri(xe, xstr, ne, st, q, nq, adaptive_quad_func,
			f, nc, tol, depth, TRUE, data) ;
  fprintf(stderr, "integration completed, t=%lg (%lg)\n",
	  g_timer_elapsed(timer, NULL), g_timer_elapsed(timer,NULL) - t) ;
  newman_tri_shape(x, &(xe[xstr*0]), &(xe[xstr*1]), &(xe[xstr*2]), NULL, 0,
		   &(g[0]), &(g[ne])) ;
  for ( i = 0 ; i < 2*ne ; i ++ ) g[i] *= -1.0/4.0/M_PI ;

  fprintf(stderr, "single layer\n") ;
  fprintf(stderr, "  adaptive: %lg %lg %lg\n", f[0], f[1], f[2]) ;
  fprintf(stderr, "  exact:    %lg %lg %lg\n", g[0], g[1], g[2]) ;
  fprintf(stderr, "  error:    %lg %lg %lg\n",
	  fabs(f[0]-g[0]), fabs(f[1]-g[1]), fabs(f[2]-g[2])) ;
  fprintf(stderr, "double layer\n") ;
  fprintf(stderr, "  adaptive: %lg %lg %lg\n", f[3], f[4], f[5]) ;
  fprintf(stderr, "  exact:    %lg %lg %lg\n", g[3], g[4], g[5]) ;
  fprintf(stderr, "  error:    %lg %lg %lg\n",
	  fabs(f[3]-g[3]), fabs(f[4]-g[4]), fabs(f[5]-g[5])) ;
  
  return 0 ;
}

gint main(gint argc, gchar **argv)

{
  gdouble xe[256], tol, rc, s0, t0, x0[3] ;
  gint ne, nq, depth, N, nx, xstr, test ;
  FILE *input ;
  gchar ch, *progname ;

  timer = g_timer_new() ;

  test = -1 ;
  depth = 0 ; tol = 1e-6 ; N = 8 ; nq = 54 ; rc = 0.05 ; nx = 33 ;
  s0 = t0 = G_MAXDOUBLE ;
  
  progname = g_strdup(g_path_get_basename(argv[0])) ;

  input = stdin ;

  while ( (ch = getopt(argc, argv, "d:e:n:N:q:r:s:T:t:w:")) != EOF ) {
    switch (ch) {
    default: g_assert_not_reached() ; break ;
    case 'd': depth = atoi(optarg) ; break ;
    case 'e': tol  = atof(optarg) ; break ;
    case 'n': nx  = atoi(optarg) ; break ;
    case 'N': N   = atoi(optarg) ; break ;
    case 'q': nq  = atoi(optarg) ; break ;
    case 'r': rc  = atof(optarg) ; break ;
    case 's': s0  = atof(optarg) ; break ;
    case 'T': test = parse_test(optarg) ; break ;      
    case 't': t0  = atof(optarg) ; break ;
    }
  }

  if ( s0 != G_MAXDOUBLE && t0 == G_MAXDOUBLE ) t0 = s0 ;
  if ( t0 != G_MAXDOUBLE && s0 == G_MAXDOUBLE ) s0 = t0 ;
  
  if ( test == 4 ) {
    koornwinder_test(N, s0, t0) ;

    return 0 ;
  }

  if ( test == 5 ) {
    koornwinder_orthogonality_test(N) ;

    return 0 ;
  }

  if ( test == 6 ) {
    blas_tests(N) ;

    return 0 ;
  }
  
  if ( test == 7 ) {
    koornwinder_interpolation_test(N) ;

    return 0 ;
  }

  read_element(input, xe, &ne, &xstr) ;
  fscanf(input, "%lg %lg %lg", &(x0[0]), &(x0[1]), &(x0[2])) ;
	 
  fprintf(stderr, "%s: %d node element\n", progname, ne) ;

  if ( test == 0 ) {
    if ( s0 != G_MAXDOUBLE ) 
      planar_triangle_shape_test(xe, xstr, ne, nq, N, NULL, s0, t0, rc,
				 depth, tol, nx) ;
    else
      planar_triangle_shape_test(xe, xstr, ne, nq, N, x0, s0, t0, rc,
				 depth, tol, nx) ;
    return 0 ;
  }

  if ( test == 1 ) {
    element_closest_point_test(xe, xstr, ne, s0, t0, x0, rc) ;
    
    return 0 ;
  }

  if ( test == 2 ) {
    self_test(xe, xstr, ne, nq, N, NULL, s0, t0, rc, depth, tol, nx) ;

    return 0 ;
  }

  if ( test == 3 ) {
    off_test(xe, xstr, ne, nq, N, NULL, s0, t0, rc, depth, tol, nx) ;

    return 0 ;
  }

  if ( test == 8 ) {
    foreach_test(xe, xstr, ne, nq, N, NULL, s0, t0, rc, depth, tol, nx) ;

    return 0 ;
  }
  
  if ( test == 9 ) {
    target_specific_test(xe, xstr, ne, nq, N, NULL, s0, t0, rc,
			 depth, tol, nx) ;

    return 0 ;
  }

  if ( test == 10 ) {
    weight_test(xe, xstr, ne, nq, N, NULL, s0, t0, rc,
		depth, tol, nx) ;

    return 0 ;
  }

  if ( test == 11 ) {
    matrix_weight_test(xe, xstr, ne, nq, N, NULL, s0, t0, rc,
		       depth, tol, nx) ;

    return 0 ;
  }

  if ( test == 12 ) {
    off_weight_test(xe, xstr, ne, nq, N, x0, depth, tol, nx) ;

    return 0 ;
  }

  if ( test == 13 ) {
    adaptive_quad_test(xe, xstr, ne, nq, N, x0, depth, tol, nx) ;

    return 0 ;
  }
  
  return 0 ;
}

