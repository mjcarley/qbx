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
#include <math.h>
#include <string.h>

#include <glib.h>

#include <blaswrap.h>

#include <qbx.h>

#include "qbx-private.h"

/*
  Target-specific QBX based on the methods of Siegel and Tornberg,
  https://doi.org/10.1016/j.jcp.2018.03.006
*/

static gint QBX_FUNCTION_NAME(laplace_quad)(QBX_REAL s, QBX_REAL t,
					    QBX_REAL w,
					    QBX_REAL *x, QBX_REAL *y,
					    QBX_REAL *n,
					    gint N,
					    gpointer data[])
{
  gint ne = *((gint *)(data[QBX_DATA_NUMBER])) ;
  QBX_REAL rc = *((QBX_REAL *)(data[QBX_DATA_RADIUS])) ;
  QBX_REAL *n0 = data[QBX_DATA_NORMAL] ;
  QBX_REAL *f  = data[QBX_DATA_WEIGHTS_S] ;  
  QBX_REAL R, L[16] ;
  QBX_REAL Cth, S2, Pnm1, Pn, dP, g1, g2 ;
  gint nn, j ;

  /*shape function on top-level element*/
  QBX_FUNCTION_NAME(qbx_element_shape_3d)(ne, s, t, L, NULL, NULL,
					  NULL, NULL, NULL) ;
  R = qbx_vector_distance(x, y) ;
  Cth = qbx_vector_diff_scalar(x, y, n0)/R ;
  S2 = 1.0/(1.0 - Cth*Cth) ;

  g1 =  qbx_vector_diff_scalar(x, y, n)/R/R ;
  g2 = -qbx_vector_scalar(n,n0)/R + Cth*g1 ;
  
  Pnm1 = 1.0 ; Pn = Cth ; dP = 0.0 ;

  w *= 0.25*M_1_PI/R ;
  nn = 0 ;
  for ( j = 0 ; j < ne ; j ++ ) {
    f[   j] += Pnm1                    *w*L[j] ;
    f[ne+j] += (g1*Pnm1*(nn+1) + g2*dP)*w*L[j] ;
  }

  /*this also works but is slower for the small values of ne in this
    function*/
  /* wt = w*Pnm1 ; */
  /* qbx_daxpy(&ne, &wt, L, qbx_1i, f, qbx_1i) ; */
  /* wt = w*(g1*Pnm1*(nn+1) + g2*dP) ; */
  /* qbx_daxpy(&ne, &wt, L, qbx_1i, &(f[ne]), qbx_1i) ; */

  
  w *= rc/R ;
  nn = 1 ;
  dP = (Pnm1 - Cth*Pn)*S2*nn ;
  for ( j = 0 ; j < ne ; j ++ ) {
    f[   j] += Pn                    *w*L[j] ;
    f[ne+j] += (g1*Pn*(nn+1) + g2*dP)*w*L[j] ;
  }
  
  for ( nn = 1 ; nn < N ; nn ++ ) {
    w *= rc/R ;
    dP = Pn ;
    Pn = ((2.0*nn+1)*Pn*Cth - (QBX_REAL)nn*Pnm1)/(nn+1) ;
    Pnm1 = dP ;
    dP = (Pnm1 - Cth*Pn)*S2*(nn+1) ;
    for ( j = 0 ; j < ne ; j ++ ) {
      f[   j] += Pn                    *w*L[j] ;
      f[ne+j] += (g1*Pn*(nn+2) + g2*dP)*w*L[j] ;
    }
  }
  
  return 0 ;
}

gint QBX_FUNCTION_NAME(qbx_laplace_ts_integrate)(QBX_REAL *xe,
						 gint xstr, gint ne,
						 QBX_REAL *q,
						 gint nq, gint oq,
						 QBX_REAL rc, gint N,
						 QBX_REAL s0, QBX_REAL t0,
						 gboolean in,
						 QBX_REAL *f,
						 gint str,
						 gint depth,
						 QBX_REAL tol,
						 QBX_REAL w)

{
  QBX_REAL xt[3], n[3], J, L[32], Ls[32], Lt[32], c[3] ;
  QBX_REAL st[] = {0.0, 0.0, 1.0, 0.0, 0.0, 1.0,
		   0.5, 0.0, 0.5, 0.5, 0.0, 0.5} ;
  gpointer data[QBX_DATA_WIDTH] ;
  
  QBX_FUNCTION_NAME(qbx_element_shape_3d)(ne, s0, t0, L, Ls, Lt,
					  NULL, NULL, NULL) ;
  QBX_FUNCTION_NAME(qbx_element_point_interp_3d)(xe, xstr, ne,
						 L, Ls, Lt, xt, n, &J,
						 NULL) ;
  if ( in ) {
    qbx_vector_shift(c,xt,n,-rc) ;
    n[0] = -n[0] ; n[1] = -n[1] ; n[2] = -n[2] ; 
  } else {
    qbx_vector_shift(c,xt,n, rc) ;
  }
  
  data[QBX_DATA_ELEMENT] = xe ; 
  data[QBX_DATA_STRIDE]  = &xstr ;
  data[QBX_DATA_NUMBER]  = &(ne) ;
  data[QBX_DATA_RADIUS]  = &(rc) ; 
  data[QBX_DATA_NORMAL]  = n ;
  data[QBX_DATA_WEIGHTS_S]  = f ;
  
#ifdef QBX_SINGLE_PRECISION
  qbx_quadrature_func_f_t func = laplace_quad_f ;
#else /*QBX_SINGLE_PRECISION*/
  qbx_quadrature_func_t func = laplace_quad ;
#endif /*QBX_SINGLE_PRECISION*/

  QBX_FUNCTION_NAME(qbx_triangle_adaptive)(func,
					   xe, xstr, ne, xe, xstr, st,
					   w, c, N, depth, q, nq, oq, tol,
					   data) ;
  return 0 ;
}

static gint QBX_FUNCTION_NAME(laplace_quad_weights)(QBX_REAL s, QBX_REAL t,
						    QBX_REAL w,
						    QBX_REAL *x, QBX_REAL *y,
						    QBX_REAL *n,
						    gint N,
						    gpointer data[])
{
  QBX_REAL rc = *((QBX_REAL *)(data[QBX_DATA_RADIUS])) ;
  QBX_REAL *n0 = data[QBX_DATA_NORMAL] ;
  QBX_REAL *Kq  = data[QBX_DATA_MATRIX] ;
  QBX_REAL *Knm = data[QBX_DATA_KNM] ;
  gint nq = *((gint *)(data[QBX_DATA_NKNM])) ;
  gint Nk = *((gint *)(data[QBX_DATA_ORDER_K])) ;
  QBX_REAL *ws  = data[QBX_DATA_WEIGHTS_S] ;  
  QBX_REAL *wd = data[QBX_DATA_WEIGHTS_D] ;  
  QBX_REAL R, wt ;
  QBX_REAL Cth, S2, Pnm1, Pn, dP, g1, g2 ;
  QBX_REAL d1 = 1.0 ;
  gint nn, j, i1 = 1 ;

  /*Koornwinder polynomials at evaluation point*/
  qbx_koornwinder_nm(Nk, s, t, 1, nq, Knm) ;
  R = qbx_vector_distance(x, y) ;
  Cth = qbx_vector_diff_scalar(x, y, n0)/R ;
  S2 = 1.0/(1.0 - Cth*Cth) ;

  g1 =  qbx_vector_diff_scalar(x, y, n)/R/R ;
  g2 = -qbx_vector_scalar(n,n0)/R + Cth*g1 ;
  
  Pnm1 = 1.0 ; Pn = Cth ; dP = 0.0 ;

  w *= 0.25*M_1_PI/R ;
  nn = 0 ;
  wt = Pnm1*w ;
  blaswrap_dgemv(TRUE, nq, nq, wt, Kq, nq, Knm, i1, d1, ws, i1) ;
  wt = (g1*Pnm1*(nn+1) + g2*dP)*w ;
  blaswrap_dgemv(TRUE, nq, nq, wt, Kq, nq, Knm, i1, d1, wd, i1) ;
  
  w *= rc/R ;
  nn = 1 ;
  dP = (Pnm1 - Cth*Pn)*S2*nn ;
  wt = Pn*w ;
  blaswrap_dgemv(TRUE, nq, nq, wt, Kq, nq, Knm, i1, d1, ws, i1) ;
  wt = (g1*Pn*(nn+1) + g2*dP)*w ;
  blaswrap_dgemv(TRUE, nq, nq, wt, Kq, nq, Knm, i1, d1, wd, i1) ;
  
  for ( nn = 1 ; nn < N ; nn ++ ) {
    w *= rc/R ;
    dP = Pn ;
    Pn = ((2.0*nn+1)*Pn*Cth - (QBX_REAL)nn*Pnm1)/(nn+1) ;
    Pnm1 = dP ;
    dP = (Pnm1 - Cth*Pn)*S2*(nn+1) ;
    wt = Pn*w ;
    blaswrap_dgemv(TRUE, nq, nq, wt, Kq, nq, Knm, i1, d1, ws, i1) ;
    wt = (g1*Pn*(nn+2) + g2*dP)*w ;
    blaswrap_dgemv(TRUE, nq, nq, wt, Kq, nq, Knm, i1, d1, wd, i1) ;
  }
  
  return 0 ;
}

gint QBX_FUNCTION_NAME(qbx_laplace_ts_self_weights)(QBX_REAL *xe,
						    gint xstr, gint ne,
						    QBX_REAL *q,
						    gint nq, gint oq,
						    QBX_REAL *Kq, gint Nk,
						    QBX_REAL rc, gint N,
						    QBX_REAL s0, QBX_REAL t0,
						    gboolean in,
						    QBX_REAL *ws, QBX_REAL *wd,
						    gint depth,
						    QBX_REAL tol,
						    QBX_REAL w)

{
  QBX_REAL xt[3], n[3], J, L[32], Ls[32], Lt[32], c[3] ;
  QBX_REAL st[] = {0.0, 0.0, 1.0, 0.0, 0.0, 1.0,
		   0.5, 0.0, 0.5, 0.5, 0.0, 0.5} ;
  QBX_REAL Knm[512] ;
  gpointer data[QBX_DATA_WIDTH] ;
  
  QBX_FUNCTION_NAME(qbx_element_shape_3d)(ne, s0, t0, L, Ls, Lt,
					  NULL, NULL, NULL) ;
  QBX_FUNCTION_NAME(qbx_element_point_interp_3d)(xe, xstr, ne,
						 L, Ls, Lt, xt, n, &J,
						 NULL) ;
  if ( in ) {
    qbx_vector_shift(c, xt, n,-rc) ;
    n[0] = -n[0] ; n[1] = -n[1] ; n[2] = -n[2] ; 
  } else {
    qbx_vector_shift(c, xt, n, rc) ;
  }
  
  data[QBX_DATA_ELEMENT]   = xe ; 
  data[QBX_DATA_STRIDE]    = &xstr ;
  data[QBX_DATA_NUMBER]    = &(ne) ;
  data[QBX_DATA_RADIUS]    = &(rc) ; 
  data[QBX_DATA_NORMAL]    = n ;
  data[QBX_DATA_MATRIX]    = Kq ;
  data[QBX_DATA_KNM]       = Knm ;
  data[QBX_DATA_NKNM]      = &nq ;
  data[QBX_DATA_ORDER_K]   = &Nk ;
  data[QBX_DATA_WEIGHTS_S] = ws ;
  data[QBX_DATA_WEIGHTS_D] = wd ;
  
#ifdef QBX_SINGLE_PRECISION
  qbx_quadrature_func_f_t func = laplace_quad_weights_f ;
#else /*QBX_SINGLE_PRECISION*/
  qbx_quadrature_func_t func = laplace_quad_weights ;
#endif /*QBX_SINGLE_PRECISION*/

  QBX_FUNCTION_NAME(qbx_triangle_adaptive)(func,
					   xe, xstr, ne, xe, xstr, st,
					   w, c, N, depth, q, nq, oq, tol,
					   data) ;
  return 0 ;
}

gint QBX_FUNCTION_NAME(qbx_triangle_laplace_ts_self_matrix)(QBX_REAL *xe,
						       gint xstr, gint ne,
						       QBX_REAL *q,
						       gint nq, gint oq,
						       QBX_REAL *Kq, gint Nk,
						       QBX_REAL rc, gint N,
						       QBX_REAL *ws,
						       QBX_REAL *wd,
						       gint depth,
						       QBX_REAL tol,
						       QBX_REAL w)

/*
  generate single and double layer matrices for element so that the
  single layer potential is evaluated by the matrix vector product

  [ws] s

  with s the source strengths at the quadrature nodes of the rule q on
  the element xe  
*/

{
  gint i, one = 1, smax = 20, Nopt, sopt ;
  QBX_REAL sc, s0, t0, J, x[3], n[3], rcopt ;

  memset(ws, 0, nq*nq*sizeof(QBX_REAL)) ;
  memset(wd, 0, nq*nq*sizeof(QBX_REAL)) ;

  for ( i = 0 ; i < nq ; i ++ ) {
    s0 = q[i*3+0] ; t0 = q[i*3+1] ;

    qbx_element_point_3d(xe, xstr, ne, s0, t0, x, n, &J, NULL) ;
    qbx_quadrature_optimal_points(w, 0.5*w/(1 << smax), 4.0*w/(1 << smax),
				  16, oq, nq, 20, smax, tol,
				  x, n, x,
				  &rcopt, &Nopt, &sopt) ;
    
    QBX_FUNCTION_NAME(qbx_laplace_ts_self_weights)(xe, xstr, ne, q, nq, oq,
						   Kq, Nk, rcopt, Nopt,
						   s0, t0, TRUE,
						   &(ws[i*nq]), &(wd[i*nq]),
						   sopt, tol, w) ;
    sc = -1.0 ;
    blaswrap_dscal(nq, sc, &(wd[i*nq]), one) ;

    QBX_FUNCTION_NAME(qbx_laplace_ts_self_weights)(xe, xstr, ne, q, nq, oq,
						   Kq, Nk, rcopt, Nopt,
						   s0, t0, FALSE,
						   &(ws[i*nq]), &(wd[i*nq]),
						   sopt, tol, w) ;
    sc = 0.5 ;
    blaswrap_dscal(nq, sc, &(ws[i*nq]), one) ;
    blaswrap_dscal(nq, sc, &(wd[i*nq]), one) ;
  }
  
  return 0 ;
}

gint QBX_FUNCTION_NAME(qbx_laplace_ts_off_weights)(QBX_REAL *xe,
						   gint xstr, gint ne,
						   QBX_REAL *q,
						   gint nq, gint oq,
						   QBX_REAL *Kq, gint Nk,
						   QBX_REAL rc, gint N,
						   QBX_REAL *xt, QBX_REAL *n,
						   gboolean in,
						   QBX_REAL *ws, QBX_REAL *wd,
						   gint depth,
						   QBX_REAL tol,
						   QBX_REAL w)

{
  QBX_REAL st[] = {0.0, 0.0, 1.0, 0.0, 0.0, 1.0,
		   0.5, 0.0, 0.5, 0.5, 0.0, 0.5} ;
  QBX_REAL Knm[512], c[3] ;
  gpointer data[QBX_DATA_WIDTH] ;

  if ( in ) {
    qbx_vector_shift(c, xt, n,-rc) ;
    n[0] = -n[0] ; n[1] = -n[1] ; n[2] = -n[2] ; 
  } else {
    qbx_vector_shift(c, xt, n, rc) ;
  }
  
  data[QBX_DATA_ELEMENT]   = xe ; 
  data[QBX_DATA_STRIDE]    = &xstr ;
  data[QBX_DATA_NUMBER]    = &(ne) ;
  data[QBX_DATA_RADIUS]    = &(rc) ; 
  data[QBX_DATA_NORMAL]    = n ;
  data[QBX_DATA_MATRIX]    = Kq ;
  data[QBX_DATA_KNM]       = Knm ;
  data[QBX_DATA_NKNM]      = &nq ;
  data[QBX_DATA_ORDER_K]   = &Nk ;
  data[QBX_DATA_WEIGHTS_S] = ws ;
  data[QBX_DATA_WEIGHTS_D] = wd ;
  
#ifdef QBX_SINGLE_PRECISION
  qbx_quadrature_func_f_t func = laplace_quad_weights_f ;
#else /*QBX_SINGLE_PRECISION*/
  qbx_quadrature_func_t func = laplace_quad_weights ;
#endif /*QBX_SINGLE_PRECISION*/

  QBX_FUNCTION_NAME(qbx_triangle_adaptive)(func,
					   xe, xstr, ne, xe, xstr, st,
					   w, c, N, depth, q, nq, oq, tol,
					   data) ;
  return 0 ;
}
