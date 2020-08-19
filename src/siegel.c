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

#include <qbx.h>

#include "qbx-private.h"

#define QBX_TS_DATA_WIDTH   8
#define QBX_TS_DATA_ELEMENT 0
#define QBX_TS_DATA_STRIDE  1
#define QBX_TS_DATA_NUMBER  2
#define QBX_TS_DATA_RADIUS  3
#define QBX_TS_DATA_NORMAL  4

/*
  Target-specific QBX based on the methods of Siegel and Tornberg,
  https://doi.org/10.1016/j.jcp.2018.03.006
*/

static gint QBX_FUNCTION_NAME(laplace_quad)(QBX_REAL s, QBX_REAL t,
					    QBX_REAL w,
					    QBX_REAL *x, QBX_REAL *y,
					    QBX_REAL *n,
					    gint N,
					    QBX_REAL *f, gint str,
					    gint nf,
					    gpointer data[])
{
  gint ne = *((gint *)(data[QBX_TS_DATA_NUMBER])) ;
  QBX_REAL rc = *((QBX_REAL *)(data[QBX_TS_DATA_RADIUS])) ;
  QBX_REAL *n0 = data[QBX_TS_DATA_NORMAL] ;
  QBX_REAL R, L[16] ;
  QBX_REAL Cth, S2, Pnm1, Pn, dP, g1, g2 ;
  gint nn, j ;

  /*shape function on top-level element*/
  QBX_FUNCTION_NAME(qbx_element_shape_3d)(ne, s, t, L, NULL, NULL,
					  NULL, NULL, NULL) ;
  R = qbx_vector_distance(x, y) ;
  Cth = qbx_vector_diff_scalar(x, y, n)/R ;
  S2 = 1.0 - Cth*Cth ;

  g1 =  qbx_vector_diff_scalar(x, y, n)/R/R ;
  g2 = -qbx_vector_scalar(n,n0)/R + Cth*g1 ;
  
  Pnm1 = 1.0 ; Pn = Cth ; dP = 0.0 ;

  w *= 0.25*M_1_PI/R ;
  nn = 0 ;
  for ( j = 0 ; j < ne ; j ++ ) {
    f[   j] += Pnm1                    *w*L[j] ;
    f[ne+j] += (g1*Pnm1*(nn+1) + g2*dP)*w*L[j] ;
  }

  w *= rc/R ;
  nn = 1 ;
  dP = (Pnm1 - Cth*Pn)/S2*nn ;
  for ( j = 0 ; j < ne ; j ++ ) {
    f[   j] += Pn                    *w*L[j] ;
    f[ne+j] += (g1*Pn*(nn+1) + g2*dP)*w*L[j] ;
  }
  
  for ( nn = 1 ; nn < N ; nn ++ ) {
    w *= rc/R ;
    dP = Pn ;
    Pn = (2.0*nn+1)/(nn+1)*Pn*Cth - (QBX_REAL)nn/(nn+1)*Pnm1 ;
    Pnm1 = dP ;
    dP = (Pnm1 - Cth*Pn)/S2*(nn+1) ;
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
						 gint nq, gint order,
						 QBX_REAL *xc,
						 QBX_REAL rc, gint N,
						 QBX_REAL s0, QBX_REAL t0,
						 QBX_REAL *f,
						 gint str,
						 gint depth,
						 QBX_REAL tol,
						 QBX_REAL w)

{
  QBX_REAL xt[3], n[3], J, L[32], Ls[32], Lt[32], c[3] ;
  QBX_REAL st[] = {0.0, 0.0, 1.0, 0.0, 0.0, 1.0,
		   0.5, 0.0, 0.5, 0.5, 0.0, 0.5} ;
  gpointer data[QBX_TS_DATA_WIDTH] ;
  
  QBX_FUNCTION_NAME(qbx_element_shape_3d)(ne, s0, t0, L, Ls, Lt,
					  NULL, NULL, NULL) ;
  QBX_FUNCTION_NAME(qbx_element_point_interp_3d)(xe, xstr, ne,
						 L, Ls, Lt, xt, n, &J,
						 NULL) ;
  qbx_vector_shift(c,xt,n,rc) ;

  data[QBX_TS_DATA_ELEMENT] = xe ; 
  data[QBX_TS_DATA_STRIDE]  = &xstr ;
  data[QBX_TS_DATA_NUMBER]  = &(ne) ;
  data[QBX_TS_DATA_RADIUS]  = &(rc) ; 
  data[QBX_TS_DATA_NORMAL]  = n ;
  
#ifdef QBX_SINGLE_PRECISION
  qbx_quadrature_func_f_t func = laplace_quad_f ;
#else /*QBX_SINGLE_PRECISION*/
  qbx_quadrature_func_t func = laplace_quad ;
#endif /*QBX_SINGLE_PRECISION*/

  QBX_FUNCTION_NAME(qbx_triangle_adaptive)(func,
					   xe, xstr, ne, xe, xstr, st,
					   w, c, N, depth, q, nq, order, tol,
					   f, str, ne, data) ;
  return 0 ;
}
