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

/*
  Target-specific QBX based on the methods of Siegel and Tornberg,
  https://doi.org/10.1016/j.jcp.2018.03.006
*/


static gint point_interp(QBX_REAL *xb, QBX_REAL *fb, gint i,
			 QBX_REAL *xe, gint xstr, gint ne,
			 QBX_REAL *fe, gint fstr, gint nf,
			 QBX_REAL s, QBX_REAL t)

{
  QBX_REAL n[3], J, L[16], dLds[16], dLdt[16] ;
  gint j, k ;
  
  QBX_FUNCTION_NAME(qbx_element_shape_3d)(ne, s, t, L, dLds, dLdt,
					  NULL, NULL, NULL) ;
  QBX_FUNCTION_NAME(qbx_element_point_interp_3d)(xe, xstr, ne,
						 L, dLds, dLdt, &(xb[xstr*i]),
						 n, &J, NULL) ;

  for ( j = 0 ; j < ne ; j ++ ) fb[fstr*i+j] = 0.0 ;
  for ( j = 0 ; j < ne ; j ++ ) {
    for ( k = 0 ; k < nf ; k ++ ) {
      fb[fstr*i+k] += L[j]*fe[j*fstr+k] ;
    }
  }
  
  return 0 ;
}
  
static gint point_copy(QBX_REAL *xb, QBX_REAL *fb, gint i,
		       QBX_REAL *xe, gint xstr, gint ne,
		       QBX_REAL *fe, gint fstr, gint nf, gint j)

{
  gint k ;
  
  xb[xstr*i+0] = xe[xstr*j+0] ;
  xb[xstr*i+1] = xe[xstr*j+1] ;
  xb[xstr*i+2] = xe[xstr*j+2] ;

  for ( k = 0 ; k < nf ; k ++ ) fb[fstr*j+k] = fe[fstr*i+k] ;

  return 0 ;
}

static gint triangle_divide_loop(QBX_REAL *xe, gint xstr, gint ne,
				 QBX_REAL *fe, gint fstr, gint nf, gint d,
				 QBX_REAL *xl, QBX_REAL *fl)

{
  g_assert(ne == 3 || ne == 6) ;

  if ( ne == 3 ) {
    switch ( d ) {
    case 0:
      point_copy(xl, fl, 0, xe, xstr, ne, fe, fstr, nf, 0) ;
      point_interp(xl, fl, 1, xe, xstr, ne, fe, fstr, nf, 0.5, 0.0) ;
      point_interp(xl, fl, 2, xe, xstr, ne, fe, fstr, nf, 0.0, 0.5) ;
      break ;
    case 1:
      point_interp(xl, fl, 0, xe, xstr, ne, fe, fstr, nf, 0.5, 0.0) ;
      point_copy(xl, fl, 1, xe, xstr, ne, fe, fstr, nf, 1) ;
      point_interp(xl, fl, 2, xe, xstr, ne, fe, fstr, nf, 0.5, 0.5) ;
      break ;
    case 2:
      point_interp(xl, fl, 0, xe, xstr, ne, fe, fstr, nf, 0.0, 0.5) ;
      point_interp(xl, fl, 1, xe, xstr, ne, fe, fstr, nf, 0.5, 0.5) ;
      point_copy(xl, fl, 2, xe, xstr, ne, fe, fstr, nf, 2) ;
      break ;
    case 3:
      point_interp(xl, fl, 0, xe, xstr, ne, fe, fstr, nf, 0.5, 0.0) ;
      point_interp(xl, fl, 1, xe, xstr, ne, fe, fstr, nf, 0.5, 0.5) ;
      point_interp(xl, fl, 2, xe, xstr, ne, fe, fstr, nf, 0.0, 0.5) ;
      break ;
    default: g_assert_not_reached() ; break ;
    }
    return 0 ;
  }
  
  switch ( d ) {
  case 0:
    point_copy(xl, fl, 0, xe, xstr, ne, fe, fstr, nf, 0) ;
    point_interp(xl, fl, 1, xe, xstr, ne, fe, fstr, nf, 0.5, 0.0) ;
    point_interp(xl, fl, 2, xe, xstr, ne, fe, fstr, nf, 0.0, 0.5) ;
    point_interp(xl, fl, 3, xe, xstr, ne, fe, fstr, nf, 0.25, 0.0 ) ;
    point_interp(xl, fl, 4, xe, xstr, ne, fe, fstr, nf, 0.25, 0.25) ;
    point_interp(xl, fl, 5, xe, xstr, ne, fe, fstr, nf, 0.0,  0.25) ;
    break ;
  case 1:
    point_interp(xl, fl, 0, xe, xstr, ne, fe, fstr, nf, 0.5, 0.0) ;
    point_copy(xl, fl, 1, xe, xstr, ne, fe, fstr, nf, 1) ;
    point_interp(xl, fl, 2, xe, xstr, ne, fe, fstr, nf, 0.5, 0.5) ;
    point_interp(xl, fl, 3, xe, xstr, ne, fe, fstr, nf, 0.75, 0.0 ) ;
    point_interp(xl, fl, 4, xe, xstr, ne, fe, fstr, nf, 0.75, 0.25) ;
    point_interp(xl, fl, 5, xe, xstr, ne, fe, fstr, nf, 0.5 , 0.25) ;
    break ;
  case 2:
    point_interp(xl, fl, 0, xe, xstr, ne, fe, fstr, nf, 0.0, 0.5) ;
    point_interp(xl, fl, 1, xe, xstr, ne, fe, fstr, nf, 0.5, 0.5) ;
    point_copy(xl, fl, 2, xe, xstr, ne, fe, fstr, nf, 2) ;
    point_interp(xl, fl, 3, xe, xstr, ne, fe, fstr, nf, 0.25, 0.5) ;
    point_interp(xl, fl, 4, xe, xstr, ne, fe, fstr, nf, 0.25, 0.75) ;
    point_interp(xl, fl, 5, xe, xstr, ne, fe, fstr, nf, 0.0, 0.75) ;
    break ;
  case 3:
    point_interp(xl, fl, 0, xe, xstr, ne, fe, fstr, nf, 0.5, 0.0) ;
    point_interp(xl, fl, 1, xe, xstr, ne, fe, fstr, nf, 0.5, 0.5) ;
    point_interp(xl, fl, 2, xe, xstr, ne, fe, fstr, nf, 0.0, 0.5) ;
    point_interp(xl, fl, 3, xe, xstr, ne, fe, fstr, nf, 0.5, 0.25) ;
    point_interp(xl, fl, 4, xe, xstr, ne, fe, fstr, nf, 0.25, 0.5) ;
    point_interp(xl, fl, 5, xe, xstr, ne, fe, fstr, nf, 0.25, 0.25) ;
    break ;
    default: g_assert_not_reached() ; break ;
  }

  return 0 ;
}

gint QBX_FUNCTION_NAME(qbx_triangle_adaptive)(
#ifdef QBX_SINGLE_PRECISION
					      qbx_quadrature_func_f_t func,
#else /*QBX_SINGLE_PRECISION*/
					      qbx_quadrature_func_t func,
#endif /*QBX_SINGLE_PRECISION*/
					      QBX_REAL *xt, gint tstr, gint ne,
					      QBX_REAL *xd, gint dstr,
					      QBX_REAL *st,
					      QBX_REAL w,
					      QBX_REAL *xc,
					      gint N,
					      gint depth,
					      QBX_REAL *q, gint nq, gint oq,
					      QBX_REAL tol,
					      QBX_REAL *f, gint fstr, gint nf,
					      gpointer data)

{
  QBX_REAL std[32], xb[32], L[32], Ls[32], Lt[32], rp, err, s, t, wt ;
  QBX_REAL J, y[3], n[3] ;
  gint i ;

  g_assert(func != NULL) ;
  
  /*check error estimate for this triangle (rough estimate but quicker
    than finding the nearest point accurately)*/
  rp = 1e6 ;
  for ( i = 0 ; i < ne ; i ++ ) {
    rp = MIN(rp, qbx_vector_distance2(xc, &(xd[i*dstr]))) ;
  }
  rp = SQRT(rp)*0.5 ;
  err = QBX_FUNCTION_NAME(qbx_quadrature_error)(N, w, oq, rp, rp) ;
  err *= 1 << depth ;
    
  if ( depth == 0 || err < tol ) {
    for ( i = 0 ; i < nq ; i ++ ) {
      s = q[3*i+0] ; t = q[3*i+1] ; wt = q[3*i+2] ;
      QBX_FUNCTION_NAME(qbx_element_shape_3d)(ne, s, t, L, Ls, Lt,
					      NULL, NULL, NULL) ;
      QBX_FUNCTION_NAME(qbx_element_point_interp_3d)(xd, dstr, ne,
						     L, Ls, Lt, y, n, &J,
						     NULL) ;
      wt *= J ;
      /*coordinates on top-level element */
      s = qbx_ddot(&ne, L, qbx_1i, st, qbx_2i) ;
      t = qbx_ddot(&ne, L, qbx_1i, &(st[1]), qbx_2i) ;
      func(s, t, wt, xc, y, n, N, f, fstr, nf, data) ;
    }

    return 0 ;
  }
  
  for ( i = 0 ; i < 4 ; i ++ ) {
    triangle_divide_loop(xd, dstr, ne, st, 2, 2, i, xb, std) ;

    qbx_triangle_adaptive(func, xt, tstr, ne, xb, dstr, std, 0.5*w,
			  xc, N, depth-1, q, nq, oq, tol, f, fstr, nf, data) ;
  }

  return 0 ;
}
