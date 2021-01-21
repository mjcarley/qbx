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

/* #define TRIANGLE_TRACE */

static void adaptive_quad_tri3(QBX_REAL *xe, gint xstr, QBX_REAL *st,
			       QBX_REAL *q, gint nq, 
#ifdef QBX_SINGLE_PRECISION
			       qbx_adaptive_func_f_t func,
#else /*QBX_SINGLE_PRECISION*/
			       qbx_adaptive_func_t func,
#endif /*QBX_SINGLE_PRECISION*/
			       QBX_REAL *quad, gint nc,
			       gpointer data)

{
  gint i ;
  QBX_REAL s, t, w, J, L[10], Ls[10], Lt[10], y[3], n[3] ;

  for ( i = 0 ; i < nq ; i ++ ) {
    s = q[3*i+0] ; t = q[3*i+1] ; w = q[3*i+2] ;
    qbx_shape_derivatives3(s, t, L, Ls, Lt) ;
    qbx_point_interp_jac3(xe, xstr, L[0], L[1], L[2],
			  Ls[0], Ls[1], Ls[2],
			  Lt[0], Lt[1], Lt[2],
			  y, n, J) ;
    w *= J ;
    s = L[0]*st[0] + L[1]*st[2] + L[2]*st[4] ;
    t = L[0]*st[1] + L[1]*st[3] + L[2]*st[5] ;
    func(s, t, w, y, n, quad, nc, data) ;
  }

  return ;
}

gint QBX_FUNCTION_NAME(qbx_adaptive_quad_tri)(QBX_REAL *xe, gint xstr, gint ne,
					      QBX_REAL *st,
					      QBX_REAL *q, gint nq,
#ifdef QBX_SINGLE_PRECISION
					      qbx_adaptive_func_f_t func,
#else /*QBX_SINGLE_PRECISION*/
					      qbx_adaptive_func_t func,
#endif /*QBX_SINGLE_PRECISION*/
					      QBX_REAL *quad, gint nc,
					      QBX_REAL tol, gint dmax,
					      gboolean init,
					      gpointer data)

{
  gint i ;
  QBX_REAL s, t, w, L[10], Ls[10], Lt[10], n[3], y[3], J ;
  QBX_REAL work[2048], *q0, *q1, *q2, *q3 ;
  QBX_REAL xe0[30], st0[20], xe1[30], st1[20], xe2[30], st2[20],
    xe3[30], st3[20] ;
  gboolean recurse ;
  
  /* fprintf(stderr, "depth %d; (%lg)\n", dmax, quad[0]) ; */

#ifdef TRIANGLE_TRACE
  fprintf(stdout, "%e %e %e %e %e %e %e %e %e\n",
	  xe[0*xstr+0], xe[0*xstr+1], xe[0*xstr+2], 
	  xe[1*xstr+0], xe[1*xstr+1], xe[1*xstr+2], 
	  xe[2*xstr+0], xe[2*xstr+1], xe[2*xstr+2]) ;
#endif
  
  if ( dmax == 0 ) return 0 ;

  g_assert(ne == 3 ) ;
  
  if ( init ) {
    memset(quad, 0, nc*sizeof(QBX_REAL)) ;
    adaptive_quad_tri3(xe, xstr, st, q, nq, func, quad, nc, data) ;
  }

  memset(work, 0, 4*nc*sizeof(QBX_REAL)) ;
  q0 = &(work[0]) ; q1 = &(q0[nc]) ; q2 = &(q1[nc]) ; q3 = &(q2[nc]) ;

  qbx_triangle_divide_loop30(xe, xstr, st, xe0, st0) ;
  adaptive_quad_tri3(xe0, xstr, st0, q, nq, func, q0, nc, data) ;
  qbx_triangle_divide_loop31(xe, xstr, st, xe1, st1) ;
  adaptive_quad_tri3(xe1, xstr, st1, q, nq, func, q1, nc, data) ;
  qbx_triangle_divide_loop32(xe, xstr, st, xe2, st2) ;
  adaptive_quad_tri3(xe2, xstr, st2, q, nq, func, q2, nc, data) ;
  qbx_triangle_divide_loop33(xe, xstr, st, xe3, st3) ;
  adaptive_quad_tri3(xe3, xstr, st3, q, nq, func, q3, nc, data) ;

  recurse = FALSE ;

  for ( i = 0 ; i < nc ; i ++ ) {
    if ( fabs(quad[i] - q0[i] - q1[i] - q2[i] - q3[i]) > tol ) {
      recurse = TRUE ; break ;
    }
    /* quad[i] = q0[i] + q1[i] + q2[i] + q3[i] ; */
  }

  if ( !recurse ) return 0 ;

  QBX_FUNCTION_NAME(qbx_adaptive_quad_tri)(xe0, xstr, ne, st0, q, nq, func,
					   q0, nc, tol, dmax-1, FALSE, data) ;
  QBX_FUNCTION_NAME(qbx_adaptive_quad_tri)(xe1, xstr, ne, st1, q, nq, func,
					   q1, nc, tol, dmax-1, FALSE, data) ;
  QBX_FUNCTION_NAME(qbx_adaptive_quad_tri)(xe2, xstr, ne, st2, q, nq, func,
					   q2, nc, tol, dmax-1, FALSE, data) ;
  QBX_FUNCTION_NAME(qbx_adaptive_quad_tri)(xe3, xstr, ne, st3, q, nq, func,
					   q3, nc, tol, dmax-1, FALSE, data) ;
  
  for ( i = 0 ; i < nc ; i ++ ) {
    quad[i] = q0[i] + q1[i] + q2[i] + q3[i] ;
  }

  return 0 ;
}
