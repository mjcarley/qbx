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


inline static void shape6(QBX_REAL s, QBX_REAL t, QBX_REAL L[])

{
  L[0] = 2.0*(1.0-s-t)*(0.5-s-t) ;
  L[1] = 2.0*s*(s-0.5) ;
  L[2] = 2.0*t*(t-0.5) ;
  L[3] = 4.0*s*(1.0-s-t) ;
  L[4] = 4.0*s*t ;
  L[5] = 4.0*t*(1.0-s-t) ;

  return ;
}

inline static void shape_derivatives6(QBX_REAL s, QBX_REAL t,
				      QBX_REAL L[],
				      QBX_REAL Ls[], QBX_REAL Lt[])

{
  L[0] = 2.0*(1.0-s-t)*(0.5-s-t) ;
  L[1] = 2.0*s*(s-0.5) ;
  L[2] = 2.0*t*(t-0.5) ;
  L[3] = 4.0*s*(1.0-s-t) ;
  L[4] = 4.0*s*t ;
  L[5] = 4.0*t*(1.0-s-t) ;

  Ls[0] = -3.0 + 4.0*s + 4.0*t ;
  Ls[1] = -1.0 + 4.0*s ;
  Ls[2] =  0.0 ;
  Ls[3] =  4.0 - 8.0*s - 4.0*t ;
  Ls[4] =  4.0*t ;
  Ls[5] = -4.0*t ;

  Lt[0] = -3.0 + 4.0*s + 4.0*t ;
  Lt[1] =  0.0 ;
  Lt[2] = -1.0 + 4.0*t ;
  Lt[3] = -4.0*s ;
  Lt[4] =  4.0*s ;
  Lt[5] =  4.0 - 4.0*s - 8.0*t ;

  return ;
}

inline static void point_interp6(QBX_REAL *xb, QBX_REAL *fb, gint i,
				 QBX_REAL *xe, gint xstr,
				 QBX_REAL *fe,
				 QBX_REAL s, QBX_REAL t)

{
  QBX_REAL L[6] ;
  
  shape6(s, t, L) ;
  xb[xstr*i+0] =
    L[0]*xe[xstr*0+0] + L[1]*xe[xstr*1+0] + L[2]*xe[xstr*2+0] +
    L[3]*xe[xstr*3+0] + L[4]*xe[xstr*4+0] + L[5]*xe[xstr*5+0] ;
  xb[xstr*i+1] =
    L[0]*xe[xstr*0+1] + L[1]*xe[xstr*1+1] + L[2]*xe[xstr*2+1] +
    L[3]*xe[xstr*3+1] + L[4]*xe[xstr*4+1] + L[5]*xe[xstr*5+1] ;
  xb[xstr*i+2] =
    L[0]*xe[xstr*0+2] + L[1]*xe[xstr*1+2] + L[2]*xe[xstr*2+2] +
    L[3]*xe[xstr*3+2] + L[4]*xe[xstr*4+2] + L[5]*xe[xstr*5+2] ;

  fb[2*i+0] =
    L[0]*fe[0] + L[1]*fe[2] + L[2]*fe[4 ] +
    L[3]*fe[6] + L[4]*fe[8] + L[5]*fe[10] ;
  fb[2*i+1] =
    L[0]*fe[1] + L[1]*fe[3] + L[2]*fe[5 ] +
    L[3]*fe[7] + L[4]*fe[9] + L[5]*fe[11] ;
  
  return ;
}

inline static void point_interp_jac6(QBX_REAL *xe, gint xstr,
				     QBX_REAL L[],
				     QBX_REAL Ls[], QBX_REAL Lt[],
				     QBX_REAL *y, QBX_REAL n[], QBX_REAL *J)

{
  QBX_REAL ys[3], yt[3] ;
  
  y[0] = 
    L[0]*xe[xstr*0+0] + L[1]*xe[xstr*1+0] + L[2]*xe[xstr*2+0] +
    L[3]*xe[xstr*3+0] + L[4]*xe[xstr*4+0] + L[5]*xe[xstr*5+0] ;
  y[1] = 
    L[0]*xe[xstr*0+1] + L[1]*xe[xstr*1+1] + L[2]*xe[xstr*2+1] +
    L[3]*xe[xstr*3+1] + L[4]*xe[xstr*4+1] + L[5]*xe[xstr*5+1] ;
  y[2] = 
    L[0]*xe[xstr*0+2] + L[1]*xe[xstr*1+2] + L[2]*xe[xstr*2+2] +
    L[3]*xe[xstr*3+2] + L[4]*xe[xstr*4+2] + L[5]*xe[xstr*5+2] ;
  ys[0] = 
    Ls[0]*xe[xstr*0+0] + Ls[1]*xe[xstr*1+0] + Ls[2]*xe[xstr*2+0] +
    Ls[3]*xe[xstr*3+0] + Ls[4]*xe[xstr*4+0] + Ls[5]*xe[xstr*5+0] ;
  ys[1] = 
    Ls[0]*xe[xstr*0+1] + Ls[1]*xe[xstr*1+1] + Ls[2]*xe[xstr*2+1] +
    Ls[3]*xe[xstr*3+1] + Ls[4]*xe[xstr*4+1] + Ls[5]*xe[xstr*5+1] ;
  ys[2] = 
    Ls[0]*xe[xstr*0+2] + Ls[1]*xe[xstr*1+2] + Ls[2]*xe[xstr*2+2] +
    Ls[3]*xe[xstr*3+2] + Ls[4]*xe[xstr*4+2] + Ls[5]*xe[xstr*5+2] ;
  yt[0] = 
    Lt[0]*xe[xstr*0+0] + Lt[1]*xe[xstr*1+0] + Lt[2]*xe[xstr*2+0] +
    Lt[3]*xe[xstr*3+0] + Lt[4]*xe[xstr*4+0] + Lt[5]*xe[xstr*5+0] ;
  yt[1] = 
    Lt[0]*xe[xstr*0+1] + Lt[1]*xe[xstr*1+1] + Lt[2]*xe[xstr*2+1] +
    Lt[3]*xe[xstr*3+1] + Lt[4]*xe[xstr*4+1] + Lt[5]*xe[xstr*5+1] ;
  yt[2] = 
    Lt[0]*xe[xstr*0+2] + Lt[1]*xe[xstr*1+2] + Lt[2]*xe[xstr*2+2] +
    Lt[3]*xe[xstr*3+2] + Lt[4]*xe[xstr*4+2] + Lt[5]*xe[xstr*5+2] ;

  qbx_vector_cross(n, ys, yt) ;
  *J = qbx_vector_length(n) ;

  n[0] /= *J ; n[1] /= *J ; n[2] /= *J ;
  
  return ;
}

static void triangle_divide_loop60(QBX_REAL *xe, gint xstr, gint ne,
				   QBX_REAL *fe, gint fstr, gint nf,
				   QBX_REAL *xl, QBX_REAL *fl)

{
  qbx_point_copy(xl, fl, 0, xe, xstr, fe, 0) ;
  point_interp6(xl, fl, 1, xe, xstr, fe, 0.5, 0.0) ;
  point_interp6(xl, fl, 2, xe, xstr, fe, 0.0, 0.5) ;
  point_interp6(xl, fl, 3, xe, xstr, fe, 0.25, 0.0 ) ;
  point_interp6(xl, fl, 4, xe, xstr, fe, 0.25, 0.25) ;
  point_interp6(xl, fl, 5, xe, xstr, fe, 0.0,  0.25) ;

  return ;
}

static void triangle_divide_loop61(QBX_REAL *xe, gint xstr, gint ne,
				   QBX_REAL *fe, gint fstr, gint nf,
				   QBX_REAL *xl, QBX_REAL *fl)
{  
  point_interp6(xl, fl, 0, xe, xstr, fe, 0.5, 0.0) ;
  qbx_point_copy(xl, fl, 1, xe, xstr, fe, 1) ;
  point_interp6(xl, fl, 2, xe, xstr, fe, 0.5, 0.5) ;
  point_interp6(xl, fl, 3, xe, xstr, fe, 0.75, 0.0 ) ;
  point_interp6(xl, fl, 4, xe, xstr, fe, 0.75, 0.25) ;
  point_interp6(xl, fl, 5, xe, xstr, fe, 0.5 , 0.25) ;

  return ;
}

static void triangle_divide_loop62(QBX_REAL *xe, gint xstr, gint ne,
				   QBX_REAL *fe, gint fstr, gint nf,
				   QBX_REAL *xl, QBX_REAL *fl)
{
  point_interp6(xl, fl, 0, xe, xstr, fe, 0.0, 0.5) ;
  point_interp6(xl, fl, 1, xe, xstr, fe, 0.5, 0.5) ;
  qbx_point_copy(xl, fl, 2, xe, xstr, fe, 2) ;
  point_interp6(xl, fl, 3, xe, xstr, fe, 0.25, 0.5) ;
  point_interp6(xl, fl, 4, xe, xstr, fe, 0.25, 0.75) ;
  point_interp6(xl, fl, 5, xe, xstr, fe, 0.0, 0.75) ;

  return ;
}

static void triangle_divide_loop63(QBX_REAL *xe, gint xstr, gint ne,
				   QBX_REAL *fe, gint fstr, gint nf,
				   QBX_REAL *xl, QBX_REAL *fl)
{
  point_interp6(xl, fl, 0, xe, xstr, fe, 0.5, 0.0) ;
  point_interp6(xl, fl, 1, xe, xstr, fe, 0.5, 0.5) ;
  point_interp6(xl, fl, 2, xe, xstr, fe, 0.0, 0.5) ;
  point_interp6(xl, fl, 3, xe, xstr, fe, 0.5, 0.25) ;
  point_interp6(xl, fl, 4, xe, xstr, fe, 0.25, 0.5) ;
  point_interp6(xl, fl, 5, xe, xstr, fe, 0.25, 0.25) ;

  return ;
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
  
  if ( ne == 3 ) {
    if ( depth == 0 || err < tol ) {
      for ( i = 0 ; i < nq ; i ++ ) {
	s = q[3*i+0] ; t = q[3*i+1] ; wt = q[3*i+2] ;
	qbx_shape_derivatives3(s, t, L, Ls, Lt) ;
	/* point_interp_jac3(xd, dstr, L, Ls, Lt, y, n, &J) ; */
	qbx_point_interp_jac3(xd, dstr,
			      L[0], L[1], L[2],
			      Ls[0], Ls[1], Ls[2],
			      Lt[0], Lt[1], Lt[2],
			      y, n, J) ;
	wt *= J ;
	s = L[0]*st[0] + L[1]*st[2] + L[2]*st[4] ;
	t = L[0]*st[1] + L[1]*st[3] + L[2]*st[5] ;
	func(s, t, wt, xc, y, n, N, data) ;
      }
      
      return 0 ;
    }
    
    qbx_triangle_divide_loop30(xd, dstr, st, xb, std) ;
    qbx_triangle_adaptive(func, xt, tstr, ne, xb, dstr, std, 0.5*w,
			  xc, N, depth-1, q, nq, oq, tol, data) ;
    qbx_triangle_divide_loop31(xd, dstr, st, xb, std) ;
    qbx_triangle_adaptive(func, xt, tstr, ne, xb, dstr, std, 0.5*w,
			  xc, N, depth-1, q, nq, oq, tol, data) ;
    qbx_triangle_divide_loop32(xd, dstr, st, xb, std) ;
    qbx_triangle_adaptive(func, xt, tstr, ne, xb, dstr, std, 0.5*w,
			  xc, N, depth-1, q, nq, oq, tol, data) ;
    qbx_triangle_divide_loop33(xd, dstr, st, xb, std) ;
    qbx_triangle_adaptive(func, xt, tstr, ne, xb, dstr, std, 0.5*w,
			  xc, N, depth-1, q, nq, oq, tol, data) ;

    return 0 ;
  }

  if ( ne == 6 ) {
    if ( depth == 0 || err < tol ) {
      for ( i = 0 ; i < nq ; i ++ ) {
	s = q[3*i+0] ; t = q[3*i+1] ; wt = q[3*i+2] ;
	shape_derivatives6(s, t, L, Ls, Lt) ;
	point_interp_jac6(xd, dstr, L, Ls, Lt, y, n, &J) ;
	wt *= J ;
	s =
	  L[0]*st[0] + L[1]*st[2] + L[2]*st[ 4] +
	  L[3]*st[6] + L[4]*st[8] + L[5]*st[10] ;
	t =
	  L[0]*st[1] + L[1]*st[3] + L[2]*st[ 5] +
	  L[3]*st[7] + L[4]*st[9] + L[5]*st[11] ;
	func(s, t, wt, xc, y, n, N, data) ;
      }

      return 0 ;
    }

    triangle_divide_loop60(xd, dstr, ne, st, 2, 2, xb, std) ;
    qbx_triangle_adaptive(func, xt, tstr, ne, xb, dstr, std, 0.5*w,
			  xc, N, depth-1, q, nq, oq, tol, data) ;
    triangle_divide_loop61(xd, dstr, ne, st, 2, 2, xb, std) ;
    qbx_triangle_adaptive(func, xt, tstr, ne, xb, dstr, std, 0.5*w,
			  xc, N, depth-1, q, nq, oq, tol, data) ;
    triangle_divide_loop62(xd, dstr, ne, st, 2, 2, xb, std) ;
    qbx_triangle_adaptive(func, xt, tstr, ne, xb, dstr, std, 0.5*w,
			  xc, N, depth-1, q, nq, oq, tol, data) ;
    triangle_divide_loop63(xd, dstr, ne, st, 2, 2, xb, std) ;
    qbx_triangle_adaptive(func, xt, tstr, ne, xb, dstr, std, 0.5*w,
			  xc, N, depth-1, q, nq, oq, tol, data) ;

    return 0 ;
  }

  g_assert_not_reached() ;
  
  return 0 ;
}
