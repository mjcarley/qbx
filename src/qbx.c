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

static gint calc_point(QBX_REAL *xe, gint xstr, gint ne,
		       QBX_REAL *L, QBX_REAL *p)

{
  gint i ;

  p[0] = p[1] = p[2] = 0.0 ;

  for ( i = 0 ; i < ne ; i ++ ) {
    p[0] += L[i]*xe[i*xstr+0] ;
    p[1] += L[i]*xe[i*xstr+1] ;
    p[2] += L[i]*xe[i*xstr+2] ;
  }
  
  return 0 ;
}

gint QBX_FUNCTION_NAME(qbx_element_shape_3d)(gint ne, QBX_REAL s, QBX_REAL t,
					     QBX_REAL *L,
					     QBX_REAL *dLds, QBX_REAL *dLdt,
					     QBX_REAL *dLdss, QBX_REAL *dLdst,
					     QBX_REAL *dLdtt)

{
  switch ( ne ) {
  default: g_assert_not_reached() ; break ;
  case 3:
    L[0] = 1.0 - s - t ; 
    L[1] =       s     ; 
    L[2] =           t ;
    if ( dLds == NULL ) break ;
    dLds[0] = -1.0 ; dLds[1] =  1.0 ; dLds[2] =  0.0 ;
    dLdt[0] = -1.0 ; dLdt[1] =  0.0 ; dLdt[2] =  1.0 ;
    if ( dLdss == NULL ) break ;
    dLdss[0] = 0.0 ; dLdss[1] =  0.0 ; dLdss[2] =  0.0 ;
    dLdst[0] = 0.0 ; dLdst[1] =  0.0 ; dLdst[2] =  0.0 ;
    dLdtt[0] = 0.0 ; dLdtt[1] =  0.0 ; dLdtt[2] =  0.0 ;
    break ;
  case 6:
    L[0] = 2.0*(1.0-s-t)*(0.5-s-t) ;
    L[1] = 2.0*s*(s-0.5);
    L[2] = 2.0*t*(t-0.5);
    L[3] = 4.0*s*(1.0-s-t);
    L[4] = 4.0*s*t;
    L[5] = 4.0*t*(1.0-s-t);

    if ( dLds == NULL ) break ;

    dLds[0] = -3.0 + 4.0*s + 4.0*t ;
    dLds[1] = -1.0 + 4.0*s;
    dLds[2] =  0.0;
    dLds[3] =  4.0 - 8.0*s - 4.0*t;
    dLds[4] =  4.0*t;
    dLds[5] = -4.0*t;

    dLdt[0] = -3.0 + 4.0*s + 4.0*t ;
    dLdt[1] =  0.0 ;
    dLdt[2] = -1.0 + 4.0*t ;
    dLdt[3] = -4.0*s ;
    dLdt[4] =  4.0*s ;
    dLdt[5] =  4.0 - 4.0*s - 8.0*t ;

    if ( dLdss == NULL ) break ;

    dLdss[0] =  4.0 ;
    dLdss[1] =  4.0 ;
    dLdss[2] =  0.0 ;
    dLdss[3] = -8.0 ;
    dLdss[4] =  0.0 ;
    dLdss[5] =  0.0 ;

    dLdst[0] =  4.0 ;
    dLdst[1] =  0.0 ;
    dLdst[2] =  0.0 ;
    dLdst[3] = -4.0 ;
    dLdst[4] =  4.0 ;
    dLdst[5] = -4.0 ;

    dLdtt[0] =  4.0 ;
    dLdtt[1] =  0.0 ;
    dLdtt[2] =  4.0 ;
    dLdtt[3] =  0.0 ;
    dLdtt[4] =  0.0 ;
    dLdtt[5] = -8.0 ;

    break ;
  }

  return 0 ;
}

gint QBX_FUNCTION_NAME(qbx_element_point_interp_3d)(QBX_REAL *xe,
						    gint xstr, gint ne,
						    QBX_REAL *L,
						    QBX_REAL *dLds,
						    QBX_REAL *dLdt,
						    QBX_REAL *y,
						    QBX_REAL *n,
						    QBX_REAL *J,
						    QBX_REAL *c)

{
  gint i ;
  QBX_REAL dyds[3]={0.0}, dydt[3]={0.0} ;
  
  y[0] = y[1] = y[2] = n[0] = n[1] = n[2] = 0.0 ;
  for ( i = 0 ; i < ne ; i ++ ) {
    y[0] += xe[i*xstr+0]*L[i] ; y[1] += xe[i*xstr+1]*L[i] ;
    y[2] += xe[i*xstr+2]*L[i] ; 
    dyds[0] += xe[i*xstr+0]*dLds[i] ; dyds[1] += xe[i*xstr+1]*dLds[i] ;
    dyds[2] += xe[i*xstr+2]*dLds[i] ; 
    dydt[0] += xe[i*xstr+0]*dLdt[i] ; dydt[1] += xe[i*xstr+1]*dLdt[i] ;
    dydt[2] += xe[i*xstr+2]*dLdt[i] ; 
  }

  vector_cross(n, dyds, dydt) ;
  
  *J = SQRT(n[0]*n[0] + n[1]*n[1] + n[2]*n[2]) ;

  n[0] /= (*J) ; n[1] /= (*J) ; n[2] /= (*J) ;

  return 0 ;
}

gint QBX_FUNCTION_NAME(qbx_element_point_3d)(QBX_REAL *xe, gint xstr, gint ne,
					     QBX_REAL s, QBX_REAL t,
					     QBX_REAL *y, QBX_REAL *n,
					     QBX_REAL *J, QBX_REAL *c)

{
  QBX_REAL L[32]={0.0}, dLds[32]={0.0}, dLdt[32]={0.0} ;
  
  QBX_FUNCTION_NAME(qbx_element_shape_3d)(ne, s, t, L, dLds, dLdt,
					  NULL, NULL, NULL) ;
  QBX_FUNCTION_NAME(qbx_element_point_interp_3d)(xe, xstr, ne, L,
						 dLds, dLdt, y, n, J, c) ;

  return 0 ;
}

QBX_REAL QBX_FUNCTION_NAME(qbx_quadrature_error)(gint N, QBX_REAL h, gint order,
						 QBX_REAL r, QBX_REAL rp)

/*
 * N:     order of multipole expansion
 * h:     panel width
 * order: order of quadrature rule, i.e. twice the length of a
 * Gaussian quadrature rule
 * r:     
 * rp:    smallest distance from panel to evaluation centre
 *
 * Siegel and Tornberg, https://doi.org/10.1016/j.jcp.2018.03.006
 */

{
  QBX_REAL eq, cn ;
  gint l ;

  eq = 0.0 ;
  cn = EXP(-4*(order/2)*rp/h)*2.0*M_PI ;
  for ( l = 0 ; l <= N ; l ++ ) {
    eq += cn ;
    cn *= (2*l+2)*(2*l+1)/(l+0.5)/(l+1)/(l+1)*(order/2*r/h) ;
  }

  eq *= h/order*2 ;
  
  return eq ;
}

QBX_REAL QBX_FUNCTION_NAME(qbx_element_area)(QBX_REAL *xe, gint xstr, gint ne,
					     QBX_REAL *qrule, gint nq)

{
  gint i ;
  QBX_REAL y[3], s, t, w, n[3], J, A ;

  A = 0.0 ;
  for ( i = 0 ; i < nq ; i ++ ) {
    s = qrule[3*i+0] ; t = qrule[3*i+1] ; w = qrule[3*i+2] ; 
    QBX_FUNCTION_NAME(qbx_element_point_3d)(xe, xstr, ne, s, t, y, n, &J,
					    NULL) ;
    A += w*J ;
  }
  
  return A ;
}

gint QBX_FUNCTION_NAME(qbx_triangle_truncation_error)(QBX_REAL *xe,
						      gint xstr, gint ne,
						      QBX_REAL *qrule, gint nq,
						      gint N,
						      QBX_REAL rc, gint depth,
						      QBX_REAL *ee)

{
  QBX_REAL Rb, ap, app1, r ;

  Rb =
    SQRT(4*QBX_FUNCTION_NAME(qbx_element_area)(xe, xstr, ne, qrule, nq)/M_PI) ;

  r = SQRT(rc*rc + Rb*Rb) ;
  if ( 2*(N/2) == N ) {
    ap = rc/r ; app1 = 1.0 ;
  } else {
    app1 = rc/r ; ap = 1.0 ;
  }
  
  ee[1] = N*app1*pow((1+SQRT(2))*rc, N+1)/pow(r, N+1) ;

  ee[0] = ap*pow((1+SQRT(2))*rc, N+1)/pow(r, N) ;

  return 0 ;
}

gint QBX_FUNCTION_NAME(qbx_truncation_optimal)(QBX_REAL Rb, QBX_REAL c,
					       QBX_REAL rp,
					       gint order,
					       gint pmax, gint smax,
					       QBX_REAL tol, gint *pq, gint *s)

{
  QBX_REAL rb, eT, b ;

  b = 1 + M_SQRT2 ;
  rb = SQRT(c*c + Rb*Rb) ;

  eT = c/rb*(b*c)*(b*c)*(b*c)/rb/rb ;

  *pq = 2 ;
  while ( eT > tol && *pq < pmax ) {
    eT *= b*c/rb ; /* pq+1, odd */
    eT *= b*c/rb ; /* *c/rb ; /\* pq+2, even *\/ */
    
    *pq += 2 ;
  }

  /* pq is now the minimum order to achieve the truncation error */
  *s = 0 ;
  while ( (QBX_FUNCTION_NAME(qbx_quadrature_error)(*pq, Rb, order, c, rp) >
	   tol ) && (*s < smax) ) {
    Rb *= 0.5 ; (*s) ++ ;
  }

  return 0 ;
}

gint QBX_FUNCTION_NAME(qbx_quadrature_optimal)(QBX_REAL Rb, QBX_REAL r0,
					       QBX_REAL r1, gint nr,
					       gint order, gint ngp,
					       gint pmax, gint smax,
					       QBX_REAL tol,
					       QBX_REAL *rc, gint *pq,
					       gint *s)

{
  gint ncalc, ncmin, pc, sc, i ;
  QBX_REAL r ;
  
  ncmin = G_MAXINT ;

  for ( i = 0 ; i <= nr ; i ++ ) {
    r = r0 + (r1-r0)*i/nr ;
    QBX_FUNCTION_NAME(qbx_truncation_optimal)(Rb, r, r, order, pmax, smax,
					      tol, &pc, &sc) ;
    ncalc = ngp*(1 << sc) ;
    if ( ncalc < ncmin ) {
      ncmin = ncalc ;
      *pq = pc ; *s = sc ;
      *rc = r ;
    }
  }
  
  return 0 ;
}

gint QBX_FUNCTION_NAME(qbx_quadrature_select)(gint nq, QBX_REAL **q,
					      gint *order)

{
#ifdef QBX_SINGLE_PRECISION
  switch ( nq ) {
  default: g_assert_not_reached() ; break ;
  case   7: *q = WANDZURA_7_F   ; *order =  5 ; break ;
  case  25: *q = WANDZURA_25_F  ; *order = 10 ; break ;
  case  54: *q = WANDZURA_54_F  ; *order = 15 ; break ;
  case  85: *q = WANDZURA_85_F  ; *order = 20 ; break ;
  case 126: *q = WANDZURA_126_F ; *order = 25 ; break ;
  case 175: *q = WANDZURA_175_F ; *order = 30 ; break ;
  case 453: *q = XIAO_GIMBUTAS_453_F ; *order = 50 ; break ;
  }
#else /*QBX_SINGLE_PRECISION*/
  switch ( nq ) {
  default: g_assert_not_reached() ; break ;
  case   7: *q = WANDZURA_7   ; *order =  5 ; break ;
  case  25: *q = WANDZURA_25  ; *order = 10 ; break ;
  case  54: *q = WANDZURA_54  ; *order = 15 ; break ;
  case  85: *q = WANDZURA_85  ; *order = 20 ; break ;
  case 126: *q = WANDZURA_126 ; *order = 25 ; break ;
  case 175: *q = WANDZURA_175 ; *order = 30 ; break ;
  case 453: *q = XIAO_GIMBUTAS_453 ; *order = 50 ; break ;
  }
#endif /*QBX_SINGLE_PRECISION*/
  
  return 0 ;
}

gint QBX_FUNCTION_NAME(qbx_quadrature_optimal_points)(QBX_REAL Rb,
						      QBX_REAL r0,
						      QBX_REAL r1,
						      gint nr,
						      gint order, gint ngp,
						      gint pmax, gint smax,
						      QBX_REAL tol,
						      QBX_REAL *x0,
						      QBX_REAL *n,
						      QBX_REAL *x,
						      QBX_REAL *rc,
						      gint *pq, gint *s)

{
  gint ncalc, ncmin, pc, sc, i ;
  QBX_REAL r, rp, c[3] ;
  
  ncmin = G_MAXINT ;

  for ( i = 0 ; i <= nr ; i ++ ) {
    r = r0 + (r1-r0)*i/nr ;
    c[0] = x0[0] + r*n[0] ;
    c[1] = x0[1] + r*n[1] ;
    c[2] = x0[2] + r*n[2] ;
    rp = SQRT((x[0] - c[0])*(x[0] - c[0]) +
	      (x[1] - c[1])*(x[1] - c[1]) +
	      (x[2] - c[2])*(x[2] - c[2])) ;
    QBX_FUNCTION_NAME(qbx_truncation_optimal)(Rb, r, rp, order, pmax, smax,
					      tol, &pc, &sc) ;
    ncalc = ngp*(1 << sc) ;
    if ( ncalc < ncmin ) {
      ncmin = ncalc ;
      *pq = pc ; *s = sc ;
      *rc = r ;
    }
  }

  return 0 ;
}

gint QBX_FUNCTION_NAME(qbx_triangle_expansion)(QBX_REAL *xe, gint xstr, gint ne,
					       QBX_REAL *fe, gint fstr, gint nf,
					       QBX_REAL *q, gint nq, gint oq,
					       QBX_REAL *x, QBX_REAL tol,
					       gint Nmax, gint smax,
					       QBX_REAL *c, QBX_REAL *C,
					       gint *N, gint *s)

{
  QBX_REAL Ae, rb ;

  Ae = QBX_FUNCTION_NAME(qbx_element_area)(xe, xstr, ne, q, nq) ;
  rb = SQRT(4.0*Ae/M_PI) ;

  
  
  return 0 ;
}
						  
gint QBX_FUNCTION_NAME(qbx_triangle_curvature)(QBX_REAL *xe, gint xstr, gint ne,
					       QBX_REAL s, QBX_REAL t,
					       QBX_REAL *kg, QBX_REAL *km)

{
  QBX_REAL p1[3], ps[3], pt[3] ;
  QBX_REAL L0[16], Ls[16], Lt[16], Lss[16], Lst[16], Ltt[16],
    ds, dt, n[3], pss[3], pst[3], ptt[3] ;
  QBX_REAL E, F, G, L, M, N ;

  QBX_FUNCTION_NAME(qbx_element_shape_3d)(ne, s, t, L0, Ls, Lt,
					  Lss, Lst, Ltt) ;
  calc_point(xe, xstr, ne, L0,  p1) ;
  calc_point(xe, xstr, ne, Ls, ps) ;
  calc_point(xe, xstr, ne, Lt, pt) ;
  calc_point(xe, xstr, ne, Lss, pss) ;
  calc_point(xe, xstr, ne, Lst, pst) ;
  calc_point(xe, xstr, ne, Ltt, ptt) ;

  vector_cross(n,ps,pt) ;
  ds = SQRT(n[0]*n[0]+n[1]*n[1]+n[2]*n[2]) ;
  n[0] /= ds ; n[1] /= ds ; n[2] /= ds ;
  
  /*first fundamental forms*/
  E = vector_scalar(ps, ps) ;
  F = vector_scalar(ps, pt) ;
  G = vector_scalar(pt, pt) ;

  /*second fundamental forms*/
  L = vector_scalar(pss, n) ;
  M = vector_scalar(pst, n) ;
  N = vector_scalar(ptt, n) ;

  *kg = (L*N - M*M)/(E*G - F*F) ;
  *km = 0.5*(E*N - 2*F*M + G*L)/(E*G - F*F) ;
  
  return 0 ;
}

gint QBX_FUNCTION_NAME(qbx_triangle_laplace_self_quad)(QBX_REAL *xe,
						       gint xstr,
						       gint ne,
						       QBX_REAL s0,
						       QBX_REAL t0,
						       QBX_REAL *q,
						       gint nq, gint oq,
						       gint Nmax, gint dmax,
						       QBX_REAL tol,
						       QBX_REAL *Is, gint istr,
						       QBX_REAL *Id, gint dstr,
						       QBX_REAL *work)

{
  QBX_REAL rp, w, x[3], c[3], n[3], J, fe[64], ee[64] ;
#ifdef QBX_SINGLE_PRECISION
  QBX_REAL *qa = WANDZURA_25_F ;
#else /*QBX_SINGLE_PRECISION*/
  QBX_REAL *qa = WANDZURA_25 ;
#endif /*QBX_SINGLE_PRECISION*/
  gint nqa = 25, N, d, i, fstr, str ;
  gboolean clear ;
  
  /*element dimension*/
  w = SQRT(4*QBX_FUNCTION_NAME(qbx_element_area)(xe, xstr, ne, qa, nqa)/M_PI) ;
  /*location of field point on element*/
  QBX_FUNCTION_NAME(qbx_element_point_3d)(xe, xstr, ne, s0, t0, x, n, &J,
					  NULL) ;
  QBX_FUNCTION_NAME(qbx_quadrature_optimal_points)(w,
						   0.125*w/(1 << dmax),
						   4.0*w/(1 << dmax),
						   16, oq, nq,
						   Nmax, dmax, tol,
						   x, n, x,
						   &rp, &N, &d) ;

  fprintf(stderr, "w = %lg; rp = %lg; N = %d; d = %d\n",
	  w, rp, N, d) ;
  
  /*expansion centre*/
  vector_shift(c,x,n,rp) ;

  /*set up nodal quantities for call to expansion generation*/
  fstr = ne ;
  memset(fe, 0, ne*fstr*sizeof(QBX_REAL)) ;
  for ( i = 0 ; i < ne ; i ++ ) {
    fe[i*fstr+i] = 1.0 ;
  }
  clear = TRUE ;

  str = ne ;

  if ( Id != NULL ) {
    QBX_FUNCTION_NAME(qbx_expansion_make_laplace_adaptive)(xe, xstr, ne,
							   fe, fstr, ne,
							   q, nq, oq,
							   c, rp, N,
							   work, str, d, tol,
							   w,
							   clear, FALSE) ;
    QBX_FUNCTION_NAME(qbx_expansion_eval_laplace)(c, N, work, str, x, Id, ee) ;
  }
  
  memset(fe, 0, ne*fstr*sizeof(QBX_REAL)) ;
  for ( i = 0 ; i < ne ; i ++ ) {
    fe[i*fstr+i] = 1.0 ;
  }

  if ( Is != NULL ) {
    QBX_FUNCTION_NAME(qbx_expansion_make_laplace_adaptive)(xe, xstr, ne,
							   fe, fstr, ne,
							   q, nq, oq,
							   c, rp, N,
							   work, str, d, tol,
							   w,
							   clear, TRUE) ;
    QBX_FUNCTION_NAME(qbx_expansion_eval_laplace)(c, N, work, str, x, Is, ee) ;
  }
  
  return 0 ;
}
					  
