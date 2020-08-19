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

  qbx_vector_cross(n, dyds, dydt) ;
  
  *J = qbx_vector_length(n) ;

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
    rp = r0 + (r1-r0)*i/nr ;
    qbx_vector_shift(c,x0,n,rp) ;
    r = qbx_vector_distance(x,c) ;
    QBX_FUNCTION_NAME(qbx_truncation_optimal)(Rb, r, rp, order, pmax, smax,
					      tol, &pc, &sc) ;
    ncalc = ngp*(1 << sc) ;
    if ( ncalc < ncmin ) {
      ncmin = ncalc ;
      *pq = pc ; *s = sc ;
      *rc = rp ;
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
  /* QBX_REAL Ae, rb ; */

  /* Ae = QBX_FUNCTION_NAME(qbx_element_area)(xe, xstr, ne, q, nq) ; */
  /* rb = SQRT(4.0*Ae/M_PI) ; */

  return 0 ;
}
						  
gint QBX_FUNCTION_NAME(qbx_triangle_curvature)(QBX_REAL *xe, gint xstr, gint ne,
					       QBX_REAL s, QBX_REAL t,
					       QBX_REAL *kg, QBX_REAL *km)

{
  QBX_REAL p1[3], ps[3], pt[3] ;
  QBX_REAL L0[16], Ls[16], Lt[16], Lss[16], Lst[16], Ltt[16],
    ds, n[3], pss[3], pst[3], ptt[3] ;
  QBX_REAL E, F, G, L, M, N ;

  QBX_FUNCTION_NAME(qbx_element_shape_3d)(ne, s, t, L0, Ls, Lt,
					  Lss, Lst, Ltt) ;
  calc_point(xe, xstr, ne, L0,  p1) ;
  calc_point(xe, xstr, ne, Ls, ps) ;
  calc_point(xe, xstr, ne, Lt, pt) ;
  calc_point(xe, xstr, ne, Lss, pss) ;
  calc_point(xe, xstr, ne, Lst, pst) ;
  calc_point(xe, xstr, ne, Ltt, ptt) ;

  qbx_vector_cross(n,ps,pt) ;
  ds = qbx_vector_length(n) ;
  n[0] /= ds ; n[1] /= ds ; n[2] /= ds ;
  
  /*first fundamental forms*/
  E = qbx_vector_scalar(ps, ps) ;
  F = qbx_vector_scalar(ps, pt) ;
  G = qbx_vector_scalar(pt, pt) ;

  /*second fundamental forms*/
  L = qbx_vector_scalar(pss, n) ;
  M = qbx_vector_scalar(pst, n) ;
  N = qbx_vector_scalar(ptt, n) ;

  *kg = (L*N - M*M)/(E*G - F*F) ;
  *km = 0.5*(E*N - 2*F*M + G*L)/(E*G - F*F) ;
  
  return 0 ;
}

static gint QBX_FUNCTION_NAME(laplace_quad_sl)(QBX_REAL s, QBX_REAL t,
					       QBX_REAL w,
					       QBX_REAL *x, QBX_REAL *y,
					       QBX_REAL *n,
					       gint N,
					       QBX_REAL *C, gint str,
					       gint nf,
					       gpointer data[])
{
  gint ne = *((gint *)(data[2])) ;
  QBX_REAL r, th, ph, L[32], sc ;
  QBX_REAL Cth, Sth, Cmph[64], Smph[64], *Pnm1, *Pn, P[128] ;
  gint nn, mm, idx, j ;
  
  QBX_FUNCTION_NAME(qbx_cartesian_to_spherical)(x, y, &r, &th, &ph) ;
  QBX_FUNCTION_NAME(qbx_element_shape_3d)(ne, s, t, L, NULL, NULL,
					  NULL, NULL, NULL) ;

  Pnm1 = &(P[0]) ; Pn = &(Pnm1[N+2]) ;
  Cth = COS(th) ; Sth = SIN(th) ;
  QBX_FUNCTION_NAME(qbx_legendre_init)(Cth, Sth, &(Pnm1[0]),
				       &(Pn[0]), &(Pn[1])) ;
  Pnm1[1] = 0.0 ;
  
  Cmph[0] = 1.0 ; Cmph[1] = COS(ph) ;
  Smph[0] = 0.0 ; Smph[1] = SIN(ph) ;
  
  w /= r ;
  
  nn = 0 ; mm = 0 ; sc = 1.0/(2*nn+1) ;
  idx = nn*nn ;
  for ( j = 0 ; j < nf ; j ++ ) C[str*idx+j] += Pnm1[mm]*w*sc*L[j] ;
  
  w /= r ;
  nn = 1 ; sc = 1.0/(2*nn+1) ;
  mm =  0 ;
  idx = nn*nn ;
  for ( j = 0 ; j < nf ; j ++ ) C[str*idx+j] += Pn[mm]*w*sc*L[j] ;
  mm =  1 ;
  idx = qbx_index_laplace_nm(nn, mm) ;
  for ( j = 0 ; j < nf ; j ++ ) {
    C[str*(idx+0)+j] += Pn[mm]*Cmph[mm]*w*sc*L[j] ;
    C[str*(idx+1)+j] -= Pn[mm]*Smph[mm]*w*sc*L[j] ;
  }
  
  for ( nn = 2 ; nn <= N ; nn ++ ) {
    QBX_FUNCTION_NAME(qbx_legendre_recursion_array)(&Pnm1, &Pn, nn-1,
						    Cth, Sth) ;
    w /= r ;
    sc = 1.0/(2*nn+1) ;      
    mm =  0 ;
    idx = nn*nn ;
    for ( j = 0 ; j < nf ; j ++ ) C[str*idx+j] += Pn[mm]*w*sc*L[j] ;
    
    Cmph[nn] = Cmph[nn-1]*Cmph[1] - Smph[nn-1]*Smph[1] ;
    Smph[nn] = Smph[nn-1]*Cmph[1] + Cmph[nn-1]*Smph[1] ;
    
    for ( mm = 1 ; mm <= nn ; mm ++ ) {
      idx = qbx_index_laplace_nm(nn, mm) ;
      for ( j = 0 ; j < nf ; j ++ ) {
	C[str*(idx+0)+j] += Pn[mm]*Cmph[mm]*w*sc*L[j] ;
	C[str*(idx+1)+j] -= Pn[mm]*Smph[mm]*w*sc*L[j] ;
      }
    }
  }

  return 0 ;
}

static gint QBX_FUNCTION_NAME(laplace_quad_dl)(QBX_REAL s, QBX_REAL t,
					       QBX_REAL w,
					       QBX_REAL *x, QBX_REAL *y,
					       QBX_REAL *n,
					       gint N,
					       QBX_REAL *C, gint str,
					       gint nf,
					       gpointer data[])
{
  gint ne = *((gint *)(data[2])) ;
  QBX_REAL r, th, ph, L[32], rn, nR, sc, fr, fi ;
  QBX_REAL Cth, Sth, Cmph[64], Smph[64], dP, nC, nph, *Pnm1, *Pn, P[128] ;
  gint nn, mm, idx, j ;
  
  QBX_FUNCTION_NAME(qbx_cartesian_to_spherical)(x, y, &r, &th, &ph) ;
  QBX_FUNCTION_NAME(qbx_element_shape_3d)(ne, s, t, L, NULL, NULL,
					  NULL, NULL, NULL) ;

  Pnm1 = &(P[0]) ; Pn = &(Pnm1[N+2]) ;
  
  Cth = COS(th) ; Sth = SIN(th) ;
  QBX_FUNCTION_NAME(qbx_legendre_init)(Cth, Sth, &(Pnm1[0]),
					 &(Pn[0]), &(Pn[1])) ;
  Pnm1[1] = Pnm1[2] = Pn[2] = 0.0 ;
  
  Cmph[0] = 1.0 ; Cmph[1] = COS(ph) ;
  Smph[0] = 0.0 ; Smph[1] = SIN(ph) ;

  nR = -(n[0]*(x[0]-y[0]) + n[1]*(x[1]-y[1]) + n[2]*(x[2]-y[2]))/r ;
  nph = (n[0]*(x[1]-y[1]) - n[1]*(x[0]-y[0]))/r/r/Sth/Sth ;
  nC = (n[2]-Cth*nR)/r ;
  
  rn = r ;
  
  nn = 0 ; mm = 0 ; sc = 1.0/(2*nn+1) ;
  dP = 0 ;
  fr = (Cmph[mm]*( dP*nC-(nn+1)*Pnm1[mm]*nR/r) -
	mm*Smph[mm]*Pnm1[mm]*nph)/rn ;
  
  idx = nn*nn ;
  for ( j = 0 ; j < nf ; j ++ ) C[str*idx+j] += fr*w*sc*L[j] ;
  
  rn *= r ;
  nn = 1 ; sc = 1.0/(2*nn+1) ;
  mm = 0 ;
  dP = SQRT((2*nn+1)/(2*nn-1)*(nn+mm)*(nn-mm))*Pnm1[mm] - nn*Cth*Pn[mm] ;
  dP /= Sth*Sth ;
  fr = (Cmph[mm]*( dP*nC-(nn+1)*Pn[mm]*nR/r) - mm*Smph[mm]*Pn[mm]*nph)/rn ;
  
  idx = nn*nn ;
  for ( j = 0 ; j < nf ; j ++ ) C[str*idx+j] += fr*w*sc*L[j] ;
  
  mm = 1 ;
  dP = SQRT((2*nn+1)/(2*nn-1)*(nn+mm)*(nn-mm))*Pnm1[mm] - nn*Cth*Pn[mm] ;
  dP /= Sth*Sth ;
  fr = (Cmph[mm]*( dP*nC-(nn+1)*Pn[mm]*nR/r) - mm*Smph[mm]*Pn[mm]*nph)/rn ;
  fi = (Smph[mm]*(-dP*nC+(nn+1)*Pn[mm]*nR/r) - mm*Cmph[mm]*Pn[mm]*nph)/rn ;
  
  idx = qbx_index_laplace_nm(nn, mm) ;
  for ( j = 0 ; j < nf ; j ++ ) {
    C[str*(idx+0)+j] += fr*w*sc*L[j] ;
    C[str*(idx+1)+j] += fi*w*sc*L[j] ;
  }

  for ( nn = 2 ; nn <= N ; nn ++ ) {
    QBX_FUNCTION_NAME(qbx_legendre_recursion_array)(&Pnm1, &Pn, nn-1,
						    Cth, Sth) ;
    rn *= r ;
    sc = 1.0/(2*nn+1) ;      
    mm =  0 ;
    dP = SQRT((2.0*nn+1)/(2*nn-1)*(nn+mm)*(nn-mm))*Pnm1[mm] - nn*Cth*Pn[mm] ;
    dP /= Sth*Sth ;
    fr = (Cmph[mm]*( dP*nC-(nn+1)*Pn[mm]*nR/r) - mm*Smph[mm]*Pn[mm]*nph)/rn ;
    fi = (Smph[mm]*(-dP*nC+(nn+1)*Pn[mm]*nR/r) - mm*Cmph[mm]*Pn[mm]*nph)/rn ;
    idx = nn*nn ;
    for ( j = 0 ; j < nf ; j ++ )  C[str*idx+j] += fr*w*sc*L[j] ;
    
    Cmph[nn] = Cmph[nn-1]*Cmph[1] - Smph[nn-1]*Smph[1] ;
    Smph[nn] = Smph[nn-1]*Cmph[1] + Cmph[nn-1]*Smph[1] ;
    
    for ( mm = 1 ; mm <= nn ; mm ++ ) {
      dP = SQRT((2.0*nn+1)/(2*nn-1)*(nn+mm)*(nn-mm))*Pnm1[mm] -
	nn*Cth*Pn[mm] ;
      dP /= Sth*Sth ;
      fr = (Cmph[mm]*( dP*nC-(nn+1)*Pn[mm]*nR/r) -
	    mm*Smph[mm]*Pn[mm]*nph)/rn ;
      fi = (Smph[mm]*(-dP*nC+(nn+1)*Pn[mm]*nR/r) -
	    mm*Cmph[mm]*Pn[mm]*nph)/rn ;
      
      idx = qbx_index_laplace_nm(nn, mm) ;
      for ( j = 0 ; j < nf ; j ++ ) {
	C[str*(idx+0)+j] += fr*w*sc*L[j] ;
	C[str*(idx+1)+j] += fi*w*sc*L[j] ;
      }
    }
  }
  
  return 0 ;
}

gint QBX_FUNCTION_NAME(select_quadrature)(QBX_REAL w, QBX_REAL tol,
					  gint nq, gint oq, gint Nmax,
					  QBX_REAL rc, QBX_REAL rp,
					  gint *N, gint *d)

{
  QBX_REAL rb, eT, b ;
  gint i ;
  
  b = 1 + M_SQRT2 ;
  rb = SQRT(rp*rp + w*w) ;
  eT = rp/rb*(b*rp)*(b*rp)*(b*rp)/rb/rb ;

  *N = 2 ; eT *= (*N) ;
  while ( eT > tol && (*N) < Nmax )  
  {
    eT *= 1.0*rp/rb ; /* pq+1, odd */
    eT *= b*rp/rb ; /* *c/rb ; /\* pq+2, even *\/ */

    *N += 2 ;
    eT *= (*N)/((*N)-2) ;
  }  ;

  (*N) += 1 ;

  (*d) ++ ;
  while ( QBX_FUNCTION_NAME(qbx_quadrature_error)((*N),
						  w/(1<<((*d))),
						  oq, rc, rp) <
  	  tol && (*d > 0) ) {
    (*d) -- ;
  }

  (*d) += 1 ;
  fprintf(stderr, "%d %d %e %e\n", *N, (*d), eT,
	 QBX_FUNCTION_NAME(qbx_quadrature_error)(*N, w/(1<<(*d)),
						 oq, rc, rp)) ;
  
  return 0 ;
}

gint QBX_FUNCTION_NAME(qbx_triangle_laplace_self_quad)(QBX_REAL *xe,
						       gint xstr,
						       gint ne,
						       QBX_REAL s0,
						       QBX_REAL t0,
						       gboolean in,
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
  
  /*element dimension*/
  w = SQRT(4*QBX_FUNCTION_NAME(qbx_element_area)(xe, xstr, ne, qa, nqa)/M_PI) ;
  /*location of field point on element*/
  QBX_FUNCTION_NAME(qbx_element_point_3d)(xe, xstr, ne, s0, t0, x, n, &J,
					  NULL) ;
  /* QBX_FUNCTION_NAME(qbx_quadrature_optimal_points)(w, */
  /* 						   0.25*w/(1 << dmax), */
  /* 						   /\* 16.0*w/(1 << dmax), *\/ */
  /* 						   0.5*w, */
  /* 						   64, oq, nq, */
  /* 						   Nmax, dmax, tol, */
  /* 						   x, n, x, */
  /* 						   &rp, &N, &d) ; */
  /* rp = 0.0625*w ; d = dmax ; */
  d = dmax ; rp = w/(1 << (d)) ;
  /* if ( rp < 1e-6 ) rp = 1e-6 ; */
  /* rp *= 2 ; */
  rp = w/64 ;
  QBX_FUNCTION_NAME(select_quadrature)(w, tol, nq, oq, 32, rp, rp, &N, &d) ;
  fprintf(stderr, "rp = %lg; N = %d; d = %d\n", rp, N, d) ;
  QBX_FUNCTION_NAME(qbx_triangle_truncation_error)(xe, xstr, ne, q, nq, N,
						   rp, d, ee) ;
  fprintf(stderr, "Truncation error: %e %e\n", ee[0], ee[1]) ;
  fprintf(stderr, "Quadrature error: %e\n",
	  QBX_FUNCTION_NAME(qbx_quadrature_error)(N, w/(1<<d), oq, rp, rp)) ;
  
  /*expansion centre*/
  if ( in ) {
    qbx_vector_shift(c,x,n,-rp) ;
  } else {
    qbx_vector_shift(c,x,n, rp) ;
  }
  
  /*set up nodal quantities for call to expansion generation*/
  fstr = ne ;
  memset(fe, 0, ne*fstr*sizeof(QBX_REAL)) ;
  for ( i = 0 ; i < ne ; i ++ ) {
    fe[i*fstr+i] = 1.0 ;
  }

  str = ne ;

  if ( Id != NULL ) {
#ifdef QBX_SINGLE_PRECISION
    qbx_quadrature_func_f_t func = laplace_quad_dl_f ;
#else /*QBX_SINGLE_PRECISION*/
    qbx_quadrature_func_t func = laplace_quad_dl ;
#endif /*QBX_SINGLE_PRECISION*/
    gpointer data[4] ;
    QBX_REAL st[] = {0.0, 0.0, 1.0, 0.0, 0.0, 1.0,
		     0.5, 0.0, 0.5, 0.5, 0.0, 0.5} ;
    
    data[0] = xe ; data[1] = &xstr ; data[2] = &ne ;
    memset(work, 0, ne*(N+1)*(N+2)/2*sizeof(QBX_REAL)) ;
    QBX_FUNCTION_NAME(qbx_triangle_adaptive)(func,
					     xe, xstr, ne, xe, xstr, st,
					     w, c, N, d, q, nq, oq, tol,
					     work, str, ne, data) ;
    QBX_FUNCTION_NAME(qbx_expansion_eval_laplace)(c, N, work, str, x, Id, ee) ;
  }

  if ( Is != NULL ) {
#ifdef QBX_SINGLE_PRECISION
    qbx_quadrature_func_f_t func = laplace_quad_sl_f ;
#else /*QBX_SINGLE_PRECISION*/
    qbx_quadrature_func_t func = laplace_quad_sl ;
#endif /*QBX_SINGLE_PRECISION*/
    gpointer data[4] ;
    QBX_REAL st[] = {0.0, 0.0, 1.0, 0.0, 0.0, 1.0,
		     0.5, 0.0, 0.5, 0.5, 0.0, 0.5} ;
    
    data[0] = xe ; data[1] = &xstr ; data[2] = &ne ;
    memset(work, 0, ne*(N+1)*(N+2)/2*sizeof(QBX_REAL)) ;
    QBX_FUNCTION_NAME(qbx_triangle_adaptive)(func,
					     xe, xstr, ne, xe, xstr, st,
					     w, c, N, d, q, nq, oq, tol,
					     work, str, ne, data) ;
    QBX_FUNCTION_NAME(qbx_expansion_eval_laplace)(c, N, work, str, x, Is, ee) ;

  }

  return 0 ;
}
					  
gint QBX_FUNCTION_NAME(qbx_triangle_laplace_quad)(QBX_REAL *xe,
						  gint xstr,
						  gint ne,
						  QBX_REAL *x,
						  QBX_REAL *q,
						  gint nq, gint oq,
						  gint Nmax, gint dmax,
						  QBX_REAL tol,
						  QBX_REAL *Is, gint istr,
						  QBX_REAL *Id, gint dstr,
						  QBX_REAL *work)

{
  QBX_REAL rp, w, c[3], x0[3], n[3], J, ee[64], s0, t0, ntol, kn, rx ;
#ifdef QBX_SINGLE_PRECISION
  QBX_REAL *qa = WANDZURA_25_F ;
#else /*QBX_SINGLE_PRECISION*/
  QBX_REAL *qa = WANDZURA_25 ;
#endif /*QBX_SINGLE_PRECISION*/
  gint nqa = 25, N, d, i, str, nimax, ni ;

  /*tolerance for location query*/
  ntol = tol ; nimax = 128 ;
  ntol = 1e-9 ;
  
  /*element dimension*/
  w = SQRT(4*QBX_FUNCTION_NAME(qbx_element_area)(xe, xstr, ne, qa, nqa)/M_PI) ;
  /*location of nearest point on element*/
  ni = QBX_FUNCTION_NAME(qbx_element_nearest_point)(xe, xstr, ne,
						    x, &s0, &t0,
						    x0, ntol, nimax, &kn) ;
  /* g_assert(ni < nimax) ; */
  if ( !(ni < nimax) ) {
    for ( i = 0 ; i < ne ; i ++ ) {
      fprintf(stderr, "%lg %lg %lg\n",
	      xe[xstr*i+0], xe[xstr*i+1], xe[xstr*i+2]) ;
    }
    fprintf(stderr, "%lg %lg %lg\n", x[0], x[1], x[2]) ;
    g_error("%s: point location failed to converge, %d iterations",
	    __FUNCTION__, ni) ;
  }

  QBX_FUNCTION_NAME(qbx_element_point_3d)(xe, xstr, ne, s0, t0, x0, n, &J,
  					  NULL) ;

  rx = qbx_vector_distance(x,x0) ;
  QBX_FUNCTION_NAME(qbx_quadrature_optimal_points)(w,
						   0, 2*rx,
						   16, oq, nq,
						   Nmax, dmax, tol,
						   x0, n, x,
						   &rp, &N, &d) ;
  
  /*expansion centre same side of element as field point*/
  if ( qbx_vector_diff_scalar(x,x0,n) < 0 ) {
    qbx_vector_shift(c, x0, n, -rp) ;
  } else {
    qbx_vector_shift(c, x0, n,  rp) ;
  }

  str = ne ;

  if ( Id != NULL ) {
#ifdef QBX_SINGLE_PRECISION
    qbx_quadrature_func_f_t func = laplace_quad_dl_f ;
#else /*QBX_SINGLE_PRECISION*/
    qbx_quadrature_func_t func = laplace_quad_dl ;
#endif /*QBX_SINGLE_PRECISION*/
    gpointer data[4] ;
    QBX_REAL st[] = {0.0, 0.0, 1.0, 0.0, 0.0, 1.0,
		     0.5, 0.0, 0.5, 0.5, 0.0, 0.5} ;
    
    data[0] = xe ; data[1] = &xstr ; data[2] = &ne ;
    memset(work, 0, ne*(N+1)*(N+2)/2*sizeof(QBX_REAL)) ;
    QBX_FUNCTION_NAME(qbx_triangle_adaptive)(func,
					     xe, xstr, ne, xe, xstr, st,
					     w, c, N, d, q, nq, oq, tol,
					     work, str, ne, data) ;
    QBX_FUNCTION_NAME(qbx_expansion_eval_laplace)(c, N, work, str, x, Id, ee) ;
  }

  if ( Is != NULL ) {
#ifdef QBX_SINGLE_PRECISION
    qbx_quadrature_func_f_t func = laplace_quad_sl_f ;
#else /*QBX_SINGLE_PRECISION*/
    qbx_quadrature_func_t func = laplace_quad_sl ;
#endif /*QBX_SINGLE_PRECISION*/
    gpointer data[4] ;
    QBX_REAL st[] = {0.0, 0.0, 1.0, 0.0, 0.0, 1.0,
		     0.5, 0.0, 0.5, 0.5, 0.0, 0.5} ;
    
    data[0] = xe ; data[1] = &xstr ; data[2] = &ne ;
    memset(work, 0, ne*(N+1)*(N+2)/2*sizeof(QBX_REAL)) ;
    QBX_FUNCTION_NAME(qbx_triangle_adaptive)(func,
					     xe, xstr, ne, xe, xstr, st,
					     w, c, N, d, q, nq, oq, tol,
					     work, str, ne, data) ;
    QBX_FUNCTION_NAME(qbx_expansion_eval_laplace)(c, N, work, str, x, Is, ee) ;

  }

  return 0 ;
}
