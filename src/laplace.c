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

static gint interp_node_func(QBX_REAL *fe, gint fstr, gint ne, gint nf,
			     QBX_REAL *L, QBX_REAL *f)

{
  gint i, j ;

  for ( i = 0 ; i < nf ; i ++ ) {
    f[i] = L[0]*fe[i] ;
    for ( j = 1 ; j < ne ; j ++ ) {
      f[i] += L[j]*fe[j*fstr+i] ;
    }
  }
  
  return 0 ;
}

gint QBX_FUNCTION_NAME(qbx_expansion_make_laplace_sl)(QBX_REAL *xe,
						      gint xstr, gint ne,
						      QBX_REAL *fe,
						      gint fstr, gint nf,
						      QBX_REAL *q, gint nq,
						      QBX_REAL *xc,
						      QBX_REAL rc, gint N,
						      QBX_REAL *C, gint str)

/*
 * xe:   element nodes
 * xstr: stride between coordinates (xe[i*xstr+0,1,2])
 * ne:   number of element nodes
 * fe:   element data
 * fstr: stride between data (fe[i*fstr+0,1,2,...])
 * nf:   number of data entries per node
 * q:    quadrature rule
 * nq:   number of quadrature points
 * xc:   expansion centre
 * rc:   expansion radius
 * N:    order of expansion
 * C:    coefficients of expansion
 * str:  coefficient stride
*/
  
{
  gint i, mm, nn, idx, j ;
  QBX_REAL s, t, w, y[3], n[3], J, r, th, ph, P[256], *Pnm1, *Pn, sc ;
  QBX_REAL Cth, Sth, Cmph[64], Smph[64], L[16], dLds[16], dLdt[16], f[16] ;

  g_assert(str > 0) ;
  g_assert(nf <= str) ;
  /* g_assert(rc > 0) ; */
  
  Pnm1 = &(P[0]) ;
  Pn = &(Pnm1[N+4]) ;

  for ( i = 0 ; i < nq ; i ++ ) {
    s = q[3*i+0] ; t = q[3*i+1] ; w = q[3*i+2] ;
    QBX_FUNCTION_NAME(qbx_element_shape_3d)(ne, s, t, L, dLds, dLdt,
					    NULL, NULL, NULL) ;
    QBX_FUNCTION_NAME(qbx_element_point_interp_3d)(xe, xstr, ne,
						   L, dLds, dLdt, y, n, &J,
						   NULL) ;

    interp_node_func(fe, fstr, ne, nf, L, f) ;
    
    QBX_FUNCTION_NAME(qbx_cartesian_to_spherical)(xc, y, &r, &th, &ph) ;
    /* if ( r < rc-1e-14 ) { */
    /*   g_error("%s: r (%lg) < rc (%lg), difference: %lg", */
    /* 	      __FUNCTION__, r, rc, fabs(r-rc)) ; */
    /* } */
    Cth = COS(th) ; Sth = SIN(th) ;
    QBX_FUNCTION_NAME(qbx_legendre_init)(Cth, Sth, &(Pnm1[0]),
					 &(Pn[0]), &(Pn[1])) ;
    Pnm1[1] = 0.0 ;
    
    Cmph[0] = 1.0 ; Cmph[1] = COS(ph) ;
    Smph[0] = 0.0 ; Smph[1] = SIN(ph) ;
    
    w *= J/r ;
    
    nn = 0 ; mm = 0 ; sc = 1.0/(2*nn+1) ;
    idx = nn*nn ;
    for ( j = 0 ; j < nf ; j ++ ) C[str*idx+j] += Pnm1[mm]*w*sc*f[j] ;

    w /= r ;
    nn = 1 ; sc = 1.0/(2*nn+1) ;
    mm =  0 ;
    idx = nn*nn ;
    for ( j = 0 ; j < nf ; j ++ ) C[str*idx+j] += Pn[mm]*w*sc*f[j] ;
    mm =  1 ;
    idx = qbx_index_laplace_nm(nn, mm) ;
    for ( j = 0 ; j < nf ; j ++ ) {
      C[str*(idx+0)+j] += Pn[mm]*Cmph[mm]*w*sc*f[j] ;
      C[str*(idx+1)+j] -= Pn[mm]*Smph[mm]*w*sc*f[j] ;
    }
    
    for ( nn = 2 ; nn <= N ; nn ++ ) {
      QBX_FUNCTION_NAME(qbx_legendre_recursion_array)(&Pnm1, &Pn, nn-1,
						      Cth, Sth) ;
      w /= r ;
      sc = 1.0/(2*nn+1) ;      
      mm =  0 ;
      idx = nn*nn ;
      for ( j = 0 ; j < nf ; j ++ ) C[str*idx+j] += Pn[mm]*w*sc*f[j] ;

      Cmph[nn] = Cmph[nn-1]*Cmph[1] - Smph[nn-1]*Smph[1] ;
      Smph[nn] = Smph[nn-1]*Cmph[1] + Cmph[nn-1]*Smph[1] ;
      
      for ( mm = 1 ; mm <= nn ; mm ++ ) {
	idx = qbx_index_laplace_nm(nn, mm) ;
	for ( j = 0 ; j < nf ; j ++ ) {
	  C[str*(idx+0)+j] += Pn[mm]*Cmph[mm]*w*sc*f[j] ;
	  C[str*(idx+1)+j] -= Pn[mm]*Smph[mm]*w*sc*f[j] ;
	}
      }
    }
  }

  return 0 ;
}

gint QBX_FUNCTION_NAME(qbx_expansion_make_laplace_dl)(QBX_REAL *xe,
						      gint xstr, gint ne,
						      QBX_REAL *fe,
						      gint fstr, gint nf,
						      QBX_REAL *q, gint nq,
						      QBX_REAL *xc, QBX_REAL rc,
						      gint N,
						      QBX_REAL *C, gint str)

{
  gint i, mm, nn, idx, j ;
  QBX_REAL s, t, w, y[3], n[3], J, r, th, ph, P[256], *Pnm1, *Pn, sc ;
  QBX_REAL Cth, Sth, Cmph[64], Smph[64], L[16], dLds[16], dLdt[16], f[16] ;
  QBX_REAL nR, nC, nph, dP, fr, fi, rn ;
    
  g_assert(str > 0) ;
  g_assert(nf <= str) ;
  
  Pnm1 = &(P[0]) ;
  Pn = &(Pnm1[N+4]) ;

  for ( i = 0 ; i < nq ; i ++ ) {
    s = q[3*i+0] ; t = q[3*i+1] ; w = q[3*i+2] ; 
    QBX_FUNCTION_NAME(qbx_element_shape_3d)(ne, s, t, L, dLds, dLdt,
					    NULL, NULL, NULL) ;
    QBX_FUNCTION_NAME(qbx_element_point_interp_3d)(xe, xstr, ne,
						   L, dLds, dLdt, y, n, &J,
						   NULL) ;

    interp_node_func(fe, fstr, ne, nf, L, f) ;

    QBX_FUNCTION_NAME(qbx_cartesian_to_spherical)(xc, y, &r, &th, &ph) ;
    Cth = COS(th) ; Sth = SIN(th) ;
    QBX_FUNCTION_NAME(qbx_legendre_init)(Cth, Sth, &(Pnm1[0]),
					 &(Pn[0]), &(Pn[1])) ;
    Pnm1[1] = Pnm1[2] = Pn[2] = 0.0 ;
    
    Cmph[0] = 1.0 ; Cmph[1] = COS(ph) ;
    Smph[0] = 0.0 ; Smph[1] = SIN(ph) ;

    nR = -(n[0]*(xc[0]-y[0]) + n[1]*(xc[1]-y[1]) + n[2]*(xc[2]-y[2]))/r ;
    nph = (n[0]*(xc[1]-y[1]) - n[1]*(xc[0]-y[0]))/r/r/Sth/Sth ;
    nC = (n[2]-Cth*nR)/r ;

    w *= J ;

    rn = r ;
    
    nn = 0 ; mm = 0 ; sc = 1.0/(2*nn+1) ;
    dP = 0 ;
    fr = (Cmph[mm]*( dP*nC-(nn+1)*Pnm1[mm]*nR/r) -
	  mm*Smph[mm]*Pnm1[mm]*nph)/rn ;
    
    idx = nn*nn ;
    for ( j = 0 ; j < nf ; j ++ ) C[str*idx+j] += fr*w*sc*f[j] ;

    rn *= r ;
    nn = 1 ; sc = 1.0/(2*nn+1) ;
    mm = 0 ;
    dP = SQRT((2*nn+1)/(2*nn-1)*(nn+mm)*(nn-mm))*Pnm1[mm] - nn*Cth*Pn[mm] ;
    dP /= Sth*Sth ;
    fr = (Cmph[mm]*( dP*nC-(nn+1)*Pn[mm]*nR/r) - mm*Smph[mm]*Pn[mm]*nph)/rn ;

    idx = nn*nn ;
    for ( j = 0 ; j < nf ; j ++ ) C[str*idx+j] += fr*w*sc*f[j] ;

    mm = 1 ;
    dP = SQRT((2*nn+1)/(2*nn-1)*(nn+mm)*(nn-mm))*Pnm1[mm] - nn*Cth*Pn[mm] ;
    dP /= Sth*Sth ;
    fr = (Cmph[mm]*( dP*nC-(nn+1)*Pn[mm]*nR/r) - mm*Smph[mm]*Pn[mm]*nph)/rn ;
    fi = (Smph[mm]*(-dP*nC+(nn+1)*Pn[mm]*nR/r) - mm*Cmph[mm]*Pn[mm]*nph)/rn ;

    idx = qbx_index_laplace_nm(nn, mm) ;
    for ( j = 0 ; j < nf ; j ++ ) {
      C[str*(idx+0)+j] += fr*w*sc*f[j] ;
      C[str*(idx+1)+j] += fi*w*sc*f[j] ;
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
      for ( j = 0 ; j < nf ; j ++ )  C[str*idx+j] += fr*w*sc*f[j] ;

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
	  C[str*(idx+0)+j] += fr*w*sc*f[j] ;
	  C[str*(idx+1)+j] += fi*w*sc*f[j] ;
	}
      }
    }
  }

  return 0 ;
}

gint QBX_FUNCTION_NAME(qbx_expansion_eval_laplace)(QBX_REAL *xc, gint N,
						   QBX_REAL *C, gint str,
						   QBX_REAL *x, QBX_REAL *f,
						   QBX_REAL *ee)

{
  QBX_REAL th, ph, r, rn, P[64]={0.0}, *Pnm1, *Pn ;
  QBX_REAL Cth, Sth, Cmph[64], Smph[64] ;
  gint mm, nn, idx, i ;

  g_assert(str > 0) ;
  /* g_assert(ne <= str) ; */
  
  Pnm1 = &(P[0]) ;
  Pn = &(Pnm1[N+1]) ;

  QBX_FUNCTION_NAME(qbx_cartesian_to_spherical)(xc, x, &r, &th, &ph) ;
  Cth = COS(th) ; Sth = SIN(th) ;
  QBX_FUNCTION_NAME(qbx_legendre_init)(Cth, Sth,
				       &(Pnm1[0]), &(Pn[0]), &(Pn[1])) ;
  Cmph[0] = 1.0 ; Cmph[1] = COS(ph) ;
  Smph[0] = 0.0 ; Smph[1] = SIN(ph) ;
  
  rn = 1.0 ;
  nn = 0 ; mm = 0 ;
  idx = nn*nn ;
  for ( i = 0 ; i < str ; i ++ ) {
    ee[i] = C[str*idx+i]*rn*Pnm1[0] ;
  }
  for ( i = 0 ; i < str ; i ++ ) f[i] += ee[i] ;

  if ( N == 0 ) return 0 ;
  
  rn *= r ;
  nn = 1 ;
  mm =  0 ;
  idx = nn*nn ;
  for ( i = 0 ; i < str ; i ++ ) {
    ee[i] = C[str*idx+i]*rn*Pn[mm] ;
  }

  mm =  1 ;
  idx = qbx_index_laplace_nm(nn, mm) ;
  for ( i = 0 ; i < str ; i ++ ) {
    ee[i] += 2.0*(C[str*(idx+0)+i]*Cmph[mm] -
		  C[str*(idx+1)+i]*Smph[mm])*rn*Pn[mm] ;
  }
  for ( i = 0 ; i < str ; i ++ ) f[i] += ee[i] ;
  
  for ( nn = 2 ; nn <= N ; nn ++ ) {
    rn *= r ;
    QBX_FUNCTION_NAME(qbx_legendre_recursion_array)(&Pnm1, &Pn, nn-1,
						    Cth, Sth) ;
    mm =  0 ;
    idx = nn*nn ;
    for ( i = 0 ; i < str ; i ++ ) {
      ee[i] = C[str*idx+i]*rn*Pn[mm] ;
    }

    Cmph[nn] = Cmph[nn-1]*Cmph[1] - Smph[nn-1]*Smph[1] ;
    Smph[nn] = Smph[nn-1]*Cmph[1] + Cmph[nn-1]*Smph[1] ;

    for ( mm =  1 ; mm <= nn ; mm ++ ) {
      idx = qbx_index_laplace_nm(nn, mm) ;
      for ( i = 0 ; i < str ; i ++ )       
	ee[i] += 2.0*(C[str*(idx+0)+i]*Cmph[mm] -
		      C[str*(idx+1)+i]*Smph[mm])*rn*Pn[mm] ;
    }
    for ( i = 0 ; i < str ; i ++ ) f[i] += ee[i] ;
  }

  return 0 ;
}

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

gint QBX_FUNCTION_NAME(qbx_expansion_make_laplace_adaptive)(QBX_REAL *xe,
							    gint xstr, gint ne,
							    QBX_REAL *fe,
							    gint fstr, gint nf,
							    QBX_REAL *q,
							    gint nq, gint order,
							    QBX_REAL *xc,
							    QBX_REAL rc, gint N,
							    QBX_REAL *L,
							    gint str,
							    gint depth,
							    QBX_REAL tol,
							    QBX_REAL w,
							    gboolean init,
							    gboolean single) 

{
  QBX_REAL xb[128], fb[128], rp, err ;
  gint i ;

  g_assert(ne == 3 || ne == 6) ;

  g_assert(ne <= str) ;
  
  if ( init ) memset(L, 0, str*(N+1)*(N+1)*sizeof(QBX_REAL)) ;

  /*check error estimate for this triangle*/
  rp = 1e6 ;
  for ( i = 0 ; i < ne ; i ++ ) {
    rp = MIN(rp, vector_distance2(xc, &(xe[i*xstr]))) ;
  }
  rp = SQRT(rp) ;
  err = QBX_FUNCTION_NAME(qbx_quadrature_error)(N, w, order, rp, rp) ;
  
  if ( depth == 0 || err < tol ) {
    if ( single ) 
      QBX_FUNCTION_NAME(qbx_expansion_make_laplace_sl)(xe, xstr, ne,
						       fe, fstr, nf, q, nq,
						       xc, rc, N, L, str) ;
    else
      QBX_FUNCTION_NAME(qbx_expansion_make_laplace_dl)(xe, xstr, ne,
						       fe, fstr, nf, q, nq,
						       xc, rc, N, L, str) ;
    return 0 ;
  }

  for ( i = 0 ; i < 4 ; i ++ ) {
    triangle_divide_loop(xe, xstr, ne, fe, fstr, nf, i, xb, fb) ;

    QBX_FUNCTION_NAME(qbx_expansion_make_laplace_adaptive)(xb, xstr, ne,
							   fb, fstr, nf,
							   q, nq, order, xc, rc,
							   N, L, str, depth-1,
							   tol, 0.5*w,
							   FALSE, single) ;
  }

  return 0 ;
}
