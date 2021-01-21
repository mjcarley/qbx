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

#include <blaswrap.h>

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

static gint QBX_FUNCTION_NAME(laplace_quad_sl)(QBX_REAL s, QBX_REAL t,
					       QBX_REAL w,
					       QBX_REAL *x, QBX_REAL *y,
					       QBX_REAL *n,
					       gint N,
					       /* QBX_REAL *C, gint str, */
					       /* gint nf, */
					       gpointer data[])
{
  gint ne = *((gint *)(data[2])) ;
  QBX_REAL *C = data[3] ;
  gint str = *((gint *)(data[4])) ;
  QBX_REAL r, th, ph, L[32], sc ;
  QBX_REAL Cth, Sth, Cmph[64], Smph[64], *Pnm1, *Pn, P[128] ;
  gint nn, mm, idx, j, nf ;

  nf = ne ;
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
					       /* QBX_REAL *C, gint str, */
					       /* gint nf, */
					       gpointer data[])
{
  gint ne = *((gint *)(data[2])) ;
  QBX_REAL *C = data[3] ;
  gint str = *((gint *)(data[4])) ;
  QBX_REAL r, th, ph, L[32], rn, nR, sc, fr, fi ;
  QBX_REAL Cth, Sth, Cmph[64], Smph[64], dP, nC, nph, *Pnm1, *Pn, P[128] ;
  gint nn, mm, idx, j, nf ;

  nf = ne ;
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
  if ( !single ) {
#ifdef QBX_SINGLE_PRECISION
    qbx_quadrature_func_f_t func = laplace_quad_dl_f ;
#else /*QBX_SINGLE_PRECISION*/
    qbx_quadrature_func_t func = laplace_quad_dl ;
#endif /*QBX_SINGLE_PRECISION*/
    gpointer data[8] ;
    QBX_REAL st[] = {0.0, 0.0, 1.0, 0.0, 0.0, 1.0,
		     0.5, 0.0, 0.5, 0.5, 0.0, 0.5} ;
    
    data[0] = xe ; data[1] = &xstr ; data[2] = &ne ;
    data[3] = L ; data[4] = &str ;
    /* memset(work, 0, ne*(N+1)*(N+2)/2*sizeof(QBX_REAL)) ; */
    QBX_FUNCTION_NAME(qbx_triangle_adaptive)(func,
					     xe, xstr, ne, xe, xstr, st,
					     w, xc, N, depth, q, nq, order, tol,
					     /* L, str, ne, */
					     data) ;
  }

  if ( single ) {
    /* fprintf(stderr, "Hello\n") ; */
#ifdef QBX_SINGLE_PRECISION
    qbx_quadrature_func_f_t func = laplace_quad_sl_f ;
#else /*QBX_SINGLE_PRECISION*/
    qbx_quadrature_func_t func = laplace_quad_sl ;
#endif /*QBX_SINGLE_PRECISION*/
    gpointer data[8] ;
    QBX_REAL st[] = {0.0, 0.0, 1.0, 0.0, 0.0, 1.0,
		     0.5, 0.0, 0.5, 0.5, 0.0, 0.5} ;
    
    data[0] = xe ; data[1] = &xstr ; data[2] = &ne ;
    data[3] = L ; data[4] = &str ;
    /* memset(work, 0, ne*(N+1)*(N+2)/2*sizeof(QBX_REAL)) ; */
    QBX_FUNCTION_NAME(qbx_triangle_adaptive)(func,
					     xe, xstr, ne, xe, xstr, st,
					     w, xc, N, depth, q, nq, order, tol,
					     data) ;
  }

  return 0 ;
}

static gint QBX_FUNCTION_NAME(laplace_quad_weights)(QBX_REAL s, QBX_REAL t,
						    QBX_REAL w,
						    QBX_REAL *x, QBX_REAL *y,
						    QBX_REAL *n,
						    gint N,
						    gpointer data[])
{
  QBX_REAL *Kq  = data[QBX_DATA_MATRIX] ;
  QBX_REAL *Knm = data[QBX_DATA_KNM] ;
  gint nq = *((gint *)(data[QBX_DATA_NKNM])) ;
  gint Nk = *((gint *)(data[QBX_DATA_ORDER_K])) ;
  QBX_REAL *ws  = data[QBX_DATA_WEIGHTS_S] ;  
  QBX_REAL *wd = data[QBX_DATA_WEIGHTS_D] ;
  QBX_REAL R, wt, G, dG ;
  QBX_REAL d1 = 1.0 ;
  gint i1 = 1 ;

  /*Koornwinder polynomials at evaluation point*/
  qbx_koornwinder_nm(Nk, s, t, 1, nq, Knm) ;
  R = qbx_vector_distance(x, y) ;

  G = 0.25*M_1_PI/R ;

  dG = qbx_vector_diff_scalar(x, y, n)/R/R*G ;

#ifndef QBX_SINGLE_PRECISION  
  wt = w*G ;
  blaswrap_dgemv(TRUE, nq, nq, wt, Kq, nq, Knm, i1, d1, ws, i1) ;
  wt = w*dG ;
  blaswrap_dgemv(TRUE, nq, nq, wt, Kq, nq, Knm, i1, d1, wd, i1) ;
#else /*QBX_SINGLE_PRECISION*/
  g_assert_not_reached() ;
#endif
  
  return 0 ;
}

gint QBX_FUNCTION_NAME(qbx_laplace_weights)(QBX_REAL *xe,
					    gint xstr, gint ne,
					    QBX_REAL *q,
					    gint nq, gint oq,
					    QBX_REAL *Kq, gint Nk,
					    QBX_REAL *xt,
					    QBX_REAL *ws, QBX_REAL *wd,
					    gint depth,
					    QBX_REAL tol,
					    QBX_REAL w)

{
  QBX_REAL st[] = {0.0, 0.0, 1.0, 0.0, 0.0, 1.0,
		   0.5, 0.0, 0.5, 0.5, 0.0, 0.5} ;
  QBX_REAL Knm[512] ;
  gpointer data[QBX_DATA_WIDTH] ;
  
  data[QBX_DATA_ELEMENT]   = xe ; 
  data[QBX_DATA_STRIDE]    = &xstr ;
  data[QBX_DATA_NUMBER]    = &(ne) ;
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
					   w, xt, 0, depth, q, nq, oq, tol,
					   data) ;
  return 0 ;
}
