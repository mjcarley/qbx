#include <stdio.h>
#include <math.h>

#include <glib.h>

#include <qbx.h>

#include "qbx-private.h"

gint QBX_FUNCTION_NAME(qbx_cartesian_to_spherical)(QBX_REAL *x0,
						   QBX_REAL *x,
						   QBX_REAL *r,
						   QBX_REAL *th,
						   QBX_REAL *ph)

{
  *r = qbx_vector_distance2(x, x0) ;
  if ( *r == 0.0 ) { *ph = *th = 0.0 ; return 0 ; }

  *r = SQRT((*r)) ;
  *ph = ATAN2(x[1]-x0[1], x[0]-x0[0]) ;

  *th = ACOS((x[2]-x0[2])/(*r)) ;

  return 0 ;
}

/*
  inputs: Pnm1, Pn, normalized Legendre functions for n-1, and n

  \bar{P}_{n}^{m} = (-1)^m 
  \sqrt((n-m)!/(n+m)!\times (2n+1)/4\pi)*P_n^|m|(\cos\theta) 

  C = \cos\theta, S = \sin\theta

  on output: arrays will contain Legendre functions for n, n+1 (note
  this means the arrays are swapped, which is why they are passed as
  pointers to pointers) so they have the same sense, but with n
  incremented

  there is no check on array bounds
*/

gint QBX_FUNCTION_NAME(qbx_legendre_recursion_array)(QBX_REAL **Pnm1,
						     QBX_REAL **Pn,
						     gint n,
						     QBX_REAL C,
						     QBX_REAL S)
  
{
  gint m ;
  QBX_REAL *pn = *Pn, *pnm1 = *Pnm1, sq2np1 ;

  /*Cheng, Crutchfield, Gimbutas, Greengard, et al. normalization*/
  /* pnm1[n+1] = -sqrt((2.0*n+3)/(2.0*n+2))*S*pn[n] ; */

  /*Gumerov and Duraiswami normalization*/
  pnm1[n+1] = SQRT((2.0*n+3)/(2.0*n+2))*S*pn[n] ;
  pn  [n+1] = 0.0 ;
  sq2np1 = SQRT(2.0*n+1) ;

  for ( m = 0 ; m <= n ; m ++ ) {
    pnm1[m] = C*sq2np1*pn[m] - SQRT((QBX_REAL)(n*n-m*m)/(2*n-1))*pnm1[m] ;
    pnm1[m] /= SQRT((QBX_REAL)((n+1)*(n+1)-m*m)/(2*n+3)) ;
  }

  /*swap the pointers*/
  *Pn = pnm1 ; *Pnm1 = pn ;

  return 0 ;
}

gint QBX_FUNCTION_NAME(qbx_legendre_init)(QBX_REAL C, QBX_REAL S, 
					  QBX_REAL *P0, QBX_REAL *P10,
					  QBX_REAL *P11)

{
  /*Cheng, Crutchfield, Gimbutas, Greengard, et al. normalization*/
  /* *P0 = 1.0 ; */
  /* *P10 = C*sqrt(3.0) ; */
  /* *P11 = -S*sqrt(1.5) ; */

  /*Gumerov and Duraiswami normalization*/
  *P0  = 1.0 ;
  *P10 = C*SQRT(3.0) ;
  *P11 = S*SQRT(1.5) ;
  *P0 /= SQRT(4*M_PI) ;
  *P10 /= SQRT(4*M_PI) ;
  *P11 /= SQRT(4*M_PI) ;

  return 0 ;
}
