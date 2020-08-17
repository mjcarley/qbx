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
					  
gint QBX_FUNCTION_NAME(qbx_expansion_rule_laplace)(QBX_REAL *xe,
						   gint xstr, gint ne,						   		   QBX_REAL *q,
						   gint nq, gint order,
						   QBX_REAL *xc,
						   QBX_REAL rc, gint N,
						   gint depth,
						   QBX_REAL w,
						   QBX_REAL tol,
						   QBX_REAL *rule,
						   gint nrmax)

{
  gint nr ;

  nr = 0 ;

  memset(rule, 0, 5*nrmax*sizeof(QBX_REAL)) ;
  
  return nr ;
}
  
