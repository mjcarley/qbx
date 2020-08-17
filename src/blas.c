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

extern gint qbx_0i[] = {0}, qbx_1i[] = {1}, qbx_2i[]  = {2} ;
extern gdouble qbx_0z[] = {0.0, 0.0} ;
extern gdouble qbx_1z[] = {1.0, 0.0} ;
extern gdouble qbx_m1z[] = {-1.0, 0.0} ;
extern gdouble qbx_0d[] = {0.0}, qbx_1d[] = {1.0}, qbx_m1d[] = {-1.0} ;

