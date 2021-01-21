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

#ifndef QBX_PRIVATE_H_INCLUDED 
#define QBX_PRIVATE_H_INCLUDED

#include <stdio.h>

#include <glib.h>

#define QBX_DATA_WIDTH     16
#define QBX_DATA_ELEMENT    0
#define QBX_DATA_STRIDE     1
#define QBX_DATA_NUMBER     2
#define QBX_DATA_RADIUS     3
#define QBX_DATA_NORMAL     4
#define QBX_DATA_MATRIX     5
#define QBX_DATA_KNM        6
#define QBX_DATA_NKNM       7
#define QBX_DATA_ORDER_K    8
#define QBX_DATA_WEIGHTS_S  9 
#define QBX_DATA_WEIGHTS_D 10 

#ifdef QBX_SINGLE_PRECISION

#define QBX_REAL gfloat

#define QBX_FUNCTION_NAME(QBX_func) QBX_func##_f

#define SQRT(QBX_x) sqrtf((QBX_x))
#define CBRT(QBX_x) cbrtf((QBX_x))
#define SIN(QBX_x) sinf((QBX_x))
#define COS(QBX_x) cosf((QBX_x))
#define ACOS(QBX_x) acosf((QBX_x))
#define ATAN(QBX_x) atanf((QBX_x))
#define ATAN2(QBX_y,QBX_x) atan2f((QBX_y),(QBX_x))
#define LOG(QBX_x) logf((QBX_x))
#define EXP(QBX_x) expf(QBX_x)

#else

#define QBX_REAL gdouble

#define QBX_FUNCTION_NAME(QBX_func) QBX_func

#define SQRT(QBX_x) sqrt((QBX_x))
#define CBRT(QBX_x) cbrt((QBX_x))
#define SIN(QBX_x) sin((QBX_x))
#define COS(QBX_x) cos((QBX_x))
#define ACOS(QBX_x) acos((QBX_x))
#define ATAN(QBX_x) atan((QBX_x))
#define ATAN2(QBX_y,QBX_x) atan2((QBX_y),(QBX_x))
#define LOG(QBX_x) log((QBX_x))
#define EXP(QBX_x) exp(QBX_x)

#endif /*QBX_SINGLE_PRECISION*/

#define SIGN(QBX_x) ((QBX_x) < 0 ? -1 : 1)

#define qbx_cos_sin_recursion(QBX_Cn,_Sn,_C,_S)	\
  do { QBX_REAL _tmp = (QBX_Cn) ;		\
  (QBX_Cn) = (QBX_Cn)*(QBX_C) - (QBX_Sn)*(QBX_S) ;		\
  (QBX_Sn) = (QBX_Sn)*(QBX_C) + (QBX_tmp)*(QBX_S) ;		\
  } while (0)

#define qbx_vector_cross(QBX_C,QBX_A,QBX_B)					\
  ((QBX_C)[0] = (QBX_A)[1]*(QBX_B)[2] - (QBX_A)[2]*(QBX_B)[1],		\
   (QBX_C)[1] = (QBX_A)[2]*(QBX_B)[0] - (QBX_A)[0]*(QBX_B)[2],		\
   (QBX_C)[2] = (QBX_A)[0]*(QBX_B)[1] - (QBX_A)[1]*(QBX_B)[0])

#define qbx_vector_scalar(QBX_A,QBX_B)  \
  (((QBX_A)[0])*((QBX_B)[0])+					\
   ((QBX_A)[1])*((QBX_B)[1])+					\
   ((QBX_A)[2])*((QBX_B)[2]))

#define qbx_vector_diff_scalar(QBX_A,QBX_B,QBX_C)				\
  (((QBX_A)[0]-(QBX_B)[0])*((QBX_C)[0]) +				\
   ((QBX_A)[1]-(QBX_B)[1])*((QBX_C)[1]) +				\
   ((QBX_A)[2]-(QBX_B)[2])*((QBX_C)[2]))

#define qbx_vector_length(QBX_A)					\
  (SQRT(((QBX_A)[0])*((QBX_A)[0])+				\
	((QBX_A)[1])*((QBX_A)[1]) +				\
	((QBX_A)[2])*((QBX_A)[2])))

#define qbx_vector_shift(QBX_A,QBX_B,QBX_C,QBX_D)	\
  { (QBX_A)[0] = (QBX_B)[0] + (QBX_C)[0]*(QBX_D) ;	\
    (QBX_A)[1] = (QBX_B)[1] + (QBX_C)[1]*(QBX_D) ;		\
    (QBX_A)[2] = (QBX_B)[2] + (QBX_C)[2]*(QBX_D) ;		\
  } while (0)

#define qbx_vector_distance2(QBX_A,QBX_B)		\
  ( ((QBX_A)[0]-(QBX_B)[0])*((QBX_A)[0]-(QBX_B)[0]) +	\
    ((QBX_A)[1]-(QBX_B)[1])*((QBX_A)[1]-(QBX_B)[1]) +	\
    ((QBX_A)[2]-(QBX_B)[2])*((QBX_A)[2]-(QBX_B)[2]) )

#define qbx_vector_distance(QBX_A,QBX_B)		\
  (SQRT((qbx_vector_distance2(QBX_A,QBX_B))))

#define qbx_point_copy(_xb,_fb,_i,_xe,_xstr,_fe,_j)	\
  do {							\
    (_xb)[(_xstr)*(_i)+0] = (_xe)[(_xstr)*(_j)+0] ;	\
    (_xb)[(_xstr)*(_i)+1] = (_xe)[(_xstr)*(_j)+1] ;	\
    (_xb)[(_xstr)*(_i)+2] = (_xe)[(_xstr)*(_j)+2] ;	\
    (_fb)[2*(_j)+0] = (_fe)[2*(_i)+0] ;			\
    (_fb)[2*(_j)+1] = (_fe)[2*(_i)+1] ;			\
  } while ( 0 ) 

#define qbx_point_interp3(_xb,_fb,_i,_xe,_xstr,_fe,_L0,_L1,_L2)		\
  do  {									\
  (_xb)[(_xstr)*(_i)+0] =						\
    (_L0)*(_xe)[(_xstr)*0+0] + (_L1)*(_xe)[(_xstr)*1+0] +		\
    (_L2)*(_xe)[(_xstr)*2+0] ;						\
  (_xb)[(_xstr)*(_i)+1] =						\
    (_L0)*(_xe)[(_xstr)*0+1] + (_L1)*(_xe)[(_xstr)*1+1] +		\
    (_L2)*(_xe)[(_xstr)*2+1] ;						\
  (_xb)[(_xstr)*(_i)+2] =						\
    (_L0)*(_xe)[(_xstr)*0+2] + (_L1)*(_xe)[(_xstr)*1+2] +		\
    (_L2)*(_xe)[(_xstr)*2+2] ;						\
  (_fb)[2*(_i)+0] = (_L0)*(_fe)[0] + (_L1)*(_fe)[2] + (_L2)*(_fe)[4] ;	\
  (_fb)[2*(_i)+1] = (_L0)*(_fe)[1] + (_L1)*(_fe)[3] + (_L2)*(_fe)[5] ;	\
} while (0)

#define qbx_point_interp_jac3(_xe,_xstr,_L0,_L1,_L2,			\
			      _Ls0,_Ls1,_Ls2,_Lt0,_Lt1,_Lt2,_y,_n,_J)	\
  {									\
  (_y)[0] = (_L0)*(_xe)[(_xstr)*0+0] + (_L1)*(_xe)[(_xstr)*1+0] +	\
    (_L2)*(_xe)[(_xstr)*2+0] ;						\
  (_y)[1] = (_L0)*(_xe)[(_xstr)*0+1] + (_L1)*(_xe)[(_xstr)*1+1] +	\
    (_L2)*(_xe)[(_xstr)*2+1] ;						\
  (_y)[2] = (_L0)*(_xe)[(_xstr)*0+2] + (_L1)*(_xe)[(_xstr)*1+2] +	\
    (_L2)*(_xe)[(_xstr)*2+2] ;						\
  (_n)[0] =								\
    ((_Ls0)*(_xe)[(_xstr)*0+1] + (_Ls1)*(_xe)[(_xstr)*1+1] +		\
     (_Ls2)*(_xe)[(_xstr)*2+1])*					\
    ((_Lt0)*(_xe)[(_xstr)*0+2] + (_Lt1)*(_xe)[(_xstr)*1+2] +		\
     (_Lt2)*(_xe)[(_xstr)*2+2]) -					\
    ((_Lt0)*(_xe)[(_xstr)*0+1] + (_Lt1)*(_xe)[(_xstr)*1+1] +		\
     (_Lt2)*(_xe)[(_xstr)*2+1])*					\
    ((_Ls0)*(_xe)[(_xstr)*0+2] + (_Ls1)*(_xe)[(_xstr)*1+2] +		\
     (_Ls2)*(_xe)[(_xstr)*2+2]) ;					\
  (_n)[1] =								\
    ((_Ls0)*(_xe)[(_xstr)*0+2] + (_Ls1)*(_xe)[(_xstr)*1+2] +		\
     (_Ls2)*(_xe)[(_xstr)*2+2])*					\
    ((_Lt0)*(_xe)[(_xstr)*0+0] + (_Lt1)*(_xe)[(_xstr)*1+0] +		\
     (_Lt2)*(_xe)[(_xstr)*2+0]) -					\
    ((_Lt0)*(_xe)[(_xstr)*0+2] + (_Lt1)*(_xe)[(_xstr)*1+2] +		\
     (_Lt2)*(_xe)[(_xstr)*2+2])*					\
    ((_Ls0)*(_xe)[(_xstr)*0+0] + (_Ls1)*(_xe)[(_xstr)*1+0] +		\
     (_Ls2)*(_xe)[(_xstr)*2+0]) ;					\
  (_n)[2] =								\
    ((_Ls0)*(_xe)[(_xstr)*0+0] + (_Ls1)*(_xe)[(_xstr)*1+0] +		\
     (_Ls2)*(_xe)[(_xstr)*2+0])*					\
    ((_Lt0)*(_xe)[(_xstr)*0+1] + (_Lt1)*(_xe)[(_xstr)*1+1] +		\
     (_Lt2)*(_xe)[(_xstr)*2+1]) -					\
    ((_Lt0)*(_xe)[(_xstr)*0+0] + (_Lt1)*(_xe)[(_xstr)*1+0] +		\
     (_Lt2)*(_xe)[(_xstr)*2+0])*					\
    ((_Ls0)*(_xe)[(_xstr)*0+1] + (_Ls1)*(_xe)[(_xstr)*1+1] +		\
     (_Ls2)*(_xe)[(_xstr)*2+1]) ;					\
  (_J) = qbx_vector_length((_n)) ;					\
  (_n)[0] /= (_J) ; (_n)[1] /= (_J) ; (_n)[2] /= (_J) ;			\
  } while (0)


#define qbx_triangle_divide_loop30(_xe,_xstr,_fe,_xl,_fl)	\
  {								\
  qbx_point_copy((_xl), (_fl), 0, (_xe), (_xstr), (_fe), 0) ;		\
  qbx_point_interp3((_xl), (_fl), 1, (_xe), (_xstr), (_fe), 0.5, 0.5, 0.0) ;\
  qbx_point_interp3((_xl), (_fl), 2, (_xe), (_xstr), (_fe), 0.5, 0.0, 0.5) ;\
} while (0)

#define qbx_triangle_divide_loop31(_xe,_xstr,_fe,_xl,_fl)	\
  {								\
  qbx_point_copy((_xl), (_fl), 1, (_xe), (_xstr), (_fe), 1) ;		\
  qbx_point_interp3((_xl), (_fl), 0, (_xe), (_xstr), (_fe), 0.5, 0.5, 0.0) ;\
  qbx_point_interp3((_xl), (_fl), 2, (_xe), (_xstr), (_fe), 0.0, 0.5, 0.5) ;\
} while (0)

#define qbx_triangle_divide_loop32(_xe,_xstr,_fe,_xl,_fl)\
  {							       \
  qbx_point_copy((_xl), (_fl), 2, (_xe), (_xstr), (_fe), 2) ;\
  qbx_point_interp3((_xl), (_fl), 0, (_xe), (_xstr), (_fe), 0.5, 0.0, 0.5) ;\
  qbx_point_interp3((_xl), (_fl), 1, (_xe), (_xstr), (_fe), 0.0, 0.5, 0.5) ;\
} while (0)

#define qbx_triangle_divide_loop33(_xe,_xstr,_fe,_xl,_fl)	\
  {									\
  qbx_point_interp3((_xl), (_fl), 0, (_xe), (_xstr), (_fe), 0.5, 0.5, 0.0) ;\
  qbx_point_interp3((_xl), (_fl), 1, (_xe), (_xstr), (_fe), 0.0, 0.5, 0.5) ;\
  qbx_point_interp3((_xl), (_fl), 2, (_xe), (_xstr), (_fe), 0.5, 0.0, 0.5) ;\
} while (0)

#define qbx_shape3(_s,_t,_L)					\
  {(_L)[0] = 1.0-(_s)-(_t); (_L)[1] = (_s); (_L)[2] = (_t);	\
  }    while (0)

#define qbx_shape_derivatives3(_s,_t,_L,_Ls,_Lt)			\
  {(_L)[0] = 1.0 - (_s) - (_t) ; (_L)[1] = (_s) ; (_L)[2] = (_t) ;	\
  (_Ls)[0] = -1.0 ; (_Ls)[1] =  1.0 ; (_Ls)[2] =  0.0 ;			\
  (_Lt)[0] = -1.0 ; (_Lt)[1] =  0.0 ; (_Lt)[2] =  1.0 ;			\
} while (0)

#define qbx_index_laplace_nm(QBX_n,QBX_m) ((QBX_n)*(QBX_n)+(2*(QBX_m))-1)

#define IS_EVEN(QBX_i) (((QBX_i)%2==0)?1:0)

#define minus_one_pow(QBX_n) ((2*((QBX_n)/2) == (QBX_n) ? 1 : -1))

#define yes_if_true(QBX_t)  ((QBX_t) == TRUE ? "yes" : "no") 

extern gdouble WANDZURA_7[], WANDZURA_25[], WANDZURA_54[], WANDZURA_85[],
  WANDZURA_126[], WANDZURA_175[], XIAO_GIMBUTAS_453[] ;
extern gfloat WANDZURA_7_F[], WANDZURA_25_F[], WANDZURA_54_F[],
  WANDZURA_85_F[], WANDZURA_126_F[], WANDZURA_175_F[],
  XIAO_GIMBUTAS_453_F[] ;

gint newman_tri(gdouble p[], gdouble x1[], gdouble x2[], gdouble x3[],
		gdouble Iq[], gdouble J[]) ;
gint newman_tri_shape(gdouble p[], gdouble x1[], gdouble x2[], gdouble x3[],
		      gdouble *Imn, gint hmax,
		      gdouble Iq[], gdouble J[]) ;

#endif /*QBX_PRIVATE_H_INCLUDED*/
