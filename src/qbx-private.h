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


/*BLAS macros*/
extern gint qbx_0i[], qbx_1i[], qbx_2i[] ;
extern gdouble qbx_0z[] ;
extern gdouble qbx_1z[] ;
extern gdouble qbx_m1z[] ;
extern gdouble qbx_0d[], qbx_1d[], qbx_m1d[] ;

extern void dgemv_(gchar *trans, gint *m, gint *n, gdouble *alpha,
		   gdouble *A, gint *lda, gdouble *v, gint *incx,
		   gdouble *beta, gdouble *y, gint *incy) ;
extern void zgemv_(gchar *trans, gint *m, gint *n, gdouble *alpha,
		   gdouble *A, gint *lda, gdouble *v, gint *incx,
		   gdouble *beta, gdouble *y, gint *incy) ;

extern void dgetri_(gint *n, gdouble *A, gint *lda, gint *ip,
		    gdouble *work, gint *lwork, gint *info) ;
extern void zgetri_(gint *n, gdouble *A, gint *lda, gint *ip,
		    gdouble *work, gint *lwork, gint *info) ;

extern void zgetrf_(gint *m, gint *n, gdouble *A, gint *lda, gint *ip,
		    gint *info) ;
extern void dgetrf_(gint *m, gint *n, gdouble *A, gint *lda, gint *ip,
		    gint *info) ;

extern void dgemm_(gchar *transa, gchar *transb,
		   gint *m, gint *n, gint *k, gdouble *alpha,
		   gdouble *A, gint *lda,
		   gdouble *B, gint *ldb,
		   gdouble *beta, gdouble *C, gint *ldc) ;
extern void zgemm_(gchar *transa, gchar *transb,
		   gint *m, gint *n, gint *k, gdouble *alpha,
		   gdouble *A, gint *lda,
		   gdouble *B, gint *ldb,
		   gdouble *beta, gdouble *C, gint *ldc) ;

extern gdouble dscal_(gint *n, gdouble *da, gdouble *dx, gint *incx) ;

extern gdouble dasum_ (gint *n, gdouble *x, gint *incx) ;
extern gdouble dzasum_(gint *n, gdouble *x, gint *incx) ;
extern gdouble dnrm2_ (gint *n, gdouble *x, gint *incx) ;
extern gdouble dznrm2_(gint *n, gdouble *x, gint *incx) ;
extern gint    idamax_(gint *n, gdouble *x, gint *incx) ;
extern gint    izamax_(gint *n, gdouble *x, gint *incx) ;
extern gdouble ddot_  (gint *n, gdouble *x, gint *incx, 
		       gdouble *y, gint *incy) ;
/* extern gsl_complex zdotu_ (gint *n, gdouble *x, gint *incx,  */
/* 			   gdouble *y, gint *incy) ; */
extern void    dcopy_(gint *n, 
		      gdouble *x, const gint *incx,
		      gdouble *y, const gint *incy) ;
extern void    daxpy_(gint *n, 
		      gdouble *a, gdouble *x, gint *incx,
		      gdouble *y, gint *incy) ;
extern void    zaxpy_(gint *n, 
		      gdouble *a, gdouble *x, gint *incx,
		      gdouble *y, gint *incy) ;

/* scale x by da */
#define qbx_dscal(_n,_al,_x,_strx) dscal_((_n),(_al),(_x),(_strx))

/* sum x[i]*y[i] */

#define qbx_ddot(_n,_x,_strx,_y,_stry)		\
  ddot_((_n), (_x), (_strx), (_y), (_stry)) 

/* y := al*A*x + bt*y */

/*trans nr nc al A lda x strx bt y stry*/

#define qbx_dgemv(_t,_m,_n,_al,_A,_lda,_x,_incx,_bt,_y,_incy)		\
  do {									\
    if ( (_t) ) {							\
      g_assert_not_reached() ;						\
    } else {								\
      dgemv_("T",(_n),(_m),(_al),(_A),(_n),(_x),(_incx),		\
	     (_bt),(_y),(_incy)) ;					\
    }									\
  } while (0)
    

#endif /*QBX_PRIVATE_H_INCLUDED*/
