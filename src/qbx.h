#ifndef QBX_H_INCLUDED
#define QBX_H_INCLUDED

#include <glib.h>

typedef gint (*qbx_quadrature_func_t)(gdouble s, gdouble t, gdouble w,
				      gdouble *x, gdouble *y, gdouble *n,
				      gint N,
				      gdouble *fq, gint fstr, gint nf,
				      gpointer data) ;
typedef gint (*qbx_quadrature_func_f_t)(gfloat s, gfloat t, gfloat w,
					gfloat *x, gfloat *y, gfloat *n,
					gint N,
					gfloat *fq, gint fstr, gint nf,
					gpointer data) ;

gint qbx_expansion_make_laplace_sl(gdouble *xe, gint xstr, gint ne,
				   gdouble *fe, gint fstr, gint nf,
				   gdouble *qrule, gint ngp,
				   gdouble *xc, gdouble rc, gint N,
				   gdouble *C, gint str) ;
gint qbx_expansion_make_laplace_dl(gdouble *xe, gint xstr,  gint ne,
				   gdouble *fe, gint fstr, gint nf,
				   gdouble *qrule, gint ngp,
				   gdouble *xc, gdouble rc,
				   gint N,
				   gdouble *L, gint str) ;
gint qbx_expansion_make_laplace_adaptive(gdouble *xe, gint xstr, gint ne,
					 gdouble *fe, gint fstr, gint nf,
					 gdouble *qrule, gint ngp, gint order,
					 gdouble *xc,
					 gdouble rc, gint N,
					 gdouble *L, gint str, gint depth,
					 gdouble tol,
					 gdouble w,
					 gboolean init, gboolean single) ;
gint qbx_cartesian_to_spherical(gdouble *x0,
				gdouble *x,
				gdouble *r,
				gdouble *th,
				gdouble *ph) ;
gint qbx_legendre_recursion_array(gdouble **Pnm1,
				  gdouble **Pn,
				  gint n,
				  gdouble C,
				  gdouble S) ;
gint qbx_legendre_init(gdouble C, gdouble S, 
		       gdouble *P0, gdouble *P10,
		       gdouble *P11) ;
gint qbx_expansion_eval_laplace(gdouble *xc, gint N, gdouble *L, gint str,
				gdouble *x, gdouble *f, gdouble *ee) ;
gint qbx_quad_3d_real(gdouble *xe, gint ne, gdouble *x,
		      gdouble *qrule, gint nq,
		      gdouble *Iq) ;
gdouble qbx_minimum_distance(gdouble *xe, gint ne, gdouble *x) ;
gdouble qbx_maximum_distance(gdouble *xe, gint ne, gdouble *x) ;
gdouble qbx_quadrature_error(gint N, gdouble h, gint ngp,
			     gdouble r, gdouble rp) ;
gint qbx_triangle_truncation_error(gdouble *xe, gint xstr, gint ne,
				   gdouble *qrule, gint nq,
				   gint N,
				   gdouble rc, gint depth, gdouble *ee) ;
gdouble qbx_element_area(gdouble *xe, gint xstr, gint ne,
			 gdouble *qrule, gint nq) ;
gint qbx_truncation_optimal(gdouble Rb, gdouble c, gdouble rp, gint order,
			    gint pmax, gint smax,
			    gdouble tol, gint *pq, gint *s) ;
gint qbx_quadrature_optimal(gdouble Rb, gdouble r0, gdouble r1, gint r,
			    gint order, gint ngp,
			    gint pmax, gint smax, gdouble tol,
			    gdouble *rc, gint *pq, gint *s) ;

gint qbx_element_shape_3d(gint ne, gdouble s, gdouble t,
			  gdouble *L, gdouble *dLds, gdouble *dLdt,
			  gdouble *dLdss, gdouble *dLdst, gdouble *dLdtt) ;
gint qbx_element_point_3d(gdouble *xe, gint xstr, gint ne,
			  gdouble s, gdouble t,
			  gdouble *y, gdouble *n,
			  gdouble *J, gdouble *c) ;
gint qbx_element_point_interp_3d(gdouble *xe, gint xstr, gint ne,
				 gdouble *L, gdouble *dLds, gdouble *dLdt,
				 gdouble *y, gdouble *n,
				 gdouble *J, gdouble *c) ;

gint qbx_quadrature_select(gint nq, gdouble **q, gint *order) ;



gint qbx_expansion_make_laplace_sl_f(gfloat *xe, gint xstr, gint ne,
				   gfloat *fe, gint fstr, gint nf,
				   gfloat *qrule, gint ngp,
				   gfloat *xc, gfloat rc, gint N,
				   gfloat *C, gint str) ;
gint qbx_expansion_make_laplace_dl_f(gfloat *xe, gint xstr,  gint ne,
				   gfloat *fe, gint fstr, gint nf,
				   gfloat *qrule, gint ngp,
				   gfloat *xc, gfloat rc,
				   gint N,
				   gfloat *L, gint str) ;
gint qbx_expansion_make_laplace_adaptive_f(gfloat *xe, gint xstr, gint ne,
					   gfloat *fe, gint fstr, gint nf,
					   gfloat *qrule, gint ngp, gint order,
					   gfloat *xc,
					   gfloat rc, gint N,
					   gfloat *L, gint str, gint depth,
					   gfloat tol,
					   gfloat w,
					   gboolean init, gboolean single) ;
gint qbx_triangle_expansion(gdouble *xe, gint xstr, gint ne,
			    gdouble *fe, gint fstr, gint nf,
			    gdouble *q, gint nq, gint oq,
			    gdouble *x, gdouble tol,
			    gint Nmax, gint smax,
			    gdouble *c, gdouble *C,
			    gint *N, gint *s) ;
gint qbx_triangle_expansion_f(gfloat *xe, gint xstr, gint ne,
			      gfloat *fe, gint fstr, gint nf,
			      gfloat *q, gint nq, gint oq,
			      gfloat *x, gfloat tol,
			      gint Nmax, gint smax,
			      gfloat *c, gfloat *C,
			      gint *N, gint *s) ;
gint qbx_triangle_laplace_self_quad(gdouble *xe, gint xstr,
				    gint ne,
				    gdouble s0, gdouble t0,
				    gboolean in,
				    gdouble *q, gint nq, gint oq,
				    gint Nmax, gint dmax,
				    gdouble tol,
				    gdouble *Is, gint istr,
				    gdouble *Id, gint dstr,
			    gdouble *work) ;
gint qbx_triangle_laplace_self_quad_f(gfloat *xe, gint xstr,
				      gint ne,
				      gfloat s0, gfloat t0,
				      gboolean in,
				      gfloat *q, gint nq, gint oq,
				      gint Nmax, gint dmax,
				      gfloat tol,
				      gfloat *Is, gint istr,
				      gfloat *Id, gint dstr,
				      gfloat *work) ;
gint qbx_triangle_laplace_quad(gdouble *xe, gint xstr, gint ne,
			       gdouble *x,
			       gdouble *q, gint nq, gint oq,
			       gint Nmax, gint dmax,
			       gdouble tol,
			       gdouble *Is, gint istr,
			       gdouble *Id, gint dstr,
			       gdouble *work) ;
gint qbx_triangle_laplace_quad_f(gfloat *xe, gint xstr, gint ne,
				 gfloat *x,
				 gfloat *q, gint nq, gint oq,
				 gint Nmax, gint dmax,
				 gfloat tol,
				 gfloat *Is, gint istr,
				 gfloat *Id, gint dstr,
				 gfloat *work) ;

gint qbx_cartesian_to_spherical_f(gfloat *x0,
				gfloat *x,
				gfloat *r,
				gfloat *th,
				gfloat *ph) ;
gint qbx_legendre_recursion_array_f(gfloat **Pnm1,
				  gfloat **Pn,
				  gint n,
				  gfloat C,
				  gfloat S) ;
gint qbx_legendre_init_f(gfloat C, gfloat S, 
		       gfloat *P0, gfloat *P10,
		       gfloat *P11) ;
gint qbx_expansion_eval_laplace_f(gfloat *xc, gint N, gfloat *L, gint str,
				gfloat *x, gfloat *f, gfloat *ee) ;
gfloat qbx_quadrature_error_f(gint N, gfloat h, gint ngp,
			     gfloat r, gfloat rp) ;
gint qbx_triangle_truncation_error_f(gfloat *xe, gint xstr, gint ne,
				   gfloat *qrule, gint nq,
				   gint N,
				   gfloat rc, gint depth, gfloat *ee) ;
gfloat qbx_element_area_f(gfloat *xe, gint xstr, gint ne,
			 gfloat *qrule, gint nq) ;
gint qbx_truncation_optimal_f(gfloat Rb, gfloat c, gfloat rp, gint order,
			      gint pmax, gint smax,
			      gfloat tol, gint *pq, gint *s) ;
gint qbx_quadrature_optimal_f(gfloat Rb, gfloat r0, gfloat r1, gint r,
			      gint order, gint ngp,
			      gint pmax, gint smax, gfloat tol,
			      gfloat *rc, gint *pq, gint *s) ;

gint qbx_element_shape_3d_f(gint ne, gfloat s, gfloat t,
			    gfloat *L, gfloat *dLds, gfloat *dLdt,
			    gfloat *dLdss, gfloat *dLdst, gfloat *dLdtt) ;
gint qbx_element_point_3d_f(gfloat *xe, gint xstr, gint ne,
			  gfloat s, gfloat t,
			  gfloat *y, gfloat *n,
			  gfloat *J, gfloat *c) ;
gint qbx_element_point_interp_3d_f(gfloat *xe, gint xstr, gint ne,
				 gfloat *L, gfloat *dLds, gfloat *dLdt,
				 gfloat *y, gfloat *n,
				 gfloat *J, gfloat *c) ;
gint qbx_triangle_curvature(gdouble *xe, gint xstr, gint ne,
			    gdouble s, gdouble t,
			    gdouble *k1, gdouble *k2) ;
gint qbx_triangle_curvature_f(gfloat *xe, gint xstr, gint ne,
			      gfloat s, gfloat t,
			      gfloat *k1, gfloat *k2) ;

gint qbx_quadrature_select_f(gint nq, gfloat **q, gint *order) ;
gint qbx_quadrature_optimal_points(gdouble Rb, gdouble r0, gdouble r1,
				   gint nr, gint order, gint ngp,
				   gint pmax, gint smax,
				   gdouble tol,
				   gdouble *x0, gdouble *n, gdouble *x,
				   gdouble *rc, gint *pq, gint *s) ;
gint qbx_quadrature_optimal_points_f(gfloat Rb, gfloat r0, gfloat r1,
				     gint nr, gint order, gint ngp,
				     gint pmax, gint smax,
				     gfloat tol,
				     gfloat *x0, gfloat *n, gfloat *x,
				     gfloat *rc, gint *pq, gint *s) ;
gint qbx_element_nearest_point(gdouble *xe, gint xstr,
			       gint ne,
			       gdouble *xf, 
			       gdouble *sn, gdouble *tn, 
			       gdouble *xn, gdouble tol,
			       gint nimax, gdouble *kn) ;
gint qbx_element_nearest_point_f(gfloat *xe, gint xstr,
				 gint ne,
				 gfloat *xf, 
				 gfloat *sn, gfloat *tn, 
				 gfloat *xn, gfloat tol,
				 gint nimax, gfloat *kn) ;
gint qbx_koornwinder_nm(gint N, gdouble u, gdouble v, gint str, gint imax,
			gdouble *Knm) ;
gint qbx_koornwinder_nm_f(gint N, gfloat u, gfloat v, gint str, gint imax,
			  gfloat *Knm) ;

gint qbx_koornwinder_interp_matrix(gdouble *q, gint nq, gdouble *A) ;
gint qbx_koornwinder_interp_matrix_f(gfloat *q, gint nq, gfloat *A) ;

gint qbx_laplace_ts_integrate(gdouble *xe, gint xstr, gint ne,
			      gdouble *q,  gint nq, gint order,
			      gdouble *xc, gdouble rc, gint N,
			      gdouble s0, gdouble t0,
			      gdouble *L,  gint str,
			      gint depth,  gdouble tol, gdouble w) ;
gint qbx_laplace_ts_integrate_f(gfloat *xe, gint xstr, gint ne,
				gfloat *q,  gint nq, gint order,
				gfloat *xc, gfloat rc, gint N,
				gfloat s0, gfloat t0,
				gfloat *L,  gint str,
				gint depth,  gfloat tol, gfloat w) ;

gint qbx_triangle_adaptive(qbx_quadrature_func_t func,
			   gdouble *xt, gint tstr, gint ne,
			   gdouble *xd, gint dstr,
			   gdouble *st,
			   gdouble w, gdouble *xc,
			   gint N,
			   gint depth,
			   gdouble *q, gint nq, gint oq,
			   gdouble tol,
			   gdouble *f, gint fstr, gint nf,
			   gpointer data) ;
gint qbx_triangle_adaptive_f(qbx_quadrature_func_f_t func,
			     gfloat *xt, gint tstr, gint ne,
			     gfloat *xd, gint dstr,
			     gfloat *st,
			     gfloat w,
			     gfloat *xc, gint N,
			     gint depth,
			     gfloat *q, gint nq, gint oq,
			     gfloat tol,
			     gfloat *f, gint fstr, gint nf,
			     gpointer data) ;


#endif /*QBX_H_INCLUDED*/
