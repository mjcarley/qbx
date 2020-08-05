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
#include <string.h>
#include <unistd.h>
#include <math.h>

#include <glib.h>

#include <qbx.h>

#include "qbx-private.h"

gchar *tests[] = {"planar_test",
		  "closest_point",
		  "self_test",
		  "off_test",
		  ""} ;

GTimer *timer ;

static gint parse_test(gchar *arg)

{
  gint i = 0 ;
  
  while ( strlen(tests[i]) != 0) {
    if ( !strcmp(tests[i], arg) ) return i ;
    i ++ ;
  }

  return -1 ;
}

gdouble distance(gdouble *x, gdouble *y)

{
  gdouble r ;

  r = sqrt((x[0] - y[0])*(x[0] - y[0]) +
	   (x[1] - y[1])*(x[1] - y[1]) +
	   (x[2] - y[2])*(x[2] - y[2]))  ;
  
  return r ;
}

static gint read_element(FILE *f, gdouble *xe, gint *ne, gint *xstr)

{
  gint i ;

  *xstr = 4 ;
  fscanf(f, "%d", ne) ;

  for ( i = 0 ; i < (*ne) ; i ++ ) {
    fscanf(f, "%lg", &(xe[i*(*xstr)+0])) ;
    fscanf(f, "%lg", &(xe[i*(*xstr)+1])) ;
    fscanf(f, "%lg", &(xe[i*(*xstr)+2])) ;
  }
  
  return 0 ;
}

static gint planar_triangle_shape_test(gdouble *xe, gint xstr, gint ne,
				       gint nq, gint N,
				       gdouble *x0, gdouble s0, gdouble t0,
				       gdouble rc, gint depth, gdouble tol,
				       gint nx)

{
  gdouble *q, C[8192], xc[3], c[3], n[3], J, fe[128], x[3], f[64], g[64] ;
  gdouble ee[16], w, rcopt, r0 ;
  gint order, str, i, Nopt, smax, sopt, fstr, nf ;
  
  g_assert(ne == 3) ;
  nf = ne ;
  
  smax = 7 ;
  
  nf = 3 ; fstr = 5 ;

  str = 2*nf ;
  qbx_quadrature_select(nq, &q, &order) ;

  if ( x0 == NULL ) {
    qbx_element_point_3d(xe, xstr, ne, s0, t0, xc, n, &J, NULL) ;
    x[0] = xc[0] ; x[1] = xc[1] ; x[2] = xc[2] ;    
  } else {
    qbx_element_point_3d(xe, xstr, ne, 0.25, 0.25, xc, n, &J, NULL) ;
    x[0] = x0[0] ; x[1] = x0[1] ; x[2] = x0[2] ;
  }

  w = sqrt(4*qbx_element_area(xe, xstr, ne, q, nq)/M_PI) ;
  /* qbx_quadrature_optimal(w, 0.5*w/(1 << smax), 4.0*w/(1 << smax), */
  /* 			 16, order, nq, 20, smax, tol, */
  /* 			 &rcopt, &Nopt, &sopt) ; */

  r0 = distance(x, xc) ;
  if ( r0 < 1e-12 ) {
    qbx_quadrature_optimal_points(w, 0.5*w/(1 << smax), 4.0*w/(1 << smax),
				  16, order, nq, 20, smax, tol,
				  xc, n, x,
				  &rcopt, &Nopt, &sopt) ;
  } else {
    qbx_quadrature_optimal_points(w, 0.125*r0, 2*r0,
				  16, order, nq, 20, smax, tol,
				  xc, n, x,
				  &rcopt, &Nopt, &sopt) ;
  }

  c[0] = xc[0] + rc*n[0] ; 
  c[1] = xc[1] + rc*n[1] ; 
  c[2] = xc[2] + rc*n[2] ; 
  
  fprintf(stderr, "planar triangle test\n") ;
  fprintf(stderr, "====================\n") ;

  fprintf(stderr, "quadrature: %d points, %dth order\n", nq, order) ;
  fprintf(stderr, "subdivision depth: %d %d\n", depth, sopt) ;
  fprintf(stderr, "N: %d %d\n", N, Nopt) ;
  fprintf(stderr, "rc: %lg %lg\n", rc, rcopt) ;
  fprintf(stderr, "xc: %lg %lg %lg\n", xc[0], xc[1], xc[2]) ;
  fprintf(stderr, "n:  %lg %lg %lg\n", n[0], n[1], n[2]) ;
  fprintf(stderr, "c:  %lg %lg %lg\n", c[0], c[1], c[2]) ;
  
  /*set nodal quantities to shape function*/
  memset(fe, 0, 3*fstr*sizeof(gdouble)) ;
  fe[0*fstr+0] = 1.0 ; fe[1*fstr+1] = 1.0 ; fe[2*fstr+2] = 1.0 ;

  fprintf(stderr, "generating expansion, t=%lg\n",
	  g_timer_elapsed(timer, NULL)) ;
  qbx_expansion_make_laplace_adaptive(xe, xstr, ne, fe, fstr, nf,
				      q, nq, order, c, rc, N,
				      C, str, depth, tol, w,
				      TRUE, TRUE) ;
  qbx_expansion_make_laplace_adaptive(xe, xstr, ne, fe, fstr, nf,
  				      q, nq, order, c, rc, N,
  				      &(C[nf]), str, depth, tol, w,
  				      FALSE, FALSE) ;
  fprintf(stderr, "expansion generated, t=%lg\n",
	  g_timer_elapsed(timer, NULL)) ;

  /*evaluate on surface and find error*/
  memset(f, 0, 2*ne*sizeof(gdouble)) ;
  qbx_expansion_eval_laplace(c, N, C, str, x, f, ee) ;
  newman_tri_shape(x, &(xe[xstr*0]), &(xe[xstr*1]), &(xe[xstr*2]), NULL, 0,
		   &(g[0]), &(g[ne])) ;
  for ( i = 0 ; i < 2*ne ; i ++ ) g[i] *= -1.0/4.0/M_PI ;
  if ( s0 != G_MAXDOUBLE ) {
    g[ne  ] = (1 - s0 - t0)/2 ; g[ne+1] = s0/2 ; g[ne+2] = t0/2 ;
  }
  
  fprintf(stderr, "expansion:") ;
  for ( i = 0 ; i < nf ; i ++ ) fprintf(stderr, " %+lg", f[i]) ;
  for ( i = 0 ; i < nf ; i ++ ) fprintf(stderr, " %+lg", f[nf+i]) ;
  fprintf(stderr, "\n") ;
  fprintf(stderr, "exact:    ") ;
  for ( i = 0 ; i < 2*ne ; i ++ ) fprintf(stderr, " %+lg", g[i]) ;
  fprintf(stderr, "\n") ;
  fprintf(stderr, "error:    ") ;
  for ( i = 0 ; i < nf ; i ++ )
    fprintf(stderr, " %+lg", fabs(f[i]-g[i])) ;
  for ( i = 0 ; i < nf ; i ++ )
    fprintf(stderr, " %+lg", fabs(f[nf+i]-g[ne+i])) ;
  fprintf(stderr, "\n") ;

  return 0 ;
}

static gint self_test(gdouble *xe, gint xstr, gint ne,
		      gint nq, gint N,
		      gdouble *x0, gdouble s0, gdouble t0,
		      gdouble rc, gint depth, gdouble tol,
		      gint nx)

{
  gdouble work[8192], n[3], J, x[3], f[64], g[64], *q, fs[8], gs[8] ;
  gint order, i, nf ;
  
  nf = ne ;

  qbx_quadrature_select(nq, &q, &order) ;

  fprintf(stderr, "triangle self-point test\n") ;
  fprintf(stderr, "========================\n") ;

  fprintf(stderr, "quadrature: %d points, %dth order\n", nq, order) ;
  fprintf(stderr, "subdivision depth: %d\n", depth) ;
  fprintf(stderr, "N: %d\n", N) ;

  fprintf(stderr, "evaluating integrals, t=%lg\n",
	  g_timer_elapsed(timer, NULL)) ;
  qbx_triangle_laplace_self_quad(xe, xstr, ne, s0, t0, q, nq, order,
				 N, depth, tol,
				 &(f[0]), 1, &(f[ne]), 1, work) ;
  fprintf(stderr, "integrals evaluated, t=%lg\n",
	  g_timer_elapsed(timer, NULL)) ;

  qbx_element_point_3d(xe, xstr, ne, s0, t0, x, n, &J, NULL) ;
  newman_tri_shape(x, &(xe[xstr*0]), &(xe[xstr*1]), &(xe[xstr*2]), NULL, 0,
		   &(g[0]), &(g[ne])) ;
  for ( i = 0 ; i < 2*ne ; i ++ ) g[i] *= -1.0/4.0/M_PI ;
  if ( s0 != G_MAXDOUBLE ) {
    qbx_element_shape_3d(ne, s0, t0, &(g[nf]), NULL, NULL, NULL, NULL, NULL) ;
    for ( i = 0 ; i < ne ; i ++ ) g[nf+i] *= 0.5 ;
  }

  memset(fs, 0, ne*sizeof(gdouble)) ;
  memset(gs, 0, ne*sizeof(gdouble)) ;

  for ( i = 0 ; i < ne ; i ++ ) {
    fs[0] += f[i]*1.0 ;
    fs[1] += f[i]*xe[i*xstr+0] ;
    fs[2] += f[i]*xe[i*xstr+1] ;
    gs[0] += g[i]*1.0 ;
    gs[1] += g[i]*xe[i*xstr+0] ;
    gs[2] += g[i]*xe[i*xstr+1] ;
  }
  
  fprintf(stderr, "single layer\n") ;
  fprintf(stderr, "expansion:") ;
  for ( i = 0 ; i < nf ; i ++ ) fprintf(stderr, " %+lg", f[i]) ;
  fprintf(stderr, "\n") ;
  fprintf(stderr, "exact:    ") ;
  for ( i = 0 ; i < nf ; i ++ ) fprintf(stderr, " %+lg", g[i]) ;
  fprintf(stderr, "\n") ;
  fprintf(stderr, "error:    ") ;
  for ( i = 0 ; i < nf ; i ++ ) fprintf(stderr, " %+lg", fabs(f[i]-g[i])) ;
  fprintf(stderr, "\n") ;
  fprintf(stderr, "source:   ") ;
  for ( i = 0 ; i < 3 ; i ++ ) fprintf(stderr, " %+lg", fs[i]) ;
  fprintf(stderr, "\n") ;
  fprintf(stderr, "          ") ;
  for ( i = 0 ; i < 3 ; i ++ ) fprintf(stderr, " %+lg", gs[i]) ;
  fprintf(stderr, "\n") ;
  fprintf(stderr, "          ") ;
  for ( i = 0 ; i < 3 ; i ++ ) fprintf(stderr, " %+lg", fabs(gs[i]-fs[i])) ;
  fprintf(stderr, "\n") ;
  
  memset(fs, 0, ne*sizeof(gdouble)) ;
  memset(gs, 0, ne*sizeof(gdouble)) ;

  for ( i = 0 ; i < ne ; i ++ ) {
    fs[0] += f[nf+i]*1.0 ;
    fs[1] += f[nf+i]*xe[i*xstr+0] ;
    fs[2] += f[nf+i]*xe[i*xstr+1] ;
    gs[0] += g[nf+i]*1.0 ;
    gs[1] += g[nf+i]*xe[i*xstr+0] ;
    gs[2] += g[nf+i]*xe[i*xstr+1] ;
  }

  fprintf(stderr, "double layer\n") ;
  fprintf(stderr, "expansion:") ;
  for ( i = 0 ; i < nf ; i ++ ) fprintf(stderr, " %+lg", f[nf+i]) ;
  fprintf(stderr, "\n") ;
  fprintf(stderr, "exact:    ") ;
  for ( i = 0 ; i < nf ; i ++ ) fprintf(stderr, " %+lg", g[nf+i]) ;
  fprintf(stderr, "\n") ;
  fprintf(stderr, "error:    ") ;
  for ( i = 0 ; i < nf ; i ++ )
    fprintf(stderr, " %+lg", fabs(f[nf+i]-g[nf+i])) ;
  fprintf(stderr, "\n") ;
  fprintf(stderr, "source:   ") ;
  for ( i = 0 ; i < 3 ; i ++ ) fprintf(stderr, " %+lg", fs[i]) ;
  fprintf(stderr, "\n") ;
  fprintf(stderr, "          ") ;
  for ( i = 0 ; i < 3 ; i ++ ) fprintf(stderr, " %+lg", gs[i]) ;
  fprintf(stderr, "\n") ;
  fprintf(stderr, "          ") ;
  for ( i = 0 ; i < 3 ; i ++ ) fprintf(stderr, " %+lg", fabs(gs[i]-fs[i])) ;
  fprintf(stderr, "\n") ;

  return 0 ;
}

static gint off_test(gdouble *xe, gint xstr, gint ne,
		     gint nq, gint N,
		     gdouble *x0, gdouble s0, gdouble t0,
		     gdouble rc, gint depth, gdouble tol,
		     gint nx)

{
  gdouble work[8192], n[3], J, x[3], f[64], g[64], *q, fs[8], gs[8] ;
  gint order, i, nf ;
  
  nf = ne ;

  qbx_quadrature_select(nq, &q, &order) ;

  fprintf(stderr, "triangle self-point test\n") ;
  fprintf(stderr, "========================\n") ;

  fprintf(stderr, "quadrature: %d points, %dth order\n", nq, order) ;
  fprintf(stderr, "subdivision depth: %d\n", depth) ;
  fprintf(stderr, "N: %d\n", N) ;

  qbx_element_point_3d(xe, xstr, ne, s0, t0, x, n, &J, NULL) ;
  qbx_vector_shift(x,x,n,rc) ;

  fprintf(stderr, "x: %lg %lg %lg\n", x[0], x[1], x[2]) ;
  
  fprintf(stderr, "evaluating integrals, t=%lg\n",
	  g_timer_elapsed(timer, NULL)) ;
  qbx_triangle_laplace_quad(xe, xstr, ne, x, q, nq, order,
			    N, depth, tol,
			    &(f[0]), 1, &(f[ne]), 1, work) ;
  fprintf(stderr, "integrals evaluated, t=%lg\n",
	  g_timer_elapsed(timer, NULL)) ;

  newman_tri_shape(x, &(xe[xstr*0]), &(xe[xstr*1]), &(xe[xstr*2]), NULL, 0,
		   &(g[0]), &(g[ne])) ;
  for ( i = 0 ; i < 2*ne ; i ++ ) g[i] *= -1.0/4.0/M_PI ;

  memset(fs, 0, ne*sizeof(gdouble)) ;
  memset(gs, 0, ne*sizeof(gdouble)) ;

  for ( i = 0 ; i < ne ; i ++ ) {
    fs[0] += f[i]*1.0 ;
    fs[1] += f[i]*xe[i*xstr+0] ;
    fs[2] += f[i]*xe[i*xstr+1] ;
    gs[0] += g[i]*1.0 ;
    gs[1] += g[i]*xe[i*xstr+0] ;
    gs[2] += g[i]*xe[i*xstr+1] ;
  }
  
  fprintf(stderr, "single layer\n") ;
  fprintf(stderr, "expansion:") ;
  for ( i = 0 ; i < nf ; i ++ ) fprintf(stderr, " %+lg", f[i]) ;
  fprintf(stderr, "\n") ;
  fprintf(stderr, "exact:    ") ;
  for ( i = 0 ; i < nf ; i ++ ) fprintf(stderr, " %+lg", g[i]) ;
  fprintf(stderr, "\n") ;
  fprintf(stderr, "error:    ") ;
  for ( i = 0 ; i < nf ; i ++ ) fprintf(stderr, " %+lg", fabs(f[i]-g[i])) ;
  fprintf(stderr, "\n") ;
  fprintf(stderr, "source:   ") ;
  for ( i = 0 ; i < 3 ; i ++ ) fprintf(stderr, " %+lg", fs[i]) ;
  fprintf(stderr, "\n") ;
  fprintf(stderr, "          ") ;
  for ( i = 0 ; i < 3 ; i ++ ) fprintf(stderr, " %+lg", gs[i]) ;
  fprintf(stderr, "\n") ;
  fprintf(stderr, "          ") ;
  for ( i = 0 ; i < 3 ; i ++ ) fprintf(stderr, " %+lg", fabs(gs[i]-fs[i])) ;
  fprintf(stderr, "\n") ;
  
  memset(fs, 0, ne*sizeof(gdouble)) ;
  memset(gs, 0, ne*sizeof(gdouble)) ;

  for ( i = 0 ; i < ne ; i ++ ) {
    fs[0] += f[nf+i]*1.0 ;
    fs[1] += f[nf+i]*xe[i*xstr+0] ;
    fs[2] += f[nf+i]*xe[i*xstr+1] ;
    gs[0] += g[nf+i]*1.0 ;
    gs[1] += g[nf+i]*xe[i*xstr+0] ;
    gs[2] += g[nf+i]*xe[i*xstr+1] ;
  }

  fprintf(stderr, "double layer\n") ;
  fprintf(stderr, "expansion:") ;
  for ( i = 0 ; i < nf ; i ++ ) fprintf(stderr, " %+lg", f[nf+i]) ;
  fprintf(stderr, "\n") ;
  fprintf(stderr, "exact:    ") ;
  for ( i = 0 ; i < nf ; i ++ ) fprintf(stderr, " %+lg", g[nf+i]) ;
  fprintf(stderr, "\n") ;
  fprintf(stderr, "error:    ") ;
  for ( i = 0 ; i < nf ; i ++ ) fprintf(stderr, " %+lg", fabs(f[nf+i]-g[nf+i])) ;
  fprintf(stderr, "\n") ;
  fprintf(stderr, "source:   ") ;
  for ( i = 0 ; i < 3 ; i ++ ) fprintf(stderr, " %+lg", fs[i]) ;
  fprintf(stderr, "\n") ;
  fprintf(stderr, "          ") ;
  for ( i = 0 ; i < 3 ; i ++ ) fprintf(stderr, " %+lg", gs[i]) ;
  fprintf(stderr, "\n") ;
  fprintf(stderr, "          ") ;
  for ( i = 0 ; i < 3 ; i ++ ) fprintf(stderr, " %+lg", fabs(gs[i]-fs[i])) ;
  fprintf(stderr, "\n") ;

  return 0 ;
}

static gint element_closest_point_test(gdouble *xe, gint xstr, gint ne,
				       gdouble s0,  gdouble t0,
				       gdouble *x0, gdouble rc)

{
  gdouble xc[3], n[3], xf[3], J, s, t, xn[3], r[3], kn, km ;
  gint ni ;
  
  fprintf(stderr, "closest point test\n") ;
  fprintf(stderr, "===================\n") ;

  fprintf(stderr, "element: %d nodes\n", ne) ;
  fprintf(stderr, "rc: %lg\n", rc) ;
  fprintf(stderr, "s,t: %lg %lg\n", s0, t0) ;

  /*closest point on element for test purposes*/
  qbx_element_point_3d(xe, xstr, ne, s0, t0, xc, n, &J, NULL) ;  
  fprintf(stderr, "xc: %lg %lg %lg\n", xc[0], xc[1], xc[2]) ;
  fprintf(stderr, "n:  %lg %lg %lg\n", n[0], n[1], n[2]) ;
  qbx_triangle_curvature(xe, xstr, ne, s0, t0, &kn, &km) ;
  fprintf(stderr, "kn: %lg\n", kn) ;
  
  xf[0] = xc[0] + n[0]*rc ;
  xf[1] = xc[1] + n[1]*rc ;
  xf[2] = xc[2] + n[2]*rc ;

  fprintf(stderr, "xf: %lg %lg %lg\n", xf[0], xf[1], xf[2]) ;

  /* s = 1/3.0 ; t = 1/3.0 ; */
  /* s = t = 0.0 ; */
  s = 1.0 ; t = 0.0 ;
  /* s = s0 ; t = t0+0.000001 ; */
  ni = qbx_element_nearest_point(xe, xstr, ne, xf, &s, &t, xn, 1e-9, 256, &kn) ;
  
  qbx_element_point_3d(xe, xstr, ne, s, t, xc, n, &J, NULL) ;  

  r[0] = xf[0] - xc[0] ; r[1] = xf[1] - xc[1] ; r[2] = xf[2] - xc[2] ;
  fprintf(stderr, "s,t: %lg %lg\n", s, t) ;
  fprintf(stderr, "ni:  %d\n", ni) ;
  fprintf(stderr, "x0:  %lg %lg %lg\n", xc[0], xc[1], xc[2]) ;
  fprintf(stderr, "kn:  %lg\n", kn) ;
  fprintf(stderr, "n:   %lg %lg %lg\n", n[0], n[1], n[2]) ;
  fprintf(stderr, "r:   %lg %lg %lg\n", r[0], r[1], r[2]) ;
  fprintf(stderr, "r.n: %lg (%lg)\n",
	  qbx_vector_scalar(r,n), qbx_vector_length(r)) ;
  
  return 0 ;
}
				       

gint main(gint argc, gchar **argv)

{
  gdouble xe[256], tol, rc, s0, t0, x0[3] ;
  gint ne, nq, depth, N, nx, xstr, test ;
  FILE *input ;
  gchar ch, *progname ;

  timer = g_timer_new() ;

  test = -1 ;
  depth = 0 ; tol = 1e-6 ; N = 8 ; nq = 54 ; rc = 0.05 ; nx = 33 ;
  s0 = t0 = G_MAXDOUBLE ;
  
  progname = g_strdup(g_path_get_basename(argv[0])) ;

  input = stdin ;

  while ( (ch = getopt(argc, argv, "d:e:n:N:q:r:s:T:t:w:")) != EOF ) {
    switch (ch) {
    default: g_assert_not_reached() ; break ;
    case 'd': depth = atoi(optarg) ; break ;
    case 'e': tol  = atof(optarg) ; break ;
    case 'n': nx  = atoi(optarg) ; break ;
    case 'N': N   = atoi(optarg) ; break ;
    case 'q': nq  = atoi(optarg) ; break ;
    case 'r': rc  = atof(optarg) ; break ;
    case 's': s0  = atof(optarg) ; break ;
    case 'T': test = parse_test(optarg) ; break ;      
    case 't': t0  = atof(optarg) ; break ;
    }
  }

  if ( s0 != G_MAXDOUBLE && t0 == G_MAXDOUBLE ) t0 = s0 ;
  if ( t0 != G_MAXDOUBLE && s0 == G_MAXDOUBLE ) s0 = t0 ;
  
  read_element(input, xe, &ne, &xstr) ;
  fscanf(input, "%lg %lg %lg", &(x0[0]), &(x0[1]), &(x0[2])) ;
	 
  fprintf(stderr, "%s: %d node element\n", progname, ne) ;

  if ( test == 0 ) {
    if ( s0 != G_MAXDOUBLE ) 
      planar_triangle_shape_test(xe, xstr, ne, nq, N, NULL, s0, t0, rc,
				 depth, tol, nx) ;
    else
      planar_triangle_shape_test(xe, xstr, ne, nq, N, x0, s0, t0, rc,
				 depth, tol, nx) ;
    return 0 ;
  }

  if ( test == 1 ) {
    element_closest_point_test(xe, xstr, ne, s0, t0, x0, rc) ;
    
    return 0 ;
  }

  if ( test == 2 ) {
    self_test(xe, xstr, ne, nq, N, NULL, s0, t0, rc, depth, tol, nx) ;

    return 0 ;
  }

  if ( test == 3 ) {
    off_test(xe, xstr, ne, nq, N, NULL, s0, t0, rc, depth, tol, nx) ;

    return 0 ;
  }
  
  return 0 ;
}

