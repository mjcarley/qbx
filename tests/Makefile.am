## Process this file with automake to produce Makefile.in

AM_CPPFLAGS = -I$(top_srcdir) -I$(top_srcdir)/src -I$(includedir)

DEPS = $(top_builddir)/src/libqbx.la

#EXTRA_DIST = \
#	qbx-tests.c.in qbx-tree-calc.c.in \
#	qbx-direct-calc.c.in \
#	qbx-laplace-tests.c.in qbx-laplace-tree-calc.c.in \
#	make-source make-source-float

bin_PROGRAMS = \
	qbx-tests

#BUILT_SOURCES = \
#	qbx-tests.c \
#	qbx-tree-calc.c \
#	qbx-direct-calc.c \
#	qbx-tests-float.c \
#	qbx-tree-calc-float.c \
#	qbx-direct-calc-float.c \
#	qbx-laplace-tests.c \
##	qbx-laplace-tests-float.c \
#	qbx-laplace-tree-calc.c \
#	qbx-laplace-tree-calc-float.c \
#	qbx-laplace-direct-calc.c \
#	qbx-laplace-direct-calc-float.c

#CLEANFILES = \
#	qbx-tests.c qbx-tests-float.c \
#	qbx-tree-calc.c qbx-tree-calc-float.c \
#	qbx-direct-calc.c qbx-direct-calc-float.c \
#	qbx-laplace-tests.c qbx-laplace-tests-float.c \
#	qbx-laplace-tree-calc.c qbx-laplace-tree-calc-float.c \
#	qbx-laplace-direct-calc.c qbx-laplace-direct-calc-float.c 

#qbx-tests.c: qbx-tests.c.in
#	./make-source qbx-tests
#
#qbx-tests-float.c: qbx-tests.c.in
#	./make-source-float qbx-tests
#
#qbx-tree-calc.c: qbx-tree-calc.c.in
#	./make-source qbx-tree-calc
#
#qbx-tree-calc-float.c: qbx-tree-calc.c.in
#	./make-source-float qbx-tree-calc
#
#qbx-direct-calc.c: qbx-direct-calc.c.in
#	./make-source qbx-direct-calc
#
#qbx-direct-calc-float.c: qbx-direct-calc.c.in
#	./make-source-float qbx-direct-calc
#
#qbx-laplace-tests.c: qbx-laplace-tests.c.in
#	./make-source qbx-laplace-tests
#
#qbx-laplace-tests-float.c: qbx-laplace-tests.c.in
#	./make-source-float qbx-laplace-tests
#
#qbx-laplace-tree-calc.c: qbx-laplace-tree-calc.c.in
#	./make-source qbx-laplace-tree-calc
#
#qbx-laplace-tree-calc-float.c: qbx-laplace-tree-calc.c.in
#	./make-source-float qbx-laplace-tree-calc
#
#qbx-laplace-direct-calc.c: qbx-laplace-direct-calc.c.in
#	./make-source qbx-laplace-direct-calc
#
#qbx-laplace-direct-calc-float.c: qbx-laplace-direct-calc.c.in
#	./make-source-float qbx-laplace-direct-calc
#

qbx_tests_LDADD = \
	$(top_builddir)/src/.libs/libqbx.a 

#
#qbx_tests_avx_SOURCES = qbx-tests.c
#qbx_tests_avx_LDADD = \
#	$(top_builddir)/src/.libs/libqbx.a \
#	$(top_builddir)/src/.libs/libqbx-base.a
#
#qbx_tests_avx_CFLAGS = -DQBX_USE_AVX $(AM_CFLAGS)
#
#qbx_tests_float_LDADD = \
#	$(top_builddir)/src/.libs/libqbx_f.a \
#	$(top_builddir)/src/.libs/libqbx-base.a
#
#qbx_tree_calc_SOURCES = qbx-tree-calc.c
#qbx_tree_calc_LDADD = \
#	$(top_builddir)/src/.libs/libqbx.a \
#	$(top_builddir)/src/.libs/libqbx-base.a
#
#qbx_tree_calc_avx_SOURCES = qbx-tree-calc.c
#qbx_tree_calc_avx_LDADD = \
#	$(top_builddir)/src/.libs/libqbx.a \
#	$(top_builddir)/src/.libs/libqbx-base.a
#
#qbx_tree_calc_avx_CFLAGS = -DQBX_USE_AVX $(AM_CFLAGS)
#
#qbx_tree_calc_float_LDADD = \
#	$(top_builddir)/src/.libs/libqbx_f.a \
#	$(top_builddir)/src/.libs/libqbx-base.a
#
#qbx_direct_calc_LDADD = \
#	$(top_builddir)/src/.libs/libqbx.a \
#	$(top_builddir)/src/.libs/libqbx-base.a
#
#qbx_direct_calc_float_LDADD = \
#	$(top_builddir)/src/.libs/libqbx_f.a \
#	$(top_builddir)/src/.libs/libqbx-base.a
#
#qbx_laplace_tests_LDADD = \
#	$(top_builddir)/src/.libs/libqbx.a \
#	$(top_builddir)/src/.libs/libqbx-base.a
#
#qbx_laplace_tests_avx_SOURCES = qbx-laplace-tests.c
#qbx_laplace_tests_avx_LDADD = \
#	$(top_builddir)/src/.libs/libqbx.a \
#	$(top_builddir)/src/.libs/libqbx-base.a
#qbx_laplace_tests_avx_CFLAGS = -DQBX_USE_AVX $(AM_CFLAGS)
#
#qbx_laplace_tests_float_LDADD = \
#	$(top_builddir)/src/.libs/libqbx_f.a \
#	$(top_builddir)/src/.libs/libqbx-base.a
#
#qbx_laplace_tree_calc_LDADD = \
#	$(top_builddir)/src/.libs/libqbx.a \
#	$(top_builddir)/src/.libs/libqbx-base.a
#
#qbx_laplace_tree_calc_avx_SOURCES = qbx-laplace-tree-calc.c
#qbx_laplace_tree_calc_avx_LDADD = \
#	$(top_builddir)/src/.libs/libqbx.a \
#	$(top_builddir)/src/.libs/libqbx-base.a
#
#qbx_laplace_tree_calc_avx_CFLAGS = -DQBX_USE_AVX $(AM_CFLAGS)
#
#qbx_laplace_tree_calc_float_LDADD = \
#	$(top_builddir)/src/.libs/libqbx_f.a \
#	$(top_builddir)/src/.libs/libqbx-base.a
#
#qbx_laplace_direct_calc_LDADD = \
#	$(top_builddir)/src/.libs/libqbx.a \
#	$(top_builddir)/src/.libs/libqbx-base.a
#
#qbx_laplace_direct_calc_float_LDADD = \
#	$(top_builddir)/src/.libs/libqbx_f.a \
#	$(top_builddir)/src/.libs/libqbx-base.a
#
