ACLOCAL_AMFLAGS = ${ACLOCAL_FLAGS} -I m4
AM_DISTCHECK_CONFIGURE_FLAGS="--enable-unit-tests"

.PHONY: all

if !BUILD_WITH_LIBDIVSUFSORT
OPT_PSUFSORT = opt/psufsort
endif

SUBDIRS = . $(OPT_PSUFSORT) libs opt src
DIST_SUBDIRS = . libs opt src opt/psufsort

