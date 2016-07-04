bin_PROGRAMS = tummer

if !BUILD_WITH_LIBDIVSUFSORT
PSUFSORT=$(top_builddir)/opt/psufsort/libpsufsort.a
# This is a hack to make sure, tummer is created with a C++ compiler
DUMMY=dummy.cxx
endif

tummer_SOURCES = tummer.c esa.c process.c sequence.c io.c global.h esa.h process.h sequence.h io.h
tummer_CPPFLAGS = $(OPENMP_CFLAGS) -I$(top_srcdir)/libs -I$(top_srcdir)/opt -std=gnu99
tummer_CFLAGS = $(OPENMP_CFLAGS) -Wall -Wextra -Wno-missing-field-initializers
tummer_CXXFLAGS = $(OPENMP_CXXFLAGS) -Wall -Wextra
tummer_LDADD = $(PSUFSORT) $(top_builddir)/libs/libpfasta.a $(top_builddir)/opt/libcompat.a
nodist_EXTRA_tummer_SOURCES = $(DUMMY)

.PHONY: perf
perf: CFLAGS+= -g -O3 -ggdb -fno-omit-frame-pointer
perf: tummer
