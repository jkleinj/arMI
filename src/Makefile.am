#===============================================================================
# Makefile.am : analytical residual column mutual information
#===============================================================================

bin_PROGRAMS = armiali armiexp

AM_CFLAGS = -Wall

armiali_SOURCES = \
armiali.c armiali.h argali.c argali.h getseqs.c getseqs.h \
matrix.h matrix.c safe.c safe.h seq.h

armiexp_SOURCES = \
armiexp.c armiexp.h argexp.c argexp.h expr.h getexprs.c getexprs.h \
matrix.h matrix.c safe.c safe.h

CLEANFILES = $(TARGET) *.o 
DISTCLEANFILES = libtool config.cache config.log

