#                                               -*- Autoconf -*-
# Process this file with autoconf to produce a configure script.

AC_PREREQ([2.59])
AC_INIT(armi, 1.2.0, Jens.Kleinjung@crick.ac.uk)
AC_CONFIG_SRCDIR([src/armi.c])
AM_INIT_AUTOMAKE([foreign -Wall -Werror])

# Checks for programs.
AC_PROG_CC_C_O
AC_PROG_INSTALL

# Checks for libraries.
AC_CHECK_LIB([m], [main], [], 
             AC_MSG_ERROR([Could not find required C math library.]))
AC_CHECK_LIB([gslcblas], main)
AC_CHECK_LIB([gsl], main)

# Checks for header files.
AC_HEADER_STDC
AC_CHECK_HEADERS([stdlib.h stdio.h])

# Checks for typedefs, structures, and compiler characteristics.
AC_TYPE_SIZE_T

# Checks for library functions.
AC_FUNC_MALLOC
AC_FUNC_REALLOC
AC_CHECK_FUNCS([pow sqrt])

# enable debug
AC_ARG_ENABLE(debug, [  --enable-debug                 turn on debugging], [ac_use_debug="yes"], [ac_use_debug="no"])
if test "$ac_use_debug" = "yes"; then
    CFLAGS="-ggdb -O0"
    FCFLAGS="-g -O0"
    AC_DEFINE(DEBUG,1,[Define to compile with DEBUG settings.])
fi

# enable profiling 
AC_ARG_ENABLE(profiling, [  --enable-profiling             turn on profiling], [ac_use_profiling="yes"], [ac_use_profiling="no"])
if test "$ac_use_profiling" = "yes"; then
    CFLAGS="-pg -O2"
    FCFLAGS="-pg -O2"
    AC_DEFINE(PROFILING,1,[Define to compile with PROFILING settings.])
fi

# enable speed
AC_ARG_ENABLE(speed, [  --enable-speed             turn on optimisation], [ac_use_speed="yes"], [ac_use_speed="no"])
if test "$ac_use_speed" = "yes"; then
    CFLAGS="-O4"
    FCFLAGS="-O4"
    AC_DEFINE(SPEED,1,[Define to compile with SPEED settings.])
fi

# MPI commpiler search with ax_mpi.m4 
# http://www.gnu.org/software/autoconf-archive/ax_mpi.html
sinclude(ax_mpi.m4)

# enable MPI
AC_ARG_ENABLE(mpi, [  --enable-mpi                   turn on MPI],
                [ac_use_mpi="yes"], [ac_use_mpi="no"])
AC_ARG_WITH(mpi, [  --with-mpi=<path>              prefix of MPI installation],
                 [  MPICC="$withval/bin/mpicc"
                    CPPFLAGS="$CPPFLAGS -I $withval/include"
                    LDFLAGS="$LDFLAGS -L $withval/lib"])

if test "$ac_use_mpi" = "yes"; then
   AX_MPI([], AC_MSG_ERROR([Could not find required MPI library.
                         Consider using --with-mpi=<my-path-to-mpi>.]))
   CC="$MPICC"
   LIBS="$MPILIBS $LIBS"
   AC_DEFINE(MPI,1,[Compiling with MPI settings.])
fi

# GSL MACROS
AH_BOTTOM([#define RETURN_IF_NULL(x) if (!x) { return ; }])

# Output files
AC_CONFIG_FILES([Makefile src/Makefile tests/Makefile])
AC_CONFIG_HEADERS([src/config.h])
AC_OUTPUT

