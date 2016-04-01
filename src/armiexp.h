/*==============================================================================
armiexp.h : analytical residual Mutual Information for expressions
(C) 2016 Jens Kleinjung
Read the COPYING file for license information.
==============================================================================*/

#ifndef ARMIEXP_H
#define ARMIEXP_H

#include <assert.h>
#include <limits.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <sys/stat.h>
#include <time.h>

#include <gsl/gsl_blas.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_permutation.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_sf_gamma.h>
#include <gsl/gsl_vector.h>

#include "argexp.h"
#include "config.h"
#include "expr.h"
#include "float.h"
#include "getexprs.h"
#include "matrix.h"
#include "safe.h"

# endif 
