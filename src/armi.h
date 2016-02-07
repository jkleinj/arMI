/*==============================================================================
armi : analytical residual Mutual Information
(C) 2015 Jens Kleinjung
Read the COPYING file for license information.
==============================================================================*/

#ifndef ARMI_H
#define ARMI_H

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

#include "arg.h"
#include "config.h"
#include "float.h"
#include "getseqs.h"
#include "matrix.h"
#include "safe.h"
#include "seq.h"

# endif 
