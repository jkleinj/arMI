/*==============================================================================
argali.h : parse command line arguments for alignment
Copyright (C) 2016 Jens Kleinjung
Read the COPYING file for license information.
==============================================================================*/

#ifndef ARGALI_H
#define ARGALI_H

#include <assert.h>
#include <float.h>
#include <getopt.h>
#include <gsl/gsl_sf_gamma.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

#include "config.h"

/*____________________________________________________________________________*/
/* structures */

/* non-parametric arguments */
typedef struct  
{
    FILE *infile;
	char *infilename;
	int random;
	int nseq; /* number of sequences */
	int lseq; /* length of alignment */
	int nsubset; /* number of sequences in subset */
	int nele; /* number of alphabet elements */
	char *prefix; /* output prefix */
} Arg;

/*____________________________________________________________________________*/
/* prototypes */
int parse_args(int argc, char **argv, Arg *arg);

#endif
