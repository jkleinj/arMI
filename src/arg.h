/*==============================================================================
arg.h : parse command line arguments
Copyright (C) 2015 Jens Kleinjung
Read the COPYING file for license information.
==============================================================================*/

#ifndef ARG_H
#define ARG_H

#include <assert.h>
#include <float.h>
#include <getopt.h>
#include <gsl/gsl_sf_gamma.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

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
	int silent; /* no stdout */
} Arg;

/*____________________________________________________________________________*/
/* prototypes */
int parse_args(int argc, char **argv, Arg *arg);

#endif
