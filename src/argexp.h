/*==============================================================================
argexp.h : parse command line arguments for expressions
Copyright (C) 2016 Jens Kleinjung
Read the COPYING file for license information.
==============================================================================*/

#ifndef ARGEXP_H
#define ARGEXP_H

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
    FILE *exprfile;
	char *exprfilename;
    FILE *rowfile;
	char *rowfilename;
    FILE *colfile;
	char *colfilename;
	int nrow; /* number of rows (genes) */
	int ncol; /* number of columns (samples) */
	int nsubset; /* number of genes in subset */
	int nele; /* number of expression levels */
	char *prefix; /* output prefix */
} Arg;

/*____________________________________________________________________________*/
/* prototypes */
int parse_args(int argc, char **argv, Arg *arg);

#endif
