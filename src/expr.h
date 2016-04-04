/*===============================================================================
expr.h : expresion data structure
Copyright (C) 2016 Jens Kleinjung
Read the COPYING file for license information.
================================================================================*/

#ifndef EXPR_H
#define EXPR_H

/*___________________________________________________________________________*/
/** data structures */

/* sequence */
typedef struct  
{
	int nrow; /* number of rows */
	int ncol; /* number of columns */
	char (*rowname)[64];
	char (*colname)[64];
	float **read; /* expression values */
	int **level; /* expression value level '(int)roundf(log2(read))' */
	int ndat; /* number of expression values */
	float min; /* minimum expression value */
	float max; /* maximum expression value */
} Expr; 

#endif

