/*===============================================================================
getexprs.h : Read expression data (normalised read counts) 
(C) 2016 Jens Kleinjung
Read the COPYING file for license information.
================================================================================*/

#ifndef GETEXPRS_H
#define GETEXPRS_H

#include <assert.h>
#include <ctype.h>
#include <float.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "config.h"
#include "expr.h"
#include "safe.h"
#include "matrix.h"

/*____________________________________________________________________________*/
/* prototypes */
int get_rownames(FILE *rownamesInFile, Expr *expr);
int get_colnames(FILE *colnamesInFile, Expr *expr);
void read_expression(FILE *exprInFile, Expr *expr);

#endif
