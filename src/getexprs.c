/*===============================================================================
getexprs.c : Read expression data (normalised read counts) 
(C) 2016 Jens Kleinjung
Read the COPYING file for license information.
================================================================================*/

#include "getexprs.h"

/*____________________________________________________________________________*/
int get_rownames(FILE *rowfile, Expr *expr)
{
	unsigned int i = 0;
	unsigned allocated = 64;
	char line[1048576];

	expr->rowname = safe_malloc(allocated * sizeof(char));

	while(fgets(line, 1048575, rowfile) != 0) {
		while (sscanf(&(line[0]), "%s", expr->rowname[i]) == 1) {
			++ i;

			if (i == allocated) {
				expr->rowname = safe_realloc(expr->rowname, (allocated += 64) * sizeof(char));
			}
		}
	}
	expr->nrow = i;

	return i;
}

/*____________________________________________________________________________*/
int get_colnames(FILE *colfile, Expr *expr)
{
	unsigned int i = 0;
	unsigned allocated = 64;
	char line[1048576];

	expr->colname = safe_malloc(allocated * sizeof(char *));

	while(fgets(line, 1048575, colfile) != 0) {
		while (sscanf(&(line[0]), "%s", expr->colname[i]) == 1) {
			++ i;

			if (i == allocated) {
				expr->colname = safe_realloc(expr->colname, (allocated += 64) * sizeof(char));
			}
		}
	}
	expr->ncol = i;

	return i;
}

/*____________________________________________________________________________*/
int read_expression(FILE *exprfile, Expr *expr)
{
	unsigned int i = 0;
	unsigned int j = 0;
	char line[4096];

	expr->read = alloc_mat2D_float(expr->read, expr->nrow, expr->ncol);

	while (fgets(line, 4096, exprfile) != 0) {	
		while (sscanf(&(line[0]), "%f", &(expr->read[i][j])) == 1) {
			++ i;
		}
		assert((i == expr->ncol) && "Expecting columns to match col file");
		i = 0;
		++ j;
	}

	assert((j == expr->nrow) && "Expecting rows to match row file");

	return 0;
}

