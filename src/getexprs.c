/*===============================================================================
getexprs.c : Read expression data (normalised read counts) 
(C) 2016 Jens Kleinjung
Read the COPYING file for license information.
================================================================================*/

#include "getexprs.h"

#define max(a,b) (((a)>(b))?(a):(b))
#define min(a,b) (((a)<(b))?(a):(b))

/*____________________________________________________________________________*/
/* read rownames of data matrix, given in a single line */
int get_rownames(FILE *rowfile, Expr *expr)
{
	unsigned int i = 0;
	unsigned allocated = 64;

	expr->rowname = safe_malloc(allocated * sizeof(char [64]));

	while(! feof(rowfile)) {
        if (fscanf(rowfile, "%s", &(expr->rowname[i][0])) == 1) {
#ifdef DEBUG
			fprintf(stderr, "%d %s\n", i, &(expr->rowname[i][0]));
#endif
			++ i; /* number of rownames */

			if (i == allocated) {
				expr->rowname = safe_realloc(expr->rowname, (allocated += 64) * sizeof(char [64]));
			}
		}
	}
	expr->nrow = i;

	return i;
}

/*____________________________________________________________________________*/
/* read colnames of data matrix, given in a single line */
int get_colnames(FILE *colfile, Expr *expr)
{
	unsigned int i = 0;
	unsigned allocated = 64;

	expr->colname = safe_malloc(allocated * sizeof(char [64]));

	while(! feof(colfile)) {
        if (fscanf(colfile, "%s", &(expr->colname[i][0])) == 1) {
#ifdef DEBUG
			fprintf(stderr, "%d %s\n", i, &(expr->colname[i][0]));
#endif
			
			++ i; /* number of colnames */

			if (i == allocated) {
				expr->colname = safe_realloc(expr->colname, (allocated += 64) * sizeof(char [64]));
			}
		}
	}
	expr->ncol = i;

	return i;
}

/*____________________________________________________________________________*/
void read_expression(FILE *exprfile, Expr *expr)
{
	unsigned int row = 0;
	unsigned int col = 0;
	unsigned int nDat = 0;

	expr->read = alloc_mat2D_float(expr->read, expr->nrow, expr->ncol);
	expr->readmin = FLT_MAX;
	expr->readmax = FLT_MIN;

	while(! feof(exprfile)) {
        /* scan expression data */
        if (fscanf(exprfile, "%f", &(expr->read[row][col])) == 1) {
			/* determine min/max values */
			expr->readmin = min(expr->readmin, expr->read[row][col]);
			expr->readmax = max(expr->readmax, expr->read[row][col]);

#ifdef DEBUG
			fprintf(stderr, "row %d, col, %d, value %f, readmin %f, readmax %f\n",
					row, col, expr->read[row][col],
					expr->readmin, expr->readmax);
#endif

			++ col;
			++ nDat;
		}
		if (col == expr->ncol) {
			++ row;
			col = 0;
		}
		if (row == expr->nrow) {
			break;
		}
	}

#ifdef DEBUG
	fprintf(stderr, "row %d\tcol %d\tnDat %d\n", row, col, nDat);
#endif

	assert((nDat == (expr->nrow * expr->ncol)) && "Dimensions of value matrix and column/row names unequal");

	expr->ndat = nDat;
}

