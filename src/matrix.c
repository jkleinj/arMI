/*==============================================================================
matrix.c : matrix routines
(C) 2015 Jens Kleinjung
Read the COPYING file for license information.
==============================================================================*/

#include "matrix.h"

/*___________________________________________________________________________*/
/** 3D double matrix */
/** allocate */
double ***alloc_mat3D_double(double ***mat3D_double, int x, int y, int z)
{
	unsigned int i, j;

	mat3D_double = (double ***)safe_malloc(x * sizeof(double **));
	for (i = 0; i < x; ++ i) {
		mat3D_double[i] = (double **)safe_malloc(y * sizeof(double *));
		for (j = 0; j < y; ++ j)
			mat3D_double[i][j] = (double *)safe_malloc(z * sizeof(double));
	}

	return mat3D_double;
}

/** initialise */
void init_mat3D_double(double ***mat3D_double, int x, int y, int z, double val)
{
    unsigned int i, j, k;

	for (i = 0; i < x; ++ i)
		for (j = 0; j < y; ++ j)
			for (k = 0; k < z; ++ k)
				mat3D_double[i][j][k] = val;
}

/** free */
void free_mat3D_double(double ***mat3D_double, int x, int y)
{
    unsigned int i, j;

    for (i = 0; i < x; ++ i) {
        for (j = 0; j < y; ++ j)
            free(mat3D_double[i][j]);
        free(mat3D_double[i]);
    }
    free(mat3D_double);
}

/*___________________________________________________________________________*/
/** 3D long double matrix */
/** allocate */
long double ***alloc_mat3D_longdouble(long double ***mat3D_longdouble, int x, int y, int z)
{
	unsigned int i, j;

	mat3D_longdouble = (long double ***)safe_malloc(x * sizeof(long double **));
	for (i = 0; i < x; ++ i) {
		mat3D_longdouble[i] = (long double **)safe_malloc(y * sizeof(long double *));
		for (j = 0; j < y; ++ j)
			mat3D_longdouble[i][j] = (long double*)safe_malloc(z * sizeof(long double));
	}

	return mat3D_longdouble;
}

/** initialise */
void init_mat3D_longdouble(long double ***mat3D_longdouble, int x, int y, int z, long double val)
{
    unsigned int i, j, k;

	for (i = 0; i < x; ++ i)
		for (j = 0; j < y; ++ j)
			for (k = 0; k < z; ++ k)
				mat3D_longdouble[i][j][k] = val;
}

/** free */
void free_mat3D_longdouble(long double ***mat3D_longdouble, int x, int y)
{
    unsigned int i, j;

    for (i = 0; i < x; ++ i) {
        for (j = 0; j < y; ++ j)
            free(mat3D_longdouble[i][j]);
        free(mat3D_longdouble[i]);
    }
    free(mat3D_longdouble);
}


