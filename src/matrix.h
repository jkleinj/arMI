/*==============================================================================
matrix.h
(C) 2015 Jens Kleinjung
Read the COPYING file for license information.
==============================================================================*/

#ifndef MATRIX_H
#define MATRIX_H

#include <stdlib.h>

#include "safe.h"

double ***alloc_mat3D_double(double ***mat3D_double, int x, int y, int z);
void init_mat3D_double(double ***mat3D_double, int x, int y, int z, double val);
void free_mat3D_double(double ***mat3D_double, int x, int y);

long double ***alloc_mat3D_longdouble(long double ***mat3D_longdouble, int x, int y, int z);
void init_mat3D_longdouble(long double ***mat3D_longdouble, int x, int y, int z, long double val);
void free_mat3D_longdouble(long double ***mat3D_longdouble, int x, int y);

#endif
