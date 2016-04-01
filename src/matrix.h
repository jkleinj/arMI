/*==============================================================================
matrix.h
(C) 2015 Jens Kleinjung
Read the COPYING file for license information.
==============================================================================*/

#ifndef MATRIX_H
#define MATRIX_H

#include <stdlib.h>

#include "safe.h"

/* prototypes */
/* double */
double ***alloc_mat3D_double(double ***mat3D_double, int x, int y, int z);
void init_mat3D_double(double ***mat3D_double, int x, int y, int z, double val);
void free_mat3D_double(double ***mat3D_double, int x, int y);

/* long double */
long double ***alloc_mat3D_longdouble(long double ***mat3D_longdouble, int x, int y, int z);
void init_mat3D_longdouble(long double ***mat3D_longdouble, int x, int y, int z, long double val);
void free_mat3D_longdouble(long double ***mat3D_longdouble, int x, int y);


/* integer */
int **alloc_mat2D_int(int **mat2D_int, int x, int y);
void init_mat2D_int(int **mat2D_int, int x, int y, int val);
void free_mat2D_int(int **mat2D_int, int x);
void print_mat2D_int(char *outFileName, int **mat2D_int, int x, int y);

int ***alloc_mat3D_int(int ***mat3D_int, int x, int y, int z);
void init_mat3D_int(int ***mat3D_int, int x, int y, int z, int val);
void free_mat3D_int(int ***mat3D_int, int x, int y);

int ****alloc_mat4D_int(int ****mat4D_int, int w, int x, int y, int z);
void init_mat4D_int(int ****mat4D_int, int w, int x, int y, int z, int val);
void free_mat4D_int(int ****mat4D_int, int w, int x, int y);

/* float */
float **alloc_mat2D_float(float **mat2D_float, int x, int y);
void init_mat2D_float(float **mat2D_float, int x, int y, float val);
void free_mat2D_float(float **mat2D_float, int x);
void print_mat2D_float(char *outFileName, float **mat2D_float, int x, int y);
void div_mat2D_float(float **mat2D_float, int x, int y, float a);

float ***alloc_mat3D_float(float ***mat3D_float, int x, int y, int z);
void init_mat3D_float(float ***mat3D_float, int x, int y, int z, float val);
void free_mat3D_float(float ***mat3D_float, int x, int y);

float ****alloc_mat4D_float(float ****mat4D_float, int w, int x, int y, int z);
void init_mat4D_float(float ****mat4D_float, int w, int x, int y, int z, float val);
void free_mat4D_float(float ****mat4D_float, int w, int x, int y);

#endif

