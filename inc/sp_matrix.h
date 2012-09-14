/* -*- Mode: C; tab-width: 2; indent-tabs-mode: nil; c-basic-offset: 2 -*- */
/*
 Copyright (C) 2011,2012 Alexey Veretennikov (alexey dot veretennikov at gmail.com)
 
 This file is part of libspmatrix.

 libspmatrix is free software: you can redistribute it and/or modify
 it under the terms of the GNU Lesser General Public License as published
 by the Free Software Foundation, either version 3 of the License, or
 (at your option) any later version.

 libspmatrix is distributed in the hope that it will be useful,
 but WITHOUT ANY WARRANTY; without even the implied warranty of
 MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 GNU Lesser General Public License for more details.

 You should have received a copy of the GNU Lesser General Public License
 along with libspmatrix.  If not, see <http://www.gnu.org/licenses/>.
*/

#ifndef __SP_MATRIX_H__
#define __SP_MATRIX_H__

#include "sp_cont.h"

typedef enum
{
  CRS = 0,                      /* Compressed Row Storage */
  CCS = 1                       /* Compressed Column Storage */
} sparse_storage_type;

/* Constants for sp_matrix.ordered field */
enum
{
  NOT_ORDERED = 0,
  ORDERED = 1
};

typedef enum
{
  PROP_GENERAL,
  PROP_SYMMETRIC,
  PROP_SYMMETRIC_PORTRAIT,
  PROP_SKEW_SYMMETRIC
} matrix_properties;

typedef enum
{
  MTX_SAME,
  MTX_EQUAL,
  MTX_SAME_PORTRAIT,
  MTX_DIFFERENT
} matrix_comparison;

/*
 * Sparse matrix row storage 
 * Internal format based on CRS or CCS
 */
typedef struct
{
  int rows_count;
  int cols_count;
  indexed_array_ptr storage;
  int ordered;                               /* if matrix was finalized */
  sparse_storage_type storage_type;          /* Storage type */
} sp_matrix;
typedef sp_matrix* sp_matrix_ptr;

/*
 * Sparse matrix CS(lower triangle)R(Skyline) format
 * used in sparse iterative solvers
 * Constructed from Sparse Matrix with assumption of the symmetric
 * matrix portrait
 */
typedef struct
{
  int rows_count;
  int cols_count;
  int nonzeros;                 /* number of nonzero elements in matrix */
  int tr_nonzeros;              /* number of nonzero elements in
                                 * upper or lower triangles */
  double *diag;                 /* rows_count elements in matrix diagonal */
  double *lower_triangle;       /* nonzero elements of the lower triangle */
  double *upper_triangle;       /* nonzero elements of the upper triangle */
  int *jptr;                    /* array of column/row indexes of the
                                 * lower/upper triangles */
  int *iptr;                    /* array of row/column offsets in jptr
                                 * for lower or upper triangles */
} sp_matrix_skyline;
typedef sp_matrix_skyline* sp_matrix_skyline_ptr;


/*
 * Sparse matrix in 3 arrays  (CRS or CCS format):
 * offsets, column/row indicies, values
 * This format is most commonly used in various solvers
 */
typedef struct
{
  sparse_storage_type storage_type;
  int rows_count;
  int cols_count;
  int nonzeros;                 /* number of nonzero elements in matrix */
  int* offsets;                 
  int* indicies;
  double* values;
} sp_matrix_yale;
typedef sp_matrix_yale* sp_matrix_yale_ptr;

/*************************************************************/
/* Sparse matrix operations                                  */

/*
 * Initializer for a sparse matrix with specified rows and columns
 * number.
 * This function doesn't allocate the memory for the matrix itself;
 * only for its structures. Matrix mtx shall be already allocated
 * bandwdith - is a start bandwidth of a matrix row
 * type - CRS or CCS sparse matrix storage types
 */
void sp_matrix_init(sp_matrix_ptr mtx,
                    int rows,
                    int cols,
                    int bandwidth,
                    sparse_storage_type type);
/*
 * Destructor for a sparse matrix
 * This function doesn't deallocate memory for the matrix itself,
 * only for its structures.
 */
void sp_matrix_free(sp_matrix_ptr mtx);

/*
 * Clear the sparse matrix.
 * Set the element values to zero keeping sparsity portrait
 */
void sp_matrix_clear(sp_matrix_ptr mtx);

/*
 * Copy sparse matrix from mtx_from to mtx_to
 * This function assumes what mtx_to is already cleared by sp_matrix_free
 * or mtx_to is a pointer to uninitialized sp_matrix structure
 */
void sp_matrix_copy(sp_matrix_ptr mtx_from,
                    sp_matrix_ptr mtx_to);

/*
 * Converts matrix storage format CRS <=> CCS
 * mtx_to shall be uninitialized sp_matrix structure
 * if the type is the same - do nothing and return 0
 * returns nonzero if conversion was successull
 */
int sp_matrix_convert(sp_matrix_ptr mtx_from,
                      sp_matrix_ptr mtx_to,
                      sparse_storage_type type);

/*
 * Converts matrix storage format CRS <=> CCS inplace
 * if the type is the same - do nothing and return 0
 * returns nonzero if conversion was successull
 */
int sp_matrix_convert_inplace(sp_matrix_ptr self,
                              sparse_storage_type type);



/*
 * Construct CSlR sparse matrix based on sp_matrix format
 * mtx - is the (reordered) sparse matrix to take data from
 * Acts as a copy-constructor
 * Limitations: works only for CRS format
 */
void sp_matrix_skyline_init(sp_matrix_skyline_ptr self,
                            sp_matrix_ptr mtx);
/*
 * Construct CSlR sparse matrix based on sp_matrix_yale format
 * mtx - is the (reordered) sparse matrix to take data from
 * Acts as a copy-constructor
 */
void sp_matrix_skyline_yale_init(sp_matrix_skyline_ptr self,
                                 sp_matrix_yale_ptr mtx);

/*
 * Destructor for a sparse matrix in CSlR format
 * This function doesn't deallocate memory for the matrix itself,
 * only for its structures.
 */
void sp_matrix_skyline_free(sp_matrix_skyline_ptr self);

/*
 * Creates the sparse matrix in Yale format
 * Acts as a copy-constructor
 * If matrix was in CCS format, produces CCS, CRS otherwise
 */
void sp_matrix_yale_init(sp_matrix_yale_ptr self,
                         sp_matrix_ptr mtx);

                          
/*
 * Copy sparse matrix from mtx_from to mtx_to
 * This function assumes what mtx_to is already cleared by sp_matrix_yale_free
 * or mtx_to is a pointer to uninitialized sp_matrix structure
 */
void sp_matrix_yale_copy(sp_matrix_yale_ptr mtx_from,
                         sp_matrix_yale_ptr mtx_to);


/*
 * Destructor for a sparse matrix in Yale format
 * This function doesn't deallocate memory for the matrix itself,
 * only for its structures.
 */
void sp_matrix_yale_free(sp_matrix_yale_ptr self);


/* getters/setters for a sparse matrix */

/* returns a pointer to the specific element
 * zero pointer if not found */
double* sp_matrix_element_ptr(sp_matrix_ptr self,int i, int j);
/* adds an element value to the matrix node (i,j) and return (i,j) */
double sp_matrix_element_add(sp_matrix_ptr self,
                             int i, int j, double value);
/*
 * Cancellation(making all zeros) of the i-th row and column
 * keeping the diagonal element the same
 * Used in FEA procedures to apply prescribed BCs
 * Returns the diagonal (i,i) value
 */
double sp_matrix_cross_cancellation(sp_matrix_ptr self, int i);

/* shortcut for adding of the matrix elements */
#define MTX(m,i,j,v) sp_matrix_element_add((m),(i),(j),(v));


/* rearrange columns of a matrix to prepare for solving SLAE */
void sp_matrix_reorder(sp_matrix_ptr self);

/* determines the matrix properties */
matrix_properties sp_matrix_properites(sp_matrix_ptr self);

/* returns the number of nonzeroes in matrix */
int sp_matrix_nonzeros(sp_matrix_ptr self);

/*
 * Matrix-vector multiplication for matrix in Yale format
 * y = A*x
 */
void sp_matrix_yale_mv(sp_matrix_yale_ptr self,double* x, double* y);

/*
 * Transposes the matrix in Yale format
 */
void sp_matrix_yale_transpose(sp_matrix_yale_ptr self,
                              sp_matrix_yale_ptr to);

/*
 * Convert sparse matrix to given format type
 * if the type is the same - do nothing
 * matrix 'to' shall be uninitialized
 */
int sp_matrix_yale_convert(sp_matrix_yale_ptr from,
                            sp_matrix_yale_ptr to,
                            sparse_storage_type type);

/*
 * Convert sparse matrix to given format type inplace
 */
int sp_matrix_yale_convert_inplace(sp_matrix_yale_ptr self,
                                   sparse_storage_type type);

/*
 * Calculates the permuted matrix C = P*A*Q
 * by given vector of inverse row permutation pinv:
 * row i becomes row k if pinv[k] = i
 * and vector of column permutation q
 * returns 0 in case of error, nonzero otherwise
 */
int sp_matrix_yale_permute(sp_matrix_yale_ptr self,
                           sp_matrix_yale_ptr permuted,
                           int* pinv,
                           int* q);

/* determines the matrix properties */
matrix_properties sp_matrix_yale_properites(sp_matrix_yale_ptr self);

/* compare 2 matricies in the same format */
matrix_comparison sp_matrix_yale_cmp(sp_matrix_yale_ptr mtx1,
                                     sp_matrix_yale_ptr mtx2);

/* Print contens of the matrix in index form to the stdout */
void sp_matrix_printf(sp_matrix_ptr self);
/* Print contents of the matrix in dense form to the stdout */
void sp_matrix_printf2(sp_matrix_ptr self);
void sp_matrix_dump(sp_matrix_ptr self, const char* filename);
/* Print contents of the matrix in array form to the stdout */
void sp_matrix_yale_printf(sp_matrix_yale_ptr self);
/* Print the matrix stats to the stdout */
void sp_matrix_yale_printf2(sp_matrix_yale_ptr self);
void sp_matrix_skyline_printf(sp_matrix_skyline_ptr self);

#endif /* __SP_MATRIX_H__ */
