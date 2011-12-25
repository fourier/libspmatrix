/* -*- Mode: C; tab-width: 2; indent-tabs-mode: nil; c-basic-offset: 2 -*- */
/*
 Copyright (C) 2011 Alexey Veretennikov (alexey dot veretennikov at gmail.com)
 
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

/*
 * Sparse matrix row/column storage array
 */
typedef struct
{
  int width;                    /* size of an array */
  int last_index;               /* last stored index, i.e. if width = 20
                                 * it will be 9 if only 10 nonzero elements
                                 * stored */
  int  *indexes;                /* array of column/row indexes */
  double *values;               /* array of values */
} indexed_array;
typedef indexed_array* indexed_array_ptr;

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
 * Sparse matrix CSLR(Skyline) format
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
 * ILU decomposition of the sparse matrix in Skyline (CSLR) format
 * ILU decomposition keeps the symmetric portrait of the sparse matrix
 */
typedef struct
{
  sp_matrix_skyline parent;
  double *ilu_diag;              /* U matrix diagonal */
  double *ilu_lowertr;           /* nonzero elements of the lower(L) matrix */
  double *ilu_uppertr;           /* nonzero elements of the upper(U) matrix */
} sp_matrix_skyline_ilu;
typedef sp_matrix_skyline_ilu* sp_matrix_skyline_ilu_ptr;

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

/* indexed_arrays operations */
/* Swap i and j elements in the indexed array.
 * Used in indexed_array_sort*/
void indexed_array_swap(indexed_array_ptr self,int i, int j);
/* Performs in-place sort of the indexed array */
void indexed_array_sort(indexed_array_ptr self, int l, int r);
/* Print contents of the indexed array to the stdout  */
void indexed_array_printf(indexed_array_ptr self);

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
 * This function assumes what mtx_to is already cleared by free_sp_matrix
 * or mtx_to is a pointer to uninitialized sp_matrix structure
 */
void sp_matrix_copy(sp_matrix_ptr mtx_from,
                    sp_matrix_ptr mtx_to);

/*
 * Converts matrix storage format CRS <=> CCS
 * mtx_to shall be uninitialized sp_matrix structure
 */
void sp_matrix_convert(sp_matrix_ptr mtx_from,
                       sp_matrix_ptr mtx_to,
                       sparse_storage_type type);


/*
 * Creates ILU decomposition of the sparse matrix 
 */
void sp_matrix_create_ilu(sp_matrix_ptr self,sp_matrix_skyline_ilu_ptr ilu);

/*
 * Construct CSLR sparse matrix based on sp_matrix format
 * mtx - is the (reordered) sparse matrix to take data from
 * Acts as a copy-constructor
 */
void sp_matrix_skyline_init(sp_matrix_skyline_ptr self,
                            sp_matrix_ptr mtx);
/*
 * Destructor for a sparse matrix in CSLR format
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

/* shortcut for adding of the matrix elements */
#define MTX(m,i,j,v) sp_matrix_element_add((m),(i),(j),(v));


/* rearrange columns of a matrix to prepare for solving SLAE */
void sp_matrix_reorder(sp_matrix_ptr self);

/* returns nonzero value if the matrix is symmetric */
int sp_matrix_issymmetric(sp_matrix_ptr self);

/* returns nonzero value if the matrix protrait is symmetric */
int sp_matrix_issymmetric_portrait(sp_matrix_ptr self);

/* returns nonzero value if the matrix is skew-symmetric */
int sp_matrix_isskew_symmetric(sp_matrix_ptr self);

/* returns the number of nonzeroes in matrix */
int sp_matrix_nonzeros(sp_matrix_ptr self);

/* Matrix-vector multiplication
 * y = A*x
 */
void sp_matrix_mv(sp_matrix_ptr self,double* x, double* y);

/*
 * Matrix-vector multiplication for matrix in Yale format
 * y = A*x
 */
void sp_matrix_yale_mv(sp_matrix_yale_ptr self,double* x, double* y);

/*
 * Constructs the elimination tree from the matrix in Yale format
 * the matrix shall be in CCS format
 * etree is the pointer to the array with rows_count elements
 * to store the elimination tree
 * returns the nonzero value if all ok
 * 0 in case of error
 */
int sp_matrix_yale_etree(sp_matrix_yale_ptr self, int* etree);

/*
 * Constructs the nonzero portrait of the kth row of the L matrix in
 * LL^T Cholesky decomposition
 * out - output array with allocated size self->rows_count
 * returns the number of first meaningfull elements in out array
 * Example: if the resulting array (out) is
 * 2 7 9 10 -1 -1 -1 -1 -1 -1 -1
 * then the return value is 4
 */
int sp_matrix_yale_ereach(sp_matrix_yale_ptr self, int* etree, int k, int* out);

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

/*
 * Solves SLAE L*x = b
 * by given L sparse matrix 
 * n - is the size of the x vector, and therefore
 * the matrix L will be used up to nth row & column.
 */
void sp_matrix_lower_solve(sp_matrix_ptr self,
                           int n,
                           double* b,
                           double* x);


/*
 * Conjugate Grade solver
 * self - matrix in Yale format
 * b - right-part vector
 * x0 - first approximation of the solution
 * max_iter - pointer to maximum number of iterations, MAX_ITER if zero;
 * will contain a number of iterations passed
 * tolerance - pointer to desired tolerance value, TOLERANCE if zero;
 * will contain norm of the residual at the end of iteration
 * x - output vector
 */
void sp_matrix_yale_solve_cg(sp_matrix_yale_ptr self,
                             double* b,
                             double* x0,
                             int* max_iter,
                             double* tolerance,
                             double* x);

/*
 * Preconditioned Conjugate Grade solver
 * Preconditioner in form of the ILU decomposition
 * self - matrix in Yale format
 * b - right-part vector
 * x0 - first approximation of the solution
 * max_iter - pointer to maximum number of iterations, MAX_ITER if zero;
 * will contain a number of iterations passed
 * tolerance - pointer to desired tolerance value, TOLERANCE if zero;
 * will contain norm of the residual at the end of iteration
 * x - output vector
 */
void sp_matrix_yale_solve_pcg_ilu(sp_matrix_yale_ptr self,
                                  sp_matrix_skyline_ilu_ptr ilu,
                                  double* b,
                                  double* x0,
                                  int* max_iter,
                                  double* tolerance,
                                  double* x);


/*
 * Create ILU decomposition of the sparse matrix in skyline format
 * lu_diag - ILU decomposition diagonal
 * lu_lowertr - lower triangle of the ILU decomposition
 * lu_uppertr - upper triangle of the ILU decomposition
 */
void sp_matrix_skyline_ilu_copy_init(sp_matrix_skyline_ilu_ptr self,
                                     sp_matrix_skyline_ptr parent);

/* Free the sparse matrix skyline & ilu decomposition structure */
void sp_matrix_skyline_ilu_free(sp_matrix_skyline_ilu_ptr self);

/*
 * by given L,U - ILU decomposition of the matrix A
 * calculates L*x = y
 */
void sp_matrix_skyline_ilu_lower_mv(sp_matrix_skyline_ilu_ptr self,
                                    double* x,
                                    double* y);
/*
 * by given L,U - ILU decomposition of the matrix A
 * calculates U*x = y
 */
void sp_matrix_skyline_ilu_upper_mv(sp_matrix_skyline_ilu_ptr self,
                                    double* x,
                                    double* y);

/*
 * by given L,U - ILU decomposition of the matrix A
 * Solves SLAE L*x = b
 * Warning! Side-Effect: modifies b
 */
void sp_matrix_skyline_ilu_lower_solve(sp_matrix_skyline_ilu_ptr self,
                                       double* b,
                                       double* x);

/*
 * by given L,U - ILU decomposition of the matrix A
 * Solves SLAE U*x = b
 * Warning! Side-Effect: modifies b 
 */
void sp_matrix_skyline_ilu_upper_solve(sp_matrix_skyline_ilu_ptr self,
                                       double* b,
                                       double* x);

/* Print contens of the matrix in index form to the stdout */
void sp_matrix_printf(sp_matrix_ptr self);
void sp_matrix_dump(sp_matrix_ptr self, const char* filename);
void sp_matrix_skyline_dump(sp_matrix_skyline_ptr self, const char* filename);

#endif /* __SP_MATRIX_H__ */
