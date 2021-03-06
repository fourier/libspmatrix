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
#ifndef _SP_DIRECT_H_
#define _SP_DIRECT_H_

#include "sp_matrix.h"

/*
 * Symbolic infromation of the sparse matrix
 * used by Cholesky decomposition
 */
typedef struct
{
  int* etree;                   /* elimination tree */
  int* post;                    /* postorder of the elimination tree */
  int* rowcounts;               /* row counts for Cholesky matrix L */
  int* colcounts;               /* column counts for Cholesky matrix L */
  int nonzeros;                 /* number of nonzeros in Cholesky matrix L */
  int* crs_indicies;            /* column indicies for CRS format of L */
  int* crs_offsets;             /* row offsets - beginning of row for CRS */
  int* ccs_indicies;            /* row indicies for CCS format of L */
  int* ccs_offsets;             /* column offsets - beginning of column */
} sp_chol_symbolic;
typedef sp_chol_symbolic* sp_chol_symbolic_ptr;

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
int sp_matrix_yale_ereach(sp_matrix_yale_ptr self,
                          int* etree,
                          int k,
                          int* out);

/*
 * Calculates row and column counts for the Cholesky decomposition
 * row counts - number of nonzero elements in rows
 * col counts - number of nonzero elements in columns
 * by given sparse matrix self and elimination tree etree
 * returns 0 in case of error, nonzero otherwise
 */
int sp_matrix_yale_chol_counts(sp_matrix_yale_ptr self,
                               int* etree,
                               int* rowcounts,
                               int* colcounts);

/*
 * Solves SLAE L*x = b
 * by given L sparse matrix with nonzero diagonal
 * Returns nonzero if successfull
 */
int sp_matrix_yale_lower_solve(sp_matrix_yale_ptr self,
                               double* b,
                               double* x);

/*
 * Solves SLAE L^T*x = b
 * by given L sparse matrix with nonzero diagonal
 * Returns nonzero if successfull
 */
int sp_matrix_yale_lower_trans_solve(sp_matrix_yale_ptr self,
                                     double* b,
                                     double* x);


/*
 * Performs the symbolic analysis used in Cholesky decomposition
 * Returns nonzero if succesfull
 */
int sp_matrix_yale_chol_symbolic(sp_matrix_yale_ptr self,
                                 sp_chol_symbolic_ptr symb);
/*
 * Deallocates Cholesky symbolic analysis structure values
 * This function doesn't deallocate memory for the struture itself,
 * only for its structures
 */
void sp_matrix_yale_symbolic_free(sp_chol_symbolic_ptr symb);


/*
 * Finds the numeric Cholesky decomposition of the given matrix.
 * Symbolic Cholesky decomposition shall already be found.
 * Returns nonzero if succesfull 
 */
int sp_matrix_yale_chol_numeric(sp_matrix_yale_ptr self,
                                sp_chol_symbolic_ptr symb,
                                sp_matrix_yale_ptr L);
/*
 * Solves the SLAE self*x=b using the Cholesky decomposition.
 * self shall be symmetric positive-definite
 * b and x sizes are the same as number of rows in self
 * Returns nonzero if successfull
 */
int sp_matrix_yale_chol_solve(sp_matrix_yale_ptr self,
                              double* b,
                              double* x);
/*
 * Solves the SLAE self*x=b by the Cholesky decomposition
 * using the preliminary calculated Cholesky symbolic decomposition
 * data.
 * self shall be symmetric positive-definite
 * b and x sizes are the same as number of rows in self
 * Returns nonzero if successfull
 */
int sp_matrix_yale_chol_symbolic_solve(sp_matrix_yale_ptr self,
                                       sp_chol_symbolic_ptr symb,
                                       double* b,
                                       double* x);


/*
 * Solves the SLAE LL'*x=b with given matrix L from the 
 * preliminary calculated Cholesky numeric decomposition
 * b and x sizes are the same as number of rows in self
 * Returns nonzero if successfull
 */
int sp_matrix_yale_chol_numeric_solve(sp_matrix_yale_ptr L,
                                      double* b,
                                      double* x);



#endif /* _SP_DIRECT_H_ */
