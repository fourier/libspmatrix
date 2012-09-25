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

#ifndef _SP_ITER_H_
#define _SP_ITER_H_

#include "sp_matrix.h"

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
 * Conjugate Gradient solver
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
 * Creates ILU decomposition of the sparse matrix 
 */
void sp_matrix_create_ilu(sp_matrix_ptr self,sp_matrix_skyline_ilu_ptr ilu);


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
/*
 * Transpose-Free Quasi-Minimal Residual solver
 * self - matrix in Yale format
 * b - right-part vector
 * x0 - first approximation of the solution
 * max_iter - pointer to maximum number of iterations, MAX_ITER if zero;
 * will contain a number of iterations passed
 * tolerance - pointer to desired tolerance value, TOLERANCE if zero;
 * will contain norm of the residual at the end of iteration
 * x - output vector
 */
void sp_matrix_yale_solve_tfqmr(sp_matrix_yale_ptr self,
                                double* b,
                                double* x0,
                                int* max_iter,
                                double* tolerance,
                                double* x);

/*
 * Conjugate Gradient Squared solver
 * self - matrix in Yale format
 * b - right-part vector
 * x0 - first approximation of the solution
 * max_iter - pointer to maximum number of iterations, MAX_ITER if zero;
 * will contain a number of iterations passed
 * tolerance - pointer to desired tolerance value, TOLERANCE if zero;
 * will contain norm of the residual at the end of iteration
 * x - output vector
 */
void sp_matrix_yale_solve_cgs(sp_matrix_yale_ptr self,
                              double* b,
                              double* x0,
                              int* max_iter,
                              double* tolerance,
                              double* x);


#endif /* _SP_ITER_H_ */
