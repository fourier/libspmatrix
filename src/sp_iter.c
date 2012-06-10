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

#include <stdlib.h>
#include <math.h>
#include <memory.h>

#include "sp_iter.h"

void sp_matrix_yale_solve_cg(sp_matrix_yale_ptr self,
                             double* b,
                             double* x0,
                             int* max_iter,
                             double* tolerance,
                             double* x)
{
  /* Conjugate Gradient Algorithm */
  /*
   * Based on the book:
   * Saad Y. Iterative methods for sparse linear systems (2ed., 2000)
   * page 178
   */
   
  /* variables */
  int i,j;
  double alpha, beta,a1,a2;
  double residn = 0;
  int size = sizeof(double)*self->rows_count;
  int msize = self->rows_count;
  int max_iterations = *max_iter;
  double tol = *tolerance;
  double* r;              /* residual */
  double* p;              /* search direction */
  double* temp;

  /* allocate memory for vectors */
  r = (double*)malloc(size);
  p = (double*)malloc(size);
  temp = (double*)malloc(size);
  /* clear vectors */
  memset(r,0,size);
  memset(p,0,size);
  memset(temp,0,size);

  /* x = x_0 */
  memcpy(x,x0,size);

  /* r_0 = b - A*x_0 */
  sp_matrix_yale_mv(self,b,r);
  for ( i = 0; i < msize; ++ i)
    r[i] = b[i] - r[i];

  /* p_0 = r_0 */
  memcpy(p,r,size);
  
  /* CG loop */
  for ( j = 0; j < max_iterations; j ++ )
  {
    /* temp = A*p_j */
    sp_matrix_yale_mv(self,p,temp);
    /* compute (r_j,r_j) and (A*p_j,p_j) */
    a1 = 0; a2 = 0;
    for (i = 0; i < msize; ++ i)
    {
      a1 += r[i]*r[i];         /* (r_j,r_j) */
      a2 += temp[i]*p[i];      /* (A*p_j,p_j) */
    }

    /*            (r_j,r_j) 
     * alpha_j = -----------
     *           (A*p_j,p_j)
     */                     
    alpha = a1/a2;              
                                
    /* x_{j+1} = x_j+alpha_j*p_j */
    for (i = 0; i < msize; ++ i)
      x[i] += alpha*p[i];
    
    /* r_{j+1} = r_j-alpha_j*A*p_j */
    for (i = 0; i < msize; ++ i)
      r[i] -= alpha*temp[i]; 

    /* check for convergence */
    residn = fabs(r[0]);
    for (i = 1; i < msize; ++ i )
      if (fabs(r[i]) > residn) residn = fabs(r[i]);
    if (residn < tol )
      break;

    /* compute (r_{j+1},r_{j+1}) */
    a2 = 0;
    for (i = 0; i < msize; ++ i)
      a2 += r[i]*r[i];

    /* b_j = (r_{j+1},r_{j+1})/(r_j,r_j) */
    beta = a2/a1;
    
    /* p_{j+1} = r_{j+1} + beta_j*p_j */
    for (i = 0; i < msize; ++ i)
      p[i] = r[i] + beta*p[i];
  }
  *max_iter = j;
  *tolerance = residn;
  
  free(r);
  free(p);
  free(temp);
}


void sp_matrix_yale_solve_pcg_ilu(sp_matrix_yale_ptr self,
                                  sp_matrix_skyline_ilu_ptr ILU,                         
                                  double* b,
                                  double* x0,
                                  int* max_iter,
                                  double* tolerance,
                                  double* x)
{
  /* Preconditioned Conjugate Gradient Algorithm */
  /*
   * Based on the book:
   * Saad Y. Iterative methods for sparse linear systems (2ed., 2000)
   * page 246
   *
   * Preconditioner: Incomplete LU decomposition (ILU)
   * M = L*U, A = M-R
   */

  /* variables */
  int i,j;
  double alpha, beta,a1,a2;
  double residn = 0;
  int size = sizeof(double)*self->rows_count;
  int msize = self->rows_count;
  int max_iterations = *max_iter;
  double tol = *tolerance;
  
  double* r;              /* residual */
  double* r1;             /* backup of the residual */
  double* p;              /* search direction */
  double* z;              /* z = M^{-1}*r */
  double* temp;

  /* allocate memory for vectors */
  r = (double*)malloc(size);
  r1 = (double*)malloc(size);
  p = (double*)malloc(size);
  z = (double*)malloc(size);
  temp = (double*)malloc(size);
  
  /* clear vectors */
  memset(r,0,size);
  memset(r1,0,size);
  memset(p,0,size);
  memset(z,0,size);
  memset(temp,0,size);

  /* x = x_0 */
  memcpy(x,x0,size);

  /* r_0 = b - A*x_0 */
  sp_matrix_yale_mv(self,b,r);
  for ( i = 0; i < msize; ++ i)
    r[i] = b[i] - r[i];
  
  /* backup residual */
  memcpy(r1,r,size);
  /* z_0 = M^{-1}*r_0 */
  /*
   * to solve system L*U*x = b
   * y = U*x, => L*y = b
   * U*x = y => x
   */ 
  sp_matrix_skyline_ilu_lower_solve(ILU,r1,temp); /* temp = L^{-1}*r */
  /* r1 now changed, temp contains solution */
  sp_matrix_skyline_ilu_upper_solve(ILU,temp,z); /* z = U^{-1}*temp */
  /* temp now changed, z contains solution*/
  
  /* p_0 = z_0 */
  memcpy(p,z,size);
  
  /* CG loop */
  for ( j = 0; j < max_iterations; j ++ )
  {
    /* temp = A*p_j */
    memset(temp,0,size);
    sp_matrix_yale_mv(self,p,temp);
    /* compute (r_j,z_j) and (A*p_j,p_j) */
    a1 = 0; a2 = 0;
    for (i = 0; i < msize; ++ i)
    {
      a1 += r[i]*z[i]; /* (r_j,z_j) */
      a2 += p[i]*temp[i];      /* (A*p_j,p_j) */
    }

    /*            (r_j,z_j) 
     * alpha_j = -----------
     *           (A*p_j,p_j)
     */                     
    alpha = a1/a2;              
                                
    /* x_{j+1} = x_j+alpha_j*p_j */
    for (i = 0; i < msize; ++ i)
      x[i] += alpha*p[i];
    
    /* r_{j+1} = r_j-alpha_j*A*p_j */
    for (i = 0; i < msize; ++ i)
      r[i] -= alpha*temp[i]; 

    /* check for convergence */
    residn = fabs(r[0]);
    for (i = 1; i < msize; ++ i )
      if (fabs(r[i]) > residn) residn = fabs(r[i]);
    if (residn < tol )
      break;

    /* z_{j+1} = M^{-1}*r_{j+1} */
    memcpy(r1,r,size);
    memset(temp,0,size);
    sp_matrix_skyline_ilu_lower_solve(ILU,r1,temp); /* temp = L^{-1}*r */
    sp_matrix_skyline_ilu_upper_solve(ILU,temp,z); /* z = U^{-1}*temp */

    
    /* compute (r_{j+1},z_{j+1}) */
    a2 = 0;
    for (i = 0; i < msize; ++ i)
      a2 += r[i]*z[i];

    /* b_j = (r_{j+1},z_{j+1})/(r_j,z_j) */
    beta = a2/a1;
    
    /* d_{j+1} = r_{j+1} + beta_j*d_j */
    for (i = 0; i < msize; ++ i)
      p[i] = z[i] + beta*p[i];
  }
  *max_iter = j;
  *tolerance = residn;
  
  /* free vectors */
  free(r);
  free(r1);
  free(z);
  free(p);
  free(temp);
}

void sp_matrix_create_ilu(sp_matrix_ptr self,sp_matrix_skyline_ilu_ptr ilu)
{
  sp_matrix_skyline A;
  /* reorder self if not already reordered */
  if (!self->ordered)
    sp_matrix_reorder(self);
  /* initialize skyline matrix for ILU decomposition */
  sp_matrix_skyline_init(&A,self);
  /*
   * create ILU decomposition of the sparse matrix in skyline format
   * taking ownership of the skyline A matrix
   */
  sp_matrix_skyline_ilu_copy_init(ilu,&A);
  /*
   * since init_copy_sp_matrix_skyline_ilu takes the ownership
   * of the A matrix it is not needed to free A matrix
   */
}


void sp_matrix_skyline_ilu_copy_init(sp_matrix_skyline_ilu_ptr self,
                                     sp_matrix_skyline_ptr parent)
{
  int i,j,k,l,q;
  double sum;
  
  /* copy parent member-wise */
  self->parent = *parent;
  /* allocate memory for ILU decomposition arrays */
  self->ilu_diag = (double*)malloc(sizeof(double)*parent->rows_count);
  self->ilu_lowertr = (double*)malloc(sizeof(double)*parent->tr_nonzeros);
  self->ilu_uppertr = (double*)malloc(sizeof(double)*parent->tr_nonzeros);
  /* clear arrays before construction of the ILU decomposition */
  memset(self->ilu_diag,0,sizeof(double)*parent->rows_count);
  memset(self->ilu_lowertr,0,sizeof(double)*parent->tr_nonzeros);
  memset(self->ilu_uppertr,0,sizeof(double)*parent->tr_nonzeros);

  for (k = 0; k < parent->rows_count; ++ k)
  {
    for ( j = parent->iptr[k]; j < parent->iptr[k+1]; ++ j)
    {
      /*
       * L_{kj} = (A_{kj} - \sum\limits_{i=1}^{j-1}L_{ki}U_{ij}/U_{jj}
       * calculate using L_{k,jptr[j]}
       */
      sum = 0;
      q = parent->jptr[j];        /* column index */
      for ( i = parent->iptr[k]; i < parent->iptr[k+1]; ++ i)
      {
        for ( l = parent->iptr[q]; l < parent->iptr[q+1]; ++ l)
        {
          /* if row and column indicies are the same */
          if ( parent->jptr[i] == parent->jptr[l] )
            sum += self->ilu_lowertr[i]*self->ilu_uppertr[l];
        }
      }
      self->ilu_lowertr[j] =
        (parent->lower_triangle[j] - sum)/self->ilu_diag[q];
    }

    /*
     * U_{kk} = A_{kk} -
     * \sum\limits_{i=1}^{k-1} L_{ki}U_{ik}
     */
    sum = 0;
    for ( i = parent->iptr[k]; i < parent->iptr[k+1]; ++ i)
      sum += self->ilu_lowertr[i]*self->ilu_uppertr[i];
    self->ilu_diag[k] = parent->diag[k] - sum;

    for (j = k; j < parent->rows_count; ++ j)
    {
      for ( q = parent->iptr[j]; q < parent->iptr[j+1]; ++ q)
        if (k == parent->jptr[q])
        {
          /*
           * U_{kj} = A_{kj} -
           * \sum\limits_{i=1}^{k-1}L_{ki}U_{ij}
           */
          sum = 0;
          /*
           * i = iptr[k]:iptr[k+1]-1 are coordinates of the
           * k-th row in lower matrix array (self->ilu_lowertr)
           * l = iptr[j]:iptr[j+1]-1 are coordinates of the
           * j-th column in upper matrix array (self->ilu_uppertr)
           */
          for ( i = parent->iptr[k]; i < parent->iptr[k+1]; ++ i)
            for ( l = parent->iptr[j]; l < parent->iptr[j+1]; ++ l)
            {
              /* if row and column indicies are the same */
              if ( parent->jptr[i] == parent->jptr[l] )
                sum += self->ilu_lowertr[i]*self->ilu_uppertr[l];
            }
          self->ilu_uppertr[q] = parent->upper_triangle[q] - sum;
        }
    }
  }
}

void sp_matrix_skyline_ilu_free(sp_matrix_skyline_ilu_ptr self)
{
  free(self->ilu_diag);
  free(self->ilu_lowertr);
  free(self->ilu_uppertr);
  sp_matrix_skyline_free(&self->parent);
}


void sp_matrix_skyline_ilu_lower_mv(sp_matrix_skyline_ilu_ptr self,
                                    double* x,
                                    double* y)
{
  int i,j;
  memset(y,0,sizeof(double)*self->parent.rows_count);
  for (i = 0; i < self->parent.rows_count; ++ i)
  {
    y[i] = x[i];
    for (j = self->parent.iptr[i]; j < self->parent.iptr[i+1]; ++ j)
      y[i] += x[self->parent.jptr[j]]*self->ilu_lowertr[j];
  }
}

void sp_matrix_skyline_ilu_upper_mv(sp_matrix_skyline_ilu_ptr self,
                                    double* x,
                                    double* y)
{
  int i,j;
  memset(y,0,sizeof(double)*self->parent.rows_count);

  for (i = 0; i < self->parent.rows_count; ++ i)
    y[i] = x[i]*self->ilu_diag[i];
  for (i = 0; i < self->parent.rows_count; ++ i)
  {
    for ( j = self->parent.iptr[i]; j < self->parent.iptr[i+1]; ++j )
      y[self->parent.jptr[j]] += x[i]*self->ilu_uppertr[j];
  }
}

void sp_matrix_skyline_ilu_lower_solve(sp_matrix_skyline_ilu_ptr self,
                                       double* b,
                                       double* x)
{
  int i,j;
  memset(x,0,sizeof(double)*self->parent.rows_count);
  
  for ( i = 0; i < self->parent.rows_count; ++ i)
  {
    for (j = self->parent.iptr[i]; j < self->parent.iptr[i+1]; ++ j)
      b[i] -= x[self->parent.jptr[j]]*self->ilu_lowertr[j];
    x[i] = b[i];
  }
}

void sp_matrix_skyline_ilu_upper_solve(sp_matrix_skyline_ilu_ptr self,
                                       double* b,
                                       double* x)
{
  int i,j;
  memset(x,0,sizeof(double)*self->parent.rows_count);

  for ( i = self->parent.rows_count-1; i >= 0; -- i)
  {
    x[i] = b[i]/self->ilu_diag[i];
    for (j = self->parent.iptr[i]; j < self->parent.iptr[i+1]; ++ j)
      b[self->parent.jptr[j]] -= x[i]*self->ilu_uppertr[j];
  }
}

