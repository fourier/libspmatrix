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

/* #include <stdlib.h> */
#include <memory.h>
#include <math.h>

#include "sp_iter.h"
#include "sp_mem.h"

/*
 * Scalar product x*y
 */ 
inline static double prod(double* x, double* y, int size)
{
  double r = 0;
  int i = 0;
  for ( ; i < size; ++ i)
    r += x[i]*y[i];
  return r;
}

/*
 * Norm 2 of the vector x: norm2(x) = sqrt(x*x)
 */ 
inline static double norm2(double* x, int size)
{
  double r = 0;
  int i = 0;
  for ( ; i < size; ++ i)
    r += x[i]*x[i];
  return sqrt(r);
}


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
  r = (double*)spalloc(size);
  p = (double*)spalloc(size);
  temp = (double*)spalloc(size);
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
    a1 = prod(r,r,msize);        /* (r_j,r_j) */
    a2 = prod(temp,p,msize);     /* (A*p_j,p_j) */

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
    residn = norm2(r,msize);
    if (residn < tol )
      break;

    /* compute (r_{j+1},r_{j+1}) */
    a2 = prod(r,r,msize);

    /* b_j = (r_{j+1},r_{j+1})/(r_j,r_j) */
    beta = a2/a1;
    
    /* p_{j+1} = r_{j+1} + beta_j*p_j */
    for (i = 0; i < msize; ++ i)
      p[i] = r[i] + beta*p[i];
  }
  *max_iter = j;
  *tolerance = residn;
  
  spfree(r);
  spfree(p);
  spfree(temp);
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
  r = (double*)spalloc(size);
  r1 = (double*)spalloc(size);
  p = (double*)spalloc(size);
  z = (double*)spalloc(size);
  temp = (double*)spalloc(size);
  
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
    a1 = prod(r,z,msize);       /* (r_j,z_j) */
    a2 = prod(p,temp,msize);    /* (A*p_j,p_j) */
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
    residn = norm2(r,msize);
    if (residn < tol )
      break;

    /* z_{j+1} = M^{-1}*r_{j+1} */
    memcpy(r1,r,size);
    memset(temp,0,size);
    sp_matrix_skyline_ilu_lower_solve(ILU,r1,temp); /* temp = L^{-1}*r */
    sp_matrix_skyline_ilu_upper_solve(ILU,temp,z); /* z = U^{-1}*temp */

    
    /* compute (r_{j+1},z_{j+1}) */
    a2 = prod(r,z,msize);

    /* b_j = (r_{j+1},z_{j+1})/(r_j,z_j) */
    beta = a2/a1;
    
    /* d_{j+1} = r_{j+1} + beta_j*d_j */
    for (i = 0; i < msize; ++ i)
      p[i] = z[i] + beta*p[i];
  }
  *max_iter = j;
  *tolerance = residn;
  
  /* free vectors */
  spfree(r);
  spfree(r1);
  spfree(z);
  spfree(p);
  spfree(temp);
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
  self->ilu_diag = (double*)spalloc(sizeof(double)*parent->rows_count);
  self->ilu_lowertr = (double*)spalloc(sizeof(double)*parent->tr_nonzeros);
  self->ilu_uppertr = (double*)spalloc(sizeof(double)*parent->tr_nonzeros);
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
  spfree(self->ilu_diag);
  spfree(self->ilu_lowertr);
  spfree(self->ilu_uppertr);
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


void sp_matrix_yale_solve_tfqmr(sp_matrix_yale_ptr self,
                                double* b,
                                double* x0,
                                int* max_iter,
                                double* tolerance,
                                double* x)
{
  /* Transpose-Free Quasi-Minimal Residual Algorithm */
  /*
   * Based on the book:
   * Saad Y. Iterative methods for sparse linear systems (2ed., 2003)
   * page 235
   */

  /* variables */
  int i,m;
  double alpha, beta;
  double residn = 0;
  double theta = 0, eta = 0, rho = 0;
  double tau;
  int size = sizeof(double)*self->rows_count;
  int msize = self->rows_count;
  int max_iterations = *max_iter;
  double tol = *tolerance;
  double c;
  double* r;              /* residual */
  double* r1;              /* r^*_0 */
  double* temp;
  double* d;
  double* v;
  double* w;
  double* u;
  
  /* allocate memory for vectors */
  r = (double*)spalloc(size);
  r1 = (double*)spalloc(size);
  temp = (double*)spalloc(size);
  d = (double*)spcalloc(msize,sizeof(double));
  v = (double*)spcalloc(msize,sizeof(double));
  w = (double*)spcalloc(msize,sizeof(double));
  u = (double*)spcalloc(msize,sizeof(double));
  /* d = 0 */
  memset(d,0,size);

  /* x = x_0 */
  memcpy(x,x0,size);

  /* r_0 = b - A*x_0 */
  sp_matrix_yale_mv(self,b,r);
  for ( i = 0; i < msize; ++ i)
    r[i] = b[i] - r[i];

  /* w_0 = r_0 */
  memcpy(w,r,size);
  /* u_0 = r_0 */
  memcpy(u,r,size);
  /* v_0 = A*u_0 */
  sp_matrix_yale_mv(self,u,v);

  /* choose r^*_0 that rho_0 = (r^*_0,r_0) != 0 */
  memcpy(r1,r,size);
  r1[0] = 1;

  rho = prod(r1,r,msize);

  tau = norm2(r,msize);
  
  /* CG loop */
  for ( m = 0; m < max_iterations; m ++ )
  {
    /* if m is even */
    if ((m&1) == 0)
    {
      /* alpha_{m+1} = alpha_m */
      /* u_{m+1} = u_m - alpha_m*v_m */
      for (i = 0; i < msize; ++ i)
        u[i] -= alpha*v[i];
    }
    else
    {
      /* alpha_m = rho_m/(v_m,r^*_0) */
      alpha = rho/prod(v,r1,msize);
    }

    /* temp = A*u_m */
    sp_matrix_yale_mv(self,u,temp);

    /* w_{m+1} = w_m - alpha_m*A*u_m */
    for (i = 0; i < msize; ++ i)
      w[i] -= alpha*temp[i];

    /* d_{m+1} = u_m + (theta^2_m/alpha_m)*eta_m*d_m */
    c = theta*theta*eta/alpha;
    for (i = 0; i < msize; ++ i)
      d[i] = u[i] + c*d[i];

    /* theta_{m+1} = norm2(w_{m+1})/tau_m */
    theta = norm2(w,msize)/tau;

    /* c_{m+1} = (1+theta^2_{m+1})^(-1/2) */
    c = 1./sqrt(1+theta*theta);

    /* tau_{m+1} = tau_m*theta_{m+1}*c_{m+1} */
    tau *= theta*c;

    /* eta_{m+1} = c^2_{m+1}*alpha_m */
    eta = c*c*alpha;
    
    /* x_{m+1} = x_m+eta_{m+1}*d_{m+1} */
    
    for (i = 0; i < msize; ++ i)
      x[i] += eta*d[i];

    /* check for convergence */
    /*
     * residual estimation:
     * norm(b-A*x_m) <= sqrt(m+1)*tau_m
     */
    if (sqrt(m+1)*tau < tol )
      break;

    /* if m is odd */
    if ( m&1 )
    {
      /* store c = rho_{m-1} */
      c = rho;
      /* rho_{m+1} = w_{m+1}*r^*_0 */
      rho = prod(w,r1,msize);
      /* beta_{m-1} = rho_{m+1}/rho_{m-1} */
      beta = rho/c;
      /* u_{m+1} = w_{m+1} + beta_{m-1}*u_m */
      for (i = 0; i < msize; ++ i)
        u[i] = w[i] + beta*u[i];
      /* v_{m+1} = A*u_{m+1} + beta_{m-1}*(A*u_m+beta_{m-1}*v_{m-1} */
      for ( i = 0; i < msize; ++ i)
        v[i] = beta*(temp[i]+beta*v[i]);
      /* temp = A*u_{m+1} */
      sp_matrix_yale_mv(self,u,temp);
      for (i = 0; i < msize; ++ i)
        v[i] += temp[i];
    }

  }
  *max_iter = m;
  *tolerance = residn;
  
  spfree(r);
  spfree(temp);
  spfree(d);
  spfree(v);
  spfree(w);
  spfree(u);
}


void sp_matrix_yale_solve_cgs(sp_matrix_yale_ptr self,
                              double* b,
                              double* x0,
                              int* max_iter,
                              double* tolerance,
                              double* x)
{
  /* Conjugate Gradient Squared Algorithm */
  /*
   * Based on the book:
   * Saad Y. Iterative methods for sparse linear systems (2ed., 2003)
   * page 229
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
  double* r1;             /* r^*_0 */
  double* p;              /* search direction */
  double* q;
  double* u;
  double* temp;

  /* allocate memory for vectors */
  r = (double*)spcalloc(msize,sizeof(double));
  r1 = (double*)spcalloc(msize,sizeof(double));
  p = (double*)spcalloc(msize,sizeof(double));
  temp = (double*)spcalloc(msize,sizeof(double));
  q = (double*)spcalloc(msize,sizeof(double));
  u = (double*)spcalloc(msize,sizeof(double));


  /* x = x_0 */
  memcpy(x,x0,size);

  /* r_0 = b - A*x_0 */
  sp_matrix_yale_mv(self,b,r);
  for ( i = 0; i < msize; ++ i)
    r[i] = b[i] - r[i];
  /* r1 - arbitrary */
  memcpy(r1,r,size);
  r1[0] = 1;
  
  /* p_0 = r_0 */
  memcpy(p,r,size);

  /* u_0 = r_0 */
  memcpy(u,r,size);
  
  /* CGS loop */
  for ( j = 0; j < max_iterations; j ++ )
  {
    /* temp = A*p_j */
    sp_matrix_yale_mv(self,p,temp);
    /* compute (r_j,r^*_0) and (A*p_j,r^*_0) */
    a1 = prod(r,r1,msize);      /* (r_j,r^*_0) */
    a2 = prod(temp,r1,msize);   /* (A*p_j,r^*_0) */
    /*                  
     *            (r_j,r^*_0) 
     * alpha_j = -----------
     *           (A*p_j,r^*_0)
     */                     
    alpha = a1/a2;              

    /*
     * q_j = u_j - alpha_j*A*p_j
     */
    for (i = 0; i < msize; ++ i)
      q[i] = u[i] - alpha*temp[i];
           
    /* x_{j+1} = x_j+alpha_j(u_j+q_j) */
    for (i = 0; i < msize; ++ i)
      x[i] += alpha*(u[i] + q[i]);

    /* temp = A*(u_j+q_j) */
    sp_matrix_yale_mvsum(self,u,q,temp);
    
    /* r_{j+1} = r_j-alpha_j*A*(u_j+q_j) */
    for (i = 0; i < msize; ++ i)
      r[i] -= alpha*temp[i]; 

    /* check for convergence */
    residn = norm2(r,msize);
    if (residn < tol )
      break;

    /* compute (r_{j+1},r^*_0) */
    a2 = prod(r,r1,msize);

    /* b_j = (r_{j+1},r^*_0)/(r_j,r^*_0) */
    beta = a2/a1;

    /* u_{j+1} = r_{j+1} + b_j*q_j */
    for (i = 0; i < msize; ++ i)
      u[i] = r[i] + beta*q[i];
    
    /* p_{j+1} = u_{j+1} + beta_j*(q_j + b_j*p_j) */
    for (i = 0; i < msize; ++ i)
      p[i] = u[i] + beta*(q[i] + beta*p[i]);
  }
  *max_iter = j;
  *tolerance = residn;
  
  spfree(r);
  spfree(p);
  spfree(temp);
  spfree(r1);
  spfree(q);
  spfree(u);
}
