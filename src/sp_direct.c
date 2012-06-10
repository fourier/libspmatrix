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

#include <assert.h>
#include <stdlib.h>
#include <float.h>
#include <math.h>
#include <stdio.h>

#include "sp_direct.h"
#include "sp_tree.h"
#include "sp_utils.h"
#include "sp_log.h"

/* trivial comparison with zero */
static int is_almost_zero(double x)
{
  return (fabs(x) < 2*FLT_EPSILON);
}


int sp_matrix_yale_etree(sp_matrix_yale_ptr self, int* tree)
{
  int k,i,j,p;
  int *parents;
  if ( !self || self->storage_type != CCS)
    return 0;

  parents = malloc(sizeof(int)*self->rows_count);

  /* initialize tree */
  for ( k = 0; k < self->rows_count; ++ k)
    tree[k] = -1;
  
  for ( k = 0; k < self->rows_count; ++ k) /* loop by Tk */
  {
    /* loop by nonzero rows in the column */
    j = 0;
    for ( p = self->offsets[k]; p < self->offsets[k+1]; ++p )
    {
      /* store the index to i for simplicity */
      i = self->indicies[p];
      if ( i < k )              /* if above the diagonal */
      {
        /* find the root for the i in the previous elimination tree */
        parents[j] = tree_find(tree,self->rows_count,i);
        /* increase the number of roots found (number of nonzeros)
         * above diagonal*/
        j++;
      }
    }
    /* now update - append found roots to the new root k */
    for (i = 0; i < j; ++ i)
      tree[parents[i]] = k;
  }
  free(parents);
  
  return 1;
}


int sp_matrix_yale_ereach(sp_matrix_yale_ptr self, int* etree, int k, int* out)
{
  int i,j,p;
  int count = 0;
  if ( self->storage_type != CCS)
    return -1;
  
  /* clear the out */
  for (j = 0; j < self->rows_count; ++ j)
    out[j] = -1;
#define OUT_MARK(i) {out[(i)] = (i); ++count; }
  /* mark k */
  OUT_MARK(k);
  /* loop by nonzero rows in the column */
  for ( p = self->offsets[k]; p < self->offsets[k+1]; ++p )
  {
    /* store the index to i for simplicity */
    i = self->indicies[p];
    if ( i < k )              /* if above the diagonal */
    {
      /* mark i */
      if (out[i] != i)
        OUT_MARK(i);
      /* find the root */
      for (j = etree[i]; j != -1 && j != etree[j]; j = etree[j])
      {
        if ( out[j] != j)
        {
          /* mark j */
          OUT_MARK(j);
        }
        else                    /* marked reached */
          break;
      }
    }
  }
#undef OUT_MARK
  /* compress the output */
#define SWAP_VALUES(x,y) { i = (x); x = (y); (y) = i; }
  k = 0;
  for ( j = 0; j < self->rows_count; ++ j)
    if ( out[j] != -1 )
    {
      SWAP_VALUES(out[j],out[k]);
      k++;
      if  (k >= count)
        break;
    }
#undef SWAP_VALUES
  return count;
}


int sp_matrix_yale_chol_counts(sp_matrix_yale_ptr self,
                               int* etree,
                               int* rowcounts,
                               int* colcounts)
{
  int result = 1;
  int n = self->rows_count;
  int i,j,k,p;
  /* array to store marked nodes. nonzero means node is marked */
  char *marked = (char*)malloc(n);
  memset(rowcounts,0,n*sizeof(int));
  memset(colcounts,0,n*sizeof(int));
#define _TREE_MARK(i) {marked[(i)] = 1; ++rowcounts[k]; }
  for (k = 0; k < n; ++ k)      /* loop by rows */
  {
    memset(marked,0,n);          /* clear all marks */
    i = 0;                       /* initialize number of row subtrees
                                  * containing k to zero */
    /* loop by nonzero rows in the column, or nonzero cols in
     * the row, depending on storage schema */
    for ( p = self->offsets[k]; p < self->offsets[k+1]; ++p )
    {
      /* store the index to i for simplicity */
      j = self->indicies[p];    /* a_kj != 0 */
      if ( j <= k )             /* if above the diagonal */
      {
        while(!marked[j]) /* marked encountered, goto next nonzero */
        {
          ++colcounts[j];
          _TREE_MARK(j);
          if (j == k)           /* root encountered, goto next nonzero */
            break;
          /* else climb up by the tree */
          j = etree[j];
        }
      }
      else                      /* assuming indicies are sorted */
        break;
    }
  }
#undef _TREE_MARK
  free(marked);
  return result;
}


int sp_matrix_yale_lower_solve(sp_matrix_yale_ptr self,
                               double* b,
                               double* x)
{
  int i,j,p;
  double value;
  int n = self->rows_count;
  assert(self);
  assert(b);
  assert(x);

  if (self->storage_type == CCS)
  {
    memcpy(x,b,n*sizeof(double));
    for ( j = 0; j < n; ++ j)
    {
      value = self->values[self->offsets[j]];
      if (is_almost_zero(value))
      {
        LOGERROR("Lower solver: diagonal element: %e",value);
        return 0;
      }
      x[j] /= value;
      for (p = self->offsets[j]+1; p < self->offsets[j+1]; ++ p)
      {
        i = self->indicies[p];
        x[i] -= self->values[p] * x[j];
      }
    }
  }
  else                          /* CRS */
  {
    memset(x,0,n*sizeof(double));
    for ( i = 0; i < n; ++ i)
    {
      /* x_i = (b_i - \sum\limits_{j=1}^{i-1} l_ij x_j)/l_ii*/
      x[i] = b[i];
      for (p = self->offsets[i];
           p < self->offsets[i+1] && (j = self->indicies[p]) < i; ++ p)
        x[i] -= self->values[p] * x[j];
      value = self->values[self->offsets[i+1]-1];
      if (is_almost_zero(value))
      {
        LOGERROR("Lower solver: diagonal element: %e",value);
        return 0;
      }
      x[i] /= value;
    }
  }
  return 1;
}

int sp_matrix_yale_lower_trans_solve(sp_matrix_yale_ptr self,
                                     double* b,
                                     double* x)
{
  int i,j,p;
  double value;
  int n = self->rows_count;  
  assert(self);
  assert(b);
  assert(x);

  if (self->storage_type == CCS)
  {
    memset(x,0,n*sizeof(double));
    for ( i = n-1; i >= 0; -- i)
    {
      /* x_i = (b_i - \sum\limits_{j=i+1}^{n} l_ji x_j)/l_ii */
      x[i] = b[i];
      for (p = self->offsets[i]+1; p < self->offsets[i+1]; ++ p)
      {
        j = self->indicies[p];
        x[i] -= self->values[p] * x[j];
      }
      value = self->values[self->offsets[i]];
      if (is_almost_zero(value))
      {
        LOGERROR("Lower solver: diagonal element: %e",value);
        return 0;
      }
      x[i] /= value;
    }
  }
  else                          /* CRS */
  {
    memcpy(x,b,n*sizeof(double));
    for ( j = n-1; j >= 0; -- j)
    {
      value = self->values[self->offsets[j+1]-1];
      if (is_almost_zero(value))
      {
        LOGERROR("Lower solver: diagonal element: %e",value);
        return 0;
      }
      x[j] /= value;
      for (p = self->offsets[j]; p < self->offsets[j+1]-1; ++ p)
      {
        i = self->indicies[p];
        x[i] -= self->values[p] * x[j];
      }
    }
  }
  return 1;
}


/*
 * Finds the symbolic Cholesky decomposition - portrait of the matrix L
 * for row/column storage type WITHOUT numeric values
 * using the given matrix and results of the premilinary symbolic analysis
 * 
 * Returns nonzero if succesfull 
 */
static
int sp_matrix_yale_chol_structure(sp_matrix_yale_ptr self,
                                  sp_chol_symbolic_ptr symb)
{
  int result = 1;
  int count,i,j,p,k;
  int *offsets, *indicies;
  offsets = calloc(self->rows_count + 1,sizeof(int));
  symb->crs_indicies = calloc(symb->nonzeros,sizeof(int));
  symb->ccs_indicies = calloc(symb->nonzeros,sizeof(int));
  indicies = calloc(symb->nonzeros,sizeof(int));  
  /* calculate offsets for CRS */
  j = 0;
  for (i = 0; i < self->rows_count; ++ i)
  {
    offsets[i] = j;
    j += symb->rowcounts[i];
  }
  offsets[i] = symb->nonzeros;
  symb->crs_offsets = memdup(offsets,(self->rows_count + 1)*sizeof(int));
  /* calculate offsets for CCS */
  j = 0;
  for (i = 0; i < self->rows_count; ++ i)
  {
    offsets[i] = j;
    j += symb->colcounts[i];
  }
  offsets[i] = symb->nonzeros;
  symb->ccs_offsets = memdup(offsets,(self->rows_count + 1)*sizeof(int));
  for (i = 0; i < self->rows_count; ++ i)
  {
    /* ereach simply defines the portrait of every i-th row */
    count = sp_matrix_yale_ereach(self,symb->etree,i,indicies);
    if (count != symb->crs_offsets[i+1]-symb->crs_offsets[i])
    {
      LOGERROR("sp_matrix_yale_chol_structure: count = %d, expected: %d."
               " Possibly not symmetric matrix",
               count,symb->crs_offsets[i+1]-symb->crs_offsets[i]);
      free(offsets);
      free(indicies);
      return 0;
    }
    /* so, a_ij != 0 where j in indicies array */
    for ( p = 0; p < count; ++ p)
    {
      j = indicies[p];
      /* a_ij != 0  */
      /* algoritm is the same as in sp_matrix_yale_transpose */
      k = offsets[j]++;
      symb->ccs_indicies[k] = i;
    }
    memcpy(symb->crs_indicies+symb->crs_offsets[i],indicies,count*sizeof(int));
  }
  free(offsets);
  free(indicies);
  return result;
}


int sp_matrix_yale_chol_symbolic(sp_matrix_yale_ptr self,
                                 sp_chol_symbolic_ptr symb)
{
#define _SYMB_VERIFY(x) if (!(x)){sp_matrix_yale_symbolic_free(symb);break;}
  int result = 0;
  int i;
  if (self && symb)
  {
    do
    {
      memset(symb,0,sizeof(sp_chol_symbolic));
      symb->etree = malloc(self->rows_count*sizeof(int));
      result = sp_matrix_yale_etree(self,symb->etree);
      _SYMB_VERIFY(result);
      symb->rowcounts = malloc(self->rows_count*sizeof(int));
      symb->colcounts = malloc(self->rows_count*sizeof(int));
      result = sp_matrix_yale_chol_counts(self,symb->etree,
                                          symb->rowcounts,
                                          symb->colcounts);
      _SYMB_VERIFY(result);
      symb->post = malloc(self->rows_count*sizeof(int));
      tree_postorder_perm(symb->etree,self->rows_count,symb->post);
      /* calculate nonzeros */
      for (i = 0; i < self->rows_count; ++ i)
      {
        symb->nonzeros += symb->rowcounts[i];
      }
      _SYMB_VERIFY((result = sp_matrix_yale_chol_structure(self,symb)));
    } while(0);
  }
#undef _SYMB_VERIFY
  return result;
}

void sp_matrix_yale_symbolic_free(sp_chol_symbolic_ptr symb)
{
  if (symb)
  {
    free(symb->etree);
    free(symb->post);
    free(symb->rowcounts);
    free(symb->colcounts);
    if (symb->crs_indicies)
      free(symb->crs_indicies);
    if (symb->crs_offsets)
      free(symb->crs_offsets);
    if (symb->ccs_indicies)
      free(symb->ccs_indicies);
    if (symb->ccs_offsets)
      free(symb->ccs_offsets);
    symb->nonzeros = 0;
    symb->etree = 0;
    symb->post = 0;
    symb->rowcounts = 0;
    symb->colcounts = 0;
    symb->crs_offsets = 0;
    symb->crs_indicies = 0;
    symb->ccs_offsets = 0;
    symb->ccs_indicies = 0;
  }
}

/*
 * Sparse Triangular solver for CCS matrix
 * n - up to n-th row.
 * x - right part, scattered. Will be overwritten on output
 * indicies - indexes of the solution
 * size - number of nonzeros in solution
 * dot - x*x
 * returns nonzero if successful
 */
static int sp_matrix_yale_sparse_lower_solve(sp_matrix_yale_ptr self,
                                             int n,
                                             double* x,
                                             int* indicies,
                                             int size,
                                             double* dot)
{
  int i,j,p,q;
  double d = 0,value = 0;
  if (self->storage_type != CCS)
    return d;
 
  /* for j in X do */
  for (p = 0; p < size; ++ p)
  {
    j = indicies[p];
    /* xj =xj/ljj */
    /* diagonal j-th element - is the first column element in L */
    value = self->values[self->offsets[j]];
    if (is_almost_zero(value))
    {
      LOGERROR("Sparse lower solver: %d diagonal element: %e",j,value);
      return 0;
    }
    x[j] = x[j]/value;
    /* for i>j  where lij != 0 do */
    for (q = self->offsets[j]+1; /* under diagonal */
         q < self->offsets[j+1] && (i = self->indicies[q]) < n;
         ++ q)
    {
      /* xi = xi - lij*xj */
      x[i] = x[i] - self->values[q]*x[j];
    } /* end for */
    
  } /* end for */
  /* calculate scalar product */
  for (p = 0; p < size; ++ p)
  {
    j = indicies[p];
    d += x[j]*x[j];
  }
  *dot = d;
  return 1;
}


int sp_matrix_yale_chol_numeric(sp_matrix_yale_ptr self,
                                sp_chol_symbolic_ptr symb,
                                sp_matrix_yale_ptr L)
{
  int result = 1;
  int i,j,p;
  int k = 0;
  int* offsets;
  double* x;
  int* rowoffsets;
  double value;
  double v,A_kk;
  if (!self || !symb || !L || self->storage_type != CCS)
    return 0;
  /* initialize L */
  L->storage_type = CCS;
  L->rows_count = self->rows_count;
  L->cols_count = self->cols_count;
  L->nonzeros = symb->nonzeros;
  L->offsets =  memdup(symb->ccs_offsets,(self->rows_count+1)*sizeof(int));
  L->indicies = memdup(symb->ccs_indicies,symb->nonzeros*sizeof(int));
  L->values = calloc(symb->nonzeros,sizeof(double));
  /* store offsets */
  offsets = memdup(symb->ccs_offsets,(self->rows_count+1)*sizeof(int));
  /* right-part vector */
  x = calloc(self->rows_count,sizeof(double));
#define _SP_CHOL_STOP {sp_matrix_yale_free(L);free(offsets);free(x);\
    LOGERROR("Cholesky decomposition: error in %d row",k);return 0;}
  /* up-looking Cholesky */
  /* loop by rows, constructing one k-th row at a time */
  for (k = 0; k < self->rows_count; ++ k)
  {
    memset(x,0,self->rows_count*sizeof(double));
    v = 0;
    /*
     * construct the scattered vector x, containing nonzeros
     * a_jk, where j < k: A(1:k-1,k)
     */
    for (p = self->offsets[k];
         p < self->offsets[k+1] && self->indicies[p] < k;
         ++p)
      x[self->indicies[p]] = self->values[p];
    A_kk = self->values[p];
    /* solve the SLAE L(1:k-1,1:k-1)*L(k,1:k-1)=A(1:k-1,k) */
    rowoffsets = symb->crs_indicies+symb->crs_offsets[k];
    if (!sp_matrix_yale_sparse_lower_solve(L,k,x,
                                           rowoffsets,
                                           symb->rowcounts[k]-1,
                                           &v))
      _SP_CHOL_STOP;
    /* store result to kth row */
    for (p = symb->crs_offsets[k];
         p < symb->crs_offsets[k+1] && (j = symb->crs_indicies[p]) < k;
         ++p)
    {
      i = offsets[j];
      L->values[i] = x[j];
      offsets[j]++;
    }
    value = A_kk-v;
    if (value < 0)
      _SP_CHOL_STOP;
    L->values[L->offsets[k]] = sqrt(value);
    offsets[k]++;
  }
  free(x);
  free(offsets);
#undef SP_CHOL_STOP
  return result;
}


int sp_matrix_yale_chol_solve(sp_matrix_yale_ptr self,
                              double* b,
                              double* x)
{
  int result = 0;
  sp_chol_symbolic symb;
  result = sp_matrix_yale_chol_symbolic(self,&symb);
  if (result)
  {
    result = sp_matrix_yale_chol_symbolic_solve(self,&symb,b,x);
    sp_matrix_yale_symbolic_free(&symb);
  }
  return result;
}


int sp_matrix_yale_chol_symbolic_solve(sp_matrix_yale_ptr self,
                                       sp_chol_symbolic_ptr symb,
                                       double* b,
                                       double* x)
{
  int result = 0;
  sp_matrix_yale L;
  result = sp_matrix_yale_chol_numeric(self,symb,&L);
  if (result)
  {
    result = sp_matrix_yale_chol_numeric_solve(&L,b,x);
    sp_matrix_yale_free(&L);
  }
  return result;
}


int sp_matrix_yale_chol_numeric_solve(sp_matrix_yale_ptr L,
                                      double* b,
                                      double* x)
{
  int result = 0;
  double *y = calloc(L->rows_count,sizeof(double));
  memset(x,0,L->rows_count*sizeof(double));
  /* SLAE: Ax=b
   * A=LL'
   * LL'x=b
   * L'x=y
   * Ly=b
   */
  /* Ly=b */
  result = sp_matrix_yale_lower_solve(L,b,y);
  if (result)
  {
    /* L'x = y  */
    result = sp_matrix_yale_lower_trans_solve(L,y,x);
  }
  free(y);
  return result;
}
