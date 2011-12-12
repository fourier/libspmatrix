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

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <math.h>
#include <assert.h>

#include "sp_matrix.h"
#include "sp_utils.h"

#include "logger.h"

#define TRUE 1
#define FALSE 0

void sp_matrix_init(sp_matrix_ptr mtx,
                    int rows,
                    int cols,
                    int bandwidth,
                    sparse_storage_type type)
{
  int i,n;
  if (mtx)
  {
    mtx->rows_count = rows;
    mtx->cols_count = cols;
    mtx->ordered = NOT_ORDERED;
    mtx->storage_type = type;
    n = type == CRS ? rows : cols;
    mtx->storage = (indexed_array*)malloc(sizeof(indexed_array)*n);
    /* create rows or cols with fixed bandwidth */
    for (i = 0; i < n; ++ i)
    {
      mtx->storage[i].width = bandwidth;
      mtx->storage[i].last_index = -1;
      mtx->storage[i].indexes = (int*)malloc(sizeof(int)*bandwidth);
      mtx->storage[i].values = (double*)malloc(sizeof(double)*bandwidth);
      memset(mtx->storage[i].indexes,0,sizeof(int)*bandwidth);
      memset(mtx->storage[i].values,0,sizeof(double)*bandwidth);
    }
  }
}


void sp_matrix_free(sp_matrix_ptr mtx)
{
  int i,n;
  if (mtx)
  {
    n = mtx->storage_type  == CRS ? mtx->rows_count : mtx->cols_count;
    for (i = 0; i < n; ++ i)
    {
      free(mtx->storage[i].indexes);
      free(mtx->storage[i].values);
    }
    free(mtx->storage);
    mtx->storage = (indexed_array*)0;
    mtx->cols_count = 0;
    mtx->rows_count = 0;
  }
}

void sp_matrix_clear(sp_matrix_ptr mtx)
{
  int i,n;
  if (mtx)
  {
    n = mtx->storage_type  == CRS ? mtx->rows_count : mtx->cols_count;
    for (i = 0; i < n; ++ i)
    {
      memset(mtx->storage[i].values,0,sizeof(double)*(mtx->storage[i].width));
    }
  }
}

void sp_matrix_copy(sp_matrix_ptr mtx_from, sp_matrix_ptr mtx_to)
{
  int i,n;
  assert(mtx_from && mtx_to);
  n = mtx_from->storage_type  == CRS ? mtx_from->rows_count :
    mtx_from->cols_count;
  mtx_to->rows_count = mtx_from->rows_count;
  mtx_to->cols_count = mtx_from->cols_count;
  mtx_to->ordered = mtx_from->ordered;
  mtx_to->storage_type = mtx_from->storage_type;
  mtx_to->storage =
    (indexed_array*)malloc(sizeof(indexed_array)*n);
  /* copy rows */
  for (i = 0; i < n; ++ i)
  {
    memset(&mtx_to->storage[i],0,sizeof(indexed_array));
    mtx_to->storage[i].width = mtx_from->storage[i].width;
    mtx_to->storage[i].last_index = mtx_from->storage[i].last_index;
    mtx_to->storage[i].indexes =
      (int*)malloc(sizeof(int)*mtx_from->storage[i].width);
    memset(mtx_to->storage[i].indexes,0,sizeof(int)*mtx_from->storage[i].width);
    mtx_to->storage[i].values =
      (double*)malloc(sizeof(double)*mtx_from->storage[i].width);
    memset(mtx_to->storage[i].values,0,sizeof(double)*mtx_from->storage[i].width);
    memcpy(mtx_to->storage[i].indexes, mtx_from->storage[i].indexes,
           sizeof(int)*mtx_from->storage[i].width);
    memcpy(mtx_to->storage[i].values, mtx_from->storage[i].values,
           sizeof(double)*mtx_from->storage[i].width);
  }
}

void sp_matrix_convert(sp_matrix_ptr mtx_from,
                       sp_matrix_ptr mtx_to,
                       sparse_storage_type type)
{
  int i,j;
  if (type == mtx_from->storage_type)
    return;
  sp_matrix_init(mtx_to,
                 mtx_from->rows_count,
                 mtx_from->cols_count,
                 mtx_from->storage[0].width,
                 type);
  if (type == CCS)              /* CRS -> CCS */
  {
    for (i = 0; i < mtx_from->rows_count; ++ i)
    {
      for (j = 0; j <= mtx_from->storage[i].last_index; ++ j)
        MTX(mtx_to,i,mtx_from->storage[i].indexes[j],
            mtx_from->storage[i].values[j]);
    }
  }
  else                          /* CCS -> CRS */
  {
    for (i = 0; i < mtx_from->cols_count; ++ i)
    {
      for (j = 0; j <= mtx_from->storage[i].last_index; ++ j)
        MTX(mtx_to,mtx_from->storage[i].indexes[j],i,
            mtx_from->storage[i].values[j]);
    }
  }
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

void sp_matrix_skyline_init(sp_matrix_skyline_ptr self,sp_matrix_ptr mtx)
{
  /*
   * Construct CSLR matrix from the sp_matrix
   * with symmetric portrait
   */
  int i,j,k,iptr,l_count,u_count,column;
  double* pvalue = 0;

  /* assert what mtx is already reordered */
  assert(mtx->ordered);
  /* currenty implemented conversion only from CRS format */
  assert(mtx->storage_type == CRS);
  /* matrix shall be with symmetric portrait */
  assert(sp_matrix_issymmetric_portrait(mtx));
  
  self->rows_count = mtx->rows_count;
  self->cols_count = mtx->cols_count;
  /*
   * get an information about number of nonzero elements
   */
  self->nonzeros = 0;
  for (i = 0; i < mtx->rows_count; ++ i)
    self->nonzeros += mtx->storage[i].last_index + 1;
  
  /* calculate number of upper-triangle elements */
  l_count = 0;
  u_count = 0;
  for (i = 0; i < mtx->rows_count; ++ i)
    for (j = 0; j <= mtx->storage[i].last_index; ++ j)
      if ( mtx->storage[i].indexes[j] > i)
        u_count ++;
      else if (mtx->storage[i].indexes[j] < i)
        l_count ++;
  /*
   * check if the number of upper triangle nonzero elements
   * is the same as number of lower triangle nonzero elements
   */
  assert(l_count == u_count);
  self->tr_nonzeros = l_count;
  
  /* allocate memory for arrays */
  self->diag = (double*)malloc(sizeof(double)*mtx->rows_count);
  self->lower_triangle = l_count ? (double*)malloc(sizeof(double)*l_count) : 0;
  self->upper_triangle = u_count ? (double*)malloc(sizeof(double)*u_count) : 0;
  self->jptr = l_count ? (int*)malloc(sizeof(int)*l_count)  : 0;
  self->iptr = (int*)malloc(sizeof(int)*(mtx->rows_count+1));

  /* fill diagonal */
  for (i = 0; i < mtx->rows_count; ++ i)
  {
    pvalue = sp_matrix_element_ptr(mtx,i,i);
    self->diag[i] = pvalue ? *pvalue : 0;
  }
  /* now fill arrays with proper values */
  u_count = 0,l_count = 0;
  for (i = 0; i < mtx->rows_count; ++ i)
  {
    iptr = -1;
    self->iptr[i] = 0;
    for (j = 0; j <= mtx->storage[i].last_index; ++ j)
    {
      if ( mtx->storage[i].indexes[j] < i)
      {
        /*
         * set a flag what we found the first nonzero element in
         * current row in lower triangle
         */
        if (iptr == -1)
          iptr  = l_count;
        /* fill lower triangle values */
        column = mtx->storage[i].indexes[j];
        self->jptr[l_count] = column;
        self->lower_triangle[l_count] = mtx->storage[i].values[j];
        /* fill upper triangle values - column-wise */
        for ( k = 0; k <= mtx->storage[column].last_index; ++ k)
          if (mtx->storage[column].indexes[k] == i)
          {
            self->upper_triangle[l_count] =
              mtx->storage[column].values[k];
            break;
          }
        l_count ++;
      }
    }
    self->iptr[i] = iptr == -1 ? l_count : iptr;
  }
  /* finalize iptr array */
  self->iptr[i] = self->tr_nonzeros;
}

void sp_matrix_skyline_free(sp_matrix_skyline_ptr self)
{
  if (self)
  {
    self->rows_count = 0;
    self->cols_count = 0;
    self->nonzeros = 0;
    self->tr_nonzeros = 0;
    free(self->diag);
    free(self->lower_triangle);
    free(self->upper_triangle);
    free(self->jptr);
    free(self->iptr);
  }
}

void sp_matrix_yale_init(sp_matrix_yale_ptr self,
                         sp_matrix_ptr mtx)
{
  int nonzeros = sp_matrix_nonzeros(mtx);
  int i,j,index;
  int n = mtx->storage_type == CRS ? mtx->rows_count : mtx->cols_count;
  /* reorder self if not already reordered */
  if (!mtx->ordered)
    sp_matrix_reorder(mtx);
  /* initialize matrix */
  memset(self,sizeof(sp_matrix_yale),0);
  self->storage_type = mtx->storage_type;
  self->rows_count = mtx->rows_count;
  self->cols_count = mtx->cols_count;
  self->nonzeros   = nonzeros;
  /* allocate memory for arrays */
  self->offsets  = calloc(n+1,      sizeof(int));
  self->indicies = calloc(nonzeros, sizeof(int));
  self->values   = calloc(nonzeros, sizeof(double));
  /* convert */
  /* loop by nonzero columns in row i */
  j = 0;
  for ( i = 0; i < n; ++ i)
  {
    self->offsets[i] = j;
    for (index = 0; index <= mtx->storage[i].last_index; ++ index)
    {
      self->indicies[j] = mtx->storage[i].indexes[index];
      self->values[j]   = mtx->storage[i].values[index];
      j++;
    }
  }
  self->offsets[i] = nonzeros;
}

void sp_matrix_yale_free(sp_matrix_yale_ptr self)
{
  free(self->offsets);
  free(self->indicies);
  free(self->values);
  memset(self,sizeof(sp_matrix_yale),0);
}


double* sp_matrix_element_ptr(sp_matrix_ptr self,int i, int j)
{
  int index;
  /* check for matrix and if i,j are proper indicies */
  assert(self && 
         (i >= 0 && i < self->rows_count ) &&
         (j >= 0 && j < self->cols_count ));
  {
    if (self->storage_type == CRS)
    {
      /* loop by nonzero columns in row i */
      for (index = 0; index <= self->storage[i].last_index; ++ index)
        if (self->storage[i].indexes[index] == j)
          return &self->storage[i].values[index];
    }
    else                        /* CCS */
    {
      /* loop by nonzero rows in column i */
      for (index = 0; index <= self->storage[j].last_index; ++ index)
        if (self->storage[j].indexes[index] == i)
          return &self->storage[j].values[index];
    }
  }
  return (double*)0;
}

double sp_matrix_element_add(sp_matrix_ptr self,int i, int j, double value)
{
  int index,new_width,I,J;
  int* indexes = (int*)0;
  double* values = (double*)0;
  /* check for matrix and if i,j are proper indicies */
  assert (self && 
          (i >= 0 && i < self->rows_count ) &&
          (j >= 0 && j < self->cols_count ));
  /* set I and J to be i and j in case of CRS or j and i otherwise */
  I = self->storage_type == CRS ? i : j;
  J = self->storage_type == CRS ? j : i;
  /* loop by nonzero columns in row/col i */
  for (index = 0; index <= self->storage[I].last_index; ++ index)
    if (self->storage[I].indexes[index] == J)
    {
      /* nonzerod element found, add to it */
      self->storage[I].values[index] += value;
      return self->storage[I].values[index];
    }
  /* needed to add a new element to the row/col */
    
  /*
   * check if bandwidth is not exceed and reallocate memory
   * if necessary
   */
  if (self->storage[I].last_index == self->storage[I].width - 1)
  {
    new_width = self->storage[I].width*2;
    if (new_width <= 0)             /* avoid crashes on bad bandwidth */
      new_width = 1;
    indexes = (int*)realloc(self->storage[I].indexes,new_width*sizeof(int));
    assert(indexes);
    self->storage[I].indexes = indexes;
    values = (double*)realloc(self->storage[I].values,new_width*sizeof(double));
    assert(values);
    self->storage[I].values = values;
    self->storage[I].width = new_width;
    self->ordered = NOT_ORDERED;
  }
  /* add an element to the row/col */
  self->storage[I].last_index++;
  self->storage[I].values[self->storage[I].last_index] = value;
  self->storage[I].indexes[self->storage[I].last_index] = J;
  return value;
}

/* Swap 2 elements of the indexed array */
void indexed_array_swap(indexed_array_ptr self,int i, int j)
{
  int tmp_idx;
  double tmp_val;
  tmp_idx = self->indexes[i];
  self->indexes[i] = self->indexes[j];
  self->indexes[j] = tmp_idx;
  tmp_val = self->values[i];
  self->values[i] = self->values[j];
  self->values[j] = tmp_val;
}

void indexed_array_sort(indexed_array_ptr self, int l, int r)
{
  /*
   * Quick sort procedure for indexed(compressed) arrays
   * for example rows for CRS sparse matrix or columns for CSC
   * sparse matrix
   */
  int pivot,i;
  int tmp_idx;

  /* boundary checks */
  if (l < r)
  {
    if ( r - l == 1)
    {
      if (self->indexes[l] > self->indexes[r])
        indexed_array_swap(self,r,l);
      return;
    }
    /* choose the pivoting element */
    pivot = (int)((r+l)/2.);
    /* in-place partition procedure - move all elements
     * lower than pivoting to the left, greater to the right */
    tmp_idx  = self->indexes[pivot];
    indexed_array_swap(self,pivot,r);
    pivot = l;
    for ( i = l; i < r; ++ i)
    {
      if (self->indexes[i] <= tmp_idx )
      {
        indexed_array_swap(self,i,pivot);
        pivot++;
      }
    }
    indexed_array_swap(self,r,pivot);
    /* repeat procedure for the left and right parts of an array */
    indexed_array_sort(self,l,pivot-1);
    indexed_array_sort(self,pivot+1,r);
  }
}

void indexed_array_printf(indexed_array_ptr self)
{
  int i;
  if (self)
  {
    printf("indexes = [");
    for (i = 0; i <= self->last_index; ++ i)
      printf("%d,\t",self->indexes[i]);
    printf("%d]\n",self->indexes[i]);
    printf("values  = [");
    for (i = 0; i <= self->last_index; ++ i)
      printf("%f,\t",self->values[i]);
    printf("%f]\n",self->values[i]);
  }
}


void sp_matrix_reorder(sp_matrix_ptr self)
{
  int i,j,n;
  int size,stored;
  n = self->storage_type == CRS ? self->rows_count : self->cols_count;
  for (i = 0; i < n; ++ i)
    indexed_array_sort(&self->storage[i],0,self->storage[i].last_index);
  
  for (i = 0; i < n; ++ i)
    for (j = 0; j < self->storage[i].last_index; ++ j)
      /* assert(self->storage[i].indexes[j] < self->storage[i].indexes[j+1]); */
      if (self->storage[i].indexes[j] >= self->storage[i].indexes[j+1])
      {
        printf("row %d\n",i);
        indexed_array_printf(&self->storage[i]);
        assert(0);
      }
  self->ordered = ORDERED;
  /* calculate fill factor of the matrix */
  stored = 0;
  size = self->rows_count*self->cols_count;
  for ( i = 0; i < n; ++ i)
    stored += self->storage[i].last_index + 1;
  /* printf("Sparse matrix compressed:\n"); */
  /* printf("- size: %dx%d\n",self->rows_count,self->cols_count); */
  /* printf("- nonzeros: %d\n",stored); */
  /* printf("- fill factor: %.2f %%\n",stored/(size/100.0)); */
  /* printf("- avergare nonzeros per row: %d\n", */
  /*        /\* (int)rint(stored/(double)self->rows_count)); *\/ */
  /*        (int)(stored/(double)self->rows_count)); */
  LOGINFO("Sparse matrix compressed:\n");
  LOGINFO("- size: %dx%d\n",self->rows_count,self->cols_count);
  LOGINFO("- nonzeros: %d\n",stored);
  LOGINFO("- fill factor: %.2f %%\n",stored/(size/100.0));
  LOGINFO("- avergare nonzeros per row: %d\n",
          /* (int)rint(stored/(double)self->rows_count)); */
          (int)(stored/(double)self->rows_count));

}


int sp_matrix_issymmetric(sp_matrix_ptr self)
{
  double *value;
  int n = self->rows_count;
  int i,j;
  if  (self->rows_count != self->cols_count )
    return FALSE;
  if (self->storage_type == CRS)
  {
    for ( i = 0; i < n; ++ i)
      for (j = 0; j < self->storage[i].last_index; ++ j)
      {
        value = sp_matrix_element_ptr(self,self->storage[i].indexes[j],i);
        if ( !value )
          return FALSE;
        if (!EQL(*value,self->storage[i].values[j]))
          return FALSE;
      }
  }
  else                          /* CCS */
  {
    for ( i = 0; i < n; ++ i)
      for (j = 0; j < self->storage[i].last_index; ++ j)
      {
        value = sp_matrix_element_ptr(self,i,self->storage[i].indexes[j]);
        if ( !value )
          return FALSE;
        if (! EQL(*value,self->storage[i].values[j]))
          return FALSE;
      }
    
  }
  return TRUE;
}


int sp_matrix_issymmetric_portrait(sp_matrix_ptr self)
{
  double *value;
  int n = self->rows_count;
  int i,j;
  if  (self->rows_count != self->cols_count )
    return FALSE;
  if (self->storage_type == CRS)
  {
    for ( i = 0; i < n; ++ i)
      for (j = 0; j < self->storage[i].last_index; ++ j)
      {
        value = sp_matrix_element_ptr(self,self->storage[i].indexes[j],i);
        if ( !value )
          return FALSE;
      }
  }
  else                          /* CCS */
  {
    for ( i = 0; i < n; ++ i)
      for (j = 0; j < self->storage[i].last_index; ++ j)
      {
        value = sp_matrix_element_ptr(self,i,self->storage[i].indexes[j]);
        if ( !value )
          return FALSE;
      }
  }
  return TRUE;
}

int sp_matrix_isskew_symmetric(sp_matrix_ptr self)
{
  double *value;
  int n = self->rows_count;
  int i,j;
  if  (self->rows_count != self->cols_count )
    return FALSE;
  if (self->storage_type == CRS)
  {
    for ( i = 0; i < n; ++ i)
      for (j = 0; j < self->storage[i].last_index; ++ j)
      {
        value = sp_matrix_element_ptr(self,self->storage[i].indexes[j],i);
        if ( !value )
          return FALSE;
        if (!EQL(*value,-self->storage[i].values[j]))
          return FALSE;
      }
  }
  else                          /* CCS */
  {
    for ( i = 0; i < n; ++ i)
      for (j = 0; j < self->storage[i].last_index; ++ j)
      {
        value = sp_matrix_element_ptr(self,i,self->storage[i].indexes[j]);
        if ( !value )
          return FALSE;
        if (! EQL(*value,-self->storage[i].values[j]))
          return FALSE;
      }
    
  }
  /* check diagonal - it shall contain zeros */
  for ( i = 0; i < n; ++ i)
  {
    if (sp_matrix_element_ptr(self,i,i))
      return FALSE;
  }
  return TRUE;
}


int sp_matrix_nonzeros(sp_matrix_ptr self)
{
  int stored = 0;
  int i;
  int n = self->storage_type == CRS ? self->rows_count : self->cols_count;
  for ( i = 0; i < n; ++ i)
    stored += self->storage[i].last_index + 1;
  return stored;
}



void sp_matrix_mv(sp_matrix_ptr self,double* x, double* y)
{
  int i,j;
  memset(y,0,sizeof(double)*self->rows_count);
  if (self->storage_type == CRS)
  {
    for ( i = 0; i < self->rows_count; ++ i)
    {
      for ( j = 0; j <= self->storage[i].last_index; ++ j)
        y[i] += self->storage[i].values[j]*x[self->storage[i].indexes[j]];
    }
  }
  else                          /* CCS */
  {
    for ( j = 0; j < self->cols_count; ++ j)
    {
      for ( i = 0; i<= self->storage[j].last_index; ++ i)
        y[self->storage[j].indexes[i]] += self->storage[j].values[i]*x[j];
    }
  }
}

void sp_matrix_yale_mv(sp_matrix_yale_ptr self,double* x, double* y)
{
  int i,j;
  memset(y,0,sizeof(double)*self->rows_count);
  if (self->storage_type == CRS)
  {
    for ( i = 0; i < self->rows_count; ++ i)
      for ( j = self->offsets[i]; j < self->offsets[i+1]; ++ j)
        y[i] += self->values[j]*x[self->indicies[j]];
  }
  else                          /* CCS */
  {
    for ( i = 0; i < self->cols_count; ++ i)
      for ( j = self->offsets[i]; j < self->offsets[i+1]; ++ j)
        y[self->indicies[j]] += self->values[j]*x[i];
  }
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
  int count;
  
  /* clear the out */
  for (j = 0; j < self->rows_count; ++ j)
    out[j] = -1;
  /* mark k */
  out[k] = k;
  count = 1;
  /* loop by nonzero rows in the column */
  for ( p = self->offsets[k]; p < self->offsets[k+1]; ++p )
  {
    /* store the index to i for simplicity */
    i = self->indicies[p];
    if ( i < k )              /* if above the diagonal */
    {
      /* find the root */
      for (j = etree[k]; j != -1 && j != etree[j]; j = etree[j])
      {
        if ( out[j] != j)
        {
          out[j] = j;
          count ++;
        }
      }
    }
  }
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



void sp_matrix_lower_solve(sp_matrix_ptr self,
                           int n,
                           double* b,
                           double* x)
{
  int i,j;
  assert(self);
  assert(b);
  assert(x);
  assert(n>0 && n <= self->rows_count);

  memset(x,0,sizeof(double)*n);
  
  if (!self->ordered)
    sp_matrix_reorder(self);
  if (self->storage_type == CCS)
  {
    for ( j = 0; j < n; ++ j)
      x[j] = b[j];
    for ( j = 0; j < n; ++ j)
    {
      x[j] /= self->storage[j].values[0]; 
      for (i = 1; i <= self->storage[j].last_index; ++ i)
        x[self->storage[j].indexes[i]] -= x[j]*self->storage[j].values[i];
    }
  }
  else                          /* CRS */
  {
    for ( i = 0; i < n; ++ i)
    {
      x[i] = b[i];
      for (j = 0; j <= self->storage[i].last_index &&
             self->storage[i].indexes[j] <= i-1; ++ j)
        x[i] -= x[self->storage[i].indexes[j]]*self->storage[i].values[j];
      x[i] /= self->storage[i].values[self->storage[i].last_index];
    }
  }
}

#if 0
void sp_matrix_solve(sp_matrix_ptr self,double* b,double* x)
{
  double tolerance = TOLERANCE;
  int max_iter = MAX_ITER;
  int i;
  double tol = 0;
  double* r = (double*)malloc(self->rows_count*sizeof(double));
  memset(r,0,self->rows_count*sizeof(double));
  /* reorder columns for to prepare to solve SLAE */
  if (!self->ordered)
    sp_matrix_reorder(self);
  /* sp_matrix_solve_pcg(self,b,b,&max_iter,&tolerance,x); */
  sp_matrix_solve_cg(self,b,b,&max_iter,&tolerance,x);
  /* Calculare residual r = A*x-b */
  sp_matrix_mv(self,x,r);
  for ( i = 0; i < self->rows_count; ++ i)
    r[i] -= b[i];
  /* Find a norm of residual vector */
  for ( i = 0; i < self->rows_count; ++ i)
    tol += r[i]*r[i];
  tol = sqrt(tol);
  /* TODO: move iter, tolerance1 and tolerance2 to the output parameters */
  printf("iter = %d, tolerance1 = %e, tolerance2 = %e\n",
         max_iter,tolerance,tol);
  free(r);
}
#endif

void sp_matrix_yale_solve_cg(sp_matrix_yale_ptr self,
                             double* b,
                             double* x0,
                             int* max_iter,
                             double* tolerance,
                             double* x)
{
  /* Conjugate Gradient Algorithm */
  /*
   * Taken from the book:
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
  for ( i = 0; i < msize; ++ i)
    x[i] = x0[i];

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
      a1 += r[i]*r[i]; /* (r_j,r_j) */
      a2 += p[i]*temp[i];      /* (A*p_j,p_j) */
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
    
    /* d_{j+1} = r_{j+1} + beta_j*d_j */
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
   * Taken from the book:
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
  for ( i = 0; i < msize; ++ i)
    x[i] = x0[i];

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

void sp_matrix_printf(sp_matrix_ptr self)
{
  int i,n;
  n = self->storage_type == CRS ? self->rows_count : self->cols_count;
  for (i = 0; i < n; ++ i)
    indexed_array_printf(&self->storage[i]);
}

void sp_matrix_dump(sp_matrix_ptr self, const char* filename)
{
  FILE* f;
  int i,j;
  double *pvalue,value;

  if ((f = fopen(filename,"w+")))
  {
    for (i = 0; i < self->rows_count; ++ i)
    {
      for (j = 0; j < self->rows_count; ++ j)
      {
        pvalue = sp_matrix_element_ptr(self,i,j);
        value = pvalue ? *pvalue : 0;
        fprintf(f,"%e ",value);
      }
      fprintf(f,"\n");
    }
    fflush(f);
    fclose(f);
  }
}


void sp_matrix_skyline_dump(sp_matrix_skyline_ptr self,const char* filename)
{
  int i;
  FILE* f;
  if ((f = fopen(filename,"w+")))
  {
    fprintf(f,"adiag = [");
    for ( i = 0; i < self->rows_count; ++ i )
      fprintf(f,"%f ",self->diag[i]);
    fprintf(f,"]\n");

    fprintf(f,"altr = [");
    for ( i = 0; i < self->tr_nonzeros; ++ i )
      fprintf(f,"%f ",self->lower_triangle[i]);
    fprintf(f,"]\n");

    fprintf(f,"autr = [");
    for ( i = 0; i < self->tr_nonzeros; ++ i )
      fprintf(f,"%f ",self->upper_triangle[i]);
    fprintf(f,"]\n");
  
    fprintf(f,"jptr = [");
    for ( i = 0; i < self->tr_nonzeros; ++ i )
      fprintf(f,"%d ",self->jptr[i]+1);
    fprintf(f,"]\n");

    fprintf(f,"iptr = [");
    for ( i = 0; i < self->rows_count; ++ i )
      fprintf(f,"%d ",self->iptr[i]+1);
    fprintf(f,"]\n");
    fclose(f);
  }
}
