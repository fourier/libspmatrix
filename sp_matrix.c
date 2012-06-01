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

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <math.h>
#include <assert.h>

#include "sp_matrix.h"
#include "sp_utils.h"
#include "sp_tree.h"
#ifdef USE_LOGGER
#include "logger.h"
#else
#define LOGINFO(...) 
#define LOG(...) printf(__VA_ARGS__);
#define LOGWARN(...) printf(__VA_ARGS__);
#define LOGERROR(...) fprintf(stderr,__VA_ARGS__);
#endif

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
      memdup(mtx_from->storage[i].indexes,
             sizeof(int)*mtx_from->storage[i].width);
    mtx_to->storage[i].values =
      memdup(mtx_from->storage[i].values,
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

void sp_matrix_yale_init2(sp_matrix_yale_ptr self,
                          sparse_storage_type type,
                          int rows_count,
                          int cols_count,
                          int nonzeros,
                          int* counts)
{
  int i,j;
  int n = type == CRS ? rows_count : cols_count;
  /* initialize matrix */
  memset(self,sizeof(sp_matrix_yale),0);
  self->storage_type = type;
  self->rows_count = rows_count;
  self->cols_count = cols_count;
  self->nonzeros   = nonzeros;
  /* allocate memory for arrays */
  self->offsets  = calloc(n+1,      sizeof(int));
  self->indicies = calloc(nonzeros, sizeof(int));
  self->values   = calloc(nonzeros, sizeof(double));
  /* calculate offsets */
  j = 0;
  for (i = 0; i < n; ++ i)
  {
    self->offsets[i] = j;
    j += counts[i];
  }
  self->offsets[i] = nonzeros;
}

void sp_matrix_yale_copy(sp_matrix_yale_ptr mtx_from,
                         sp_matrix_yale_ptr mtx_to)
{
  int n = (mtx_from->storage_type == CRS) ?
    mtx_from->rows_count : mtx_from->cols_count;
  /* initialize matrix */
  memset(mtx_to,sizeof(sp_matrix_yale),0);
  mtx_to->storage_type = mtx_from->storage_type;
  mtx_to->rows_count   = mtx_from->rows_count;
  mtx_to->cols_count   = mtx_from->cols_count;
  mtx_to->nonzeros     = mtx_from->nonzeros;
  /* copy data */
  mtx_to->offsets  = memdup(mtx_from->offsets,(n+1)*sizeof(int));
  mtx_to->indicies = memdup(mtx_from->indicies,mtx_from->nonzeros*sizeof(int));
  mtx_to->values   = memdup(mtx_from->values,mtx_from->nonzeros*sizeof(double));
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


void sp_matrix_printf2(sp_matrix_ptr self)
{
  double *p;
  int i,j;
  if (self)
  {
    for (i = 0; i < self->rows_count; ++ i)
    {
      for (j = 0; j < self->cols_count; ++ j)
      {
        p = sp_matrix_element_ptr(self,i,j);
        printf("%.2f ", p ? *p : 0);
      }
      printf("\n");
    }
  }
}


void sp_matrix_yale_transpose(sp_matrix_yale_ptr self,
                              sp_matrix_yale_ptr to)
{
  int n = self->storage_type == CRS ? self->rows_count : self->cols_count;
  int i,j,k,p;
  int* offsets = calloc(n+1,sizeof(int));
  /* 1. row/column counts shifted by 1 */
  for ( i = 0; i < self->nonzeros; ++i)
  {
    offsets[self->indicies[i]+1]++;
  }
  /* 2. initialize an empty matrix */
  sp_matrix_yale_init2(to,self->storage_type,
                       self->rows_count,self->cols_count,
                       self->nonzeros,offsets+1);

  /* 3. offsets - partial sums of counts of row/columns */
  memcpy(offsets,to->offsets,(n+1)*sizeof(int));
  for ( i = 0; i < n; ++i)
  {
    for ( p = self->offsets[i]; p < self->offsets[i+1]; ++p )
    {
      j = self->indicies[p];
      /* a_ij != 0 */
      /* in new matrix a'_ji = a_ij 
       * offsets[j] - shall point to the beginning of the jth column
       * therefore offsets[j]++ - next nonzero row in this column,
       * which is the element a_ji, and therefore its row index = i
       */ 
      k = offsets[j]++;
      to->indicies[k] = i;
      to->values[k]   = self->values[p];
    }
  }
  free(offsets);
}

void sp_matrix_yale_convert(sp_matrix_yale_ptr from,
                            sp_matrix_yale_ptr to,
                            sparse_storage_type type)
{
  if (from->storage_type == type)
    return;
  sp_matrix_yale_transpose(from,to);
  to->storage_type = type;
}



int sp_matrix_yale_permute(sp_matrix_yale_ptr self,
                           sp_matrix_yale_ptr permuted,
                           int* p,
                           int* q)
{
  int i,j,n;
  int result = 1;
  sp_matrix mtx;
  /* init matrix with average bandwidth */
  n = self->storage_type == CRS ? self->rows_count : self->cols_count;
  sp_matrix_init(&mtx,
                 self->rows_count,self->cols_count,
                 (int)(n/(double)self->nonzeros),
                 self->storage_type);
  
  for (i = 0; i < n; ++ i)
  {
    for (j = self->offsets[i]; j < self->offsets[i+1]; ++ j)
    {
      if ( self->storage_type == CRS )
      {
        MTX(&mtx,p[i],q[self->indicies[j]],self->values[j]);
      }
      else
      {
        MTX(&mtx,p[self->indicies[j]],q[i],self->values[j]);
      }
    }
  }
  sp_matrix_yale_init(permuted,&mtx);
  sp_matrix_free(&mtx);
  return result;
}




void sp_matrix_printf(sp_matrix_ptr self)
{
  int i,n;
  n = self->storage_type == CRS ? self->rows_count : self->cols_count;
  for (i = 0; i < n; ++ i)
    indexed_array_printf(&self->storage[i]);
}

void sp_matrix_yale_printf(sp_matrix_yale_ptr self)
{
  int i;
  int n = self->storage_type == CRS ? self->rows_count : self->cols_count;
  printf("offsets: [ ");
  for ( i = 0; i < n; ++ i)
    printf("%d ", self->offsets[i]+1);
  printf("%d ]\n", self->offsets[i]);
  printf("indicies: [ ");
  for ( i = 0; i < self->nonzeros; ++ i)
    printf("%d ", self->indicies[i]+1);
  printf("]\n");
  printf("values: [ ");
  for ( i = 0; i < self->nonzeros; ++ i)
    printf("%.1f ", self->values[i]);
  printf("]\n");
}

void sp_matrix_yale_printf2(sp_matrix_yale_ptr self)
{
  printf("Storage type: %s\n", self->storage_type == CRS ? "CRS" : "CCS");
  printf("Size: %dx%d\n", self->rows_count,self->cols_count);
  printf("Nonzeros: %d\n", self->nonzeros);
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
