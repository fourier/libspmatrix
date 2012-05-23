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
#include <string.h>

#include "sp_direct.h"
#include "sp_tree.h"


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
  return result;
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

int sp_matrix_yale_symbolic_init(sp_matrix_yale_ptr self,
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
    symb->nonzeros = 0;
    symb->etree = 0;
    symb->rowcounts = 0;
    symb->colcounts = 0;
  }
}
