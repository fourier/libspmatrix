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

#include "sp_algo.h"



disjoint_set_union_ptr dsu_alloc(long size)
{
  int i;
  disjoint_set_union_ptr self = malloc(sizeof(disjoint_set_union));
  self->size = size;
  self->values = malloc(sizeof(dsu_value)*size);
  for (i = 0; i < size; ++ i)
    self->values[i] = DSU_DEFAULT_VALUE;
  return self;
}

disjoint_set_union_ptr dsu_free(disjoint_set_union_ptr self)
{
  if (self)
  {
    free(self->values);
    free(self);
    self = 0;
  }
  return self;
}

disjoint_set_union_ptr dsu_make_set(disjoint_set_union_ptr self,
                                    dsu_value value)
{
  if (self)
  {
    if ( value >= 0 && value < self->size) 
    self->values[value] = value;
    else
      fprintf(stderr,
              "dsu_make_set: value %d out of range %ld", value, self->size);
  }
  return self;
}

dsu_value dsu_find(disjoint_set_union_ptr self, dsu_value value)
{
  if ( !(value >= 0 && value < self->size))
  {
    fprintf(stderr, "dsu_find: value %d out of range %ld", value, self->size);
    return DSU_DEFAULT_VALUE;
  }
    
  if (!self)
    return DSU_DEFAULT_VALUE;
  if (self->values[value] == DSU_DEFAULT_VALUE)
    return DSU_DEFAULT_VALUE;
  if (self->values[value] == value)
    return value;
  /* path compression heuristic  */
  return self->values[value] = dsu_find(self,self->values[value]);
}

disjoint_set_union_ptr dsu_union(disjoint_set_union_ptr self,
                                 dsu_value x,
                                 dsu_value y)
{
  if (self)
  {
    x = dsu_find(self,x);
    y = dsu_find(self,y);
    if (x == y)
      return self;

    /* randomly select to which tree append */
    /* if (rand() & 1) */
    /*   self->values[x] = y; */
    /* else */
    /*   self->values[y] = x; */

    /* append y to x */
    self->values[y] = x;
  }
  return self;
}
