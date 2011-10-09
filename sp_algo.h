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

#ifndef _SP_ALGO_H_
#define _SP_ALGO_H_

#define DSU_DEFAULT_VALUE -1

typedef int dsu_value;

typedef struct 
{
  dsu_value* values;
  long size;
} disjoint_set_union;
typedef disjoint_set_union* disjoint_set_union_ptr;

disjoint_set_union_ptr dsu_alloc(long size);
disjoint_set_union_ptr dsu_free(disjoint_set_union_ptr self);

disjoint_set_union_ptr dsu_make_set(disjoint_set_union_ptr self,
                                    dsu_value value);
dsu_value dsu_find(disjoint_set_union_ptr self, dsu_value value);
disjoint_set_union_ptr dsu_union(disjoint_set_union_ptr self,
                                 dsu_value x,
                                 dsu_value y);


#endif /* _SP_ALGO_H_ */
