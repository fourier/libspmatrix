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

#ifndef _SP_TREE_H_
#define _SP_TREE_H_

/*
 * Dynamic array data structure
 */

typedef struct
{
  int allocated;
  int step_size;
  int* items;
  
} int_array;
typedef int_array* int_array_ptr;

/*
 * Stack data structure
 */

typedef struct
{
  int top;
  int_array data;
} int_stack;
typedef int_stack* int_stack_ptr;

/*
 * Queue data structure
 */


/*
 * Dynamic array functions
 */
void int_array_init(int_array_ptr self, int initial_size, int step_size);
void int_array_extend(int_array_ptr self);
void int_array_free(int_array_ptr self);


/*
 * Stack functions
 */

int_stack_ptr int_stack_alloc(int initial_size, int step_size);
int_stack_ptr int_stack_free(int_stack_ptr self);

void int_stack_push(int_stack_ptr self, int value);
void int_stack_pop(int_stack_ptr self);

int int_stack_isempty(int_stack_ptr self);
int int_stack_top(int_stack_ptr self);

/*
 * Queue functions
 */


#endif /* _SP_TREE_H_ */
