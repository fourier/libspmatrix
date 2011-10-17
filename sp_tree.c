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

#include <stdlib.h>
#include <limits.h>


#include "sp_tree.h"

void int_array_init(int_array_ptr self, int initial_size, int step_size)
{
  self->allocated = initial_size;
  self->step_size = step_size;
  self->items = malloc(sizeof(int)*self->allocated);
}

void int_array_extend(int_array_ptr self)
{
    self->items = realloc(self->items,
                          sizeof(int)*(self->allocated + self->step_size));
    self->allocated += self->step_size;
}

void int_array_free(int_array_ptr self)
{
  free(self->items);
}

int_stack_ptr int_stack_alloc(int initial_size, int step_size)
{
  int_stack_ptr stack = malloc(sizeof(int_stack));
  stack->top = -1;
  int_array_init(&stack->data,initial_size,step_size);
  return stack;
}


int_stack_ptr int_stack_free(int_stack_ptr self)
{
  if (self)
  {
    int_array_free(&self->data);
    free(self);
    self = 0;
  }
  return self;
}


void int_stack_push(int_stack_ptr self, int value)
{
  if (self->top == self->data.allocated-1)
  {
    int_array_extend(&self->data);
  }
  self->data.items[++self->top] = value;
}


void int_stack_pop(int_stack_ptr self)
{
  if (self->top >= 0)
    self->top--;
}


int int_stack_isempty(int_stack_ptr self)
{
  return self->top == -1;
}


int int_stack_top(int_stack_ptr self)
{
  if (self->top >= 0)
    return self->data.items[self->top];
  return INT_MAX;
}

