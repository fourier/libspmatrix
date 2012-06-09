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
#include <limits.h>
#include <string.h>

#include "sp_cont.h"

void int_array_init(int_array_ptr self, int initial_size, int step_size)
{
  self->allocated = initial_size;
  self->step_size = step_size;
  self->items = calloc(self->allocated,sizeof(int));
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


static
int_queue_element* int_queue_element_alloc(int value, int_queue_element* next)
{
  int_queue_element* element = malloc(sizeof(int_queue_element));
  element->value = value;
  element->next  = next;
  return element;
}

static
int_queue_element* int_queue_element_free(int_queue_element* self)
{
  free(self);
  return 0;
}

int_queue_ptr int_queue_alloc()
{
  int_queue_ptr self = malloc(sizeof(int_queue));
  memset(self,0, sizeof(int_queue));
  return self;
}


int_queue_ptr int_queue_free(int_queue_ptr self)
{
  if (self)
  {
    while(!int_queue_isempty(self))
      int_queue_pop(self);
    free(self);
    self = 0;
  }
  return self;
}


int int_queue_isempty(int_queue_ptr self)
{
  return self->size == 0;
}


int int_queue_front(int_queue_ptr self)
{
  return self->size ? self->first->value : INT_MAX;
}


void int_queue_pop(int_queue_ptr self)
{
  if ( self->size )
  {
    int_queue_element* first = self->first->next;
    (void)int_queue_element_free(self->first);
    self->first = first;
    self->size--;
    if ( !self->size )
      self->last = 0;
  }
}


void int_queue_push(int_queue_ptr self, int value)
{
  if (self->size)
  {
    self->last->next = int_queue_element_alloc(value,0);
    self->last = self->last->next;
  }
  else
  {
    self->first = int_queue_element_alloc(value,0);
    self->last = self->first;
  }
  self->size++;
}


