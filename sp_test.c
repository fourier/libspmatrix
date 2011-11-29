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

#include "sp_test.h"
#include <memory.h>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <limits.h>
#include <assert.h>
#include <setjmp.h>

/*
 * Structure encapsulating test
 */
typedef struct
{
  const char* test_name;
  test_func_t test_func;
} sp_test_object;


/*
 * Queue data structure
 */
typedef struct sp_test_queue_element_tag
{
  sp_test_object value;
  struct sp_test_queue_element_tag* next;
} sp_test_queue_element;

typedef struct 
{
  int size;
  sp_test_queue_element* first;
  sp_test_queue_element* last;
} sp_test_queue;
typedef sp_test_queue* sp_test_queue_ptr;

static
sp_test_queue_element* sp_test_queue_element_alloc(sp_test_object value, sp_test_queue_element* next)
{
  sp_test_queue_element* element = malloc(sizeof(sp_test_queue_element));
  element->value = value;
  element->next  = next;
  return element;
}

static
sp_test_queue_element* sp_test_queue_element_free(sp_test_queue_element* self)
{
  free(self);
  return 0;
}

static sp_test_queue_ptr sp_test_queue_alloc()
{
  sp_test_queue_ptr self = malloc(sizeof(sp_test_queue));
  memset(self,0, sizeof(sp_test_queue));
  return self;
}

static sp_test_object sp_test_queue_front(sp_test_queue_ptr self)
{
  assert(self->size);
  return self->first->value;
}

static int sp_test_queue_isempty(sp_test_queue_ptr self)
{
  return self->size == 0;
}

static void sp_test_queue_pop(sp_test_queue_ptr self)
{
  if ( self->size )
  {
    sp_test_queue_element* first = self->first->next;
    (void)sp_test_queue_element_free(self->first);
    self->first = first;
    self->size--;
    if ( !self->size )
      self->last = 0;
  }
}

static sp_test_queue_ptr sp_test_queue_free(sp_test_queue_ptr self)
{
  if (self)
  {
    while(!sp_test_queue_isempty(self))
      sp_test_queue_pop(self);
    free(self);
    self = 0;
  }
  return self;
}

static void sp_test_queue_push(sp_test_queue_ptr self, sp_test_object value)
{
  if (self->size)
  {
    self->last->next = sp_test_queue_element_alloc(value,0);
    self->last = self->last->next;
  }
  else
  {
    self->first = sp_test_queue_element_alloc(value,0);
    self->last = self->first;
  }
  self->size++;
}

sp_test_queue_ptr g_test_queue;
jmp_buf g_test_restart_jmp_point;
/* int setjmp(jmp_buf env); */
/* void longjmp(jmp_buf env, int value); */

void sp_add_test(test_func_t func, const char* name)
{
  sp_test_object test = {name,func};
  if (!g_test_queue)
    g_test_queue = sp_test_queue_alloc();
  
  sp_test_queue_push(g_test_queue, test);
}

void sp_run_tests()
{
  if (g_test_queue)
    while(!sp_test_queue_isempty(g_test_queue))
    {
      sp_test_object current_test = sp_test_queue_front(g_test_queue);
      if (!setjmp(g_test_restart_jmp_point))
      {
        current_test.test_func();
        printf("test %s *passed*\n",
               current_test.test_name);
      }
      else
      {
        printf("test %s *failed*\n",
               current_test.test_name);
      }
      sp_test_queue_pop(g_test_queue);
    }
  g_test_queue = sp_test_queue_free(g_test_queue);
}

void sp_assertion_failed(const char* file, int line, const char* condition)
{
  printf("%s:%d: error: Assertion failed: %s\n",file,line,condition);
  longjmp(g_test_restart_jmp_point,0);
}

void sp_expectation_failed(const char* file, int line, const char* condition)
{
  printf("%s:%d: warning: Expectation failed: %s\n",file,line,condition);
}
