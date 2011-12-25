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
 * Tests Queue element
 */
typedef struct sp_test_queue_element_tag
{
  sp_test_object value;
  struct sp_test_queue_element_tag* next;
} sp_test_queue_element;

/*
 * Tests Queue
 */
typedef struct 
{
  int size;
  sp_test_queue_element* first;
  sp_test_queue_element* last;
} sp_test_queue;
typedef sp_test_queue* sp_test_queue_ptr;

/*
 * Tests suites list 
 */
typedef struct sp_test_suite_list_tag
{
  struct sp_test_suite_list_tag* next;
  sp_test_queue_ptr tests;
  sp_test_suite_ptr suite;
} sp_test_suite_list;
typedef sp_test_suite_list* sp_test_suite_list_ptr;


/*
 * Tests queue implementation
 */

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
sp_test_suite_list_ptr g_test_suite;

jmp_buf g_test_restart_jmp_point;



void sp_add_test(test_func_t func, const char* name)
{
  sp_test_object test = {name,func};
  if (!g_test_queue)
    g_test_queue = sp_test_queue_alloc();
  
  sp_test_queue_push(g_test_queue, test);
}

static int sp_is_test(sp_test_object test, const char* name)
{
  return strcmp(test.test_name,name) == 0;
}

static int sp_test_exist(sp_test_queue_ptr test_queue, const char* name)
{
  sp_test_queue_element* first = test_queue->first;
  sp_test_queue_element* last  = test_queue->last;
  for ( ; first != last; first = first->next)
    if (sp_is_test(first->value,name) )
      return 1;
  return sp_is_test(first->value,name);
}


static void sp_perform_tests(sp_test_queue_ptr test_queue,
                             int argc,
                             const char* argv[])
{
  sp_test_object current_test;
  const char* test_name = 0;
  if (argc > 1)
    test_name = argv[1];
  while(!sp_test_queue_isempty(test_queue))
  {
    current_test = sp_test_queue_front(test_queue);
    if (!test_name || (test_name && sp_is_test(current_test,test_name)))
    {
      if (!setjmp(g_test_restart_jmp_point))
      {
        current_test.test_func();
        printf("test %s [PASSED]\n",
               current_test.test_name);
      }
      else
      {
        printf("test %s [FAILED]\n",
               current_test.test_name);
      }
    }
    sp_test_queue_pop(test_queue);
  }
}

static void sp_perform_suites(sp_test_suite_list_ptr suite_head,
                              int argc,
                              const char* argv[])
{
  sp_test_suite_list_ptr suite;
  const char* test_name = 0;
  if (argc > 1)
    test_name = argv[1];

  while(suite_head)
  {
    suite = suite_head->next;
    if (!sp_test_queue_isempty(suite_head->tests))
    {
      if (!test_name ||
          (test_name && sp_test_exist(suite_head->tests, test_name)))
      {
        printf("Suite %s start:\n",suite_head->suite->test_suite_name);
        if (suite_head->suite->test_suite_init)
          suite_head->suite->test_suite_init();
        sp_perform_tests(suite_head->tests,argc,argv);
        if (suite_head->suite->test_suite_fini)
          suite_head->suite->test_suite_fini();
        printf("Suite %s end\n",suite_head->suite->test_suite_name);
      }
    }
    sp_test_queue_free(suite_head->tests);
    free(suite_head->suite);
    free(suite_head);
    suite_head = suite;
  }
}

void sp_run_tests(int argc, const char* argv[])
{
  /*
   * TODO:
   * 1) fix memory leaks
   * 2) possibility to list available tests with --list cmdline option
   * 3) remove 'Running..' lines when no tests/suites available(or found)
   * 4) add tests results backend - log file, html report etc.
   */
  if (g_test_queue)
  {
    printf("Running individual tests...\n");
    sp_perform_tests(g_test_queue, argc, argv);
    g_test_queue = sp_test_queue_free(g_test_queue);
  }
  if (g_test_suite)
  {
    printf("Running suites...\n");
    sp_perform_suites(g_test_suite, argc, argv);
    g_test_suite = 0;
  }
}

sp_test_suite_ptr sp_add_suite(const char* name,
                               void(*test_suite_init)(),
                               void(*test_suite_fini)())
{
  /* create global suite if not created yet */
  sp_test_suite_list_ptr suite,next;
  if (!g_test_suite)
  {
    g_test_suite = calloc(1,sizeof(sp_test_suite_list));
    suite = g_test_suite;
  }
  else
  {
    /*
     * otherwise find the last suite and set the next pointer
     * to the newly created suite
     */
    suite = calloc(1,sizeof(sp_test_suite_list));
    next = g_test_suite;
    while(next->next)
      next = next->next;
    next->next = suite;
  }
  /* create suite */
  suite->suite = calloc(1,sizeof(sp_test_suite));
  suite->tests = sp_test_queue_alloc();
  suite->suite->test_suite_name = name;
  suite->suite->test_suite_init = test_suite_init;
  suite->suite->test_suite_fini = test_suite_fini;
  return suite->suite;
}

void sp_add_suite_test(sp_test_suite_ptr suite,
                       test_func_t func,
                       const char* name)
{
  sp_test_suite_list_ptr next;
  sp_test_object test = {name,func};
  /* find appropriate suite */
  next = g_test_suite;
  while (next->suite != suite && next->next)
    next = next->next;
  if (next && next->suite == suite)
  {
    sp_test_queue_push(next->tests, test);
  }
  else
    fprintf(stderr,"ERROR: Suite %s not found, cannot add test %s\n",
            suite->test_suite_name,
            name);
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


