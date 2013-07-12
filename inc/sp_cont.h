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

#ifndef _SP_CONT_H_
#define _SP_CONT_H_

/*
 * Sparse matrix row/column storage array
 */
typedef struct
{
  int width;                    /* size of an array */
  int last_index;               /* last stored index, i.e. if width = 20
                                 * it will be 9 if only 10 nonzero elements
                                 * stored */
  int  *indexes;                /* array of column/row indexes */
  double *values;               /* array of values */
} indexed_array;
typedef indexed_array* indexed_array_ptr;


/*
 * Dynamic array data structure
 */
typedef struct
{
  int allocated;
  int step_size;
  int* items;
  int current;
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
typedef struct int_queue_element_tag
{
  int value;
  struct int_queue_element_tag* next;
} int_queue_element;

typedef struct 
{
  int size;
  int_queue_element* first;
  int_queue_element* last;
} int_queue;
typedef int_queue* int_queue_ptr;

/*
 * indexed_arrays operations
 */
/* Performs in-place sort of the indexed array */
void indexed_array_sort(indexed_array_ptr self, int l, int r);
/* Print contents of the indexed array to the stdout  */
void indexed_array_printf(indexed_array_ptr self);



/*
 * Dynamic array functions
 */
void int_array_init(int_array_ptr self, int initial_size, int step_size);
void int_array_extend(int_array_ptr self);
void int_array_free(int_array_ptr self);
void int_array_add(int_array_ptr self, int value);


/*
 * Stack functions
 */

/* allocate stack with specified initial size and increment step */
int_stack_ptr int_stack_alloc(int initial_size, int step_size);
/* deallocate stack and all its data */
int_stack_ptr int_stack_free(int_stack_ptr self);

/* push element to the top of the stack */
void int_stack_push(int_stack_ptr self, int value);
/* pop element. Do nothing if the stack is empty */
void int_stack_pop(int_stack_ptr self);

/* returns nonzero if stack is empty */
int int_stack_isempty(int_stack_ptr self);
/* returns value of the element on top of the stack */
int int_stack_top(int_stack_ptr self);

/*
 * Queue functions
 */

/* allocate queue and all its structures */
int_queue_ptr int_queue_alloc();
/* deallocate queue and all its structures */
int_queue_ptr int_queue_free(int_queue_ptr self);

/* Returns nonzero if the queue is empty, and 0 otherwise. */
int int_queue_isempty(int_queue_ptr self);

/* Returns a value at the front of a non-empty queue. */
int int_queue_front(int_queue_ptr self);

/* Removes the item at the front of a non-empty queue. */
void int_queue_pop(int_queue_ptr self);

/* Inserts the argument value at the back of the queue. */
void int_queue_push(int_queue_ptr self, int value);

#endif /* _SP_CONT_H_ */
