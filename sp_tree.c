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
#include <string.h>

#include "sp_tree.h"
#include "sp_cont.h"


void tree_dfs(int* tree, int size, traverse_func_t func, void* arg)
{
  int i,j,current;
  /* initialize container */
  int_stack_ptr stack = int_stack_alloc(size,1);
  /* loop by all roots in the forest */
  for ( i = 0; i < size; ++ i)
  {
    if (tree[i] == -1)
    {
      int_stack_push(stack,i);
      /* start to walk through the tree starting with roots */
      while (!int_stack_isempty(stack))
      {
        current = int_stack_top(stack);
        int_stack_pop(stack);
        if (func(current,arg))      /* if predicate func is true, break */
          break;
        /* find all children and push them into container */
        for ( j = 0; j < size; ++ j)
          if (tree[j] == current)   /* current is the parent of j */
            int_stack_push(stack,j);
      }
    }
  }
  int_stack_free(stack);
}

void tree_bfs(int* tree, int size, traverse_func_t func, void* arg)
{
  int i,j,current;
  /* initialize container */
  int_queue_ptr queue = int_queue_alloc();
  /* loop by all roots in the forest */
  for ( i = 0; i < size; ++ i)
  {
    if (tree[i] == -1)
    {
      int_queue_push(queue,i);
      /* start to walk through the tree starting with roots */
      while (!int_queue_isempty(queue))
      {
        current = int_queue_front(queue);
        int_queue_pop(queue);
        if (func(current,arg))      /* if predicate func is true, break */
          break;
        /* find all children and push them into container */
        for ( j = 0; j < size; ++ j)
          if (tree[j] == current)   /* current is the parent of j */
            int_queue_push(queue,j);
      }
    }
  }
  int_queue_free(queue);
}

void tree_postorder_perm(int* tree, int size, int* postorder)
{
  /* Algorithm:
   * 
   * function postorder(T)
   * k = 0
   * for each root node j in T do
   *   dfstree j
   *
   * function dfstree(j)
   * for each child i of j do
   *   dfstree(i)
   * post[k] = j
   * k = k + 1
   */
  int counter = 0;
  int i,child,current;
  int_stack_ptr* children;
  /* initialize stack for DFS procedure */
  int_stack_ptr stack = int_stack_alloc(size,1);
  /* allocate memory for childen
   * children of every node stored in the stack */
  children = malloc(sizeof(int_stack_ptr)*size);
  for ( i = 0; i < size; ++ i)
    children[i] = int_stack_alloc(1,2);
  /* now fill the stacks of children with the actual children
   * in reverse order to keep children in the stack in ascending order
   * i.e. ->[1 2 3] */
  for ( i = size -1 ; i >= 0; -- i)
    if (tree[i] != -1)
      int_stack_push(children[tree[i]],i);

  /* loop by all roots in the forest */
  for ( i = 0; i < size; ++ i)
  {
    /* root node in the forest. if forest is the 1 tree, only one
     * root will exist */
    if (tree[i] == -1)
    {
      int_stack_push(stack,i);
      /* start to walk through the tree starting with roots */
      while (!int_stack_isempty(stack))
      {
        current = int_stack_top(stack);
        /* if exist not visited children */
        if (!int_stack_isempty(children[current]))
        {
          /* take the not visited child */
          child = int_stack_top(children[current]);
          /* and put it to the stack */
          int_stack_push(stack,child);
          /* remove child from the list of not visited */
          int_stack_pop(children[current]);
        }
        else                    /* no more not visited children */
        {
          /* remove current from the stack to traverse up by the tree  */
          int_stack_pop(stack);
          /* write the last visited node to the place
           * specified by the counter  */
          postorder[counter] = current;
          /* increase the counter for postordering */
          counter++;
        }
      }
    }
  }
  int_stack_free(stack);
}

void tree_postorder(int* tree, int size, int* postorder, int* result)
{
  int i;
  for ( i = 0; i < size; ++ i)
  {
    result[i] = tree[postorder[i]];
  }
}
