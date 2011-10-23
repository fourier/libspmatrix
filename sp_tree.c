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

typedef struct
{
  int* tree;
  int* postordered;
  int_stack_ptr* children;
  int counter;
} tree_postorder_data;

static int tree_postorder_helper(int n, void* arg)
{
  printf("visiting %d\n",n);
  tree_postorder_data* data = (tree_postorder_data*)arg;
  int_stack_ptr parents_childen;
  int parent;
  /* check if current node has no children */
  if (int_stack_isempty(data->children[n]))
  {
    parent = data->tree[n];
    /* check if current node has parent */
    if (parent != -1)
    {
      /* has parent */
      parents_childen = data->children[parent];
      /* remove one child(ourself) from parent's children list */
      if (!int_stack_isempty(parents_childen))
        int_stack_pop(parents_childen);
      /* and finally increase the counter */
      data->counter ++;
      printf("node %d: %d\n",n,data->counter);
      /* data->postordered[n] = data->counter; */
    }
  }
  return 0;
}

void tree_postorder(int* tree, int size, int* postordered)
{
  int k = 0;
  int i;
  tree_postorder_data data;
  /* create the data for postoredering */
  data.tree = tree;
  data.postordered = postordered;
  data.counter = 0;
  /* allocate memory for childen
   * children of every node stored in the stack */
  data.children = malloc(sizeof(int_stack_ptr)*size);
  for ( i = 0; i < size; ++ i)
    data.children[i] = int_stack_alloc(1,2);
  /* now fill the stacks of children with the actual children */
  for ( i = 0; i < size; ++ i)
    if (tree[i] != -1)
      int_stack_push(data.children[tree[i]],i);
  
  tree_dfs(tree,size,tree_postorder_helper,&data);

  for ( i = 0; i < size; ++ i)
  {
    printf("%d: ",i+1);
    while(!int_stack_isempty(data.children[i]))
    {
      printf("%d ", int_stack_top(data.children[i])+1);
      int_stack_pop(data.children[i]);
    }
    printf("\n");
  }
}
