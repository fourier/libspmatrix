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

#include <stdio.h>

#include "sp_tree.h"
#include "sp_mem.h"
#include "sp_cont.h"


void forest_dfs_preorder(int* tree, int size, traverse_func_t func, void* arg)
{
  int i,j,current;
  /* initialize container */
  int_stack_ptr stack = int_stack_alloc(size,1);
  /* loop by all roots in the forest */
  for ( i = 0; i < size; ++ i)
  {
    if (tree[i] == EMPTY_PARENT)
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

void tree_dfs_preorder(int* tree, int size, int root, traverse_func_t func, void* arg)
{
  int j,current;
  /* initialize container */
  int_stack_ptr stack = int_stack_alloc(size,1);
  int_stack_push(stack,root);
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
  int_stack_free(stack);
}

static int dfs_postorder_helper(int n, void* arg)
{
  int_stack_push((int_stack_ptr)arg,n);
  return 0;
}

void forest_dfs_postorder(int* tree, int size, traverse_func_t func, void* arg)
{
  int_stack_ptr stack = int_stack_alloc(size,1);
  forest_dfs_preorder(tree, size, dfs_postorder_helper, stack);
  while(!int_stack_isempty(stack))
  {
    func(int_stack_top(stack), arg);
    int_stack_pop(stack);
  }
  int_stack_free(stack);
}

void tree_dfs_postorder(int* tree, int size, int root, traverse_func_t func, void* arg)
{
  int_stack_ptr stack = int_stack_alloc(size,1);
  tree_dfs_preorder(tree, size, root, dfs_postorder_helper, stack);
  while(!int_stack_isempty(stack))
  {
    func(int_stack_top(stack), arg);
    int_stack_pop(stack);
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
    if (tree[i] == EMPTY_PARENT)
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

int_stack_ptr tree_children(int* tree, int size, int k)
{
  int i;
  /* initialize stack */
  int_stack_ptr children = int_stack_alloc(1,1);
  for (i = 0; i < size; ++ i)
    if (tree[i] == k)
      int_stack_push(children,i);
  return children;
}

int_stack_ptr tree_roots(int* tree, int size)
{
  return tree_children(tree, size, EMPTY_PARENT);
}

struct postorder_perm_helper
{
  int k;
  int* post;
};

static int tree_postorder_perm_helper(int n, void* arg)
{
  struct postorder_perm_helper* helper = (struct postorder_perm_helper*)arg;
  helper->post[helper->k] = n;
  helper->k++;
  return 0;
}

/* recursive version: a working code for reference only, not actual use */

static int dfs_recursive			/* return the new value of k */
(
    int p,		/* start a DFS at node p */
    int k,		/* start the node numbering at k */
    int Post [ ],	/* Post ordering, modified on output */
    int Head [ ],	/* Head [p] = youngest child of p; EMPTY on output */
    int Next [ ]	/* Next [j] = sibling of j; unmodified */
)
{
    int j ;
    /* start a DFS at each child of node p */
    for (j = Head [p] ; j != EMPTY_PARENT ; j = Next [j])
    {
	/* start a DFS at child node j */
	k = dfs_recursive (j, k, Post, Head, Next) ;
    }
    Post [k++] = p ;	/* order node p as the kth node */
    Head [p] = EMPTY_PARENT ;	/* link list p no longer needed */
    return (k) ;	/* the next node will be numbered k */
}

static int dfs_iterative		/* return the new value of k */
(
    int p,		/* start the DFS at a root node p */
    int k,		/* start the node numbering at k */
    int Post [ ],	/* Post ordering, modified on output */
    int Head [ ],	/* Head [p] = youngest child of p; EMPTY on output */
    int Next [ ],	/* Next [j] = sibling of j; unmodified */
    int Pstack [ ]	/* workspace of size n, undefined on input or output */
)
{
  int j, phead ;

  /* put the root node on the stack */
  Pstack [0] = p ;
  phead = 0 ;
  /* while the stack is not empty, do: */
  while (phead >= 0)
  {
    /* grab the node p from top of the stack and get its youngest child j */
    p = Pstack [phead] ;
    j = Head [p] ;
    if (j == EMPTY_PARENT)
    {
	    /* all children of p ordered.  remove p from stack and order it */
	    phead-- ;
	    Post [k++] = p ;	/* order node p as the kth node */
    }
    else
    {
	    /* leave p on the stack.  Start a DFS at child node j by putting
	     * j on the stack and removing j from the list of children of p. */
	    Head [p] = Next [j] ;
	    Pstack [++phead] = j ;
    }
  }
  return (k) ;	/* the next node will be numbered k */
}

static void tree_postorder_perm3(int* tree, int size, int* postorder, int iterative)
{
  int j,p;
  int* Next = spalloc(sizeof(int)*size);
  int* Head = spalloc(sizeof(int)*size);
  int* Pstack = spalloc(sizeof(int)*size);

  for (j = 0; j < size; ++ j) Next[j] = Head[j] = EMPTY_PARENT;
  	/* in reverse order so children are in ascending order in each list */
	for (j = size-1 ; j >= 0 ; j--)
	{
	    p = tree [j] ;
	    if (p >= 0 && p < size)
	    {
        /* add j to the list of children for node p */
        Next [j] = Head [p] ;
        Head [p] = j ;
	    }
	}
	/* Head [p] = j if j is the youngest (least-numbered) child of p */
	/* Next [j1] = j2 if j2 is the next-oldest sibling of j1 */

  int k = 0;
  for ( j = 0 ; j < size; ++ j)
    if ( tree[j] == EMPTY_PARENT )
    {
        if (iterative)
            dfs_iterative(j, k, postorder, Head, Next, Pstack);
        else
            dfs_recursive(j, k, postorder, Head, Next);
    }

  
}


static void tree_postorder_perm2(int* tree, int size, int* postorder)
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

  int i;
  struct postorder_perm_helper helper;
  helper.k = 0;
  helper.post = postorder;
  for ( i = 0 ; i < size; ++ i)
    if ( tree[i] == EMPTY_PARENT )
      tree_dfs_postorder(tree, size, i, tree_postorder_perm_helper, &helper);
}

/* depth-first search and postorder of a tree rooted at node j */
int cs_tdfs (int j, int k, int *head, const int *next, int *post, int *stack)
{
    int i, p, top = 0 ;
    if (!head || !next || !post || !stack) return (-1) ;    /* check inputs */
    stack [0] = j ;                 /* place j on the stack */
    while (top >= 0)                /* while (stack is not empty) */
    {
        p = stack [top] ;           /* p = top of stack */
        i = head [p] ;              /* i = youngest child of p */
        if (i == -1)
        {
            top-- ;                 /* p has no unordered children left */
            post [k++] = p ;        /* node p is the kth postordered node */
        }
        else
        {
            head [p] = next [i] ;   /* remove i from children of p */
            stack [++top] = i ;     /* start dfs on child node i */
        }
    }
    return (k) ;
}

/* post order a forest */
static void tree_postorder_perm4(int* tree, int size, int* postorder)
{
    int j, k = 0, *w, *head, *next, *stack ;
    w = spalloc (3*size*sizeof (int)) ;                 /* get workspace */
    head = w ; next = w + size ; stack = w + 2*size ;
    for (j = 0 ; j < size ; j++) head [j] = -1 ;           /* empty linked lists */
    for (j = size-1 ; j >= 0 ; j--)            /* traverse nodes in reverse order*/
    {
        if (tree [j] == -1) continue ;    /* j is a root */
        next [j] = head [tree [j]] ;      /* add j to list of its parent */
        head [tree [j]] = j ;
    }
    for (j = 0 ; j < size ; j++)
    {
        if (tree [j] != -1) continue ;    /* skip j if it is not a root */
        k = cs_tdfs (j, k, head, next, postorder, stack) ;
    }
    spfree(w);
}


static void tree_postorder_perm1(int* tree, int size, int* postorder)
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
  children = spalloc(sizeof(int_stack_ptr)*size);
  for ( i = 0; i < size; ++ i)
    children[i] = int_stack_alloc(1,2);
  /* now fill the stacks of children with the actual children
   * in reverse order to keep children in the stack in ascending order
   * i.e. ->[1 2 3] */
  for ( i = size -1 ; i >= 0; -- i)
    if (tree[i] != EMPTY_PARENT)          /* if i has parent (tree[i]) */
      int_stack_push(children[tree[i]],i);

  /* loop by all roots in the forest */
  for ( i = 0; i < size; ++ i)
  {
    /* root node in the forest. if forest is the 1 tree, only one
     * root will exist */
    if (tree[i] == EMPTY_PARENT)
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
  for ( i = 0; i < size; ++ i)
    int_stack_free(children[i]);
  int_stack_free(stack);
  spfree(children);
}


void tree_postorder_perm(int* tree, int size, int* postorder)
{
    /*tree_postorder_perm3(tree, size, postorder, 1);*/
    tree_postorder_perm4(tree, size, postorder);
}

void tree_node_levels(int* tree, int size, int* level)
{
  int i,j;
  for ( i = 0; i < size; ++ i)
  {
    level[i] = 0;
    /* calculate level of i-th node */
    j = tree[i];
    for ( ; j != EMPTY_PARENT; j = tree[j])
      level[i] ++;
  }
}

void tree_first_descendant(int* tree, int size, int* postorder, int* first)
{
  /* Algorithm:
   * in postordered tree every node k has proper descendants numbered
   * from k-d to k-1.
   * Therefore we take every 'postordered' node and traverse up the tree
   * marking setting all nodes in first array(ancestors) to have the
   * first descendant value as the node
   */
  int i,j,k;

  for ( i = 0; i < size; ++ i)
    first[i] = EMPTY_PARENT;
  for ( i = 0; i < size; ++ i)
  {
    k = postorder[i];        /* k is the node i in postordered tree */
    /* printf("node %d, postordered: %d\n",k+1,i+1); */
    /* traverse up the tree */
    for ( j = k; j != EMPTY_PARENT && first[j] == EMPTY_PARENT; j = tree[j]) 
    {
      first[j] = i;
      /* printf("traverse up to node %d\n",j+1); */
    }
  }
}


int tree_find(int* tree, int size, int value)
{
  int i;
  int temp;
  for (i = 0; i < size; ++ i)
  {
    temp = tree[value];
    if ( temp == EMPTY_PARENT)
      return value;
    value = temp;
  }
  return value;
}

int tree_is_etree(int* tree, int size)
{
  int result = 1, i = 0;
  for ( ; i < size; ++ i)
    if (tree[i] <= i && tree[i] != EMPTY_PARENT)
    {
      result = 0;
      break;
    }
  return result;
}

void tree_dot_printf(int* tree, int size)
{
  int i = 0;
  if (tree)
  {
    printf("digraph etree {\n  rankdir = BT;\n");
    printf("  d2tgraphstyle=\"scale=0.4\"\n  node [shape=circle];\n");
    for (; i < size; ++ i)
      if (tree[i] != EMPTY_PARENT)
        printf("  %d->%d\n",i+1,tree[i]+1);
    printf("}\n");
  }
}
