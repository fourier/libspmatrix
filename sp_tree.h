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

typedef int (*traverse_func_t)(int n, void*);

/*
 * Depth First Search of the tree(or forest), applying the function
 * func of type traverse_func_t with the argument arg to every
 * traversed node
 */
void tree_dfs(int* tree, int size, traverse_func_t func, void* arg);

/*
 * Breadth First Search of the tree(or forest), applying the function
 * func of type traverse_func_t with the argument arg to every
 * traversed node
 */
void tree_bfs(int* tree, int size, traverse_func_t func, void* arg);

/*
 * Postorder the tree. A postordered tree is the tree where for every
 * node k all its children <= k-1
 * This function returns the permutation required
 * to apply to the source tree to get the postordered tree.
 * For example postorder[k] = i means what node i of the source tree is
 * the node k of postordered tree.
 * Example(1-based):
 * given tree [4,3,-1,3]:
 *
 *    3
 *   / \
 *  2   4
 *      |
 *      1
 *
 * postordered tree([4,3,4,-1]):
 *
 *    4
 *   / \
 *  1   3
 *      |
 *      2
 *
 * keepting the order: if c1 < c2 < ... < c_n - children of some node, to
 * postorder(c1) < postorder(c2) < ... < postorder(c_n).
 * Postorder returned by this function is [2,1,4,3].
 * TODO: how to get the postordered tree in this array form ?
 */
void tree_postorder_perm(int* tree, int size, int* postorder);

/*
 * Calculates the level of every node and store it in appropriate position
 * in the level array
 * Example: given the tree
 *    4
 *   / \
 *  1   3
 *      |
 *      2
 *
 * Levels arrays shall look like: [1;2;1;0]
 */
void tree_node_levels(int* tree, int size, int* level);

/*
 * Calculates first descendant nodes for every node in array first
 * First descendant of the node j is the smallest postordered node
 * of any descendants of j.
 * See T.A.Davis, Direct Methods for Sparse Linear Systems(2006) p.47
 */
void tree_first_descendant(int* tree, int size, int* postorder, int* first);

/*
 * Find the root of the node 'value' in the tree
 * represented by array 'tree' of size 'size'
 * there the position number is the tree node, and
 * the value is the parent value. 
 * -1 means that the node has no parents
 * Example: tree = [5,2,-1,5,-1,6,-1,-1]
 * tree_fine shall return 6 for 0,4,5; 2 for 1, and
 * the same values for others (like 2 for 2 etc.)
 * This algorithm is used for calculation of the
 * Elimination tree
 * */
int tree_find(int* tree, int size, int value);

#endif /* _SP_TREE_H_ */
