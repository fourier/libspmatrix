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

#ifndef _DEMO_FEM2D_H_
#define _DEMO_FEM2D_H_

typedef enum
{
  FIXED_X = 0,
  FIXED_Y,
  FIXED_XY
} prescribed_node_type;

typedef struct
{
  double x;
  double y;
} point_2d;

typedef int triangle_3p[3];

typedef struct
{
  int horizontal_blocks_count;    /* number of horizontal blocks */
  int vertical_blocks_count;      /* number of vertical blocks */
  int triangles_count;            /* number of triangles */
  int points_count;               /* number of points */
  point_2d* points;
  triangle_3p* triangles;
} geometry_2d;

typedef struct
{
  int point_index;
  point_2d point;
  prescribed_node_type type;
} prescribed_point;

typedef struct
{
  int prescribed_count;
  prescribed_point* points;
} prescr_boundary_2d;

typedef struct
{
  int rows;
  int cols;
  double **A;
} dense_mtx;


void generate_brick_mesh(int N,int M,
                         double x, double y,
                         double dx, double dy,
                         geometry_2d* g,
                         prescr_boundary_2d* b);

/* size of the triangular 2d - element */
double element_size(geometry_2d* g, int element_no);

/*
 * Construct matrix of derivatives of shape functions - B
 * for given element and in given point
 */
void create_b_matrix(geometry_2d* g, int element_no, dense_mtx* B);

/*
 * Construct local stiffness matrix for element element_no
 */
void create_local_mtx(geometry_2d* g, int element_no, dense_mtx* K);

/*
 * Operations on dense matricies
 */

/* initialize dense matrix rows x cols filled with zeros */
void dense_mtx_init(dense_mtx* m, int rows, int cols);
/* free the memory for the dense matrix */
void dense_mtx_free(dense_mtx* m);
/* multiply dense matrix M = A*B, initializing the M matrix */
void dense_mtx_mul_a_b(dense_mtx* A, dense_mtx* B, dense_mtx* M);
/* multiply dense matrix M = A'*B, initializing the M matrix */
void dense_mtx_mul_at_b(dense_mtx* A, dense_mtx* B, dense_mtx* M);

#endif /* _DEMO_FEM2D_H_ */
