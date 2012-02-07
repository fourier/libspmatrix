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
#include <stdlib.h>
#include "demo_fem2d.h"

/*
  1------2------3------4------5
  |\     |\     |\     |\     |
  |  \   |  \   |  \   |  \   |
  |    \ |    \ |    \ |    \ |
  6------7------8------9-----10
  |\     |\     |\     |\     |
  |  \   |  \   |  \   |  \   |
  |    \ |    \ |    \ |    \ |
  11----12-----13-----14-----15
  |\     |\     |\     |\     |
  |  \   |  \   |  \   |  \   |
  |    \ |    \ |    \ |    \ |
  16----17-----16-----19-----20

  triangles:
  [1,6,7],   [1,7,2]
  [2,7,9],   [2,8,3]
  [3,8,9],   [3,9,4]
  [4,9,10],  [4,10,5]
  ...

*/


void generate_brick_mesh(int N,int M,
                         double x, double y,
                         double dx, double dy,
                         geometry_2d* g,
                         prescr_boundary_2d* b)
{
  int i,j,k,base;
  g->horizontal_blocks_count = N;
  g->vertical_blocks_count = M;

  g->points_count = (N+1)*(M+1);
  g->triangles_count = N*M*2;

  g->points = calloc(g->points_count,sizeof(point_2d));
  g->triangles = calloc(g->triangles_count, sizeof(triangle_3p));

  /* create points */
  k = 0;
  for (i = 0; i < N+1; ++ i)
  {
    for (j = 0; j < M+1; ++ j)
    {
      g->points[k].x = x + i*dx;
      g->points[k].y = y + j*dy;
      k++;
    }
  }
  /* create triangels */
  k = 0;
  for (i = 0; i < N; ++ i)
  {
    for (j = 0; j < M; ++ j)
    {
      base = j*(N+1)+i;
      g->triangles[k].p1 = base;
      g->triangles[k].p2 = base + (N+1);
      g->triangles[k].p3 = base + (N+1) + 1;
      k++;
      g->triangles[k].p1 = base;
      g->triangles[k].p2 = base + (N+1) + 1;
      g->triangles[k].p3 = base + 1;
      k++;
    }
  }
  /* create boundary conditions in form of prescribed
   * displacements 
   * considering offset = 1% of the length, it is normal for
   * small deformations */

  /* leftmost points - fixed, rightmost - with prescribed offset */
  b->prescribed_count = (M+1)*2; 
  b->points = calloc(b->prescribed_count, sizeof(prescribed_point));
  for ( i = 0; i < M+1; ++ i)
  {
    k = i*(N+1);
    b->points[i].type = FIXED_XY;
    b->points[i].point_index = k;
    b->points[i].point.x = 0;
    b->points[i].point.y = 0;

    b->points[i+M+1].type = FIXED_X;
    b->points[i+M+1].point_index = N+k;
    b->points[i+M+1].point.x = (dx*N-x)/100.0;
    b->points[i+M+1].point.y = 0;
  }
  
}
