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
#include <math.h>
#include <string.h>

#include "demo_fem2d.h"
#include "sp_matrix.h"
#include "sp_file.h"

static void apply_bc(sp_matrix_ptr m, int idx)
{
  int i,msize;
  double *pvalue;
  msize = m->rows_count;
  for (i = 0; i < msize; ++ i)
  {
    pvalue = sp_matrix_element_ptr(m,idx,i);
    if (pvalue)
      *pvalue = 0;
    pvalue = sp_matrix_element_ptr(m,i,idx);
    if (pvalue)
      *pvalue = 0;
  }
  pvalue = sp_matrix_element_ptr(m,idx,idx);
  if (pvalue)
    *pvalue = 1;
}

static void usage(const char* progname)
{
  printf("Usage: %s N M export_filename.mtx\n",progname);
  printf("where:\n");
  printf(" N - number of finite element \"rows\"(3 on a picture below)\n");
  printf(" M - number of finite element \"columns\"(4 on a picture below)\n");
  printf(" export_filename.mtx - name of the sparse matrix for output\n");
  printf("This demo generates the positive-definite sparse matrix\n");
  printf("from the following plane stress task:\n\n");
  printf(
    " \\ \n"
    " -1-----2-----3-----4-----5\n"
    " /|\\    |\\    |\\    |\\    |\n"
    "  |  \\  |  \\  |  \\  |  \\  |\n"
    " \\|    \\|    \\|    \\|    \\|\n"
    " -6-----7-----8-----9----10\n"
    " /|\\    |\\    |\\    |\\    |\n"
    "  |  \\  |  \\  |  \\  |  \\  |\n"
    " \\|    \\|    \\|    \\|    \\|\n"
    " -11---12----13----14----15\n"
    " /|\\    |\\    |\\    |\\    |\n"
    "  |  \\  |  \\  |  \\  |  \\  |\n"
    " \\|    \\|    \\|    \\|    \\|\n"
    " -16---17----18----19----20\n"
    " /\n\n");
  printf("Left side of the plate is fixed\n");
  printf("Matrial constants E = 1e9, nu = 0.3\n");
}


int main(int argc, const char* argv[])
{
  int N = 4;           /* number of vertical blocks */
  int M = 3;           /* number of horizontal blocks */
  int i,j,k,I,J;
  int msize;
  const double x = 1.0, y=1.0;                 /* upper-left point */
  const double dx = 1.0,dy = 1.0;              /* size of the block */
  geometry_2d g;
  prescr_boundary_2d b;
  dense_mtx K;
  sp_matrix m;
  sp_matrix_yale yale;
  const char* ptr = 0;
  if (argc < 4)
  {
    usage(argv[0]);
    return 1;
  }
  M = atoi(argv[1]);
  N = atoi(argv[2]);
  if (M <= 0 || N <= 0)
  {
    usage(argv[0]);
    return 1;
  }
  ptr = argv[3];
    
  generate_brick_mesh(N,M,x,y,dx,dy,&g,&b);
#if 0
  for (i = 0; i < g.triangles_count; ++ i)
    printf("[%d, %d, %d]\n",
           g.triangles[i][0]+1,
           g.triangles[i][1]+1,
           g.triangles[i][2]+1);
  printf("boundary:\n");
  for (i = 0; i < b.prescribed_count; ++ i)
    printf("node %d, type %s, dx = %.2f, dy = %.2f\n",
           b.points[i].point_index + 1,
           b.points[i].type == FIXED_XY ? "XY" :
           (b.points[i].type == FIXED_X ? "X" : "Y"),
           b.points[i].point.x,
           b.points[i].point.y);
#endif
  msize = g.points_count*2;
  /* initialize global matrix */
  sp_matrix_init(&m,msize,msize,
                 sqrt(msize/2.),CCS);
  for (k = 0; k < g.triangles_count; ++ k)
  {
    /* create local stiffness */
    create_local_mtx(&g,k,&K);
    /* global dof */
    for (i = 0; i < 3; ++ i)
      for (j = 0; j < 3; ++ j)
      {
        /* i: l*2 */
        I = g.triangles[k][i]*2;
        J = g.triangles[k][j]*2;
        sp_matrix_element_add(&m,I,J,K.A[i*2][j*2]);
        sp_matrix_element_add(&m,I+1,J,K.A[i*2+1][j*2]);
        sp_matrix_element_add(&m,I,J+1,K.A[i*2][j*2+1]);
        sp_matrix_element_add(&m,I+1,J+1,K.A[i*2+1][j*2+1]);
      }

    /* distribute in global matrix */
    dense_mtx_free(&K);
  }
  /* apply bc */
  for (i = 0; i < b.prescribed_count; ++ i)
  {
    switch (b.points[i].type)
    {
    case FIXED_X:
      /* apply_bc(&m,b.points[i].point_index*2); */
      break;
    case FIXED_Y:
      /* apply_bc(&m,b.points[i].point_index*2+1); */
      break;
    case FIXED_XY:
      apply_bc(&m,b.points[i].point_index*2);
      apply_bc(&m,b.points[i].point_index*2+1);
      break;
    default:
      break;
    };
  }
  sp_matrix_yale_init(&yale,&m);
  sp_matrix_free(&m);
  
  if (sp_matrix_yale_save_file(&yale,ptr))
  {
    printf("Sparse matrix saved to %s\n", ptr);
    sp_matrix_yale_printf2(&yale);
  }
  
  sp_matrix_yale_free(&yale);

  return 0;
}

