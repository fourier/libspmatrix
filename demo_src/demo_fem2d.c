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
#include <string.h>

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
  for (i = 0; i < M+1; ++ i)    /* row */
  {
    for (j = 0; j < N+1; ++ j)  /* column */
    {
      g->points[k].x = x + j*dx;
      g->points[k].y = y + i*dy;
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
      g->triangles[k][0] = base;
      g->triangles[k][1] = base + (N+1) + 1;
      g->triangles[k][2] = base + (N+1);
      k++;
      g->triangles[k][0] = base;
      g->triangles[k][1] = base + 1;
      g->triangles[k][2] = base + (N+1) + 1;
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

static
double det3x3(double (*m)[3])
{
  double result;
  result = m[0][0]*(m[1][1]*m[2][2]-m[1][2]*m[2][1]) - 
    m[0][1]*(m[1][0]*m[2][2]-m[1][2]*m[2][0]) + 
    m[0][2]*(m[1][0]*m[2][1]-m[1][1]*m[2][0]);
  return result;  
}
  

double element_size(geometry_2d* g, int element_no)
{
  double el[3][3];
  el[0][2] = 1;
  el[1][2] = 1;
  el[2][2] = 1;

  el[0][0] = g->points[g->triangles[element_no][0]].x;
  el[0][1] = g->points[g->triangles[element_no][0]].y;
  el[1][0] = g->points[g->triangles[element_no][1]].x;
  el[1][1] = g->points[g->triangles[element_no][1]].y;
  el[2][0] = g->points[g->triangles[element_no][2]].x;
  el[2][1] = g->points[g->triangles[element_no][2]].y;
  return det3x3(el)/2.0;
}

double local(geometry_2d* g, int element_no,int l, double x, double y)
{
  double a = 0,b = 0,c = 0;
  const double S = element_size(g,element_no);
  const double x1 = g->points[g->triangles[element_no][0]].x;
  const double y1 = g->points[g->triangles[element_no][0]].y;
  const double x2 = g->points[g->triangles[element_no][1]].x;
  const double y2 = g->points[g->triangles[element_no][1]].y;
  const double x3 = g->points[g->triangles[element_no][2]].x;
  const double y3 = g->points[g->triangles[element_no][2]].y;

  switch (l)
  {
  case 0:
    a = x2*y3-x3*y2;
    b = y2-y3;
    c = x3-x2;
    break;
  case 1:
    a = x3*y1-x1*y3;
    b = y3-y1;
    c = x1-x3;
    break;
  case 2:
    a = x1*y2-x2*y1;
    b = y1-y2;
    c = x2-x1;
    break;
  default:
    break;
  }
  return 0.5*(a+b*x+c*y)/S;
}

void create_b_matrix(geometry_2d* g, int element_no, dense_mtx* B)
{
  const double x1 = g->points[g->triangles[element_no][0]].x;
  const double y1 = g->points[g->triangles[element_no][0]].y;
  const double x2 = g->points[g->triangles[element_no][1]].x;
  const double y2 = g->points[g->triangles[element_no][1]].y;
  const double x3 = g->points[g->triangles[element_no][2]].x;
  const double y3 = g->points[g->triangles[element_no][2]].y;
  int i,j;
  double x21,x13,x32,y12,y31,y23;
  double A;
  dense_mtx_init(B,3,6);

  /* get derivatives */
  y23 = y2-y3;
	y31 = y3-y1;
	y12 = y1-y2;
	x32 = x3-x2;
	x13 = x1-x3;
	x21 = x2-x1;

  /* fill the B matrix ... */
  /* first row of B matrix */
  B->A[0][0] = y23; B->A[0][2] = y31; B->A[0][4] = y12;
	/* second row of B matrix */
	B->A[1][1] = x32; B->A[1][3] = x13; B->A[1][5] = x21;
	/* third row of B matrix */
	B->A[2][0] = x32; B->A[2][1] = y23;
  B->A[2][2] = x13; B->A[2][3] = y31;
  B->A[2][4] = x21; B->A[2][5] = y12;
  /* ... and divide by 2*A */
  A = element_size(g,element_no)*2;
  for (i = 0; i < B->rows; ++ i)
    for (j = 0; j < B->cols; ++ j)
      B->A[i][j] /= A;
}


static dense_mtx* elasticity_matrix()
{
  static dense_mtx* result = 0;
  static dense_mtx A;
  /*
   * Elasticity matrix:
   * D * mul
   * mul = E/(1-nu^2)
   * D = [ 1 nu        0 ]      
   *     [ nu 1        0 ]
   *     [ 0  0 (1-nu)/2 ]
   */
  const double E = 1e9;
  const double nu = 0.3;
  const double mul = E/(1-nu*nu);
  if (!result)
  {
    dense_mtx_init(&A,3,3);
    A.A[0][0] = mul;
    A.A[1][1] = mul;
    A.A[2][2] = 0.5*(1-nu)*mul;
    A.A[0][1] = nu*mul;
    A.A[1][0] = nu*mul;

    result = &A;
  }
  return result;
}

void create_local_mtx(geometry_2d* g, int element_no, dense_mtx* K)
{
  dense_mtx tmp;
  dense_mtx B;
  double A = element_size(g,element_no);
  int i,j;
  create_b_matrix(g,element_no,&B);
  /* K = B'*A*B */
  /* 1. tmp = B'*A */
  dense_mtx_mul_at_b(&B,elasticity_matrix(),&tmp);
  /* 2. K = tmp*B */
  dense_mtx_mul_a_b(&tmp,&B,K);

  for (i = 0; i < K->rows; ++ i)
    for (j = 0; j < K->cols; ++ j)
      K->A[i][j] *= A;
  
  dense_mtx_free(&tmp);
  dense_mtx_free(&B);
}


void dense_mtx_init(dense_mtx* m, int rows, int cols)
{
  int i;
  memset(m,0,sizeof(dense_mtx));
  m->rows = rows;
  m->cols = cols;
  m->A = calloc(rows,sizeof(double*));
  for (i = 0; i < rows; ++ i)
    m->A[i] = calloc(cols,sizeof(double));
}

void dense_mtx_eye_init(dense_mtx* m, int N)
{
  int i;
  dense_mtx_init(m,N,N);
  for (i = 0; i < N; ++ i)
    m->A[i][i] = 1;
}

void dense_mtx_free(dense_mtx* m)
{
  int i;
  for (i = 0; i < m->rows; ++ i)
    free(m->A[i]);
  free(m->A);
  memset(m,0,sizeof(dense_mtx));
}

void dense_mtx_mul_a_b(dense_mtx* A, dense_mtx* B, dense_mtx* M)
{
  int i,j,k;
  dense_mtx_init(M,A->rows,B->cols);
  for (i = 0; i < A->rows; ++ i)
    for (j = 0; j < B->cols; ++ j)
      for (k = 0; k < A->cols; ++ k)
        M->A[i][j] += A->A[i][k]*B->A[k][j];
}


void dense_mtx_mul_at_b(dense_mtx* A, dense_mtx* B, dense_mtx* M)
{
  int i,j,k;
  dense_mtx_init(M,A->cols,B->cols);
  for (i = 0; i < A->cols; ++ i)
    for (j = 0; j < B->cols; ++ j)
      for (k = 0; k < A->rows; ++ k)
        M->A[i][j] += A->A[k][i]*B->A[k][j];
}


void dense_mtx_printf(dense_mtx* A)
{
  int i,j;
  printf("%dx%d:\n",A->rows,A->cols);
  for (i = 0; i < A->rows; ++ i)
  {
    for (j = 0; j < A->cols; ++j )
      printf("%f ",A->A[i][j]);
    printf("\n");
  }
}
