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

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <memory.h>

#include "sp_matrix.h"
#include "sp_utils.h"
#include "sp_file.h"
#include "sp_cont.h"
#include "sp_tree.h"

#include "logger.h"

static int test_sp_matrix()
{
  int result = 1;
  sp_matrix mtx,mtx2,mtx3;
  double b[] = {1, 2, 3, 4, 3, 2, 1};
  double expected[] = {25, 34, 40, 45, 42, 16, 23};
  int size = sizeof(expected)/sizeof(double);
  int i;
  double x[] = {0,0,0,0,0,0,0};
  /*
   * Sparse matrix
   * 9  0  0  3  1  0  1
   * 0  11 2  1  0  0  2
   * 0  1  10 2  0  0  0
   * 0  0  2  9  1  0  0
   * 1  0  0  1  12 0  1
   * 0  0  0  0  0  8  0
   * 2  2  0  0  3  0  8
   */

  sp_matrix_init(&mtx,7,7,5,CRS);

  MTX(&mtx,0,0,9);MTX(&mtx,0,3,3);MTX(&mtx,0,4,1);MTX(&mtx,0,6,1);
  MTX(&mtx,1,1,11);MTX(&mtx,1,2,2);MTX(&mtx,1,3,1);MTX(&mtx,1,6,2);
  MTX(&mtx,2,1,1);MTX(&mtx,2,2,10);MTX(&mtx,2,3,2);
  MTX(&mtx,3,2,2);MTX(&mtx,3,3,9);MTX(&mtx,3,4,1);
  MTX(&mtx,4,0,1);MTX(&mtx,4,3,1);MTX(&mtx,4,4,12);MTX(&mtx,4,6,1);
  MTX(&mtx,5,5,8);
  MTX(&mtx,6,0,2);MTX(&mtx,6,1,2);MTX(&mtx,6,4,3);MTX(&mtx,6,6,8);

  sp_matrix_reorder(&mtx);
  
  /* 1st test - matrix-vector multiplication */
  sp_matrix_mv(&mtx,b,x);
  for (i = 0; i < size; ++ i)
    result &= EQL(x[i],expected[i]);
  
  /* 2nd test - conversion btw different storage types */
  if (result)
  {
    sp_matrix_convert(&mtx,&mtx2,CCS);
    sp_matrix_mv(&mtx2,b,x);
    for (i = 0; i < size; ++ i)
    {
      result &= EQL(x[i],expected[i]);
    }
  }
  if (result)
  {
    sp_matrix_convert(&mtx2,&mtx3,CRS);
    sp_matrix_mv(&mtx3,b,x);
    for (i = 0; i < size; ++ i)
    {
      result &= EQL(x[i],expected[i]);
    }
    sp_matrix_free(&mtx2);
    sp_matrix_free(&mtx3);
  }
  
  sp_matrix_free(&mtx);
  printf("test_sp_matrix result:\t*%s*\n",result ? "pass" : "fail");
  return result;
}

static int test_yale()
{
  int result = 0;
  int i;
  sp_matrix mtx;
  sp_matrix_yale yale;
  /* 1. CRS */
  /* matrix:
   *
   * [ 1 2 0 0 ]
   * [ 0 3 9 0 ]
   * [ 0 1 4 0 ]
   */
  /* expected result for CRS: */
  int offsets1[] = {0,2,4,6};
  int indicies1[] = {0,1,1,2,1,2};
  double values1[] = {1,2,3,9,1,4};

  double x[] = {1,2,3,4};
  double y[3] = {0};
  double b[] = {5,33,14};
  
  /* expected result for CCS: */
  int offsets2[] = {0,3,5,7,9,11};
  int indicies2[] = {0,2,4,0,3,1,4,0,3,1,4};
  double values2[] = {1,2,5,-3,4,-2,-5,-1,-4,3,6};

  sp_matrix_init(&mtx,3,4,2,CRS);
  MTX(&mtx,0,0,1);MTX(&mtx,0,1,2);
  MTX(&mtx,1,1,3);MTX(&mtx,1,2,9);
  MTX(&mtx,2,1,1);MTX(&mtx,2,2,4);
  
  sp_matrix_yale_init(&yale,&mtx);
  result = memcmp(yale.offsets,offsets1,4*sizeof(int)) == 0;
  result &= memcmp(yale.indicies,indicies1,6*sizeof(int)) == 0;
  result &= memcmp(yale.values,values1,6*sizeof(double)) == 0;

  sp_matrix_yale_mv(&yale,x,y);
  for (i = 0; i < 3; ++ i)
    result &= EQL(y[i],b[i]);  

  sp_matrix_free(&mtx);
  sp_matrix_yale_free(&yale);
  /* 2. CCS */
  /* matrix: */
  /*
   * 1. -3.  0. -1.  0.
   * 0.  0. -2.  0.  3.
   * 2.  0.  0.  0.  0.
   * 0.  4.  0. -4.  0.
   * 5.  0. -5.  0.  6.
   */
  sp_matrix_init(&mtx,5,5,3,CCS);
  /* 1. -3.  0. -1.  0. */
  MTX(&mtx,0,0,1);MTX(&mtx,0,1,-3);MTX(&mtx,0,3,-1);
  /* 0.  0. -2.  0.  3. */
  MTX(&mtx,1,2,-2);MTX(&mtx,1,4,3);
  /* 2.  0.  0.  0.  0. */
  MTX(&mtx,2,0,2);
  /* 0.  4.  0. -4.  0. */
  MTX(&mtx,3,1,4);MTX(&mtx,3,3,-4);
  /* 5.  0. -5.  0.  6. */
  MTX(&mtx,4,0,5);MTX(&mtx,4,2,-5);MTX(&mtx,4,4,6);
  
  sp_matrix_yale_init(&yale,&mtx);
  result = memcmp(yale.offsets,offsets2,6*sizeof(int)) == 0;
  result &= memcmp(yale.indicies,indicies2,11*sizeof(int)) == 0;
  result &= memcmp(yale.values,values2,11*sizeof(double)) == 0;
  
  sp_matrix_free(&mtx);
  sp_matrix_yale_free(&yale);
  
  printf("test_yale:\t*%s*\n",result ? "pass" : "fail");
  return result;
}

static int test_triangle_solver()
{
  int result = 1;
  int i;
  sp_matrix mtx,mtx2;
  double x[5] = {0};
  double x_expected[] = {1,2,-3,5,-7};
  double b[] = {-1, 5, -10, 40, -71};

  /*
   * |-1  0  0  0  0 |   | 1 |   |-1 |
   * | 1  2  0  0  0 |   | 2 |   | 5 |
   * |-1  0  3  0  0 | x |-3 | = |-10|
   * | 0  5  0  6  0 |   | 5 |   | 40|
   * | 0  0 -2  0 11 |   |-7 |   |-71|
   */
  sp_matrix_init(&mtx,5,5,3,CRS);
  MTX(&mtx,0,0,-1);
  MTX(&mtx,1,0,1);MTX(&mtx,1,1,2);
  MTX(&mtx,2,0,-1);MTX(&mtx,2,2,3);
  MTX(&mtx,3,1,5);MTX(&mtx,3,3,6);
  MTX(&mtx,4,2,-2);MTX(&mtx,4,4,11);

  sp_matrix_lower_solve(&mtx,5,b,x);
  for (i = 0; i < 5; ++ i)
    result &= EQL(x_expected[i],x[i]);
  
  if (result)
  {
    sp_matrix_convert(&mtx,&mtx2,CCS);
    memset(x,0,sizeof(double)*5);
    sp_matrix_lower_solve(&mtx,5,b,x);
    for (i = 0; i < 5; ++ i)
      result &= EQL(x_expected[i],x[i]);
    sp_matrix_free(&mtx2);
  }
  
  sp_matrix_free(&mtx);
  
  printf("test_triangle_solver result:\t*%s*\n",result ? "pass" : "fail");
  return result;
}

static int test_cg_solver()
{
  int result = 1;
  sp_matrix mtx;
  sp_matrix_yale yale;
  double v[3] = {0}, x[3] = {0}, z[3] = {0};
  int max_iter = 20000;
  const double desired_tolearance = 1e-15;
  double tolerance = desired_tolearance;

  /* matrix solver test  */

  /* Test 1: */
  /*
   * | 1 0 -2 |   | 1 |   |-5 |
   * | 0 1  0 | x | 2 | = | 2 | 
   * |-2 0  5 |   | 3 |   |13 |
   */
  memset(x,0,3);
  v[0] = -5;
  v[1] = 2;
  v[2] = 13;
  sp_matrix_init(&mtx,3,3,2,CRS);
  
  MTX(&mtx,0,0,1);MTX(&mtx,0,2,-2);
  MTX(&mtx,1,1,1);
  MTX(&mtx,2,0,-2);MTX(&mtx,2,2,5);

  sp_matrix_reorder(&mtx);
  sp_matrix_yale_init(&yale,&mtx);
  sp_matrix_yale_solve_cg(&yale,v,v,&max_iter,&tolerance,x);

  /* check for convergence */
  sp_matrix_mv(&mtx,x,z);
  
  tolerance = sqrt(pow(z[0]-v[0],2)+pow(z[1]-v[1],2)+pow(z[2]-v[2],2));
  result = tolerance < desired_tolearance*10;
  
  sp_matrix_free(&mtx);
  sp_matrix_yale_free(&yale);
  printf("test_cg_solver result:\t*%s*\n",result ? "pass" : "fail");
  return result;
}

static int test_ilu()
{
  int result = 1;
  sp_matrix mtx;
  sp_matrix_skyline m;
  sp_matrix_skyline_ilu ILU;
  double x_exact[] = {1,2,3,0,3,2,1};
  double x[7];
  double b[7];
  int i;
  /* test data for ILU decomposition test */
  double lu_diag_expected[] = {9.000000,
                               11.000000,
                               9.818182,
                               7.888889,
                               11.823161,
                               8.000000,
                               7.205303};
  double lu_lowertr_expected[] = {0.090909,
                                  0.222222,
                                  0.090909,
                                  0.185185,
                                  0.111111,
                                  0.084507,
                                  0.222222,
                                  0.181818,
                                  0.234944};
  double lu_uppertr_expected[] = {2.000000,
                                  3.000000,
                                  1.000000,
                                  1.909091,
                                  1.000000,
                                  0.777778,
                                  1.000000,
                                  2.000000,
                                  0.888889};

  memset(b,0,sizeof(b));
  memset(x,0,sizeof(x));
  
  /* Sparse matrix from Balandin
   * 9  0  0  3  1  0  1
   * 0  11 2  1  0  0  2
   * 0  1  10 2  0  0  0
   * 2  1  2  9  1  0  0
   * 1  0  0  1  12 0  1
   * 0  0  0  0  0  8  0
   * 2  2  0  0  3  0  8
   *
   * Test for
   * 1) skyline format
   * 2) ILU decomposition
   * 3) LU - solvers for ILU decomposition
   */

  sp_matrix_init(&mtx,7,7,5,CRS);

  MTX(&mtx,0,0,9);MTX(&mtx,0,3,3);MTX(&mtx,0,4,1);MTX(&mtx,0,6,1);
  MTX(&mtx,1,1,11);MTX(&mtx,1,2,2);MTX(&mtx,1,3,1);MTX(&mtx,1,6,2);
  MTX(&mtx,2,1,1);MTX(&mtx,2,2,10);MTX(&mtx,2,3,2);
  MTX(&mtx,3,0,2);MTX(&mtx,3,1,1);MTX(&mtx,3,2,2);MTX(&mtx,3,3,9);
  MTX(&mtx,3,4,1);
  MTX(&mtx,4,0,1);MTX(&mtx,4,3,1);MTX(&mtx,4,4,12);MTX(&mtx,4,6,1);
  MTX(&mtx,5,5,8);
  MTX(&mtx,6,0,2);MTX(&mtx,6,1,2);MTX(&mtx,6,4,3);MTX(&mtx,6,6,8);

  sp_matrix_reorder(&mtx);
  sp_matrix_skyline_init(&m,&mtx);
  sp_matrix_skyline_ilu_copy_init(&ILU,&m);

  for (i = 0; i <  m.rows_count; ++ i)
    result &= fabs(ILU.ilu_diag[i] - lu_diag_expected[i]) < 1e-5;
  
  if (result)
  {
    for (i = 0; i <  m.tr_nonzeros; ++ i)
      result &= fabs(ILU.ilu_lowertr[i] - lu_lowertr_expected[i]) < 1e-5;
  }
  
  if (result)
  {       
    for (i = 0; i <  m.tr_nonzeros; ++ i)
      result &= fabs(ILU.ilu_uppertr[i] - lu_uppertr_expected[i]) < 1e-5;
  }

  /*
   * test for solving Lx=b
   */
  if (result)
  {
    /* prepare a right-part vector */
    sp_matrix_skyline_ilu_lower_mv(&ILU,x_exact,b);
    
    /* solve for x */
    sp_matrix_skyline_ilu_lower_solve(&ILU,b,x);
    /* test result */
    for ( i = 0; i < m.rows_count; ++ i)
      result &= EQL(x[i],x_exact[i]);
  }
  
  /*
   * test for solving Ux=b
   */
  if (result)
  {
    memset(b,0,sizeof(b));
    memset(x,0,sizeof(x));
    /* prepare a right-part vector */
    sp_matrix_skyline_ilu_upper_mv(&ILU,x_exact,b);
    
    /* solve for x */
    sp_matrix_skyline_ilu_upper_solve(&ILU,b,x);
    /* test result */
    for ( i = 0; i < m.rows_count; ++ i)
      result &= EQL(x[i],x_exact[i]);

    sp_matrix_free(&mtx);
    sp_matrix_skyline_ilu_free(&ILU);
  }
  
  printf("test_ilu result:\t*%s*\n",result ? "pass" : "fail");
  return result;
}

static int test_pcg_ilu_solver()
{
  int result = 1;
  sp_matrix mtx;
  sp_matrix_yale yale;
  sp_matrix_skyline_ilu ilu;
  double v[3] = {0}, x[3] = {0}, z[3] = {0};
  int max_iter = 20000;
  const double desired_tolearance = 1e-15;
  double tolerance = desired_tolearance;

  /* matrix solver test  */

  /* Test 1: */
  /*
   * | 1 0 -2 |   | 1 |   |-5 |
   * | 0 1  0 | x | 2 | = | 2 | 
   * |-2 0  5 |   | 3 |   |13 |
   */
  memset(x,0,3);
  v[0] = -5;
  v[1] = 2;
  v[2] = 13;
  sp_matrix_init(&mtx,3,3,2,CRS);
  
  MTX(&mtx,0,0,1);MTX(&mtx,0,2,-2);
  MTX(&mtx,1,1,1);
  MTX(&mtx,2,0,-2);MTX(&mtx,2,2,5);
  sp_matrix_reorder(&mtx);


  sp_matrix_yale_init(&yale,&mtx);
  sp_matrix_create_ilu(&mtx,&ilu);
  sp_matrix_yale_solve_pcg_ilu(&yale,&ilu,v,v,&max_iter,&tolerance,x);

  /* check for convergence */
  sp_matrix_mv(&mtx,x,z);

  tolerance = sqrt(pow(z[0]-v[0],2)+pow(z[1]-v[1],2)+pow(z[2]-v[2],2));
  result = tolerance < desired_tolearance*10;
  
  sp_matrix_skyline_ilu_free(&ilu);
  sp_matrix_free(&mtx);
  sp_matrix_yale_free(&yale);
  printf("test_pcg_ilu_solver result:\t*%s*\n",result ? "pass" : "fail");
  return result;
}

static int test_cholesky()
{
  int result = 0;
  /* initial matrix */
  /* {90, 6, 4, 46, 29, 0, 26}, */
  /* {6, 127, 34, 22, 7, 0, 38}, */
  /* {4, 34, 108, 40, 2, 0, 4}, */
  /* {46, 22, 40, 96, 24, 0, 6}, */
  /* {29, 7, 2, 24, 155, 0, 37}, */
  /* {0, 0, 0, 0, 0, 64, 0}, */
  /* {26, 38, 4, 6, 37, 0, 70} */
  sp_matrix mtx;
  sp_matrix_skyline m;

  /* expected decomposition */
  double cholesky_expected[7][7] = 
    {{9.48683, 0.632456, 0.421637, 4.84883, 3.05687, 0., 2.74064},
     {0.,11.2517, 2.99807, 1.68271, 0.450304, 0., 3.22323},
     {0., 0., 9.94152,3.31043, -0.0642691, 0., -0.685914},
     {0., 0., 0., 7.66149, 1.12678,0., -1.36292},
     {0., 0., 0., 0., 12.0075, 0., 2.38705},
     {0., 0., 0.,0., 0., 8., 0.},
     {0., 0., 0., 0., 0., 0., 6.6388}};

  /* fill initial matrix */
  sp_matrix_init(&mtx,7,7,5,CRS);

/* {90, 6, 4, 46, 29, 0, 26}, */
  MTX(&mtx,0,0,90);MTX(&mtx,0,1,6);MTX(&mtx,0,2,4);MTX(&mtx,0,3,46);
  MTX(&mtx,0,4,29);MTX(&mtx,0,6,26);
/* {6, 127, 34, 22, 7, 0, 38}, */
  MTX(&mtx,1,0,6);MTX(&mtx,1,1,127);MTX(&mtx,1,2,34);MTX(&mtx,1,3,22);
  MTX(&mtx,1,4,7);MTX(&mtx,1,6,38);
/* {4, 34, 108, 40, 2, 0, 4}, */
  MTX(&mtx,2,0,4);MTX(&mtx,2,1,34);MTX(&mtx,2,2,108);MTX(&mtx,2,3,40);
  MTX(&mtx,2,4,2);MTX(&mtx,2,6,4);
/* {46, 22, 40, 96, 24, 0, 6}, */
  MTX(&mtx,3,0,46);MTX(&mtx,3,1,22);MTX(&mtx,3,2,40);MTX(&mtx,3,3,96);
  MTX(&mtx,3,4,24);MTX(&mtx,3,6,6);
/* {29, 7, 2, 24, 155, 0, 37}, */
  MTX(&mtx,4,0,29);MTX(&mtx,4,1,7);MTX(&mtx,4,2,2);MTX(&mtx,4,3,24);
  MTX(&mtx,4,4,155);MTX(&mtx,4,6,37);
/* {0, 0, 0, 0, 0, 64, 0}, */
  MTX(&mtx,5,5,64);
/* {26, 38, 4, 6, 37, 0, 70} */
  MTX(&mtx,6,0,26);MTX(&mtx,6,1,38);MTX(&mtx,6,2,4);MTX(&mtx,6,3,6);
  MTX(&mtx,6,4,37);MTX(&mtx,6,6,70);

  /* prepare initial matrix for conversion to Skyline format */
  sp_matrix_reorder(&mtx);
  /* initialize skyline format from given CRS format */
  sp_matrix_skyline_init(&m,&mtx);

  /* clear matrix */
  sp_matrix_free(&mtx);
  sp_matrix_skyline_free(&m);
  
  printf("test_cholesky result:\t*%s*\n",result ? "pass" : "fail");
  return result;
}

static int test_load()
{
  int result = 0;
  sp_matrix_yale mtx;
  /* result = sp_matrix_load_file("bcsstk09.rsa",CCS); */
  /* result = sp_matrix_load_file("af23560.rua",CCS); */
  /* result = sp_matrix_load_file("kershaw_rua.hb",CCS); */
  result = sp_matrix_yale_load_file(&mtx, "5by5_rua.hb");
  if (result)
  {
    /* sp_matrix_save_file(result,"export.mtx"); */
    sp_matrix_yale_free(&mtx);
  }
  printf("test_load:\t*%s*\n",result ? "pass" : "fail");
  return result;
}

static int test_etree()
{
  int result = 0;
  int i;
  sp_matrix mtx;
  sp_matrix_yale yale;
  int etree_expected[] = {6,3,8,6,8,7,9,10,10,11,0};
  int* etree = 0;
  int* ereach = 0;
  /* fill initial matrix */
  sp_matrix_init(&mtx,11,11,5,CCS);

  /* Octave representation:
     m = [1,0,0,0,0,1,1,0,0,0,0;
     0,1,1,0,0,0,0,1,0,0,0;
     0,1,1,0,0,0,0,0,0,1,1;
     0,0,0,1,0,1,0,0,0,1,0;
     0,0,0,0,1,0,0,1,0,0,1;
     1,0,0,1,0,1,0,0,1,1,0;
     1,0,0,0,0,0,1,0,0,0,1;
     0,1,0,0,1,0,0,1,0,1,1;
     0,0,0,0,0,1,0,0,1,0,0;
     0,0,1,1,0,1,0,1,0,1,1;
     0,0,1,0,1,0,1,1,0,1,1];
  */
  
  /* 1 row: 1,0,0,0,0,1,1,0,0,0,0; */
  MTX(&mtx,0,0,1);MTX(&mtx,0,5,1);MTX(&mtx,0,6,1);
  /* 2 row: 0,1,1,0,0,0,0,1,0,0,0; */
  MTX(&mtx,1,1,1);MTX(&mtx,1,2,1);MTX(&mtx,1,7,1);
  /* 3 row: 0,1,1,0,0,0,0,0,0,1,1; */
  MTX(&mtx,2,1,1);MTX(&mtx,2,2,1);MTX(&mtx,2,9,1);MTX(&mtx,2,10,1);
  /* 4 row: 0,0,0,1,0,1,0,0,0,1,0; */
  MTX(&mtx,3,3,1);MTX(&mtx,3,5,1);MTX(&mtx,3,9,1);
  /* 5 row: 0,0,0,0,1,0,0,1,0,0,1; */
  MTX(&mtx,4,4,1);MTX(&mtx,4,7,1);MTX(&mtx,4,10,1);
  /* 6 row: 1,0,0,1,0,1,0,0,1,1,0; */
  MTX(&mtx,5,0,1);MTX(&mtx,5,3,1);MTX(&mtx,5,5,1);MTX(&mtx,5,8,1);
  MTX(&mtx,5,9,1);
  /* 7 row: 1,0,0,0,0,0,1,0,0,0,1; */
  MTX(&mtx,6,0,1);MTX(&mtx,6,6,1);MTX(&mtx,6,10,1);
  /* 8 row: 0,1,0,0,1,0,0,1,0,1,1;*/
  MTX(&mtx,7,1,1);MTX(&mtx,7,4,1);MTX(&mtx,7,7,1);MTX(&mtx,7,9,1);
  MTX(&mtx,7,10,1);
  /* 9 row: 0,0,0,0,0,1,0,0,1,0,0; */
  MTX(&mtx,8,5,1);MTX(&mtx,8,8,1);
  /* 10 row: 0,0,1,1,0,1,0,1,0,1,1; */
  MTX(&mtx,9,2,1);MTX(&mtx,9,3,1);MTX(&mtx,9,5,1);MTX(&mtx,9,7,1);
  MTX(&mtx,9,9,1);MTX(&mtx,9,10,1);
  /* 11 row: 0,0,1,0,1,0,1,1,0,1,1 */
  MTX(&mtx,10,2,1);MTX(&mtx,10,4,1);MTX(&mtx,10,6,1);MTX(&mtx,10,7,1);
  MTX(&mtx,10,9,1);MTX(&mtx,10,10,1);

  sp_matrix_yale_init(&yale,&mtx);

  etree = sp_matrix_yale_etree(&yale);
  result = etree[0] == etree_expected[0] - 1;
  do
  {
    for ( i = 1; i < 11; ++ i)
      result &= etree[i] == etree_expected[i] - 1;
    if (!result)
      break;

    ereach = malloc(11*sizeof(int));
    for ( i = 0; i < 11; ++ i)
    {
      sp_matrix_yale_ereach(&yale,etree,i,ereach);
      /* TODO: add test here */
      /* printf("L_%d:\t",i); */
      /* for (j = 0; j < 11; ++ j) */
      /*   if (ereach[j] != -1) */
      /*     printf("%d ",ereach[j]); */
      /* printf("\n"); */
    }
  } while(0);
  
  free(etree);
  free(ereach);
  sp_matrix_yale_free(&yale);
  sp_matrix_free(&mtx);
  printf("test_etree:\t*%s*\n",result ? "pass" : "fail");
  return result;
}

static int test_stack()
{
  int result = 0;
  int i;
  int_stack_ptr stack;
  const int count = 10;
  stack = int_stack_alloc(5,2);
  do
  {
    if (!int_stack_isempty(stack))
      break;
    
    int_stack_push(stack,10);
    int_stack_push(stack,20);
  
    if (int_stack_top(stack) != 20 )
      break;
    int_stack_pop(stack);
    if (int_stack_top(stack) != 10 )
      break;
    int_stack_pop(stack);

    if (!int_stack_isempty(stack))
      break;

    for ( i = 0; i < count; ++ i)
      int_stack_push(stack,i);
    for ( i = 0; i < count; ++ i)
    {
      if (int_stack_top(stack) != count - 1 - i)
        break;
      int_stack_pop(stack);
    }
    if ( i != count)
      break;
    
    if (!int_stack_isempty(stack))
      break;
    result = 1;
  } while(0);
  stack = int_stack_free(stack);
  
  printf("test_stack:\t*%s*\n",result ? "pass" : "fail");
  return result;
}

static int test_queue()
{
  int result = 0;
  int i;
  int_queue_ptr queue;
  const int count = 10;
  queue = int_queue_alloc();
  do
  {
    if (!int_queue_isempty(queue))
      break;
    
    int_queue_push(queue,10);
    int_queue_push(queue,20);
  
    if (int_queue_front(queue) != 10 )
      break;
    int_queue_pop(queue);
    if (int_queue_front(queue) != 20 )
      break;
    int_queue_pop(queue);

    if (!int_queue_isempty(queue))
      break;

    for ( i = 0; i < count; ++ i)
      int_queue_push(queue,i);
    for ( i = 0; i < count; ++ i)
    {
      if (int_queue_front(queue) != i)
        break;
      int_queue_pop(queue);
    }
    if ( i != count)
      break;
    
    if (!int_queue_isempty(queue))
      break;
    result = 1;
  } while(0);
  queue = int_queue_free(queue);
  
  printf("test_queue:\t*%s*\n",result ? "pass" : "fail");
  return result;
}

typedef struct
{
  int current;
  int* result;
} search_result;

static int mark_elt(int n, void* arg)
{
  search_result* search = (search_result*)arg;
  search->result[search->current] = n;
  search->current++;
  return 0;
}


static int test_tree_search()
{
  int result = 0;
  int i;
  int tree[] = {5,2,7,5,7,6,8,9,9,10,-1};
  int bfs[]  = {10,9,7,8,2,4,6,1,5,0,3};
  int dfs[]  = {10,9,8,6,5,3,0,7,4,2,1};
  search_result search;

  search.current = 0;
  search.result  = malloc(11*sizeof(int));
  do
  {
    /* test Breadth First Search */
    tree_bfs(tree,11,mark_elt,&search);
    result = search.result[0] == bfs[0];
    for ( i = 1; i < 11; ++ i)
      result &= search.result[i] == bfs[i];
    if (!result)
      break;

    /* test Deep First Search */
    search.current = 0;
    tree_dfs(tree,11,mark_elt,&search);
    result = search.result[0] == dfs[0];
    for ( i = 1; i < 11; ++ i)
      result &= search.result[i] == dfs[i];
  } while(0);
  free(search.result);
  printf("test_tree_search:\t*%s*\n",result ? "pass" : "fail");
  return result;
}


int main(/* int argc, char *argv[] */)
{
  /* logger */
  logger_parameters params;
  memset(&params,0,sizeof(params));
  params.log_file_path = "spmatrix.log";
  params.log_level = LOG_LEVEL_ALL;
  logger_init_with_params(&params);

  /* tests */
  test_sp_matrix();
  test_yale();
  test_triangle_solver();
  test_cg_solver();
  test_ilu();
  test_pcg_ilu_solver();
  test_cholesky();
  test_load();
  test_etree();
  test_stack();
  test_queue();
  test_tree_search();  
 
  /* finalize logger */
  logger_fini();
  return 0;
}
