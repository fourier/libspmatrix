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

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <memory.h>

#include "sp_matrix.h"
#include "sp_utils.h"
#include "sp_file.h"
#include "sp_cont.h"
#include "sp_tree.h"
#include "sp_perm.h"
#include "sp_test.h"
#include "logger.h"

static void sp_matrix_create_convert_mv()
{
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
    ASSERT_TRUE(EQL(x[i],expected[i]));
  
  /* 2nd test - conversion btw different storage types */
  sp_matrix_convert(&mtx,&mtx2,CCS);
  sp_matrix_mv(&mtx2,b,x);
  for (i = 0; i < size; ++ i)
    ASSERT_TRUE(EQL(x[i],expected[i]));
    
  sp_matrix_convert(&mtx2,&mtx3,CRS);
  sp_matrix_mv(&mtx3,b,x);
  for (i = 0; i < size; ++ i)
    ASSERT_TRUE(EQL(x[i],expected[i]));
  sp_matrix_free(&mtx2);
  sp_matrix_free(&mtx3);
  
  sp_matrix_free(&mtx);
}

static void yale_format()
{
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
  ASSERT_TRUE(memcmp(yale.offsets,offsets1,4*sizeof(int)) == 0);
  ASSERT_TRUE(memcmp(yale.indicies,indicies1,6*sizeof(int)) == 0);
  ASSERT_TRUE(memcmp(yale.values,values1,6*sizeof(double)) == 0);

  sp_matrix_yale_mv(&yale,x,y);
  for (i = 0; i < 3; ++ i)
    ASSERT_TRUE(EQL(y[i],b[i]));
  
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
  ASSERT_TRUE(memcmp(yale.offsets,offsets2,6*sizeof(int)) == 0);
  ASSERT_TRUE(memcmp(yale.indicies,indicies2,11*sizeof(int)) == 0);
  ASSERT_TRUE(memcmp(yale.values,values2,11*sizeof(double)) == 0);
  
  sp_matrix_free(&mtx);
  sp_matrix_yale_free(&yale);

}

static void permutations()
{
  int i;
  int p[] = {0,2,3,1};
  int pinv[] = {0,3,1,2};
  int result[4];
  sp_perm_inverse(p,4,result);
  for ( i = 0; i < 4; ++ i)
    ASSERT_TRUE(pinv[i] == result[i]);
}

static void sparse_permutations()
{
  /* Matrix(octave format)
     m = [1,2,0,-1;
          0,-2,0,0;
          2,1,0,-5;
          0,0,2,3]
  */
  /*
     Permutations:
     row permutation
     p = [1,3,4,2]
     column permutation
     q = [2,1,3,4]
     inverse row permutation
     p = [1,4,2,3]
  */
  /* int p[4] =    {0,2,3,1}; */
  int pinv[4] = {0,3,1,2};
  int q[4] =    {1,0,2,3};
  /*
   * Expected result: m(p,q):
   * 
   *   2   1   0  -1
   *   1   2   0  -5
   *   0   0   2   3
   *  -2   0   0   0
   *
   *  Expected result in CRS form:
   *  offsets  = [1        4        7     9 9]
   *  indicies = [1  2  4  1  2  4  3  4  1]
   *  values   = [2  1 -1  1  2 -5  2  3 -2]
   *  
   *  Expected result in CCS form:
   *  offsets  = [1        4     6  7       9]
   *  indicies = [1  2  4  1  2  3  1  2  3]
   *  values   = [2  1 -2  1  2  2 -1 -5  3]
   */
  int offsets_expected_crs[]  = {1, 4, 7, 9, 9};
  int indicies_expected_crs[] = {1, 2, 4, 1, 2, 4, 3, 4, 1};
  int values_expected_crs[]   = {2, 1,-1, 1, 2,-5, 2, 3,-2};

  int offsets_expected_ccs[]  = {1, 4, 6, 7, 9};
  int indicies_expected_ccs[] = {1, 2, 4, 1, 2, 3, 1, 2, 3};
  int values_expected_ccs[]   = {2, 1,-2, 1, 2, 2,-1,-5, 3};
  
  sp_matrix mtx_crs, mtx_ccs;
  sp_matrix_yale yale_crs, yale_ccs;
  sp_matrix_yale permuted_crs, permuted_ccs;
  int i;
  sp_matrix_init(&mtx_crs,4,4,2,CRS);
  MTX(&mtx_crs,0,0,1);MTX(&mtx_crs,0,1,2);MTX(&mtx_crs,0,3,-1);
  MTX(&mtx_crs,1,1,-2);
  MTX(&mtx_crs,2,0,2);MTX(&mtx_crs,2,1,1);MTX(&mtx_crs,2,3,-5);
  MTX(&mtx_crs,3,2,2);MTX(&mtx_crs,3,3,3);

  sp_matrix_yale_init(&yale_crs,&mtx_crs);
  sp_matrix_convert(&mtx_crs,&mtx_ccs,CCS);
  sp_matrix_yale_init(&yale_ccs,&mtx_ccs);
  
  ASSERT_TRUE(sp_matrix_yale_permute(&yale_crs,&permuted_crs,pinv,q));

  for (i = 0; i < 4; ++ i)
    ASSERT_TRUE(EQL(permuted_crs.offsets[i]+1,offsets_expected_crs[i]));
  ASSERT_TRUE(EQL(permuted_crs.offsets[i],offsets_expected_crs[i]));  
  for (i = 0; i < 9; ++ i)
    ASSERT_TRUE(EQL(permuted_crs.indicies[i]+1,indicies_expected_crs[i]));
  for (i = 0; i < 9; ++ i)
    ASSERT_TRUE(EQL(permuted_crs.values[i],values_expected_crs[i]));

  ASSERT_TRUE(sp_matrix_yale_permute(&yale_ccs,&permuted_ccs,pinv,q));
  for (i = 0; i < 4; ++ i)
    ASSERT_TRUE(EQL(permuted_ccs.offsets[i]+1,offsets_expected_ccs[i]));
  ASSERT_TRUE(EQL(permuted_ccs.offsets[i],offsets_expected_ccs[i]));
  for (i = 0; i < 9; ++ i)
    ASSERT_TRUE(EQL(permuted_ccs.indicies[i]+1,indicies_expected_ccs[i]));
  for (i = 0; i < 9; ++ i)
    ASSERT_TRUE(EQL(permuted_ccs.values[i],values_expected_ccs[i]));
  
  
  sp_matrix_yale_free(&permuted_crs);
  sp_matrix_free(&mtx_crs);
  sp_matrix_yale_free(&yale_crs);
  sp_matrix_yale_free(&permuted_ccs);
  sp_matrix_free(&mtx_ccs);
  sp_matrix_yale_free(&yale_ccs);

}

static void triangle_solver()
{
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
    ASSERT_TRUE(EQL(x_expected[i],x[i]));
  
  sp_matrix_convert(&mtx,&mtx2,CCS);
  memset(x,0,sizeof(double)*5);
  sp_matrix_lower_solve(&mtx,5,b,x);
  for (i = 0; i < 5; ++ i)
    ASSERT_TRUE(EQL(x_expected[i],x[i]));

  sp_matrix_free(&mtx2);
  sp_matrix_free(&mtx);
}

static void cg_solver()
{
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
  ASSERT_TRUE(tolerance < desired_tolearance*10);
  
  sp_matrix_free(&mtx);
  sp_matrix_yale_free(&yale);
}

static void ilu_and_skyline()
{
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
    ASSERT_TRUE(fabs(ILU.ilu_diag[i] - lu_diag_expected[i]) < 1e-5);
  
  for (i = 0; i <  m.tr_nonzeros; ++ i)
    ASSERT_TRUE(fabs(ILU.ilu_lowertr[i] - lu_lowertr_expected[i]) < 1e-5);

  for (i = 0; i <  m.tr_nonzeros; ++ i)
    ASSERT_TRUE(fabs(ILU.ilu_uppertr[i] - lu_uppertr_expected[i]) < 1e-5);
  
  /*
   * test for solving Lx=b
   */
  /* prepare a right-part vector */
  sp_matrix_skyline_ilu_lower_mv(&ILU,x_exact,b);
    
  /* solve for x */
  sp_matrix_skyline_ilu_lower_solve(&ILU,b,x);
  /* test result */
  for ( i = 0; i < m.rows_count; ++ i)
    ASSERT_TRUE(EQL(x[i],x_exact[i]));
    
  
  /*
   * test for solving Ux=b
   */
  memset(b,0,sizeof(b));
  memset(x,0,sizeof(x));
  /* prepare a right-part vector */
  sp_matrix_skyline_ilu_upper_mv(&ILU,x_exact,b);
    
  /* solve for x */
  sp_matrix_skyline_ilu_upper_solve(&ILU,b,x);
  /* test result */
  for ( i = 0; i < m.rows_count; ++ i)
    ASSERT_TRUE(EQL(x[i],x_exact[i]));
  
  sp_matrix_free(&mtx);
  sp_matrix_skyline_ilu_free(&ILU);
}

static void pcg_ilu_solver()
{
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
  ASSERT_TRUE(tolerance < desired_tolearance*10);
  
  sp_matrix_skyline_ilu_free(&ilu);
  sp_matrix_free(&mtx);
  sp_matrix_yale_free(&yale);
}

static void cholesky()
{
  /* initial matrix(octave representation):
     m = [90, 6, 4, 46, 29, 0, 26;
     6, 127, 34, 22, 7, 0, 38;
     4, 34, 108, 40, 2, 0, 4;
     46, 22, 40, 96, 24, 0, 6;
     29, 7, 2, 24, 155, 0, 37;
     0, 0, 0, 0, 0, 64, 0;
     26, 38, 4, 6, 37, 0, 70];
  */
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
}


static void load_from_files()
{
  sp_matrix_yale mtx;
  int result;
  EXPECT_TRUE((result = sp_matrix_yale_load_file(&mtx,"bcsstk09.rsa",CRS)));
  if (result) sp_matrix_yale_free(&mtx);
  EXPECT_TRUE((result = sp_matrix_yale_load_file(&mtx,"af23560.rua",CRS)));
  if (result) sp_matrix_yale_free(&mtx);  
  EXPECT_TRUE((result = sp_matrix_yale_load_file(&mtx,"kershaw_rua.hb",CRS)));
  if (result) sp_matrix_yale_free(&mtx);
  EXPECT_TRUE((result = sp_matrix_yale_load_file(&mtx, "5by5_rua.hb",CRS)));

  /* sp_matrix_save_file(result,"export.mtx"); */
  if (result) sp_matrix_yale_free(&mtx);  
}


sp_matrix mtx;
sp_matrix_yale yale;

static void test_etree_init()
{
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
  /* sp_matrix_save_file(&mtx, "davis.txt"); */
}

static void test_etree_fini()
{
  sp_matrix_yale_free(&yale);
  sp_matrix_free(&mtx);
}

static void etree_create_etree()
{
  int i;
  /*
   * octave representation of the expected result:
   * etree(sparse(m))
   * Image:
   *        11
   *        |
   *        10
   *       /  \
   *      /    9
   *     /     |
   *    8      7
   *   / \     |
   *  3   \    6
   *  |   |   / \
   *  2   5  1   4
   *
   */
  int etree_expected[] = {6,3,8,6,8,7,9,10,10,11,0};
  int etree[11];
  ASSERT_TRUE(sp_matrix_yale_etree(&yale,etree));
  ASSERT_TRUE(etree[0] == etree_expected[0] - 1);
  
  for ( i = 1; i < 11; ++ i)
    ASSERT_TRUE(etree[i] == etree_expected[i] - 1);
  /* tree_dot_printf(etree,11); */
}

static void etree_postorder()
{
  int i,j,k,l;
  sp_matrix_yale yale_post;
  /*
   * octave representation of the expected result of postordering:
   * p = [2, 3, 5, 8, 1, 4, 6, 7, 9, 10, 11];
   * Image:
   *
   *         11
   *         |
   *         10
   *        /  \
   *       /    9
   *      /     |
   *     4      8
   *    / \     |
   *   2   \    7
   *   |   |   / \
   *   1   3  5   6
   *
   */
  int postorder_expected[] = {2, 3, 5, 8, 1, 4, 6, 7, 9, 10, 11};
  int etree[11];
  int postorder[11];
  int postinv[11];
  /*
   * octave representation of the postordered matrix:
   * m(p,p)
   */
  int postordered_matrix[][11] = {{1,1,0,1,0,0,0,0,0,0,0},
                                  {1,1,0,0,0,0,0,0,0,1,1},
                                  {0,0,1,1,0,0,0,0,0,0,1},
                                  {1,0,1,1,0,0,0,0,0,1,1},
                                  {0,0,0,0,1,0,1,1,0,0,0},
                                  {0,0,0,0,0,1,1,0,0,1,0},
                                  {0,0,0,0,1,1,1,0,1,1,0},
                                  {0,0,0,0,1,0,0,1,0,0,1},
                                  {0,0,0,0,0,0,1,0,1,0,0},
                                  {0,1,0,1,0,1,1,0,0,1,1},
                                  {0,1,1,1,0,0,0,1,0,1,1}};
  /* for ( j = 0; j < 11; ++ j) */
  /* { */
  /*   for ( i = 0; i < 11; ++ i) */
  /*     if ( postordered_matrix[j][i] ) */
  /*       printf("A_%d,%d ",j+1,i+1); */
  /*   printf("\n"); */
  /* } */

  ASSERT_TRUE(sp_matrix_yale_etree(&yale,etree));
  
  tree_postorder_perm(etree,11,postorder);
  for ( i = 0; i < 11; ++ i)
    ASSERT_TRUE(postorder[i] == postorder_expected[i]-1);
  
  /* for ( i = 0; i < 11; ++ i) */
  /*   printf("post[%d] = %d,\tnode %d in source tree\ */
  /* is node %d in postordered \n", */
  /*          i+1, postorder[i]+1, postorder[i]+1,i+1 ); */

  for ( i = 0; i < 11; ++ i)
  {
    k = postorder[i];
    for ( j = 0; j < 11; ++ j)
    {
      l = postorder[j];
      if (sp_matrix_element_ptr(&mtx,k,l))
        ASSERT_TRUE(postordered_matrix[i][j]);
      if (postordered_matrix[i][j])
        ASSERT_TRUE(sp_matrix_element_ptr(&mtx,k,l));
    }
  }

  sp_perm_inverse(postorder,11,postinv);
  sp_matrix_yale_permute(&yale,&yale_post,postinv,postinv);
  for ( i = 0; i < 11; ++ i)
  {
    for ( j = yale_post.offsets[i]; j < yale_post.offsets[i+1]; ++ j)
    {
      if (yale_post.storage_type == CRS)
      {
        k = i;
        l = yale_post.indicies[j];
      }
      else
      {
        k = yale_post.indicies[j];
        l = i;
      }
      ASSERT_TRUE(postordered_matrix[k][l]);
    }
  }
      
}

static void etree_ereach()
{
  int i,j,k,count;
  int etree[11];
  int ereach[11];
  /* Cholesky factor L of M */
  int chol_portrait[][11] =
    {{1,0,0,0,0,0,0,0,0,0,0},
     {0,1,0,0,0,0,0,0,0,0,0},
     {0,1,1,0,0,0,0,0,0,0,0},
     {0,0,0,1,0,0,0,0,0,0,0},
     {0,0,0,0,1,0,0,0,0,0,0},
     {1,0,0,1,0,1,0,0,0,0,0},
     {1,0,0,0,0,1,1,0,0,0,0},
     {0,1,1,0,1,0,0,1,0,0,0},
     {0,0,0,0,0,1,1,0,1,0,0},
     {0,0,1,1,0,1,1,1,1,1,0},
     {0,0,1,0,1,0,1,1,1,1,1}};
    
  ASSERT_TRUE(sp_matrix_yale_etree(&yale,etree));
  for ( i = 0; i < 11; ++ i)
  {
    ASSERT_TRUE((count = sp_matrix_yale_ereach(&yale,etree,i,ereach)) > 0);
    for (j = 0; j < count; ++ j)
      ASSERT_TRUE(chol_portrait[i][ereach[j]]);
    k = 0;
    for (j = 0; j < 11; ++ j)
      if (chol_portrait[i][j]) k++;
    ASSERT_TRUE(k == count);
  }
}

static void etree_rowcolcounts()
{
  int i,j,k,count;
  int etree[11];
  int rowcounts[11];
  int colcounts[11];
  int rowcounts_expected[11] = {1,1,2,1,1,3,3,4,3,7,7};
  int colcounts_expected[11] = {3,3,4,3,3,4,4,3,3,2,1};
  
  ASSERT_TRUE(sp_matrix_yale_etree(&yale,etree));
  ASSERT_TRUE(sp_matrix_yale_chol_counts(&yale,etree,rowcounts,colcounts));
  for ( i = 0; i < 11; ++ i)
  {
    ASSERT_TRUE(rowcounts[i] == rowcounts_expected[i]);
    ASSERT_TRUE(colcounts[i] == colcounts_expected[i]);
  }
}

static void etree_rowcount()
{
  /* values received from the T.A.Davis's firstdesc function, but 1-based */
  int first_expected[] = {5,1,1,6,3,5,5,1,5,1,1};
  /* values calculated by calculation of the maxfirst with
   * forcefully postordered matrix, 1-based */
  int maxfirst_expected[] = {-1, 1, -1, 3, -1, -1, 6, 5, 5, 6, 5};
  int i,j,p,k;
  int etree[11];
  int postorder[11];
  int postinv[11];
  int first[11];
  int level[11];
  int maxfirst[11];
  int maxfirst1[11];
  int skeleton[11][11] = {
    {1,0,0,0,0,0,0,0,0,0,0},
    {0,1,0,0,0,0,0,0,0,0,0},
    {0,0,1,0,0,0,0,0,0,0,0},
    {0,0,0,1,0,0,0,0,0,0,0},
    {0,0,0,0,1,0,0,0,0,0,0},
    {0,0,0,0,0,1,0,0,0,0,0},
    {0,0,0,0,0,0,1,0,0,0,0},
    {0,0,0,0,0,0,0,1,0,0,0},
    {0,0,0,0,0,0,0,0,1,0,0},
    {0,0,0,0,0,0,0,0,0,1,0},
    {0,0,0,0,0,0,0,0,0,0,1}};
  sp_matrix_yale yale_post;
  ASSERT_TRUE(sp_matrix_yale_etree(&yale,etree));
  tree_postorder_perm(etree,11,postorder);
  tree_node_levels(etree,11,level);
  /* for ( i = 0; i < 11; ++ i) */
  /*   printf("%d ",level[i]); */
  /* printf("\n"); */
  printf("Postordered elimination tree:\n");
  printf("         11\n");
  printf("         |\n");
  printf("         10\n");
  printf("        /  \\\n");
  printf("       /    9\n");
  printf("      /     |\n");
  printf("     4      8\n");
  printf("    / \\     |\n");
  printf("   2   \\    7\n");
  printf("   |   |   / \\\n");
  printf("   1   3  5   6\n");

  tree_first_descendant(etree,11,postorder,first);
  for ( i = 0; i < 11; ++ i)
    ASSERT_TRUE(first[i]+1 == first_expected[i]);

  /* let's calulate the skeleton matrix */
  /* matrix and etree shall be postordered */
  sp_perm_inverse(postorder,11,postinv);
  sp_matrix_yale_permute(&yale,&yale_post,postinv,postinv);
   
  for (j = 0; j < 11; ++ j)
    maxfirst[j] = -1, maxfirst1[j] = -1;

  for (j = 0; j < 11; ++ j)
  {
    /* j is k in postordered tree */
    k = postorder[j]; /* we didn't permute the 'first' array,
                         so it's time to do it here */

    for ( p = yale_post.offsets[j]; p < yale_post.offsets[j+1]; ++p )
    {
      i = yale_post.indicies[p];
      if ( i > j)
      {
        printf("A_%d,%d ",i+1,j+1);
        if (first[k] > maxfirst[i])
        {
          printf("node %d is a leaf in the %d-th subtree\n",j+1,i+1);
          maxfirst[i] = first[k];
          skeleton[j][i] = skeleton[i][j] = 1;
        }
      }
    }
    printf("\n");
  }
  for (i = 0; i < 11; ++ i)
  {
    for (j = 0; j < 11; ++ j)
      printf("%d ",skeleton[i][j]);
    printf("\n");
  }
  printf("maxfirst = [");
  for (i = 0; i < 10; ++ i)
    printf("%d, ", maxfirst[i] == -1 ? -1 : maxfirst[i] + 1);
  printf("%d]\n", maxfirst[i] == -1 ? -1 : maxfirst[i] + 1);
  
  /* the next step is to calculate skeleton matrix w/o construction
     of the temporary yale_post - postordered matrix */

  for (j = 0; j < 11; ++ j)
    maxfirst1[j] = -1;
  int r,s;
  for (j = 0; j < 11; ++ j)
  {
    /* j - column number */
    /* j is k in postordered tree */
    k = postorder[j];
    s = postinv[j];           /* new col number */
    printf("j = %d, postorder[j] = %d, postinv[j] = %d\n",
           j,postorder[j],postinv[j]);
    for ( p = yale.offsets[k]; p < yale.offsets[k+1]; ++p )
    {
      i = yale.indicies[p];
      /* i - row number */
      r = postinv[i];           /* new row number */
      if (r > s)
      {
        printf("A_%d,%d ",r+1,s+1);
        if (first[k] > maxfirst1[r])
          maxfirst1[r] = first[k];
      }
    }
    printf("\n");
  }
  printf("maxfirst1 = [");
  for (i = 0; i < 10; ++ i)
    printf("%d, ", maxfirst1[i] == -1 ? -1 : maxfirst1[i] + 1);
  printf("%d]\n", maxfirst1[i] == -1 ? -1 : maxfirst1[i] + 1);
  for (i = 0; i < 11; ++ i)
    ASSERT_TRUE(maxfirst[i] == maxfirst1[i]);
  
}



static void stack_container()
{
  int i;
  int_stack_ptr stack;
  const int count = 10;
  stack = int_stack_alloc(5,2);
  ASSERT_TRUE(int_stack_isempty(stack));

  int_stack_push(stack,10);
  int_stack_push(stack,20);
  ASSERT_TRUE(int_stack_top(stack) == 20 );

  int_stack_pop(stack);
  ASSERT_TRUE(int_stack_top(stack) == 10 );

  int_stack_pop(stack);

  ASSERT_TRUE(int_stack_isempty(stack));

  for ( i = 0; i < count; ++ i)
    int_stack_push(stack,i);
  for ( i = 0; i < count; ++ i)
  {
    ASSERT_TRUE(int_stack_top(stack) == count - 1 - i);
    int_stack_pop(stack);
  }
  ASSERT_TRUE(i == count);
    
  ASSERT_TRUE(int_stack_isempty(stack));
  stack = int_stack_free(stack);
}

static void queue_container()
{
  int i;
  int_queue_ptr queue;
  const int count = 10;
  queue = int_queue_alloc();

  ASSERT_TRUE(int_queue_isempty(queue));
    
  int_queue_push(queue,10);
  int_queue_push(queue,20);
  
  ASSERT_TRUE(int_queue_front(queue) == 10 );

  int_queue_pop(queue);
  ASSERT_TRUE(int_queue_front(queue) == 20 );

  int_queue_pop(queue);
  ASSERT_TRUE(int_queue_isempty(queue));

  for ( i = 0; i < count; ++ i)
    int_queue_push(queue,i);
  for ( i = 0; i < count; ++ i)
  {
    ASSERT_TRUE(int_queue_front(queue) == i);
    int_queue_pop(queue);
  }
  ASSERT_TRUE( i == count);
  ASSERT_TRUE(int_queue_isempty(queue));

  queue = int_queue_free(queue);
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


static void tree_search()
{
  int i;
  int tree[] = {5,2,7,5,7,6,8,9,9,10,-1};
  int bfs[]  = {10,9,7,8,2,4,6,1,5,0,3};
  int dfs[]  = {10,9,8,6,5,3,0,7,4,2,1};
  search_result search;

  search.current = 0;
  search.result  = malloc(11*sizeof(int));

  /* test Breadth First Search */
  tree_bfs(tree,11,mark_elt,&search);
  ASSERT_TRUE(search.result[0] == bfs[0]);
    
  for ( i = 1; i < 11; ++ i)
    ASSERT_TRUE(search.result[i] == bfs[i]);

  /* test Deep First Search */
  search.current = 0;
  tree_dfs(tree,11,mark_elt,&search);
  ASSERT_TRUE(search.result[0] == dfs[0]);
  for ( i = 1; i < 11; ++ i)
    ASSERT_TRUE(search.result[i] == dfs[i]);

  free(search.result);
}

static void save_vector(int* v, int size, const char* fname)
{
  int i = 0;
  FILE *f = fopen(fname,"wt+");
  if (f)
  {
    for (; i < size; ++ i)
      fprintf(f,"%d\n", v[i]);
    fclose(f);
  }
}

static void big_matrix()
{
  int result;
  sp_matrix_yale yale;
  int *etree, *rowcounts, *colcounts;
  int i,j;
  EXPECT_TRUE((result = sp_matrix_yale_load_file(&yale,"bcsstk11.mtx",CCS)));
  if (result)
  {
    /* in order to test one need to have bcsstk11.mtx
     * 1. download mmread.m from matrixmarket site and put it into the
     *    ../octave-mm directory
     * 2. start octave in ../octave-mm directory
     * 3. copy bcsstk11.mtx to ../octave-mm directory
     * 4. run this test
     * 5. Perform in octave:
     * ------------------------------
     *    A = mmread("bcsstk11.mtx");
     *    e = etree(A);
     *    for i = 1:rows
     *      if e(i) == 0
     *        j = j+1;
     *      end
     *    end
     * ------------------------------
     *  j shall be equal to 9
     * ------------------------------
     * [COUNT, H, PARENT, POST, R] = symbfact(A);
     * e1 = load("-ascii","etree.txt");
     * max(abs(e-e1'))
     * cols = load("-ascii","cols.txt");
     * max(abs(COUNT-cols))
     * ------------------------------
     * last to operations(max and max) shall give 0
     */
    sp_matrix_yale_printf2(&yale);
    etree = (int*)malloc(yale.rows_count*sizeof(int));
    rowcounts = (int*)malloc(yale.rows_count*sizeof(int));
    colcounts = (int*)malloc(yale.rows_count*sizeof(int));

    ASSERT_TRUE(sp_matrix_yale_etree(&yale,etree));
    j = 0;
    for (i = 0; i < yale.rows_count; ++ i)
      if (etree[i] == -1)
        j++;
    
    printf("Number of roots in the Elimination tree = %d\n",j);
    ASSERT_TRUE(sp_matrix_yale_chol_counts(&yale,etree,rowcounts,colcounts));
    for (i = 0; i < yale.rows_count; ++ i)
      etree[i]++;
    save_vector(etree,yale.rows_count,"../octave-mm/etree.txt");
    save_vector(rowcounts,yale.rows_count,"../octave-mm/rows.txt");
    save_vector(colcounts,yale.rows_count,"../octave-mm/cols.txt");
    
    free(etree);
    free(rowcounts);
    free(colcounts);
    
    sp_matrix_yale_free(&yale);
  }
}

/* #define DEF_TEST(name) static void name() { printf( #name "\n"); } */

/* DEF_TEST(init1) */
/* DEF_TEST(fini1) */
/* DEF_TEST(test1_1) */
/* DEF_TEST(test1_2) */

/* DEF_TEST(init2) */
/* DEF_TEST(fini2) */
/* DEF_TEST(test2_1) */
/* DEF_TEST(test2_2) */


int main(int argc, const char *argv[])
{
  sp_test_suite *suite1;
  /* logger */
  logger_parameters params;
  memset(&params,0,sizeof(params));
  params.log_file_path = "spmatrix.log";
  params.log_level = LOG_LEVEL_ALL;
  logger_init_with_params(&params);

  /* tests */
  SP_ADD_TEST(sp_matrix_create_convert_mv);
  SP_ADD_TEST(yale_format);
  SP_ADD_TEST(permutations);
  SP_ADD_TEST(sparse_permutations);
  SP_ADD_TEST(triangle_solver);
  SP_ADD_TEST(cg_solver);
  SP_ADD_TEST(ilu_and_skyline);
  SP_ADD_TEST(pcg_ilu_solver);
  SP_ADD_TEST(cholesky);
  SP_ADD_TEST(load_from_files);
  SP_ADD_TEST(stack_container);
  SP_ADD_TEST(queue_container);
  SP_ADD_TEST(tree_search);
  
  suite1 = sp_add_suite("etree",test_etree_init,test_etree_fini);
  SP_ADD_SUITE_TEST(suite1,etree_create_etree);
  SP_ADD_SUITE_TEST(suite1,etree_postorder);
  SP_ADD_SUITE_TEST(suite1,etree_ereach);
  SP_ADD_SUITE_TEST(suite1,etree_rowcolcounts);
  /* SP_ADD_SUITE_TEST(suite1,etree_rowcount); */
  SP_ADD_TEST(big_matrix);
  
  sp_run_tests(argc,argv);

  /* finalize logger */
  logger_fini();

  return 0;
}
