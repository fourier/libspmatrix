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

#include <time.h>
#include <sys/time.h>


#include "sp_matrix.h"
#include "sp_direct.h"
#include "sp_iter.h"
#include "sp_file.h"


#ifdef __MACH__
#include <mach/clock.h>
#include <mach/mach.h>
#endif

static void portable_gettime(struct timespec *ts)
{
  /* OS X does not have clock_gettime, use clock_get_time */
#ifdef __MACH__ 
  clock_serv_t cclock;
  mach_timespec_t mts;
  host_get_clock_service(mach_host_self(), CALENDAR_CLOCK, &cclock);
  clock_get_time(cclock, &mts);
  mach_port_deallocate(mach_task_self(), cclock);
  ts->tv_sec = mts.tv_sec;
  ts->tv_nsec = mts.tv_nsec;
#else
  clock_gettime(CLOCK_REALTIME, ts);
#endif
}

static void print_time_difference(struct timespec* t1,struct timespec* t2)
{
  long diff_sec, diff_micsec;
  diff_sec = t2->tv_sec - t1->tv_sec;
  diff_micsec = (t2->tv_nsec - t1->tv_nsec)/1000;
  if ( diff_micsec < 0)
  {
    diff_sec--;
    diff_micsec = 1000000+diff_micsec;
  }
  if ( diff_sec)
  {
    if (diff_micsec > 1000)
      printf("seconds: %ld, milliseconds: %.4f",
             diff_sec, diff_micsec/1000.);
    else
      printf("seconds: %ld, microseconds: %ld",
             diff_sec, diff_micsec);
  }
  else
  {
    if (diff_micsec > 1000)     /* time in milliseconds */
      printf("milliseconds: %.4f",
             diff_micsec/1000.);
    else
      printf("microseconds: %ld",
             diff_micsec);
  }
  printf("\n");
}


static void usage(const char* prog)
{
  printf("Usage: %s matrix-file\n",prog);
  exit(0);
}

static void print_error(const double* x, const double* x0, int size)
{
  double err = 0;
  int i = 0;
  for (; i < size; ++ i)
    if (fabs(x[i]-x0[i]) > err)
      err = fabs(x[i]-x0[i]);
  printf("%e\n",err);
}

int main(int argc, char *argv[])
{
  int i;
  sp_matrix_yale mtx,L;
  sp_chol_symbolic symb;
  struct timespec t1,t2,t3;
  double* x, *b, *x0;
  double desired_tolerance[3] = {1e-7,1e-12,1e-15};
  double tolerance;
  const int max_iter = 20000;
  int iter = max_iter;
  if(argc < 2)
    usage(argv[0]);
  if (sp_matrix_yale_load_file(&mtx,argv[1],CCS))
  {
    printf("Matrix %s statistics:\n",argv[1]);
    sp_matrix_yale_printf2(&mtx);
    portable_gettime(&t1);
    if (!sp_matrix_yale_chol_symbolic(&mtx,&symb))
      printf("Unable to create symbolic Cholesky decomposition\n");
    else
    {
      portable_gettime(&t2);
      if(!sp_matrix_yale_chol_numeric(&mtx,&symb,&L))
        printf("Unable to create numeric Cholesky decomposition\n");
      else
      {
        portable_gettime(&t3);
        printf("Cholesky decomposition statistics:\n");
        sp_matrix_yale_printf2(&L);
        printf("Nonzeros size increase:");
        printf("from %d to %d is %.2f %% size increase\n",
               mtx.nonzeros,L.nonzeros,L.nonzeros/(mtx.nonzeros/100.)-100.0);
        printf("Cholesky symblic decomposition calculation time: ");
        print_time_difference(&t1,&t2);
        printf("Cholesky numeric decomposition calculation time: ");
        print_time_difference(&t2,&t3);
        /* verify SLAE */
        x0 = calloc(mtx.rows_count,sizeof(double));
        x = calloc(mtx.rows_count,sizeof(double));
        b = calloc(mtx.rows_count,sizeof(double));
        for (i = 0; i < mtx.rows_count; ++ i)
          x0[i] = pow(-1,i%3) * (i % 10);
        sp_matrix_yale_mv(&mtx,x0,b);
        /* right part b and exact solution x0 constructed */
        portable_gettime(&t1);
        sp_matrix_yale_chol_numeric_solve(&L,b,x);
        portable_gettime(&t2);
        printf("Solving SLAE using Cholesky decomposition time: ");
        print_time_difference(&t1,&t2);
        printf("SLAE using Cholesky decomposition max error: ");
        print_error(x0,x,mtx.rows_count);
        /* sp_matrix_yale_convert_inplace(&mtx,CRS); */
        for (i = 0; i < 3; ++ i)
        {
          tolerance = desired_tolerance[i];
          iter = max_iter;
          /* memset(x,0,mtx.rows_count*sizeof(double)); */
          /* memcpy(x,b,mtx.rows_count*sizeof(double)); */
          portable_gettime(&t1);
          sp_matrix_yale_solve_cg(&mtx,b,b,&iter,&tolerance,x);
          portable_gettime(&t2);
          printf("Solving SLAE using Conjugate Gradient method");
          printf(" with tolerance %e(iterations: %d) time: ",
                 tolerance,iter);
          print_time_difference(&t1,&t2);
          printf("SLAE using Conjugate Gradient with tolerance");
          printf(" %e(iterations: %d) max error: ",desired_tolerance[i],iter);
          print_error(x0,x,mtx.rows_count);
        }
        sp_matrix_yale_free(&L);
        free(x0);
        free(x);
        free(b);
      }
      sp_matrix_yale_symbolic_free(&symb);
    }
    sp_matrix_yale_free(&mtx);
  }
  else
    printf("Unable to load file %s\n",argv[1]);
  return 0;
}
