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

int main(int argc, char *argv[])
{
  const int N = 4;           /* number of vertical blocks */
  const int M = 3;           /* number of horizontal blocks */
  int i,j;
  const double x = 1.0, y=1.0;                 /* upper-left point */
  const double dx = 1.0,dy = 1.0;              /* size of the block */
  geometry_2d g;
  prescr_boundary_2d b;
  generate_brick_mesh(N,M,x,y,dx,dy,&g,&b);

  printf("  1-----2-----3-----4-----5\n"
         "  |\\    |\\    |\\    |\\    |\n"
         "  |  \\  |  \\  |  \\  |  \\  |\n"
         "  |    \\|    \\|    \\|    \\|\n"
         "  6-----7-----8-----9----10\n"
         "  |\\    |\\    |\\    |\\    |\n"
         "  |  \\  |  \\  |  \\  |  \\  |\n"
         "  |    \\|    \\|    \\|    \\|\n"
         "  11---12----13----14----15\n"
         "  |\\    |\\    |\\    |\\    |\n"
         "  |  \\  |  \\  |  \\  |  \\  |\n"
         "  |    \\|    \\|    \\|    \\|\n"
         "  16---17----18----19----20\n");
  for (i = 0; i < g.triangles_count; ++ i)
    printf("[%d, %d, %d]\n",
           g.triangles[i].p1+1,
           g.triangles[i].p2+1,
           g.triangles[i].p3+1);
    

  return 0;
}

