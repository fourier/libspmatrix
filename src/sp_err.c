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

#include <stdlib.h>
#include <execinfo.h>

#include "sp_log.h"
#include "sp_err.h"

/* Print the backtrace to LOGERROR. */
static void log_backtrace (void)
{
  void *array[30];
  size_t size;
  char **strings;
  size_t i;

  /* number of stack frames obtained */
  size = backtrace (array, 30);
  strings = backtrace_symbols (array, size);
    
  for (i = 0; i < size; i++)
    LOGERROR("#%d  %s", i, strings[i]);
     
  free (strings);
}


void sp_error()
{
  LOGERROR("Unrecoverable error occured. Backtrace");
  log_backtrace();
  exit(1);
}
