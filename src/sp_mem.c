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
#include <string.h>

#include "sp_mem.h"
#include "sp_err.h"
#include "sp_log.h"

#define CHUNK_HEAD(x) (void*)((size_t*)x-1)
#define CHUNK_SIZE(x) *((size_t*)x-1)
#define CHUNK_PTR(x) ((size_t*)x+1)

size_t allocated = 0;

size_t spallocated()
{
  return allocated;
}

void* spalloc(size_t size)
{
  void* chunk = malloc(size+sizeof(size_t));
  if (!chunk)
  {
    LOGERROR("spalloc: cannot allocate memory block of size %zd",size);
    sp_error();
  }
  *(size_t*)(chunk) = size;
  allocated += size;
  return CHUNK_PTR(chunk);
}

void* spcalloc(size_t nmemb, size_t size)
{
  void* chunk = spalloc(nmemb*size);
  if (!chunk)
  {
    LOGERROR("spcalloc: cannot allocate memory block of size %zd",
             size*nmemb);    
    sp_error();
  }
  memset(chunk,0,size*nmemb);
  return chunk;
}

void* sprealloc(void* ptr, size_t size)
{
  size_t old_sz = CHUNK_SIZE(ptr);
  void* chunk = realloc(CHUNK_HEAD(ptr), size+sizeof(size_t));
  if (!chunk)
  {
    LOGERROR("sprealloc: cannot reallocate memory block of size %zd",size);
    sp_error();
  }
  *(size_t*)(chunk) = size;
  allocated += (size - old_sz); 
  
  return CHUNK_PTR(chunk);
}

void spfree(void* ptr)
{
  size_t sz = CHUNK_SIZE(ptr);
  allocated -= sz;
  free(CHUNK_HEAD(ptr));
}

void* memdup(const void* src, int bytes)
{
  void* result = spalloc(bytes);
  if (result)
  {
    memcpy(result,src,bytes);
  }
  return result;
}
