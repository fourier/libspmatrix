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

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "sp_file.h"
#include "sp_utils.h"

/*
 * Data desrtiption for Matrix Market file format
 */

const char* MM_HEADER_STRING = "%%MatrixMarket";

typedef enum
{
  MM_MATRIX,
  MM_OTHER
} object_type;

typedef enum
{
  MM_COORDINATE,
  MM_ARRAY
} mm_storage_type;

typedef enum
{
  MM_REAL,
  MM_INTEGER,
  MM_COMPLEX,
  MM_PATTERN
} mm_elements_type;

typedef enum
{
  MM_GENERAL,
  MM_SYMMETRIC,
  MM_SKEW_SYMMETRIC,
  MM_HERMITIAN
} mm_matrix_portrait_type;

typedef struct
{
  int object;
  int storage;
  int elements;
  int portrait;
} mm_header;

static const char* mm_read_header(const char* contents, mm_header* header)
{
  const char* ptr = contents;
  const char* header_ptr = MM_HEADER_STRING;
  const char* word;
  /* read the necessary format header */
  while (*header_ptr)
  {
    if ( *header_ptr != *ptr )
    {
      fprintf(stderr, "Error reading MatrixMarket header!");
      return contents;
    }
    header_ptr++, ptr++;
  }
  /* extract object type */
  ptr = sp_extract_next_word(ptr,&word);
  printf("object = %s\n", word);
  header->object = !sp_istrcmp(word,"matrix") ? MM_MATRIX : MM_OTHER;
  free((char*)word);

  /* extract storage type */
  ptr = sp_extract_next_word(ptr,&word);
  printf("storage = %s\n", word);
  if (!sp_istrcmp(word,"coordinate"))
    header->storage = MM_COORDINATE;
  else if (!sp_istrcmp(word,"array"))
    header->storage = MM_ARRAY;
  else
    fprintf(stderr, "Unkown storage type: %s", word);
  free((char*)word);

  /* extract elements type */
  ptr = sp_extract_next_word(ptr,&word);
  printf("elements = %s\n", word);
  if (!sp_istrcmp(word,"real"))
    header->elements = MM_REAL;
  else if (!sp_istrcmp(word,"integer"))
    header->elements = MM_INTEGER;
  else if (!sp_istrcmp(word,"complex"))
    header->elements = MM_COMPLEX;
  else if (!sp_istrcmp(word,"pattern"))
    header->elements = MM_PATTERN;
  else
    fprintf(stderr, "Unkown elements type: %s", word);
  free((char*)word);

  /* extract portrait */
  ptr = sp_extract_next_word(ptr,&word);
  printf("portrait = %s\n", word);
  if (!sp_istrcmp(word,"general"))
    header->portrait = MM_GENERAL;
  else if (!sp_istrcmp(word,"symmetric"))
    header->portrait = MM_SYMMETRIC;
  else if (!sp_istrcmp(word,"skew-symmetric"))
    header->portrait = MM_SKEW_SYMMETRIC;
  else if (!sp_istrcmp(word,"Hermitian"))
    header->portrait = MM_HERMITIAN;
  else
    fprintf(stderr, "Unkown portrait type: %s", word);
  free((char*)word);

  return ptr;
}

static sp_matrix_ptr sp_matrix_load_file_mm(const char* filename)
{
  sp_matrix_ptr self = 0;
  mm_header header;
  FILE* file;
  const int block_size = 1024;
  const char* ptr;
  /* contents buffer */
  char* contents = calloc(block_size+1,1);
  /* auxulary counters */
  int read_chunk = 0,read = 0;
  file = fopen(filename,"rt");
  if (!file)
  {
    fprintf(stderr,"Cannot read file %s", filename);
    return self;
  }
  /* read file contents  */
  while(!feof(file))
  {
    read_chunk = fread(contents+read,1,block_size, file);
    read += read_chunk;
    
    contents = (char*)realloc(contents,read + block_size);
  }
  contents[read] = '\0';
  /* close the file */
  fclose(file);

  /* TODO: implementation */

  ptr = mm_read_header(contents, &header);
  
  /* free resources */
  free(contents);
  return self;
}

static sp_matrix_ptr sp_matrix_load_file_hb(const char* filename)
{
  sp_matrix_ptr self = 0;

  return self;
}


sp_matrix_ptr sp_matrix_load_file(const char* filename)
{
  sp_matrix_ptr self = 0;
  /* determine file extension */
  const char* ext = sp_parse_file_extension(filename);
  if (!ext)
  {
    fprintf(stderr,"File type is not supported: %s\n", filename);
    return self;
  }
  if ( !sp_istrcmp(ext,"mtx") )
    return sp_matrix_load_file_mm(filename);
  else
    fprintf(stderr,"File type is not supported: .%s\n", ext);
  
  return self;
}
