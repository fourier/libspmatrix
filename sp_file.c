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
#include <assert.h>

#include "sp_file.h"
#include "sp_utils.h"

/*
 * Data desrtiption for Matrix Market file format
 */

const char* MM_HEADER_STRING = "%%MatrixMarket";
const char* MM_HEADER_OUTPUT_ABOUT = "% Created by libspmatrix (c) "
  "Alexey Veretennikov, https://github.com/fourier/libspmatrix\n";

/*
 * length of the HB file line is 80, but sometimes it can be more
 * possibly because some bugs in export software
 * Counter-example: bcsstk09.rsa
 */
#define HB_LINE_SIZE 100

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
      fprintf(stderr, "Error reading MatrixMarket header!\n");
      return contents;
    }
    header_ptr++, ptr++;
  }
  /* extract object type */
  ptr = sp_extract_next_word(ptr,&word);
  header->object = !sp_istrcmp(word,"matrix") ? MM_MATRIX : MM_OTHER;
  free((char*)word);

  /* extract storage type */
  ptr = sp_extract_next_word(ptr,&word);
  if (!sp_istrcmp(word,"coordinate"))
    header->storage = MM_COORDINATE;
  else if (!sp_istrcmp(word,"array"))
    header->storage = MM_ARRAY;
  else
  {
    fprintf(stderr, "Unkown storage type: %s\n", word);
    ptr = contents;
  }
  free((char*)word);

  /* extract elements type */
  ptr = sp_extract_next_word(ptr,&word);
  if (!sp_istrcmp(word,"real"))
    header->elements = MM_REAL;
  else if (!sp_istrcmp(word,"integer"))
    header->elements = MM_INTEGER;
  else if (!sp_istrcmp(word,"complex"))
    header->elements = MM_COMPLEX;
  else if (!sp_istrcmp(word,"pattern"))
    header->elements = MM_PATTERN;
  else
  {
    fprintf(stderr, "Unkown elements type: %s\n", word);
    ptr = contents;
  }
  free((char*)word);

  /* extract portrait */
  ptr = sp_extract_next_word(ptr,&word);
  if (!sp_istrcmp(word,"general"))
    header->portrait = MM_GENERAL;
  else if (!sp_istrcmp(word,"symmetric"))
    header->portrait = MM_SYMMETRIC;
  else if (!sp_istrcmp(word,"skew-symmetric"))
    header->portrait = MM_SKEW_SYMMETRIC;
  else if (!sp_istrcmp(word,"Hermitian"))
    header->portrait = MM_HERMITIAN;
  else
  {
    fprintf(stderr, "Unkown portrait type: %s\n", word);
    ptr = contents;
  }
  free((char*)word);

  return ptr;
}

/* return 1 if mm-matrix format is supported, 0 otherwise  */
static int mm_validate_header(mm_header* header)
{
  if (header->object != MM_MATRIX)
  {
    fprintf(stderr,"MM objects other than matrix are not supported\n");
    return 0;
  }
  if (header->storage != MM_COORDINATE)
  {
    fprintf(stderr,"Currently only coordinate format of MM files supported\n");
    return 0;
  }
  if (header->elements == MM_COMPLEX)
  {
    fprintf(stderr,"Complex matricies not supported\n");
    return 0;
  }
  if (header->portrait == MM_HERMITIAN)
  {
    fprintf(stderr,"Hermitian matricies not supported\n");
    return 0;
  }
  return 1;
}

/*
 * Load matrix in MatrixMarket format
 */
static sp_matrix_ptr sp_matrix_load_file_mm(const char* filename, int storage_type)
{
  sp_matrix_ptr self = 0;
  mm_header header;
  const char* ptr;
  const char* line;
  int block_number = 0;
  int rows, cols, nonzeros;
  int element_number = -1;
  int i,j;
  double value;
  /* read the file contents  */
  char* contents = sp_read_text_file(filename);
  if (!contents)
    return self;

  /* split by lines */
  line = strtok(contents,"\n");
  while (line)
  {
    if (block_number == 0)      /* header */
    {
      ptr = mm_read_header(line, &header);
      if (ptr == line)          /* error */
        break;
      if (!mm_validate_header(&header))
      {
        printf("Not supported matrix type\n");
        break;
      }
      /* goto next step - parse sizes */
      block_number++;
    }
    else if (block_number == 1) /* sizes */
    {
      ptr = sp_skip_whitespaces(line);
      if (*ptr && *ptr != '%')  /* not empty line and not comment */
      {
        /* read sizes */
        if (sscanf(line, "%d %d %d", &rows, &cols, &nonzeros) != 3)
        {
          fprintf(stderr, "Unable to parse sizes\n");
          break;
        }
        self = calloc(1,sizeof(sp_matrix));
        sp_matrix_init(self,rows,cols,nonzeros/rows + 1,storage_type);
        element_number = 0;     /* start with first element */
        block_number++;
      }
    }
    else if (block_number == 2) /* data */
    {
      ptr = sp_skip_whitespaces(line);
      if (*ptr && *ptr != '%')  /* not empty line and not comment */
      {
        if ( header.elements == MM_PATTERN)
        {
          if (sscanf(line, "%d %d", &i, &j) != 2)
          {
            fprintf(stderr, "Unable to parse line %s\n",line);
            break;
          }
          value = 1;
        }
        else
        {
          if (sscanf(line, "%d %d %lg", &i, &j, &value) != 3)
          {
            fprintf(stderr, "Unable to parse line %s\n",line);
            break;
          }

        }
        sp_matrix_element_add(self,i-1,j-1,value);
        /* handle symmetry property */
        if ( i != j )
        {
          if ( header.portrait == MM_SYMMETRIC )
            sp_matrix_element_add(self,j-1,i-1,value);
          else if (header.portrait == MM_SKEW_SYMMETRIC)
            sp_matrix_element_add(self,j-1,i-1,-value);
        }
        ++ element_number;
      }
    }
    line = strtok(NULL,"\n");
  }

  /* free resources */
  free(contents);
  /* check for error */
  if ( element_number != nonzeros)
  {
    fprintf(stderr,"Error loading matrix, expected %d nonzeros, parsed %d\n",
            nonzeros, element_number);
    sp_matrix_free(self);
    free(self);
    self = 0;
  }
  else
  {
    sp_matrix_compress(self);
  }
  return self;
}

static int hb_extract_positional_format(const char* from, int size,
                                        fortran_io_format* fmt)
{
  int result = 0;
  char* buf = malloc(size+1);
  char* ptr = buf;
  int i = 0;
  while (i < size && *from)
  {
    *ptr++ = *from++;
    ++i;
  }
  buf[size] = '\0';
  result = sp_parse_fortran_format(buf,fmt);
  free(buf);
  return result;
}

/*
 * Load matrix in Harwell Boeing format
 */
static sp_matrix_ptr sp_matrix_load_file_hb(const char* filename,
                                            int storage_type)
{
  sp_matrix_ptr self = 0;
  int i,j,n;
  /* 1 for '\n' */

  /* HB format line limitation 80 chars */
  char buf[HB_LINE_SIZE+1];
  char* ptr;
  /* constants from HB format */
  /* for line 2 */
  int totcrd, ptrcrd, indcrd, valcrd, rhscrd;
  /* for line 3 */
  int nrow, ncol, nnzero;
  /* for line 4 */
  fortran_io_format ptrfmt, indfmt, valfmt, rhsfmt;
  /* temporary exctracted values */
  fortran_number* fortran_numbers;
  /* data in column-wise triplet form */
  int* colptr    = 0;                /* location of first entry */
  int* rowind    = 0;                /* row indicies */
  double* values = 0;                /* numerical valus */
  /* example: */
  /*
   * 1. -3.  0. -1.  0.
   * 0.  0. -2.  0.  3.
   * 2.  0.  0.  0.  0.
   * 0.  4.  0. -4.  0.
   * 5.  0. -5.  0.  6.
   *
   * array values, for 1-based indicies:
   * subscripts | 1   2   3   4   5   6   7   8   9   10  11
   * --------------------------------------------------------
   * colptr     | 1   4   6   8  10  12
   * rowind     | 1   3   5   1   4   2   5   1   4    2   5
   * values     | 1.  2.  5. -3.  4. -2. -5. -1. -4.   3.  6.
   */
  
  /*
   * temporary variables - counter of extracted by sp_extract_fortran_numbers
   * numbers and  number of column indicies (colptr array)
   */
  int extracted;

  FILE* file = fopen(filename,"rt");
  if (!file)
  {
    fprintf(stderr, "Unable to open file %s for reading\n",filename);
    return 0;
  }

  /*
   * Line 1.
   * TITLE, (72 characters)
   * KEY, (8 characters)
   */
  fgets(buf,HB_LINE_SIZE,file);
  /* skip them */

  /*
   * Line 2.
   * TOTCRD, integer, total number of data lines, (14 characters)
   * PTRCRD, integer, number of data lines for pointers, (14 characters)
   * INDCRD, integer, number of data lines for row or variable indices,
   *   (14 characters)
   * VALCRD, integer, number of data lines for numerical values of matrix
   *    entries, (14 characters)
   * RHSCRD, integer, number of data lines for right hand side vectors,
   *    starting guesses, and solutions, (14 characters)
   */
  fgets(buf,HB_LINE_SIZE,file);
  ptr = buf;
  totcrd = sp_extract_positional_int(ptr,14);
  ptr += 14;
  ptrcrd = sp_extract_positional_int(ptr,14);
  ptr += 14;
  indcrd = sp_extract_positional_int(ptr,14);
  ptr += 14;
  valcrd = sp_extract_positional_int(ptr,14);
  ptr += 14;
  rhscrd = sp_extract_positional_int(ptr,14);
  if (rhscrd)
    fprintf(stderr, "Right-part vector is not supported, ignoring\n");
  /*
   * Line 3.
   * MXTYPE, matrix type (see table), (3 characters)
   * blank space, (11 characters)
   * NROW, integer, number of rows or variables, (14 characters)
   * NCOL, integer, number of columns or elements, (14 characters)
   * NNZERO, integer, number of row or variable indices. For "assembled"
   *   matrices, this is just the number of nonzero entries. (14 characters)
   * NELTVL, integer, number of elemental matrix entries. For "assembled"
   *   matrices, this is 0. (14 characters)
   */
  fgets(buf,HB_LINE_SIZE,file);
  /* we support only real matrix */
  if (buf[0] != 'R')
  {
    fprintf(stderr,"Complex or Pattern matrix not supported\n");
    fclose(file);
    return 0;
  }
  if (buf[1] == 'H')
  {
    fprintf(stderr,"Complex Hermitian matrix not supported\n");
    fclose(file);
    return 0;
  }
  if (buf[2] == 'E')
  {
    fprintf(stderr,"Elemental matrix not supported\n");
    fclose(file);
    return 0;
  }
  ptr = buf + 14;
  nrow = sp_extract_positional_int(ptr,14);
  ptr += 14;
  ncol = sp_extract_positional_int(ptr,14);
  ptr += 14;
  nnzero = sp_extract_positional_int(ptr,14);
  /*
   * Line 4.
   * PTRFMT, FORTRAN I/O format for pointers, (16 characters)
   * INDFMT, FORTRAN I/O format for row or variable indices, (16 characters)
   * VALFMT, FORTRAN I/O format for matrix entries, (20 characters)
   * RHSFMT, FORTRAN I/O format for right hand sides, initial guesses, and
   *   solutions, (20 characters)
   */
  fgets(buf,HB_LINE_SIZE,file);
  ptr = buf;
  if (!hb_extract_positional_format(ptr,16,&ptrfmt))
  {
    fprintf(stderr,"Unknown format: %s",ptr);
    fclose(file);
    return 0;
  }
  ptr += 16;
  if (!hb_extract_positional_format(ptr,16,&indfmt))
  {
    fprintf(stderr,"Unknown format: %s",ptr);
    fclose(file);
    return 0;
  }
  ptr += 16;
  if (!hb_extract_positional_format(ptr,20,&valfmt))
  {
    fprintf(stderr,"Unknown format: %s",ptr);
    fclose(file);
    return 0;
  }
  if (rhscrd)
  {
    ptr += 20;
    if (!hb_extract_positional_format(ptr,20,&rhsfmt))
    {
      fprintf(stderr,"Unknown format: %s",ptr);
      fclose(file);
      return 0;
    }
  }
  /*
   * Line 5: (only present if 0 <RHSCRD!)
   * RHSTYP, describes the right hand side information, (3 characters)
   * blank space, (11 characters)
   * NRHS, integer, the number of right hand sides, (14 characters)
   * NRHSIX, integer, number of row indices, (14 characters)
   */
  if ( rhscrd )
    fgets(buf,HB_LINE_SIZE,file);

  /* header parsing done, parsing the data */
  
  /* Section 1. pointers */
  fortran_numbers = calloc(ptrfmt.repeat,sizeof(fortran_number));
  n = 0;
  while (ptrcrd --)
  {
    fgets(buf,HB_LINE_SIZE,file);
    ptr = (char*)sp_extract_fortran_numbers(buf,
                                            &ptrfmt,
                                            fortran_numbers,
                                            &extracted);
    if (ptr == &buf[0])
    {
      fprintf(stderr,"Unable to parse pointers: %s\n", buf);
      fclose(file);
      free(fortran_numbers);
      free(colptr);
      return 0;
    }
    colptr = colptr ? realloc(colptr, (n+extracted)*sizeof(int)) :
      calloc(n+extracted,sizeof(int));
    for (i = n; i < n + extracted; ++ i) 
      colptr[i] = fortran_numbers[i-n].integer;
    n += extracted;
  }
  free(fortran_numbers);
  if ( n != nrow + 1)
  {
    fprintf(stderr,"Unable to parse row indicies: parsed = %d != "
            "%d nonzeros\n", n, nnzero);
    fclose(file);
    free(colptr);
    free(rowind);
    return 0;
  }


  /* Section 2. rows */
  fortran_numbers  = calloc(indfmt.repeat, sizeof(fortran_number));
  n = 0;
  while (indcrd --)
  {
    fgets(buf,HB_LINE_SIZE,file);
    ptr = (char*)sp_extract_fortran_numbers(buf,
                                            &indfmt,
                                            fortran_numbers,
                                            &extracted);
    if (ptr == &buf[0])
    {
      fprintf(stderr,"Unable to parse row indicies: %s\n", buf);
      fclose(file);
      free(fortran_numbers);
      free(colptr);
      free(rowind);
      return 0;
    }
    
    rowind = rowind ? realloc(rowind, (n+extracted)*sizeof(int)) :
      calloc(n+extracted, sizeof(int));
    for (i = n; i < n + extracted; ++ i) 
      rowind[i] = fortran_numbers[i-n].integer;
    n += extracted;
  }
  free(fortran_numbers);
  if ( n != nnzero)
  {
    fprintf(stderr,"Unable to parse row indicies: parsed = %d != "
            "%d nonzeros\n", n, nnzero);
    fclose(file);
    free(colptr);
    free(rowind);
    return 0;
  }

  /* Section 3. values */
  fortran_numbers  = calloc(valfmt.repeat, sizeof(fortran_number));
  n = 0;
  while (valcrd --)
  {
    fgets(buf,HB_LINE_SIZE,file);
    ptr = (char*)sp_extract_fortran_numbers(buf,
                                            &valfmt,
                                            fortran_numbers,
                                            &extracted);
    if (ptr == &buf[0])
    {
      fprintf(stderr,"Unable to parse values: %s\n", buf);
      fclose(file);
      free(fortran_numbers);
      free(colptr);
      free(rowind);
      free(values);
      return 0;
    }
    
    values = values ? realloc(values, (n+extracted)*sizeof(double)) :
      calloc(n+extracted, sizeof(double));
    for (i = n; i < n + extracted; ++ i) 
      values[i] = fortran_numbers[i-n].real;
    n += extracted;
  }
  free(fortran_numbers);
  if ( n != nnzero)
  {
    fprintf(stderr,"Unable to parse values: parsed = %d != "
            "%d nonzeros\n", n, nnzero);
    fclose(file);
    free(colptr);
    free(rowind);
    free(values);
    return 0;
  }
  printf("So far HB file %s parsed successfully\n",filename);
  
  return self;
}

sp_matrix_ptr sp_matrix_load_file(const char* filename, int storage_type)
{
  sp_matrix_ptr self = 0;
  /* determine file extension */
  const char* ext = sp_parse_file_extension(filename);
  if (!ext)
  {
    fprintf(stderr,"File type is not supported: %s\n", filename);
    return self;
  }
  if (storage_type != CRS && storage_type != CCS)
  {
    fprintf(stderr, "Unknown storage type %d\n", storage_type);
    return self;
  }
  if ( !sp_istrcmp(ext,"mtx") )
    return sp_matrix_load_file_mm(filename, storage_type);
  else if (!sp_istrcmp(ext,"hb") ||
           !sp_istrcmp(ext,"rua") ||
           !sp_istrcmp(ext,"rsa") ||
           !sp_istrcmp(ext,"rza") ||
           !sp_istrcmp(ext,"rra"))
    return sp_matrix_load_file_hb(filename, storage_type);
  else
    fprintf(stderr,"File type is not supported: .%s\n", ext);

  return self;
}



static int sp_matrix_save_file_mm(sp_matrix_ptr self, const char* filename)
{
  int result = 1;
  FILE* file = fopen(filename,"wt+");
  int matrix_type;
  int i,j,n,nonzeros;
  int size;
  double value;
  /* MM format limitation for the line is 1024 characters */
  char buf[1024+1];
  if (!file)
  {
    fprintf(stderr,"Error opening file %s for writing",filename);
    return 0;
  }
  if ( !self->ordered )         /* order matrix */
    sp_matrix_compress(self);

  if (sp_matrix_issymmetric(self))
    matrix_type = MM_SYMMETRIC;
  else if (sp_matrix_isskew_symmetric(self))
    matrix_type = MM_SKEW_SYMMETRIC;
  else
    matrix_type = MM_GENERAL;

  size = sprintf(buf,"%s matrix coordinate real %s\n",
                 MM_HEADER_STRING,
                 matrix_type == MM_GENERAL ? "general" :
                 (matrix_type == MM_SYMMETRIC ? "symmetric" :
                  "skew-symmetric"));
  fwrite(buf, 1, size, file);
  fwrite(MM_HEADER_OUTPUT_ABOUT,1,strlen(MM_HEADER_OUTPUT_ABOUT),file);

  /* calculate nonzeros */
  nonzeros = 0;
  n = self->storage_type == CRS ? self->rows_count : self->cols_count;
  for ( i = 0; i < n; ++ i)
    nonzeros += self->storage[i].last_index + 1;

  /*
   * in skew symmetric matrix nzeros/2 elements in lower triangle
   * (diagonal always empty)
   */
  if ( matrix_type == MM_SKEW_SYMMETRIC)
    nonzeros  /= 2;
  /*
   * in symmetric matrix (nzeros-diagonal)/2 elements in lower triangle
   */
  if ( matrix_type == MM_SYMMETRIC)
  {
    j = 0;
    for ( i = 0; i < n; ++ i)
      if (sp_matrix_element_ptr(self,i,i))
        j++;
    nonzeros = (nonzeros + j)/2;
  }
  /* write sizes */
  size = sprintf(buf,"%d %d %d\n",self->rows_count,self->cols_count,nonzeros);
  fwrite(buf,1,size,file);

  /* write elements */
  for ( i = 0; i < n; ++ i)
  {
    for ( j = 0; j <= self->storage[i].last_index; ++ j)
    {
      if (matrix_type != MM_GENERAL && self->storage[i].indexes[j] > i)
        break;
      if (self->storage_type == CRS)
        size = sprintf(buf,"%d %d %.16e\n",i+1,self->storage[i].indexes[j]+1,
                       self->storage[i].values[j]);
      else                      /* CCS */
      {
        value = self->storage[i].values[j];
        /* lower triangle = -upper triangle for skew symmetic matix */
        if (matrix_type == MM_SKEW_SYMMETRIC)
          value = -value;
        if (matrix_type == MM_GENERAL)
          size = sprintf(buf,"%d %d %.16e\n",self->storage[i].indexes[j]+1,
                       i+1,
                       value);
        else                    /* transposed */
          size = sprintf(buf,"%d %d %.16e\n",i+1,self->storage[i].indexes[j]+1,
                         value);
      }
      fwrite(buf,1,size,file);
    }
  }

  fflush(file);
  fclose(file);
  return result;
}

/*
 * Save the sparse martix from the file.
 * File format guessed from the extension
 * Currently supported formats:
 * MM (matrix market), file extension .mtx
 * see http://math.nist.gov/MatrixMarket/formats.html#MMformat
 * Returns 0 if not possible to write(or unknown file format)
 */
int sp_matrix_save_file(sp_matrix_ptr self, const char* filename)
{
  /* determine file extension */
  const char* ext = sp_parse_file_extension(filename);
  if (!ext)
  {
    fprintf(stderr,"File type is not supported: %s\n", filename);
    return 0;
  }
  if ( !sp_istrcmp(ext,"mtx") )
    return sp_matrix_save_file_mm(self,filename);
  return 0;
}
