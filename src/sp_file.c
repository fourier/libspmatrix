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
#include <assert.h>

#include "sp_file.h"
#include "sp_mem.h"
#include "sp_utils.h"
#include "sp_log.h"
#include "sp_cont.h"

/*
 * Enum - supported export formats
 */
typedef enum
{
  FMT_MM,
  FMT_TXT,
  FMT_DAT,
  FMT_UNSUPPORTED
} supported_format;
/*
 * Data desrtiption for Matrix Market file format
 */

const char* MM_HEADER_STRING = "%%MatrixMarket";
const char* MM_HEADER_OUTPUT_ABOUT = "% Created by libspmatrix (c) "
  "Alexey Veretennikov, https://github.com/fourier/libspmatrix\n";
const char* DAT_HEADER_OUTPUT_ABOUT = "# Created by libspmatrix (c) "
"Alexey Veretennikov, https://github.com/fourier/libspmatrix\n";
const char* DAT_HEADER_FORMAT = "# name: %s\n"
"# type: sparse matrix\n"
"# nnz: %d\n"
"# rows: %d\n"
"# columns: %d\n";

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
      LOGERROR("Error reading MatrixMarket header!");
      return contents;
    }
    header_ptr++, ptr++;
  }
  /* extract object type */
  ptr = sp_extract_next_word(ptr,&word);
  header->object = !sp_istrcmp(word,"matrix") ? MM_MATRIX : MM_OTHER;
  spfree((char*)word);

  /* extract storage type */
  ptr = sp_extract_next_word(ptr,&word);
  if (!sp_istrcmp(word,"coordinate"))
    header->storage = MM_COORDINATE;
  else if (!sp_istrcmp(word,"array"))
    header->storage = MM_ARRAY;
  else
  {
    LOGERROR("Unkown storage type: %s", word);
    ptr = contents;
  }
  spfree((char*)word);

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
    LOGERROR("Unkown elements type: %s", word);
    ptr = contents;
  }
  spfree((char*)word);

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
    LOGERROR("Unkown portrait type: %s", word);
    ptr = contents;
  }
  spfree((char*)word);

  return ptr;
}

/* return 1 if mm-matrix format is supported, 0 otherwise  */
static int mm_validate_header(mm_header* header)
{
  if (header->object != MM_MATRIX)
  {
    LOGERROR("MM objects other than matrix are not supported");
    return 0;
  }
  if (header->storage != MM_COORDINATE)
  {
    LOGERROR("Currently only coordinate format of MM files supported");
    return 0;
  }
  if (header->elements == MM_COMPLEX)
  {
    LOGERROR("Complex matricies not supported");
    return 0;
  }
  if (header->portrait == MM_HERMITIAN)
  {
    LOGERROR("Hermitian matricies not supported");
    return 0;
  }
  return 1;
}

/*
 * Load matrix in MatrixMarket format
 */
static int sp_matrix_yale_load_file_mm(sp_matrix_yale_ptr self,
                                       const char* filename,
                                       sparse_storage_type type)
{
  int result = 0;
  sp_matrix mtx;
  mm_header header;
  const char* ptr;
  const char* line;
  int block_number = 0;
  int rows, cols, nonzeros = 0;
  int element_number = -1;
  int i,j;
  double value;
  char* contents;
  /* clear the structures */
  memset(&mtx,0,sizeof(mtx));
  memset(&header,0,sizeof(header));
  /* read the file contents  */
  contents = sp_read_text_file(filename);
  if (!contents)
    return 0;

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
        LOGERROR("Not supported matrix type");
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
          LOGERROR("Unable to parse sizes");
          break;
        }
        sp_matrix_init(&mtx,rows,cols,nonzeros/rows + 1,type);
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
            LOGERROR("Unable to parse line %s",line);
            break;
          }
          value = 1;
        }
        else
        {
          if (sscanf(line, "%d %d %lg", &i, &j, &value) != 3)
          {
            LOGERROR("Unable to parse line %s",line);
            break;
          }

        }
        sp_matrix_element_add(&mtx,i-1,j-1,value);
        /* handle symmetry property */
        if ( i != j )
        {
          if ( header.portrait == MM_SYMMETRIC )
            sp_matrix_element_add(&mtx,j-1,i-1,value);
          else if (header.portrait == MM_SKEW_SYMMETRIC)
            sp_matrix_element_add(&mtx,j-1,i-1,-value);
        }
        ++ element_number;
      }
    }
    line = strtok(NULL,"\n");
  }

  /* free resources */
  spfree(contents);
  /* check for error */
  if ( element_number != nonzeros)
  {
    LOGERROR("Error loading matrix, expected %d nonzeros, parsed %d",
             nonzeros, element_number);
    sp_matrix_free(&mtx);
  }
  else
  {
    sp_matrix_yale_init(self,&mtx);
    sp_matrix_free(&mtx);
    result = 1;
  }
  return result;
}

static int hb_extract_positional_format(const char* from, int size,
                                        fortran_io_format* fmt)
{
  int result = 0;
  char* buf = spalloc(size+1);
  char* ptr = buf;
  int i = 0;
  while (i < size && *from)
  {
    *ptr++ = *from++;
    ++i;
  }
  buf[size] = '\0';
  result = sp_parse_fortran_format(buf,fmt);
  spfree(buf);
  return result;
}

/*
 * Load matrix in Harwell Boeing format
 */
static int sp_matrix_yale_load_file_hb(sp_matrix_yale_ptr self,
                                       const char* filename)
{
  int i,j,p,n;
  sp_matrix mtx;
  matrix_properties props = PROP_GENERAL;
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
    LOGERROR("Unable to open file %s for reading",filename);
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
  if (totcrd != ptrcrd + indcrd + valcrd + rhscrd)
  {
    LOGERROR("Load file in HB format: total number of data lines %d not equal to %d+%d+%d+%d = %d",
             totcrd, ptrcrd, indcrd, valcrd, rhscrd, ptrcrd + indcrd + valcrd + rhscrd);
    return 0;
  }
  if (rhscrd)
    LOGWARN("Right-part vector is not supported, ignoring");
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
    LOGERROR("Complex or Pattern matrix not supported");
    fclose(file);
    return 0;
  }
  if (buf[1] == 'H')
  {
    LOGERROR("Complex Hermitian matrix not supported");
    fclose(file);
    return 0;
  }
  if (buf[2] == 'E')
  {
    LOGERROR("Elemental matrix not supported");
    fclose(file);
    return 0;
  }
  switch(buf[1])
  {
  case 'S': props = PROP_SYMMETRIC; break;
  case 'Z': props = PROP_SKEW_SYMMETRIC; break;
  default:  props = PROP_GENERAL; break;
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
    LOGERROR("Unknown format: %s",ptr);
    fclose(file);
    return 0;
  }
  ptr += 16;
  if (!hb_extract_positional_format(ptr,16,&indfmt))
  {
    LOGERROR("Unknown format: %s",ptr);
    fclose(file);
    return 0;
  }
  ptr += 16;
  if (!hb_extract_positional_format(ptr,20,&valfmt))
  {
    LOGERROR("Unknown format: %s",ptr);
    fclose(file);
    return 0;
  }
  if (rhscrd)
  {
    ptr += 20;
    if (!hb_extract_positional_format(ptr,20,&rhsfmt))
    {
      LOGERROR("Unknown format: %s",ptr);
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
  fortran_numbers = spcalloc(ptrfmt.repeat,sizeof(fortran_number));
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
      LOGERROR("Unable to parse pointers: %s", buf);
      fclose(file);
      spfree(fortran_numbers);
      spfree(colptr);
      return 0;
    }
    colptr = colptr ? sprealloc(colptr, (n+extracted)*sizeof(int)) :
      spcalloc(n+extracted,sizeof(int));
    for (i = n; i < n + extracted; ++ i) 
      colptr[i] = fortran_numbers[i-n].integer;
    n += extracted;
  }
  spfree(fortran_numbers);
  if ( n != nrow + 1 || !colptr)
  {
    LOGERROR("Unable to parse row indicies: parsed = %d != "
             "%d nonzeros", n, nnzero);
    fclose(file);
    spfree(colptr);
    spfree(rowind);
    return 0;
  }

  if (colptr[nrow] != nnzero+1)
  {
    LOGERROR("Unable to parse row indicies: last index = %d != "
             "%d nonzeros", colptr[nrow], nnzero);
    fclose(file);
    spfree(colptr);
    spfree(rowind);
    return 0;
  }


  /* Section 2. rows */
  fortran_numbers  = spcalloc(indfmt.repeat, sizeof(fortran_number));
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
      LOGERROR("Unable to parse row indicies: %s", buf);
      fclose(file);
      spfree(fortran_numbers);
      spfree(colptr);
      spfree(rowind);
      return 0;
    }
    
    rowind = rowind ? sprealloc(rowind, (n+extracted)*sizeof(int)) :
      spcalloc(n+extracted, sizeof(int));
    for (i = n; i < n + extracted; ++ i) 
      rowind[i] = fortran_numbers[i-n].integer;
    n += extracted;
  }
  spfree(fortran_numbers);
  if ( n != nnzero)
  {
    LOGERROR("Unable to parse row indicies: parsed = %d != "
             "%d nonzeros", n, nnzero);
    fclose(file);
    spfree(colptr);
    spfree(rowind);
    return 0;
  }

  /* Section 3. values */
  fortran_numbers  = spcalloc(valfmt.repeat, sizeof(fortran_number));
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
      LOGERROR("Unable to parse values: %s", buf);
      fclose(file);
      spfree(fortran_numbers);
      spfree(colptr);
      spfree(rowind);
      spfree(values);
      return 0;
    }
    
    values = values ? sprealloc(values, (n+extracted)*sizeof(double)) :
      spcalloc(n+extracted, sizeof(double));
    for (i = n; i < n + extracted; ++ i) 
      values[i] = fortran_numbers[i-n].real;
    n += extracted;
  }
  spfree(fortran_numbers);
  fclose(file);
  if ( n != nnzero)
  {
    LOGERROR("Unable to parse values: parsed = %d != "
            "%d nonzeros", n, nnzero);
    spfree(colptr);
    spfree(rowind);
    spfree(values);
    return 0;
  }
  LOGINFO("So far HB file %s parsed successfully",filename);
  /* all indicies are 1 based in HB format */
  for ( i = 0; i < ncol+1; ++ i)
    colptr[i]--;
  for ( i = 0; i < nnzero; ++ i)
  {
    rowind[i]--;
    assert(rowind[i] < nrow);
  }

  if (props == PROP_GENERAL)
  {
    self->rows_count = nrow;
    self->cols_count = ncol;
    self->nonzeros   = nnzero;
    self->storage_type = CCS;
    self->offsets = colptr;
    self->indicies = rowind;
    self->values = values;
  }
  else                          /* symmetric and skew-symmetic */
  {
    /* easiest way is to construct the new matrix, since for
     * the symmetric and skew-symmetric matricies only one half of
     * elements stored
     */
    /* bandwidth */
    n = nnzero/ncol + 1;
    sp_matrix_init(&mtx,nrow,ncol,n+1,CCS);
    for (i = 0; i < ncol; ++ i)
    {
      for (p = colptr[i]; p < colptr[i+1]; ++ p)
      {
        j = rowind[p];
        assert(j < nrow);
        MTX(&mtx,j,i,values[p]);
        if (i != j)
          MTX(&mtx,i,j,(props == PROP_SYMMETRIC ? values[p] : -values[p]));
      }
    }
    sp_matrix_yale_init(self,&mtx);
    sp_matrix_free(&mtx);
    spfree(colptr);
    spfree(rowind);
    spfree(values);
  }
  /*  */
  return 1;
}

int sp_matrix_yale_load_file(sp_matrix_yale_ptr self,
                             const char* filename,
                             sparse_storage_type type)
{
  /* determine file extension */
  const char* ext = sp_parse_file_extension(filename);
  if (!ext)
  {
    LOGERROR("File type is not supported: %s", filename);
    return 0;
  }
  if ( !sp_istrcmp(ext,"mtx") )
    return sp_matrix_yale_load_file_mm(self, filename,type);
  else if (!sp_istrcmp(ext,"hb") ||
           !sp_istrcmp(ext,"rua") ||
           !sp_istrcmp(ext,"rsa") ||
           !sp_istrcmp(ext,"rza") ||
           !sp_istrcmp(ext,"rra"))
    return sp_matrix_yale_load_file_hb(self, filename);
  else
    LOGERROR("File type is not supported: *.%s", ext);

  return 0;
}


static int sp_matrix_save_file_triplet(sp_matrix_ptr self,
                                       FILE* file,
                                       int matrix_type,
                                       int base)
{
  int result = 1;
  int i,j,n,size;
  double value;
  char buf[64];

  n = self->storage_type == CRS ? self->rows_count : self->cols_count;
  
  /* write elements */
  for ( i = 0; i < n; ++ i)
  {
    for ( j = 0; j <= self->storage[i].last_index; ++ j)
    {
      if (matrix_type != MM_GENERAL && self->storage[i].indexes[j] > i)
      {
        break;
      }
      if (self->storage_type == CRS)
        size = sprintf(buf,"%d %d %.16e\n",i+base,
                       self->storage[i].indexes[j]+base,
                       self->storage[i].values[j]);
      else                      /* CCS */
      {
        value = self->storage[i].values[j];
        /* lower triangle = -upper triangle for skew symmetic matix */
        if (matrix_type == MM_SKEW_SYMMETRIC)
          value = -value;
        if (matrix_type == MM_GENERAL)
          size = sprintf(buf,"%d %d %.16e\n",self->storage[i].indexes[j]+base,
                       i+base,
                       value);
        else                    /* transposed */
          size = sprintf(buf,"%d %d %.16e\n",i+base,
                         self->storage[i].indexes[j]+base,
                         value);
      }
      fwrite(buf,1,size,file);
    }
  }
  return result;
}

static int sp_matrix_type_mm(sp_matrix_ptr self)
{
  int matrix_type = MM_GENERAL;
  matrix_properties props;
  if ( !self->ordered )         /* order matrix */
    sp_matrix_reorder(self);
  props = sp_matrix_properites(self);
  switch(props)
  {
  case PROP_SYMMETRIC:
    matrix_type = MM_SYMMETRIC; break;
  case PROP_SKEW_SYMMETRIC:
    matrix_type = MM_SKEW_SYMMETRIC; break;
  case PROP_GENERAL:
  case PROP_SYMMETRIC_PORTRAIT:
  default:
    break;
  }
  return matrix_type;
}

static int sp_matrix_save_file_mm(sp_matrix_ptr self, const char* filename)
{
  int result = 1;
  FILE* file = fopen(filename,"wt+");
  int matrix_type;
  int i,j,n,nonzeros;
  int size;
  /* MM format limitation for the line is 1024 characters */
  char buf[1024+1];
  if (!file)
  {
    LOGERROR("Error opening file %s for writing",filename);
    return 0;
  }
  matrix_type = sp_matrix_type_mm(self);
  
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

  if (!sp_matrix_save_file_triplet(self,file,matrix_type,1))
  {
    LOGERROR("Cannot save file!");
    result = 0;
  }
  
  fflush(file);
  fclose(file);
  return result;
}


static int sp_matrix_save_file_txt(sp_matrix_ptr self, const char* filename)
{
  int result = 1;
  int matrix_type;
  FILE* file = fopen(filename,"wt+");
  if (!file)
  {
    LOGERROR("Error opening file %s for writing",filename);
    return 0;
  }
  matrix_type = MM_GENERAL; /* sp_matrix_type_mm(self); */
  if (!sp_matrix_save_file_triplet(self,file,matrix_type,0))
  {
    LOGERROR("Cannot save file!");
    result = 0;
  }
  fflush(file);
  fclose(file);
  return result;
}


static supported_format guess_export_format(const char* filename)
{
  /* determine file extension */
  const char* ext = sp_parse_file_extension(filename);
  if (!ext)
  {
    LOGERROR("File type is not supported: %s", filename);
    return 0;
  }
  if ( !sp_istrcmp(ext,"mtx") )
    return FMT_MM;
  else if ( !sp_istrcmp(ext,"txt") )
    return FMT_TXT;
  else if ( !sp_istrcmp(ext, "dat") )
    return FMT_DAT;
  return FMT_UNSUPPORTED;
}


int sp_matrix_save_file(sp_matrix_ptr self, const char* filename)
{
  supported_format fmt = guess_export_format(filename);
  switch(fmt)
  {
  case FMT_MM:  return sp_matrix_save_file_mm(self,filename);
  case FMT_TXT: return sp_matrix_save_file_txt(self,filename);
  case FMT_UNSUPPORTED:
  default:
    break;
  }
  return 0;
}

static int sp_matrix_yale_save_file_triplet(sp_matrix_yale_ptr self,
                                            FILE* file,
                                            int matrix_type,
                                            int base)
{
  int result = 1;
  int i,p,n,size;
  double value;
  char buf[64];
  
  n = self->storage_type == CRS ? self->rows_count : self->cols_count;
  
  /* write elements */
  for ( i = 0; i < n; ++ i)
  {
    for ( p = self->offsets[i]; p < self->offsets[i+1]; ++ p )
    {
      if (matrix_type != MM_GENERAL && self->indicies[p] > i)
      {
        break;
      }
      if (self->storage_type == CRS)
        size = sprintf(buf,"%d %d %.16e\n",i+base,
                       self->indicies[p]+base,
                       self->values[p]);
      else                      /* CCS */
      {
        value = self->values[p];
        /* lower triangle = -upper triangle for skew symmetic matix */
        if (matrix_type == MM_SKEW_SYMMETRIC)
          value = -value;
        if (matrix_type == MM_GENERAL)
          size = sprintf(buf,"%d %d %.16e\n",self->indicies[p]+base,
                         i+base,
                         value);
        else                    /* transposed */
          size = sprintf(buf,"%d %d %.16e\n",i+base,
                         self->indicies[p]+base,
                         value);
      }
      fwrite(buf,1,size,file);
    }
  }
  return result;
}

static
int sp_matrix_yale_save_file_mm(sp_matrix_yale_ptr self,const char* filename)
{
  int result = 1;
  FILE* file = fopen(filename,"wt+");
  int matrix_type;
  matrix_properties props;
  int i,j,p,n,nonzeros;
  int size;
  /* MM format limitation for the line is 1024 characters */
  char buf[1024+1];
  if (!file)
  {
    fprintf(stderr,"Error opening file %s for writing",filename);
    return 0;
  }

  n = self->storage_type == CRS ? self->rows_count : self->cols_count;
  
  props = sp_matrix_yale_properites(self);
  switch (props)
  {
  case PROP_SYMMETRIC:
    matrix_type = MM_SYMMETRIC; break;
  case PROP_SKEW_SYMMETRIC:
    matrix_type = MM_SKEW_SYMMETRIC; break;
  case PROP_SYMMETRIC_PORTRAIT:
  case PROP_GENERAL:
  default:
    matrix_type = MM_GENERAL;
  }
  
  size = sprintf(buf,"%s matrix coordinate real %s\n",
                 MM_HEADER_STRING,
                 matrix_type == MM_GENERAL ? "general" :
                 (matrix_type == MM_SYMMETRIC ? "symmetric" :
                  "skew-symmetric"));
  fwrite(buf, 1, size, file);
  fwrite(MM_HEADER_OUTPUT_ABOUT,1,strlen(MM_HEADER_OUTPUT_ABOUT),file);

  /* calculate nonzeros */
  nonzeros = self->nonzeros;
  /*
   * in skew symmetric matrix nzeros/2 elements in lower triangle
   * (diagonal always empty)
   */
  if ( matrix_type == MM_SKEW_SYMMETRIC)
    nonzeros  /= 2;
  /*
   * in symmetric matrix (nzeros+diagonal)/2 elements in lower triangle
   */
  if ( matrix_type == MM_SYMMETRIC)
  {
    j = 0;
    for ( i = 0; i < n; ++ i)
    {
      for (p = self->offsets[i]; p < self->offsets[i+1]; ++ p)
      {
        if (self->indicies[p] > i)
          break;
        if (self->indicies[p] == i)
        {
          j++;
          break;
        }
      }
    }
    nonzeros = (nonzeros + j)/2;
  }
  
  /* write sizes */
  size = sprintf(buf,"%d %d %d\n",self->rows_count,self->cols_count,nonzeros);
  fwrite(buf,1,size,file);

  if (!sp_matrix_yale_save_file_triplet(self,file,matrix_type,1))
  {
    LOGERROR("Cannot save file!");
    result = 0;
  }
  
  fflush(file);
  fclose(file);
  return result;
}

static
int sp_matrix_yale_save_file_txt(sp_matrix_yale_ptr self,const char* filename)
{
  int result = 1;
  int matrix_type;
  FILE* file = fopen(filename,"wt+");
  if (!file)
  {
    LOGERROR("Error opening file %s for writing",filename);
    return 0;
  }
  matrix_type = MM_GENERAL; /* sp_matrix_type_mm(self); */
  if (!sp_matrix_yale_save_file_triplet(self,file,matrix_type,0))
  {
    LOGERROR("Cannot save file!");
    result = 0;
  }
  fflush(file);
  fclose(file);
  return result;
}

static
int sp_matrix_yale_save_file_dat(sp_matrix_yale_ptr self,const char* filename)
{
  int result = 1;
  int matrix_type;
  char buf[1024];
  char* name = spalloc(strlen(filename)+1);
  const char* ptr = filename;
  int size;
  FILE* file = fopen(filename,"wt+");
  if (!file)
  {
    LOGERROR("Error opening file %s for writing",filename);
    spfree(name);
    return 0;
  }
  fwrite(DAT_HEADER_OUTPUT_ABOUT,1,strlen(DAT_HEADER_OUTPUT_ABOUT),file);
  ptr = sp_parse_file_name(filename);
  sp_parse_file_basename(ptr, name);
  
  size = sprintf(buf, DAT_HEADER_FORMAT,name,
                 self->nonzeros,self->rows_count,self->cols_count);
  fwrite(buf, 1, size, file);
  spfree(name);
  matrix_type = MM_GENERAL;
  if (!sp_matrix_yale_save_file_triplet(self,file,matrix_type,0))
  {
    LOGERROR("Cannot save file!");
    result = 0;
  }
  fflush(file);
  fclose(file);
  return result;
}

int sp_matrix_yale_save_file(sp_matrix_yale_ptr self, const char* filename)
{
  supported_format fmt = guess_export_format(filename);
  switch(fmt)
  {
  case FMT_MM:  return sp_matrix_yale_save_file_mm(self,filename);
  case FMT_TXT: return sp_matrix_yale_save_file_txt(self,filename);
  case FMT_DAT: return sp_matrix_yale_save_file_dat(self,filename);
  case FMT_UNSUPPORTED:
  default:
    break;
  }
  return 0;
}

void sp_save_int_vector(int* v, int size, const char* fname)
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

void sp_save_double_vector(double* v, int size, const char* fname)
{
  int i = 0;
  FILE *f = fopen(fname,"wt+");
  if (f)
  {
    for (; i < size; ++ i)
      fprintf(f,"%e\n", v[i]);
    fclose(f);
  }
}

int sp_load_int_vector(int** v, int* size, const char* fname)
{
  int result = 0;
  int i = 0;
  int x;
  FILE *f = fopen(fname,"rt");
  int_array arr;
  if (f)
  {
    int_array_init(&arr, 10, 10);
    while (fscanf(f, "%d\n", &x) == 1)
    {
      int_array_add(&arr, x);
      i++;
    }
    *v = arr.items;
    *size = i;
    fclose(f);
    result = 1;
  }
  return result;
}


int sp_load_double_vector(double** v, int* size, const char* fname)
{
  return 0;
}