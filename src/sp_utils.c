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
#include <ctype.h>

#include "sp_utils.h"
#include "sp_mem.h"
#include "sp_log.h"


static int is_from(const char c, const char* from)
{
  if (from && c)
  {
    while (*from)
    {
      if (c == *from)
        return 1;
      from++;
    }
  }
  return 0;
}

/*
 * Case-insensitive ASCII 7bit string comparison
 */
int sp_istrcmp ( const char * str1, const char * str2 )
{
  int result = 0;
  while (*str1 && *str2 )
  {
    result = toupper(*str1) - toupper (*str2);
    if (result)
      return result;
    str1++,str2++;
  } 
  return toupper(*str1) - toupper (*str2);
}

/*
 * Returns the copy of n characters of the string
 * Caller shall free the returned string
 */
char* sp_strndup(const char *s, size_t n)
{
  char* str = spalloc(n+1);
  memcpy(str,s,n);
  str[n] = '\0';
  return str;
}

static const char* ptr_after_char(const char* str,char c)
{
  const char* ptr = strrchr(str, c);
  return ptr ? ptr + 1 : ptr;
}


/*
 * Extract file extension from the filename
 */
const char* sp_parse_file_extension(const char* filename)
{
  return ptr_after_char(filename,'.');
}

/*
 * Extract file name from the path
 */
const char* sp_parse_file_name(const char* filename)
{
  const char* ptr = ptr_after_char(filename,'/');
  return ptr ? ptr : filename;
}

/*
 * Extract file name without extension. Naive implementation
 * cannot handle full paths like /the/path.to/file
 * where the result will be /the/path
 */
void sp_parse_file_basename(const char* filename, char* basename)
{
  const char* ptr;
  size_t count, size;
  count = size = strlen(filename);
  ptr = filename + count;
  while (count)
  {
    ptr--,count--;
    if ( *ptr == '.' )
      break;
  }
  if (!count) 
    count = size;
  
  memcpy(basename,filename,count);
  basename[count] = 0;
  
}

/*
 * Skip whitespaces
 */
const char* sp_skip_whitespaces(const char* line)
{
  while (*line && (*line == ' ' || *line == '\t'))
    line++;
  return line;
}

/*
 * Skip alphanumeric characters and characters from chars array
 */
const char* sp_skip_alnum(const char* line, const char* chars)
{
  while (*line && (isalnum(*line) || is_from(*line,chars)))
    line++;
  return line;
}

/*
 * Extract word from the string
 * return pointer to the word into the word argument
 * Caller shall free the word
 */
const char* sp_extract_next_word(const char* line, const char** word)
{
  const char* end;
  /* skip whitespaces to the object type */
  line = sp_skip_whitespaces(line);
  /* extract word */
  end = sp_skip_alnum(line,"-");
  *word = sp_strndup(line,end - line);
  return end;
}


/*
 * Read the text file and return buffer with the file's contens
 * returns 0 if unable to read
 */
char* sp_read_text_file(const char* filename)
{
  FILE* file;
  const int block_size = 1024;
  /* contents buffer */
  char* contents = spcalloc(block_size+1,1);
  /* auxulary counters */
  size_t read_chunk = 0,read = 0;

  file = fopen(filename,"rt");
  if (!file)
  {
    LOGERROR("Cannot read file %s", filename);
    spfree(contents);
    return 0;
  }
  /* read file contents  */
  while(!feof(file))
  {
    read_chunk = fread(contents+read,1,block_size, file);
    read += read_chunk;
    
    contents = (char*)sprealloc(contents,read + block_size);
  }
  contents[read] = '\0';
  /* close the file */
  fclose(file);
  return contents;
}


/* Extracts the integer of size bytes from the buffer from */
int sp_extract_positional_int(const char* from, size_t size)
{
  int result;
  char* buf = spalloc(size+1);
  char* ptr = buf;
  size_t i = 0;
  while (i < size && *from)
  {
    *ptr++ = *from++;
    ++i;
  }
  buf[size] = '\0';
  result = atoi(buf);
  spfree(buf);
  return result;
}

/* Extracts the float of size bytes from the buffer from */
double sp_extract_positional_float(const char* from, int size)
{
  double result;
  char* buf = spcalloc(size+1,1);
  char* ptr = buf;
  int i = 0;
  while (i < size && *from)
  {
    *ptr++ = *from++;
    ++i;
  }
  buf[size] = '\0';
  ptr = buf;
  /* try to parse strings like 0.156211824D+02 */
  while (*ptr) { if(*ptr == 'D') *ptr = 'E'; ++ptr; }
  result = atof(buf);
  spfree(buf);
  return result;  
}



/*
 * Very simple Parser for FORTRAN IO Format specifiers
 * Subset of Fortran IO format grammar BNF specification
 * (taken from http://bebop.cs.berkeley.edu/smc/hb.html)
 * <format> --> \(<format-string>\)
 * <format-string> --> <count>?<rest>
 * <count> --> \d+
 * <rest> --> <fixedid>|<intid>|<fltid>|<doubleid>|<generalid>
 * <fixedid> --> F<field-width>\.<digits-after-decimal-point>
 * <intid> --> I<field-width>(\.<min-num-digits>)?
 * <fltid> --> E<field-width>\.<decimal-significand-length>(E<num-digits-in-exponent>)?
 * <doubleid> --> D<field-width>\.<decimal-significand-length>(E<num-digits-in-exponent>)?
 * <generalid> --> G<field-width>\.<decimal-significand-length>(E<num-digits-in-exponent>)?
 * <field-width> --> \d+
 * <min-num-digits> --> \d+  (defaults to zero)
 * <num-digits-in-exponent> --> \d+
 * <decimal-significand-length> --> \d+
 * <digits-after-decimal-point> --> \d+
 */
int sp_parse_fortran_format(const char* string, fortran_io_format* format)
{
  const char* ptr = sp_skip_whitespaces(string);
  const char* end;
  const char* ptr1;
  memset(format,0,sizeof(fortran_io_format));
  format->repeat = 1;
  /* <format> --> \(<format-string>\) */
  /* find start \( */
  if (!*ptr || *ptr != '(')
  {
    LOGERROR("Incorrect FORTRAN IO format: %s, no '('",string);
    return 0;
  }
  /* find end \) */
  end = ptr;
  while(*end && *end != ')') ++end;
  if (!*end  || *end != ')')
  {
    LOGERROR("Incorrect FORTRAN IO format: %s, no ')'",string);
    return 0;
  }
  ptr++;
  /* <format-string> --> <count>?<rest> */
  ptr1 = ptr;
  /* <count> --> \d+ */
  while (isdigit(*ptr1) && ptr1 != end) ++ptr1;
  format->repeat = ptr1 > ptr ? sp_extract_positional_int(ptr,ptr1-ptr) : 0;
  /* <rest> --> <fixedid>|<intid>|<fltid>|<doubleid>|<generalid> */
  ptr = ptr1;
  if (!is_from(*ptr1,"FIEDGfiedg"))
  {
    LOGERROR("Incorrect FORTRAN IO format: %s, %c != [FIEDG]",
            string,*ptr1);
    return 0;
  }
  format->type = toupper(*ptr);
  ptr1 = ++ptr;
  /* [FIEDG]\d+ */
  while(isdigit(*ptr1) && ptr1 != end) ++ptr1;
  if (ptr1 == ptr)
  {
    LOGERROR("Incorrect FORTRAN IO format: %s, field-width = ''",
             string);
    return 0;
  }
  format->width = sp_extract_positional_int(ptr,ptr1-ptr);
  
  /* <intid> --> I<field-width>(\.<min-num-digits>)? */
  if (format->type == 'I')
  {
    if (ptr1 == end)
      return 1;
    if (*ptr1 != '.')
    {
      LOGERROR("Incorrect FORTRAN IO format: %s, type = 'I'",
               string);
      return 0;
    }
    ptr1++;
    ptr = ptr1;
    while(isdigit(*ptr1) && ptr1 != end) ++ptr1;
    if (ptr1 == ptr)
    {
      LOGERROR("Incorrect FORTRAN IO format: %s, type = 'I'",
               string);
      return 0;
    }
    format->num1 = sp_extract_positional_int(ptr,ptr1-ptr);
  }
  /* <fixedid> --> F<field-width>\.<digits-after-decimal-point> */
  else if (format->type == 'F')
  {
    if (*ptr1 != '.')
    {
      LOGERROR("Incorrect FORTRAN IO format: %s, type = 'F'",
               string);
      return 0;
    }
    ptr1++;
    ptr = ptr1;
    while(isdigit(*ptr1) && ptr1 != end) ++ptr1;
    if (ptr1 == ptr)
    {
      LOGERROR("Incorrect FORTRAN IO format: %s, type = 'F'",
               string);
      return 0;
    }
    format->num1 = sp_extract_positional_int(ptr,ptr1-ptr);
  }
  else                          /* E,D,G */
  {
    if (*ptr1 != '.')
    {
      LOGERROR("Incorrect FORTRAN IO format: %s, type = '%c'",
               string,format->type);
      return 0;
    }
    ptr1++;
    ptr = ptr1;
    while(isdigit(*ptr1) && ptr1 != end) ++ptr1;
    if (ptr1 == ptr)
    {
      LOGERROR("Incorrect FORTRAN IO format: %s, type = '%c'",
               string,format->type);
      return 0;
    }
    format->num1 = sp_extract_positional_int(ptr,ptr1-ptr);
    /* (E<num-digits-in-exponent>)? */
    if (ptr1 != end)
    {
      if (*ptr1 != 'E')
      {
        LOGERROR("Incorrect FORTRAN IO format: %s, type = '%c'",
                 string,format->type);
        return 0;
      }
      ptr1++;
      ptr = ptr1;
      while(isdigit(*ptr1) && ptr1 != end) ++ptr1;
      if (ptr1 == ptr)
      {
        LOGERROR("Incorrect FORTRAN IO format: %s, type = '%c'",
                 string,format->type);
        return 0;
      }
      format->num2 = sp_extract_positional_int(ptr,ptr1-ptr);
    }
  }
  if ( ptr1 != end)
  {
    LOGERROR("Incorrect FORTRAN IO format: %s, type = '%c'",
             string,format->type);
    return 0;
  }
  return 1;
}


static const char* sp_extract_fortran_numbers_I(const char* string, 
                                                const fortran_io_format* format,
                                                /* output */
                                                fortran_number* number,
                                                int* extracted)
{
  const char* result = string;
  int size = 0;
  const char* ptr = string, *ptr1;
  int i;
  /* width */
  while (*ptr && *ptr != '\n')
    size++,ptr++;
  /* if (size < format->repeat * format->width) */
  /*   return string; */
  ptr = string;
  *extracted = 0;
  for ( i = 0; i < format->repeat; ++ i)
  {
    ptr1 = sp_skip_whitespaces(ptr);
    if (*ptr1 && *ptr1 != '\n' && ptr1-string < size &&
        ptr1-ptr < format->width)
    {
      number[*extracted].integer =
        sp_extract_positional_int(ptr,format->width);
      number[*extracted].fortran_type = FORTRAN_TYPE_INTEGER;
      (*extracted) ++;
    }
    else
      break;
    ptr += format->width;
  }
  result = string + size;
  return result;
}


static
const char* sp_extract_fortran_numbers_F(const char* string, 
                                         const fortran_io_format* format,
                                         /* output */
                                         fortran_number* number,
                                         int* extracted)
{
  const char* result = string;
  int size = 0;
  const char* ptr = string, *ptr1;
  int i;
  /* width */
  while (*ptr && *ptr != '\n')
    size++,ptr++;
  /* if (size < format->repeat * format->width) */
  /*   return string; */
  ptr = string;
  *extracted = 0;
  for ( i = 0; i < format->repeat; ++ i)
  {
    ptr1 = sp_skip_whitespaces(ptr);
    if (*ptr1 && *ptr1 != '\n' && ptr1-string < size &&
        ptr1-ptr < format->width)
    {
      number[*extracted].real =
        sp_extract_positional_float(ptr,format->width);
      number[*extracted].fortran_type = FORTRAN_TYPE_DOUBLE;
      (*extracted) ++;
    }
    else
      break;
    ptr += format->width;
  }
  result = string + size;
  return result;
}


static
const char* sp_extract_fortran_numbers_ED(const char* string, 
                                          const fortran_io_format* format,
                                          /* output */
                                          fortran_number* number,
                                          int* extracted)
{
  return sp_extract_fortran_numbers_F(string,format,number,extracted);
}



static
const char* sp_extract_fortran_numbers_G(const char* string, 
                                         const fortran_io_format* format,
                                         /* output */
                                         fortran_number* number,
                                         int* extracted)
{
  return sp_extract_fortran_numbers_F(string,format,number,extracted);
}


/*
 * Extract number specified by fortran format
 * Returns the same string pointer in case of error
 * See for example
 * http://cpan.uwinnipeg.ca/htdocs/Fortran-Format/Fortran/Format.pm.html#I_i_w_i
 */
const char* sp_extract_fortran_numbers(const char* string, 
                                       const fortran_io_format* format,
                                       /* output */
                                       fortran_number* number,
                                       int* extracted)
{
  const char* result = string;
  if (format->type == 'I')
    result = sp_extract_fortran_numbers_I(string,format,number,extracted);
  else if (format->type == 'F')
    result = sp_extract_fortran_numbers_F(string,format,number,extracted);
  else if (format->type == 'E' || format->type == 'D')
    result = sp_extract_fortran_numbers_ED(string,format,number,extracted);
  else if (format->type == 'G')
    result = sp_extract_fortran_numbers_G(string,format,number,extracted);
  return result;
}

