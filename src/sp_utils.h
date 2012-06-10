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

#ifndef _SP_UTILS_H_
#define _SP_UTILS_H_

#include <string.h>
#include <math.h>

/*
 * Simple equals macro for double values. Assuming values are always > e-16
 * For additional information, see
 * http://www.rsdn.ru/forum/cpp/2640596.1.aspx
 */
#ifndef DBL_EPSILON
#define DBL_EPSILON 2.2204460492503131e-16
#endif
#define FMAX(x,y) ((x) > (y) ? (x) : (y))
#define EQL(x,y) ((fabs((x)-(y))<= DBL_EPSILON*2) ? 1:0)

enum
{
  FORTRAN_TYPE_INTEGER,
  FORTRAN_TYPE_DOUBLE
};

typedef struct 
{
  char type;
  int repeat;                   /* repeat field */
  int width;                    /* field-width field */
  int num1;                     /* min-num-digits or decimal-significand
                                 * or digits-after-decimal-point
                                 */
  int num2;                     /* num-digits-in-exponent */
} fortran_io_format;

typedef struct
{
  int fortran_type;
  int integer;
  double real;
} fortran_number;


/*
 * Case-insensitive ASCII 7bit string comparison
 */
int sp_istrcmp ( const char * str1, const char * str2 );

/*
 * Returns the copy of n characters of the string
 * Caller shall free the returned string
 */
char* sp_strndup(const char *s, size_t n);

/*
 * Extract file extension from the filename
 */
const char* sp_parse_file_extension(const char* filename);

/*
 * Skip whitespaces
 */
const char* sp_skip_whitespaces(const char* line);

/*
 * Skip alphanumeric characters and characters from chars array
 */
const char* sp_skip_alnum(const char* line, const char* chars);

/*
 * Extract word from the string
 * return pointer to the word into the word argument
 * Caller shall free the word
 */
const char* sp_extract_next_word(const char* line, const char** word);

/*
 * Read the text file and return buffer with the file's contens
 * returns 0 if unable to read
 */
char* sp_read_text_file(const char* filename);


/* Extracts the integer of size bytes from the buffer from */
int sp_extract_positional_int(const char* from, int size);

/* Extracts the float of size bytes from the buffer from */
double sp_extract_positional_float(const char* from, int size);


/*
 * Very simple Parser for FORTRAN IO Format specifiers
 */
int sp_parse_fortran_format(const char* string, fortran_io_format* format);

/*
 * Extract numbers specified by fortran format
 * numbers is a pointer to an array of fortran number(used an output)
 * number of elements in the array is format->repeat value
 * extracted = number of extracted values. <= format->repeat
 * Returns the same string pointer in case of error
 */
const char* sp_extract_fortran_numbers(const char* string, 
                                       const fortran_io_format* format,
                                       /* output */
                                       fortran_number* numbers,
                                       int* extracted);

#endif /* _SP_UTILS_H_ */
