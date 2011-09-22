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

#include <stdlib.h>
#include <string.h>
#include <ctype.h>

#include "sp_utils.h"


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
  char* str = malloc(n+1);
  memcpy(str,s,n);
  str[n] = '\0';
  return str;
}



/*
 * Extract file extension from the filename
 */
const char* sp_parse_file_extension(const char* filename)
{
  const char* ptr;
  int len = strlen(filename);
  ptr = filename + len;
  while (len)
  {
    ptr--,len--;
    if ( *ptr == '.' )
      break;
  }
  return len ? ptr + 1 : 0;
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
 * Skip alphanumeric characters
 */
const char* sp_skip_alnum(const char* line)
{
  while (*line && isalnum(*line))
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
  end = sp_skip_alnum(line);
  *word = sp_strndup(line,end - line);
  return end;
}
