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
#ifndef _SP_FILE_H_
#define _SP_FILE_H_

#include "sp_matrix.h"


/*
 * Load the sparse martix from the file.
 * Currently supported formats:
 * MM (matrix market), see http://math.nist.gov/MatrixMarket/formats.html#MMformat
 */
sp_matrix_ptr sp_matrix_load_file(const char* filename);



#endif /* _SP_FILE_H_ */
