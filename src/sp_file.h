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
#ifndef _SP_FILE_H_
#define _SP_FILE_H_

#include "sp_matrix.h"

/*
 * Load the sparse martix from the file.
 * type = CRS or CCS - preferred storage type
 * Currently supported formats:
 * MM (Matrix Market) (*.mtx)
 * Harwell-Boeing format (*.hb, *.r[su]a)
 * if file is in Harwell-Boeing format, storage type is ignored
 */
int sp_matrix_yale_load_file(sp_matrix_yale_ptr self,
                             const char* filename,
                             sparse_storage_type type);

/*
 * Save the sparse martix from the file.
 * File format guessed from the extension
 * Currently supported formats:
 * MM (matrix market), file extension .mtx
 * see http://math.nist.gov/MatrixMarket/formats.html#MMformat
 * txt, 0-based triplet format, each line is 0-based tripet:
 * row, col, value
 * Returns 0 if not possible to write(or unknown file format)
 * Side-effect: matrix gets ordered
 */
int sp_matrix_save_file(sp_matrix_ptr self, const char* filename);
int sp_matrix_yale_save_file(sp_matrix_yale_ptr self, const char* filename);



#endif /* _SP_FILE_H_ */
