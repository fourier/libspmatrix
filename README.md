libspmatrix Library Description
===============================

Introduction
------------
libspmatrix is a simple library for working with sparse matrices. It was designed to provide the generic sparse matrix techniques for the finite element analysis.
Therefore all methods to solve linear algebraic equation systems implemented in the library requires matrices to be positive linear-definite. It supports different storage schemes as well as different file formats. 

Features
--------
 * Support for different sparse matrix formats: CRS/CCS(3-arrays) or Yale format, CS(lower triangle)R(Skyline) format, internal dynamic-arrays based format
 * Solvers: 
   * Conjugate Gradient method
   * Preconditioned Conjugate Gradient(with ILU preconditioner)
   * Conjugate Gradient **Squared** method (for non-symmetric matrices)
   * Transpose-Free QMR method (for non-symmetric matrices)
   * Sparse Cholesky Solver (direct solver)
 * Sparse matrix file input formats: Matrix Market, Harwell-Boeing
 * Sparse matrix file output formats: Matrix Market, txt 0-based triplet format (each line is 0-based triplet: row, column, value)
 * Sparse Cholesky Solver based on book (T.Davis Direct Solvers for Sparse Lineer systems), therefore supported operations like elimination tree construction, symbolic Cholesky decomposition, numeric Cholesky decomposition 
 * Test suite based on own tiny unit test framework
 * Command line tool to compare performance of different SLAE solvers for given file with a sparse matrix
 * Small 2D Finite element application used to generate typical matrices arising in FEA applications  
 
Typical usage
-------------
The typical usage is the following(dictated by the typical FEM-procedure):
 * Create the matrix of specified size in format allowing fast element insertion
 * Fill it with values (following the typical for the FEM schema: construct local stiffness matrix and distribute its values in global stiffness matrix by **adding** values to existing ones)
 * Compress (reorder) columns/rows in order to prepare to conversion to the format more appropriate for the computation purposes
 * Convert matrix to the appropriate format for faster SLAE solving
 * (If necessary create appropriate symbolic/numeric decompositions)
 * (Optionally store matrix to the file, in order to analyze with specific software)
 * Solve SLAE

Directory structure
-------------------
<pre>
.
|-- bin - Contains binaries for utilities after the build
|-- demo_src - Contains demo application implementing simple FEA task for 2D case
|-- examples - Contains example matrices used for test purposes. All matrices are taken from MatrixMarket
|-- inc - library include files 
|-- lib - Output directory for the library binary
|-- obj - generated objects
|-- solver_src - Contains utility used to compare effectiveness of different SLAE solvers
|-- src - library sources
`-- test_src - Unit tests for the library
</pre>

Getting started
---------------
To use the library, build it (on Linux and MacOS X simply run **make**). The library itself will be generated in *lib/* directory. Link your sources against the library and add *inc/* as an include path to your project. 

Matrix formats
==============
The format used to operate with sparse matrices is Yale format(http://en.wikipedia.org/wiki/Sparse_matrix#Yale_format). The structure describing matrix in this format is **sp_matrix_yale**. 'Generic' matrix format used to create the matrix manually represented by the structure **sp_matrix**.
Both declaration could be found in *sp_matrix.h* in *inc/* directory. In order to initialize and free structures of these formats the appropriate init/free functions used.

How to create matrix
====================
There are 3 main ways to construct the matrix in this format:
 * Create it in the internal format and convert to the Yale format
 * Load it from the file
 * Operate with the Yale format directly. Use it on your risk.

Manual construction of the sparse matrix
========================================
Since Yale matrix format was designed to operate with matrices, not to create them, it is easier to create a matrix in any other format and convert it to the Yale format. The internal sparse matrix format is List of lists(LIL)-based format with dynamic arrays as an appropriate row/column storage(instead of lists).
Lets introduce the matrix creation/conversion process by example.
Suppose we want to create and convert to Yale format the matrix
<pre>
  | 0  1  0 |
  | 0 -1 -1 |
  | 1  0  2 |
</pre>
First we need to declare structures
```c
sp_matrix mtx;
sp_matrix_yale yale;
```
Then we need to initialize the first structure used to create matrix:
```c
sp_matrix_init(&mtx,3,3,2,CRS);
```
Here second and third arguments are number of rows and columns, forth - bandwidth(approximate number of nonzero elements in row), fifth - storage format (stored row-wise)
Now we can fill the matrix with values.
Since by default matrix considered to be zero matrix, we can 'add' values to appropriate matrix position - the way sparse matrix constructed in FEM. We can use the function *sp_matrix_element_add* for this purpose. However there is already exists a macro **MTX** for this operation:
 1. With the sp_matrix_element_add function call:
```c
sp_matrix_element_add(&mtx,0,1,1);
sp_matrix_element_add(&mtx,1,1,-1); sp_matrix_element_add(&mtx,1,2,-1);
sp_matrix_element_add(&mtx,2,0,1); sp_matrix_element_add(&mtx,2,2,2);
```
 2. With the **MTX** macro
```c
MTX(&mtx,0,1,1);
MTX(&mtx,1,1,-1); MTX(&mtx,1,2,-1);
MTX(&mtx,2,0,1);MTX(&mtx,2,2,2);
```

In order to modify matrix element one can get a *pointer* to this element:
```c
double *pvalue = sp_matrix_element_ptr(&mtx,0,1);
```
This function will return 0 if an element is not exist.



Loading sparse matrices from file
=================================



