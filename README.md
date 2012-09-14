libspmatrix Library Description
===============================

Introduction
------------
libspmatrix is a simple library for working with sparse matrices. It was designed to provide the generic sparse matrix techniques for the finite element analysis.
Therefore all methods to solve linear algebraic equation systems implemented in the library requires matrices to be positive linear-definite. It supports different storage schemes as well as different file formats. 

Features
--------
 * Support for different sparse matrix formats: CRS/CCS(3-arrays) or Yale format, CS(lower triangle)R(Skyline) format, internal dynamic-arrays based format
 * Solvers: Conjugate Gradedient method, Preconditioned Conjugate Gradient(with ILU preconditioner), Sparse Cholesky Solver (direct solver)
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
 * Fill it with values (following the typical for the FEM schema: construct local stiffness matrix and distribute its values in global stiffness matrix by *adding* values to existing ones)
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




