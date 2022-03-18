README for MINPACK
==================

Minpack includes software for solving nonlinear equations and
nonlinear least squares problems.  Five algorithmic paths each include
a core subroutine and an easy-to-use driver.  The algorithms proceed
either from an analytic specification of the Jacobian matrix or
directly from the problem functions.  The paths include facilities for
systems of equations with a banded Jacobian matrix, for least squares
problems with a large amount of data, and for checking the consistency
of the Jacobian matrix with the functions.

Jorge Moré, Burt Garbow, and Ken Hillstrom at Argonne National Laboratory.


References
==========

* M. J. D. Powell, *A Hybrid Method for Nonlinear Equations*.
  Numerical Methods for Nonlinear Algebraic Equations, P. Rabinowitz, editor. Gordon and Breach, 1970.

* Jorge J. Moré, *The Levenberg-Marquardt Algorithm, Implementation and Theory*. 
  Numerical Analysis, G. A. Watson, editor. Lecture Notes in Mathematics 630, Springer-Verlag, 1977.
  https://doi.org/10.1007/BFb0067700

* J. J. Moré, B. S. Garbow, and K. E. Hillstrom, *User Guide for MINPACK-1*, 
  Argonne National Laboratory Report ANL-80-74, Argonne, Ill., 1980.

* J. J. Moré, D. C. Sorensen, K. E. Hillstrom, and B. S. Garbow, *The MINPACK Project*, 
  in Sources and Development of Mathematical Software, W. J. Cowell, ed., Prentice-Hall, pages 88-111, 1984.

* J.J. Moré, B.S. Garbow and K.E. Hillstrom, *Testing unconstrained optimization software*. 
  ACM Trans. Math. Soft. 7, 1 (March 1981), 17-41. https://doi.org/10.1145/355934.355943


MINPACK versions
================

* The original MINPACK-1 Fortran version from 1980. https://www.netlib.org/minpack/
  
  Extensive original documentation. Extensive battery of tests (described in
  Algorithm 566 from ACM TOMS).


* MINPACK-2. https://ftp.mcs.anl.gov/pub/MINPACK-2/

  A MINPACK-2 project is announced. Contains some code and some new test problems,
  but is not backward compatible with the original MINPACK-1.


* CMINPACK from Manolis Lourakis.  http://www.netlib.org/minpack/cminpack.tar
  
  In 2002, Lourakis (lourakis at ics forth gr) published a C version derived
  from the fortran code using `f2c` and limited manual editing. This does not
  offer much more than directly linking against the original Fortran library.


* C/C++ MINPACK from Frédéric Devernay. http://devernay.free.fr/hacks/cminpack/
  
  Starting from 2007 and until 2021(at least), Devernay published much
  improved C versions (readable C code, friendlier C interface, no GOTOs,
  some tests, man pages, debian package, some CUDA support, CMAKE support,
  some LAPACK calls). A notable downside is that this does not run the
  original full test suite.
	
* Modernized Fortran code. https://github.com/fortran-lang/minpack

  Several people (including John Burkardt) have produced modernized Fortran
  versions of the original code (Fortran 90, free form, etc.).  Notably,
  there is an ongoing effort (as of 2022) to maintain a modern Fortran
  version of the original code.

Testing Minpack
===============

Testing MINPACK is quite difficult.  The original test code processed all test
problems and printed a summary that had to be manually inspected.  This is
unsuitable for continuous integration.  However, this has the merit of
existing already.  A basic approach consists in checking that the test
programs generate the same output as the orginal MINPACK Fortran code. 

It must be noted that this include cases where MINPACK actually fails to
converge to a solution --- these are "known failures".


Basic approach
--------------
The main problem is that checking that the test programs produce the exact
same output as the original fortran code is only a partial solution.  It
allows to check whether some code produces the exact same result as the
original. This was used to check that manually editing the C code did not
introduce bugs.

This detects bizare situations.  For instance, automatically translating the
fortran code to C (using f2c) and compiling it yields the exact same results
on x86 CPUs, but different results on Power8 CPUs (this happens both with gcc
and clang).

Also, compiling the C code using the "-O1" and "-O3" flags yields the same
results on x86 CPUs, but different results on Power8 CPUs (this happens only
with gcc, not clang).

Are these bugs in f2c? Or bugs in gcc and clang on Power8? Or is it something
else (badly defined language semantics)?

However, very minor modifications to the code yield different numetrical
values. This is unsuitable to test them.  Also, using BLAS and LAPACK make
the whole thing even more complicated, because then numerical values depend
on the actual BLAS used, and on specificities of the hardware.


Refined approach
----------------

The test problems have (mostly) known solutions (as documented in ACM TOMS
algorithm 566). A much better solution consists in checking that the test
programs generate the "expected" outputs up to a certain precision --- or
fail to converge when the original MINPACK also failed.

This has the advantage of being tolerant to small deviations in the numerical
results.

However, this is still not bullet-proof, as discussed below.


In this repository
==================

This repository contains several versions of MINPACK-1 converted to in C and
potentially altered by Charles Bouillaguet(charles.bouillaguet@lip6.fr).  All
errors are his own!

All versions should be drop-in replacement for the original Fortran code.

Tests
-----

The original 32 test cases have been ported to C, and modified to detect
success or failure. This uses the "Test Anything Protocol" that originates
from perl (using the `prove` test harness which is part of perl).  The new
set of tests is thus suitable for Continuous Integration.

When the original Fortran code fails to converge to the expected solution, the
tests are maked as "known failure".

Some limited benchmarking code has been added.


`fortran` version
-----------------

This is the original fortran code, along with the C tests.


`f2ced` version
-----------------

This is the result of converting the original Fortran code to C using f2c, 
along with the C tests.


`base` version
--------------

The code in the `base/` folder has been obtained by converting the original
fortran code to C (using `f2c`), restructuring the code to avoid most GOTOs,
and some manual editing.

The `lmstr` (non-linear least squares with limited storage) function has not
been included.

The main difference with the Fortran code is that the `dpmpar` function has
been removed, and replaced by constants from the `float.h` standard header.

There are tiny differences between the constants hardcoded in `dpmpar` and the
actual values (`2.22044604926e-16` vs `2.2204460492503130808e-16` for the
actual machine epsilon of an IEE 754 `double`).  *These changes alter the
behavior of the code*.

Using a different compiler alters the behavior of the code (e.g. tests fail
with `gcc`, succeed with `clang`). Changing compiler options may make some
test fails (e.g. pass with `-O1` but fail with `-O2`; actually observed with
`gcc` on POWER8 processors).


`lapackified` version
---------------------

The code in the `lapackified/` folder differs from the `base/` code as follows:

- It uses LAPACK to compute QR factorizations and related operations. This
  reduces the amount of code (completely removes the `qrfac` function) and
  makes it much faster on large instances.

- It uses some BLAS functions when applicable, notably to reduce code size.

- The `mldif` and `lmder` functions are very similar; they have been merged.

- The `hybrid` and `hybrj` functions are very similar; they have been merged.

- All functions have been converted to zero-based indexing.

All-in-all, the `lapackified` code is about 25% smaller and much faster than
the `base` code.

However, the actual numerical results will depend on the actual BLAS
implementation.  Most BLAS libraries adjust to the underlying hardware.
Therefore, numerical results may vary from machine to machine.

The reference BLAS yields different numerical values than OpenBLAS/ATLAS. OpenBLAS
is multi-threaded; changing the number of threads alters the numerical
values. Running in `valgrind` using OpenBLAS or the ATLAS BLAS changes the
numerical values. At least one test fails under valgrind (does not converge). 
This does not happen with the reference BLAS. 
