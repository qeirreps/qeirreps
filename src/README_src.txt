# Copyright (c) 2020 Akishi Matsugatani, Seishiro Ono, Yusuke Nomura, Haruki Watanabe

======================================

Here, qeirreps/src is the source directory.

qeirreps.f90:
The Fortran 90 source code of qeirreps.

Makefile:
The setting file for compile qeirreps.

======================================

How to compile qeirreps

0. Requirement
   - Fortran compiler
   - BLAS library
   - LAPACK library (intel MKL)

1. Edit "Makefile" for your computational situation.
   Put the link for Fortran compiler and MKL library.

2. Compile qeirreps. Type "make".
   The executable binary "qeirreps.x" will be produced.

======================================
