# Copyright (c) 2020 Akishi Matsugatani, Seishiro Ono, Yusuke Nomura, Haruki Watanabe

======================================

Here, qeirreps is the root directory.

1. See "Preparation" section in this document to install Quantum Espresso and qe2respack.
2. Go to "src" to compile qeirreps. See src/README_src.txt for detailed information.
3. See "example", "example/README_ex.txt", and each example directory to know how to use qeirreps.
4. See "reference" and "reference/README_ref.txt" to confirm the results.

======================================

# Preparation

In order to use qeirreps, you need Quantum Espresso (QE) and qe2respack.
Here, we introduce how to get and install these applications.

1. Install Quantum Espresso (QE) for DFT calculation of the target material.
   See "https://www.quantum-espresso.org/".
   
   When you compile QE, use the following option.
   
   ./configure --enable-xml=no
   
   This option is necessary to construct the input file of "qeirreps" through "qe2respack".
   It enables qe2respack to use iotk toolkits.

2. Install qe2respack to prepare input files of qeirreps from the output of QE.
   See "https://github.com/mnmpdadish/respackDev/tree/master/util/qe2respack" to get qe2respack.
   qe2respack is available in "respackDev-master/util/qe2respack/".

   See "README.txt" in that directory, edit "Makefile" to specify the location of QE,
   and type "make" to compile qe2respack.

   Our program qeirreps currently uses qe2respack to read the output files of QE.


======================================


# Contents

src:
The source directory. The Fortran 90 program "qeirreps.f90" is available here.
See src/README_src.txt for compilation.

example:
The directory to store the example files.
See example/README_ex.txt for detailed information.

reference:
The directory to store the reference files.
See reference/README_ref.txt for detailed information.

pseudo:
The directory to store the files of pseudo potential.
See pseudo/README_psd.txt for detailed information.


======================================
