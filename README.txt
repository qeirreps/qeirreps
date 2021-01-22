# Copyright (c) 2020 Akishi Matsugatani, Seishiro Ono, Yusuke Nomura, Haruki Watanabe

======================================

Here, qeirreps is the root directory.

1. See the section "Preparation" in this document to install Quantum Espresso and qe2respack.
2. Go to "src" to compile qeirreps. See src/README_src.txt for detailed information.
3. See the section "Usage" in this document to know how to use qeirreps.
4. See "example," "example/README_ex.txt," and each example directory. They work as the tutorials.
5. See "reference" and "reference/README_ref.txt" to confirm the results.

======================================

# Contents

src:
The source directory. The Fortran 90 program "qeirreps.f90" is here.
See src/README_src.txt for compilation.

example:
The directory to store the example files.
There are some input files for calculations by QE, qe2respack, and qeirreps.
See example/README_ex.txt for detailed information.

reference:
The directory to store the reference files. 
One can check the results in example directories.
See reference/README_ref.txt for detailed information.

pseudo:
The directory to store the pseudopotential files.
See pseudo/README_psd.txt for detailed information.


======================================

# Preparation

In order to use qeirreps, you need Quantum Espresso (QE) and qe2respack.
Here, we introduce how to get and install these applications.

1. Install Quantum Espresso (QE) for DFT calculation of the target material.
   See "https://www.quantum-espresso.org/."
   We confirm that our program works with the latest version of QE, Quantum ESPRESSO v. 6.7.
   The example files in "reference" are calculated with QE v. 6.6.

2. Download qe2respack.py to prepare input files of qeirreps from the output of QE.
   See the branch of respack "maxime2" in GitHub repository "https://github.com/mnmpdadish/respackDev/tree/maxime2" to get qe2respack.py.
   qe2respack.py is in "util/qe2respack/."
   
   REMARK: This program qe2respack belongs to the program package RESPACK (https://sites.google.com/view/kazuma7k6r).
           For the latest version of RESPACK, qe2respack does not support DFT calculation with spin-orbit coupling.
           For the usage of qeirreps, please use the specific version obtained as above.


======================================


# Usage

Here, we introduce how to use qeirreps.
The example directory is also helpful as tutorials. See "example/README_ex.txt" and each example directory for more information.


1. Preparing input files of qeirreps

   qeirreps works based on the output of QE. To prepare the input files of qeirreps, the following three steps are needed to be done one by one:

1-1. Self-consistent first-principles (scf) calculation of a target material.
1-2. Non-self-consistent first-principles (nscf) calculation of the material for each high-symmetry momentum.
1-3. Data conversion from QE output files to qeirreps input files.

1-1 and 1-2: carried out by the original functions of QE.
   The wavefunction data will be produced in the directory OUTDIR/PREFIX.save, which is specified in the input files of QE.
   
   REMARK:
   Norm-conserving calculations are necessary.
   Set the option "wf_collect = .TRUE."
   Use the pseudopotentials optimized for norm-conserving calculations. See "pseudo/README_psd.txt" for detailed information.

1-3: carried out by qe2respack.
   The output files in OUTDIR/PREFIX.save should be converted by qe2respack into the form of input files of qeirreps.
   Go to the work directory (referred to as "DIRECTRY_NAME" in the following).
   Type as follows at the directory "DIRECTORY_NAME."

   $ python PATH_OF_qe2respack/qe2respack.py OUTDIR/PREFIX.save

   "PATH_OF_qe2respack" is the directory which has the file qe2respack.py.
   "OUTDIR/PREFIX.save" is the directory produced by QE in step 1-2.

   Nine files including "dat.wfn" will be generated in the directory "dir-wfn." qeirreps reads these files in step 2.


2. Running qeirreps 

   Go to the work directory (referred to as "DIRECTRY_NAME" in 1-3).
   Create a directory named "output" if it does not exist. 
   Run qeirreps by typing as 

   $ PATH_OF_qeirreps/qeirreps.x DIRECTORY_NAME

   "PATH_OF_qeirreps" is the location of the qeirreps executable file.
   "DIRECTORY_NAME" is the directory referred to in 1-3, which contains "dir-wfn" and "output."

   Some text files will be exported in "output," for example, "character.dat."
   Check the document and files in reference directory for more information.

   For materials with inversion symmetry, qeirreps can also evaluate the Z4 index.
   An option of filling should be added to the command as 

   $ PATH_OF_qeirreps/qeirreps.x DIRECTORY_NAME FILLING z4

   "FILLING" is the number of electrons per unit cell of the target material.
   This filling is shown in the standard output of scf calculation by QE as "number of electrons = FILLING."
   "z4" is the option to calculate Z4 index.

   Here, the option "z4" can be omitted and the command

   $ PATH_OF_qeirreps/qeirreps.x DIRECTORY_NAME FILLING

   also produces the same result.

   Z4 index will be exported as "z4.dat" in the directory "output."

   REMARK:
   There are some different definitions of Z4 index.
   Our program qeirreps calculates the Z4 index as shown in equation (4) in the following paper:
   E. Khalaf, H. C. Po, A. Vishwanath, and H. Watanabe, Phys. Rev. X 8, 031070 (2018) 
   ( https://journals.aps.org/prx/abstract/10.1103/PhysRevX.8.031070 ) 
 
   qeirreps can produce the input file for CheckTopologicalMat ( https://www.cryst.ehu.es/cgi-bin/cryst/programs/topological.pl )
   One can obtain the additional file "Bilbao.txt." by using the following command

   $ PATH_OF_qeirreps/qeirreps.x DIRECTORY_NAME FILLING ctm 
   where "ctm" is the option to export an input file for CheckTopologicalMat.   
   
   This file can be used as an input file of CheckTopologicalMat. 
   This tool tells us topological properties of the computed material such as the value of symmetry indicators and gaplessness.

   To use this tool, one should run nscf calculation on the high-symmetry k-points of each space group.
   To get the list of high-symmetry k-points for each space group, one should access the CheckTopologialMat and download the program named vasp2trace.
   As described in README.pdf of vasp2trace, the k-points are listed in the file max_KPOINTS_VASP.
   For the nscf calculation, one should specify the k-points as listed in max_KPOINTS_VASP.    

======================================

