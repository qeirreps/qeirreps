# Copyright (c) 2020 Akishi Matsugatani, Seishiro Ono, Yusuke Nomura, Haruki Watanabe

======================================

Here, qeirreps/reference/Si_nonsoc is the directory which stores the output files for DFT calculation of Si without spin-orbit coupling.

See also qeirreps/example/Si_nonsoc to test the calculation processes.

======================================

# Contents

Si.scf.in, Si.band.in, Si.bands.in, Si.rep.in:
The input files of QE
Check detailed information in example/Si_nonsoc

Si.scf.out:
The log file of QE for self-consistent first-principles (scf) calculation

Si.band.out:
The log file of QE for non-self-consistent first-principles (nscf) calculation to obtain the band structure

Si.bands.out:
The log file of QE to obtain the plot of the band structure

Si.rep.out:
The log file of QE for non-self-consistent first-principles (nscf) calculation to obtain the irreducible representation

Si.band, Si.band.rap, Si.band.gnu:
The result files of nscf calculation with input files "Si.band.in" and "Si.bands.in"
One can plot the band structure from "Si.band.gnu" by using gnuplot.

dir-wfn:
The directory to store the input files for "qeirreps".
"qe2respack" writes the result files here.

output:
The directory to store the result files of "qeirreps".
"qeirreps" writes the result files here.
See the section "Output of qeirreps" in this document for detailed information.

exsample_Si.nb:
The Mathematica notebook to read the results of "qeirreps".
See the section "Output of qeirreps" in this document for detailed information.

qeirreps.log:
The log text of standard output of "qeirreps.x".
This contains following information.
  - the name of the reading directory
    (if you specify the filling) the number of filling 
  - ncomp: the number of components in the Bloch wave function.
           For calculations Without spin-orbit coupling, ncomp=1.
           For calculations With spin-orbit coupling, ncomp=2.
  - the primitive/reciprocal lattice vectors
  - the symmetry operation types for each space group element
  - NBAND: the number of bands considered in the calculation 
  - (if you specify the filling) Z4 index

======================================

# Output of qeirreps

The output files of qeirreps are in "output".

There are two types of files, "*.dat" and "*_import.txt".
The files named "*.dat" are ordinary output files.
The files named "*_import.txt" are constructed for Mathematica usage.


1. "*.dat"


1-1. pg.dat:
     The point group part of symmetry operators in the cartesian coordinate.
     
     The indices of symmetry operators and the 3 by 3 matrices are described as     

                1
        1.00000000000000       5.551115123125783E-017 -5.551115123125783E-017
       5.551115123125783E-017   1.00000000000000      -5.551115123125783E-017
       1.110223024625157E-016  2.775557561562891E-016   1.00000000000000     


1-2. tg.dat:
     The translation part of symmetry operators in the crystalline coordinate.

     The indices of symmetry operators and the 3 dimensional vectors are described as

                2
       0/  2   0/  2  -1/  2
 
1-3. srg.dat:
     The spin rotation part of symmetry operators in the cartesian coordinate.

     The indices of symmetry operators and the 2 by 2 matrices are described as

                1
       1.000000  0.000000  0.000000  0.000000
       0.000000  0.000000  1.000000  0.000000

     where 1st(2nd) and 3rd(4th) column represent real(imaginary) part of the matrix element, respectively.

1-4. character.dat:
     The characters for G_k of each band.

     k-points are described in the the crystalline coordinate like

     k=( 0.50000, 0.00000, 0.00000),

     For each k-point, the symmetry operations of G_k are described by the indices like

      symmetry operator indices in G_k:
        1    5    9   14   18   23   25   29   33   38   42   47 

     These indices correspond to the numbers produced in standard output, pg.dat, tg.dat, and srg.dat. Check them to confirm the symmetry operations.

1-5. factor_system_spin.dat:
     The factor system of G_k on each k-point.

     k-points are described in the crystalline coordinate like

     k=( 0.50000, 0.00000, 0.00000),   

     The factor system on each k-point is described as a matrix.
     Each row and column represents the index of symmetry operators produced in standard output, pg.dat, tg.dat, and srg.dat. Check them to confirm the symmetry operations.
     If the symmetry operator does not belong to G_k, the value of factor system is described as zero.


2. For Mathematica usage, "*_import.txt"

   See and edit the Mathematica notebook file "example_Si.nb" to confirm the information of each file.

2-1. pg_import.txt: imported to "pgList".
     The rotation part of symmetry operators.

     pgList[[i]]: The 3 by 3 matrix of the rotation part of the i-th symmetry operator in cartesian coordinate.

2-2. tg_import.txt: imported to "tgList".
     The translation part of symmetry operators.

     tgList[[i]]: The 3 dimensional vector of the translation part of the i-th symmetry operator in crystalline coordinate.

2-3. srg_import.txt: imported to "srgList".
     The spin rotation part of symmetry operators in cartesian coordinate.
     You can ignore this when you run calculation without spin-orbit coupling.

     srgList[[i]]: The 2 by 2 matrix of the spin rotation part of the i-th symmetry operator.

2-4. character_import.txt: imported to "BandList".
     The characters for G_k of each band.

     BandList[[i,1]]: The 3 dimensional k-vector of the i-th k-point in the crystalline coordinate.
     BandList[[i,2,n,j]]: The character of the n-th band for the j-th symmetry operator. 
   		          If the j-th symmetry is not an element of G_k on the i-th k-point, this value becomes "Null".
   
2-5. factor_system_spin_import.txt: imported to "factorSystemSpin".
     The factor system of G_k on each k-point.

     factorSystemSpin[[i,1]]: The 3 dimensional k-vector of the i-th k-point in the crystalline coordinate.
     factorSystemSpin[[i,2,m,n]]: The factor system between m-th and n-th symmetry operators.
	                          If the m-th or n-th symmetry is not an element of G_k on the i-th k-point, this value becomes "Null".
						      
======================================
