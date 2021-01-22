# Copyright (c) 2020 Akishi Matsugatani, Seishiro Ono, Yusuke Nomura, Haruki Watanabe

======================================

Here, qeirreps/reference/Bi_soc is the directory which stores the output files for DFT calculation of bismuth with spin-orbit coupling.

See also qeirreps/example/Bi_soc to test the calculation processes.

======================================

# Contents

Bi.scf.in, Bi.band.in, Bi.bands.in, Bi.rep.in:
The input files of QE
Check detailed information in example/Bi_soc

Bi.scf.out:
The log file of QE for self-consistent first-principles (scf) calculation

Bi.band.out:
The log file of QE for non-self-consistent first-principles (nscf) calculation to obtain the band structure

Bi.bands.out:
The log file of QE to obtain the plot of the band structure

Bi.rep.out:
The log file of QE for non-self-consistent first-principles (nscf) calculation to obtain the irreducible representation

Bi.band, Bi.band.rap, Bi.band.gnu:
The result files of nscf calculation with input files "Bi.band.in" and "Bi.bands.in"
One can plot the band structure from "Bi.band.gnu" by using gnuplot.

dir-wfn:
The directory to store the input files for "qeirreps."
"qe2respack" writes the result files here.

output:
The directory to store the result files of "qeirreps."
"qeirreps" writes the result files here.
See the section "Output of qeirreps" in this document for detailed information.

example_Bi.nb:
The Mathematica notebook to read the results of "qeirreps."
See the section "Output of qeirreps" in this document for detailed information.

qeirreps.log:
The log text of standard output of "qeirreps.x."
This contains the following information.
  - the name of the reading directory
    (if you specify the filling) the number of filling 
  - ncomp: the number of components in the Bloch wave function.
           For calculations without spin-orbit coupling, ncomp=1.
           For calculations with spin-orbit coupling, ncomp=2.
  - the primitive/reciprocal lattice vectors
  - the symmetry operation types for each space group element
  - NBAND: the number of bands considered in the calculation 
  - (if you specify the filling) Z4 index

for_CTM:
The directory to store the output files with the "ctm" option.
This contains Bi.rep.in, Bi.rep.out, dir-wfn, output, and qeirreps.log.
Bi.rep.in has a different choice of k-points to export an input file of CheckTopologicalMat. Bi.rep.out is the log file of nscf calculations.
dir-wfn also has a different choice of k-points.
output has the additional file named Bilbao.txt, which is the input file for CheckTopologicalMat.
qeirreps.log is the log file of qeirreps with the "ctm" option.

======================================

# Output of qeirreps

The output files of qeirreps are in "output."

There are three types of files, "*.dat," "*_import.txt," and "Bilbao.txt."
The files named "*.dat" are ordinary output files.
The files named "*_import.txt" are constructed for Mathematica usage.
The files named "Bilbao.txt" is the input file for CheckTopologicalMat.


1. "*.dat"


1-1. pg.dat:
     The point group part of symmetry operators in the cartesian coordinate.
     
     The indices of symmetry operators and the 3 by 3 matrices are described as     

                1
        1.00000000000000       0.000000000000000E+000  0.000000000000000E+000
      -1.568992671209136E-017   1.00000000000000      -5.551115123125783E-017
       2.353374252904200E-017  0.000000000000000E+000   1.00000000000000     


1-2. tg.dat:
     The translation part of symmetry operators in the crystalline coordinate.

     The indices of symmetry operators and the 3 dimensional vectors are described as

                1
       0/  1   0/  1   0/  1
 
1-3. srg.dat:
     The spin rotation part of symmetry operators in the cartesian coordinate.

     The indices of symmetry operators and the 2 by 2 matrices are described as

                1
       1.000000  0.000000  0.000000  0.000000
       0.000000  0.000000  1.000000  0.000000

     where 1st(2nd) and 3rd(4th) columns represent the real(imaginary) part of the matrix element.

1-4. character.dat:
     The characters of G_k of each band.

     k-points are described in the crystalline coordinate like

     k=( 0.50000, 0.00000,-0.00000),

     For each k-point, the symmetry operations of G_k are described by the indices like

      symmetry operator indices in G_k:
         1    5    7   11 

     These indices correspond to the numbers produced in standard output, pg.dat, tg.dat, and srg.dat. Check them to confirm the symmetry operations.

1-5. factor_system_spin.dat:
     The factor system of G_k on each k-point.

     k-points are described in the crystalline coordinate like

     k=( 0.50000, 0.00000,-0.00000),   

     The factor system on each k-point is described as a matrix.
     Each row and column represents the index of symmetry operators produced in standard output, pg.dat, tg.dat, and srg.dat. Check them to confirm the symmetry operations.
     If the symmetry operator does not belong to G_k, the value of factor system is described as zero.

1-6. irreps_list.dat:
     The character tables of irreducible representations of G_k.

     k-points are described in the crystalline coordinate like

     k=( 0.50000, 0.00000,-0.00000),

     For each k-point, the symmetry operations of G_k are described by the indices like

      symmetry operator indices in G_k:
         1    5    7   11 

     These indices correspond to the numbers produced in standard output, pg.dat, tg.dat, and srg.dat. Check them to confirm the symmetry operations.

1-7. irreps_number.dat:
     The number of irreducible representations of G_k in the n-th energy level.

     k-points are described in the crystalline coordinate like

     k=( 0.50000, 0.00000,-0.00000),

     For each energy level, the integers show the number of irreducible representations 1, 2, 3... from left to right.
     The indices of irreducible representations are described in irreps_list.dat. Check it to confirm the character tables of the irreducible representations.

1-8. z4.dat:
     This files is stored in ./output/, not in ./for_CTM/output/.
     The value of Z4 index defined from the sum of inversion parities of occupied bands at all TRIMs.
     

2. For Mathematica usage, "*_import.txt"

   See and edit the Mathematica notebook file "example_Bi.nb" to confirm the information of each file.

2-1. pg_import.txt: imported to "pgList."
     The rotation part of symmetry operators.

     pgList[[i]]: The 3 by 3 matrix of the rotation part of the i-th symmetry operator in cartesian coordinate.

2-2. tg_import.txt: imported to "tgList."
     The translation part of symmetry operators.

     tgList[[i]]: The 3 dimensional vector of the translation part of the i-th symmetry operator in crystalline coordinate.

2-3. srg_import.txt: imported to "srgList."
     The spin rotation part of symmetry operators in cartesian coordinate.
     You can ignore this when you run calculation without spin-orbit coupling.

     srgList[[i]]: The 2 by 2 matrix of the spin rotation part of the i-th symmetry operator.

2-4. character_import.txt: imported to "BandList."
     The characters of G_k of each band.

     BandList[[i,1]]: The 3 dimensional k-vector of the i-th k-point in the crystalline coordinate.
     BandList[[i,2,n,j]]: The character of the n-th band for the j-th symmetry operator. 
   		          If the j-th symmetry is not an element of G_k on the i-th k-point, this value becomes "Null."
   
2-5. factor_system_spin_import.txt: imported to "factorSystemSpin."
     The factor system of G_k on each k-point.

     factorSystemSpin[[i,1]]: The 3 dimensional k-vector of the i-th k-point in the crystalline coordinate.
     factorSystemSpin[[i,2,m,n]]: The factor system between m-th and n-th symmetry operators.
	                          If the m-th or n-th symmetry is not an element of G_k on the i-th k-point, this value becomes "Null."


3. Bilbao.txt:
   This file is stored in ./for_CTM/output/, not in ./output/.
   The input file for CheckTopologicalMat ( https://www.cryst.ehu.es/cgi-bin/cryst/programs/topological.pl )
   CheckTopologicalMat tells us topological properties of computed material.

						      
======================================
