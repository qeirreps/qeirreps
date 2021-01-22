# Copyright (c) 2020 Akishi Matsugatani, Seishiro Ono, Yusuke Nomura, Haruki Watanabe

======================================

Here, qeirreps/reference/NaCdAs_nonsoc is the directory which stores the output files for DFT calculation of NaCdAs without spin-orbit coupling.

See also qeirreps/example/NaCdAs_nonsoc to test the calculation processes.

======================================

# Contents

NaCdAs.scf.in, NaCdAs.band.in, NaCdAs.bands.in, NaCdAs.rep.in:
The input files of QE
Check detailed information in qeirreps/example/NaCdAs_nonsoc

NaCdAs.scf.out:
The log file of QE for self-consistent first-principles (scf) calculation

NaCdAs.band.out:
The log file of QE for non-self-consistent first-principles (nscf) calculation to obtain the band structure

NaCdAs.bands.out:
The log file of QE to obtain the plot of the band structure

NaCdAs.rep.out:
The log file of QE for non-self-consistent first-principles (nscf) calculation to obtain the irreducible representation

NaCdAs.band, NaCdAs.band.rap, NaCdAs.band.gnu:
The result files of nscf calculation with input files "NaCdAs.band.in" and "NaCdAs.bands.in"
One can plot the band structure from "NaCdAs.band.gnu" by using gnuplot.

dir-wfn:
The directory to store the input files for "qeirreps."
"qe2respack" writes the result files here.
REMARK: In the case of NaCdAs, dat.wfn is not stored because the file size is too large.

output:
The directory to store the result files of "qeirreps."
"qeirreps" writes the result files here.
See the section "Output of qeirreps" in this document for detailed information.

example_NaCdAs.nb:
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
REMARK:
One can find a warning in the standard output as

	Warning: number of irreps is not an integer
	 at k=  0.500000000000000       0.500000000000000       0.000000000000000E+000 
	 and level=          27 :
	 (0.458579004940740,-4.264861192059216E-018)

This is the warning in the processes to count the irreducible representations in each energy level.
This is the warning for the highest energy level in the unoccupied part. This indicates that the highest energy level has degeneracy, but only a part of degenerate bands are computed in the nscf calculation.
There is no problem with the other energy levels, such as the occupied bands.


======================================

# Output of qeirreps

The output files of qeirreps are in "output."

There are two types of files, "*.dat" and "*_import.txt."
The files named "*.dat" are ordinary output files.
The files named "*_import.txt" are constructed for Mathematica usage.


1. "*.dat"


1-1. pg.dat:
     The point group part of symmetry operators in the cartesian coordinate.
     
     The indices of symmetry operators and the 3 by 3 matrices are described as     

                1
        1.00000000000000       0.000000000000000E+000  0.000000000000000E+000
       0.000000000000000E+000   1.00000000000000       0.000000000000000E+000
       0.000000000000000E+000  0.000000000000000E+000   1.00000000000000     


1-2. tg.dat:
     The translation part of symmetry operators in the crystalline coordinate.

     The indices of symmetry operators and the 3 dimensional vectors are described as

                1
       0/  2   0/  2   0/  2

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

     k=( 0.50000, 0.00000, 0.00000),

     For each k-point, the symmetry operations of G_k are described by the indices like

      symmetry operator indices in G_k:
        1    2    3    4    5    6    7    8 

     These indices correspond to the numbers produced in standard output, pg.dat, tg.dat, and srg.dat. Check them to confirm the symmetry operations.

1-5. factor_system_spin.dat:
     The factor system of G_k on each k-point.

     k-points are described in the crystalline coordinate like

     k=( 0.50000, 0.00000, 0.00000),

     The factor system on each k-point is described as a matrix.
     Each row and column represents the index of symmetry operators produced in standard output, pg.dat, tg.dat, and srg.dat. Check them to confirm the symmetry operations.
     If the symmetry operator does not belong to G_k, the value of factor system is described as zero.

1-6. irreps_list.dat:
     The character tables of irreducible representations of G_k.

     k-points are described in the crystalline coordinate like

     k=( 0.50000, 0.00000, 0.00000),

     For each k-point, the symmetry operations of G_k are described by the indices like

      symmetry operator indices in G_k:
        1    2    3    4    5    6    7    8 

     These indices correspond to the numbers produced in standard output, pg.dat, tg.dat, and srg.dat. Check them to confirm the symmetry operations.

1-7. irreps_number.dat:
     The number of irreducible representations of G_k in the n-th energy level.

     k-points are described in the crystalline coordinate like

     k=( 0.50000, 0.00000, 0.00000),

     For each energy level, the integers show the number of irreducible representations 1, 2, 3... from left to right.
     The indices of irreducible representations are described in irreps_list.dat. Check it to confirm the character tables of the irreducible representations.


2. For Mathematica usage, "*_import.txt"

   See and edit the Mathematica notebook file "example_NaCdAs.nb" to confirm the information of each file.

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
						      
======================================
