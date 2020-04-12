# Copyright (c) 2020 Akishi Matsugatani, Seishiro Ono, Yusuke Nomura, Haruki Watanabe

======================================

Here, qeirreps/reference/Si_soc is the directory 
which stores the output files for DFT calculation of silicon
without spin-orbit coupling.

See also qeirreps/example/Si_soc to confirm the detail of input files.

======================================

# Files and directories

Si.scf.in, Si.band.in, Si.bands.in, Si.rep.in:
The input files of QE
Check detailed information in qeirreps/example/Si_soc

Si.scf.out:
The log file of QE for self-consistent first-principle (scf) calculation

Si.band.out:
The log file of QE for non-self-consistent first-principle (nscf) calculation
to obtain the band structure

Si.bands.out:
The log file of QE to obtain the plot of the band structure

Si.rep.out:
The log file of QE for non-self-consistent first-principle (nscf) calculation
to obtain the irreducible representation

Si.band, Si.band.rap, Si.band.gnu:
The result files of nscf calculation with input "Si.band.in" and "Si.bands.in"
You can plot the band structure by gnuplot from "Si.band.gnu".

dir-wfn:
The directory to store the input files for "qeirreps".
"qe2respack" writes the result files of it here.

output:
The directory to store the result files of "qeirreps".
"qeirreps" writes the result files of it here.

exsample_Si.nb:
The Mathematica notebook to read the results of "qeirreps".
See "Output of qeirreps" in the next part for more information.

qeirreps.log:
The log text of standard output when you run "qeirreps.x".
This contains following information.

the name of the reading directory
(if you specify the filling) the number of filling 
ncomp: the number of components in the Bloch wave function.
       For calculations Without spin-orbit coupling, ncomp=1.
       For calculations With spin-orbit coupling, ncomp=2.
the lattice/reciprocal lattice vectors
the symmetry operation types for each space group element
the detailed information of symmetry operators
NBAND: the number of bands considered in the calculation 
(if you specify the filling) the summation of parities for time-reversal-invariant-momenta

======================================

# Output of qeirreps

The output files of qeirreps are available in "output".

There are two types of files, "*.dat" and "*_import.txt".
The files named "*.dat" are constructed for brief check.
The files named "*_import.txt" are constructed for Mathematica usage.


1. For brief check, "*.dat"

1-1. character.dat:
     The characters for G_k of each band.
     (REMARK: These are not characters for G_k/T but for G_k 
	      where T is the lattice translational group.
              Then these characters contain phase factors
              as Exp[-i \vec{k}\cdot\vec{t}_g].)

     k-points are described in the crystalline coordinate like
     "k=( 0.50000, 0.00000,-0.00000),".

     For each k-point, the symmetry operations of G_k are described by the indices like
     "symmetry operator indices in G_k:
         1    2    3    4    5    6    7    8    9   10   11   12 "
     These indices correspond to the numbers produced in standard output.
     Check "qeirreps.log" to confirm the symmetry operations.

1-2. factor_system_nonsymmorphic.dat:
     The nonsymmorphic part of the factor system of G_k/T on each k-point.
     This value is defined as follows.

     \omega(g,g') = Exp[i \vec{k}\cdot(\vec{t}_{g'} - p_g\vec{t}_{g'})]
   
     k-points are described in the crystalline coordinate like
     "k=( 0.50000, 0.00000,-0.00000),".   

     The factor system on each k-point is described as a matrix.
     Each row and column represents the index of symmetry operators described in standard output.
     Check "qeirreps.log" to confirm the indices of symmetry operations.

     If the symmetry operator does not exist in G_k/T, the value of factor system is described as zero.

1-3. factor_system_spin.dat:
     The spin part of the factor system of G_k/T on each k-point.

     k-points are described in the crystalline coordinate like
     "k=( 0.50000, 0.00000,-0.00000),".   

     The factor system on each k-point is described as a matrix.
     Each row and column represents the index of symmetry operators described in standard output.
     Check "qeirreps.log" to confirm the indices of symmetry operations.

     If the symmetry operator does not exist in G_k/T, the value of factor system is described as zero.


2. For Mathematica usage, "*_import.txt"

   See and edit the Mathematica notebook file (i.e., "example/Si_soc/Si.nb")
   to confirm the information of each file.

2-1. rg_import.txt: imported to "pgList".
     The rotation part of symmetry operators.

     pgList[[i]]: The 3 by 3 matrix of the rotation part of the i-th symmetry operator in coordinate space.

2-2. tg_import.txt: imported to "tgList".
     The translation part of symmetry operators.

     tgList[[i]]: The 3 dimensional vector of the translation part of the i-th symmetry operator.

2-3. srg_import.txt: imported to "srgList".
     The spin rotation part of symmetry operators.
     You can ignore this when you run calculation without spin-orbit coupling.

     srgList[[i]]: The 2 by 2 matrix of the spin rotation part of the i-th symmetry operator.

2-4. kg_import.txt: imported to "kOpeList".
     The symmetry operator indices which are elements of G_k on each k-point.
   
     kOpeList[[i,1]]: The 3 dimensional k-vector of the i-th k-point.
     kOpeList[[i,2]]: The symmetry operator index which belongs to G_k of the i-th k-point.

2-5. character_import.txt: imported to "BandList".
     The characters for G_k of each band.
     (REMARK: These are not characters for G_k/T but for G_k 
	      where T is the lattice translational group.
              Then these characters contain phase factors
              as Exp[-i \vec{k}\cdot\vec{t}_g].)

     BandList[[i,1]]: The 3 dimensional k-vector of the i-th k-point.
     BandList[[i,2,n,j]]: The character of the n-th band for the j-th symmetry operator. 
   		        If the j-th symmetry is not an element of G_k on the i-th k-point,
		        this value becomes "Null".

2-6. factor_system_nonsymmorphic_import.txt: imported to "factorSystemNS".
     The nonsymmorphic part of the factor system of G_k/T on each k-point.
     This value is defined as follows.

     \omega(g,g') = Exp[i \vec{k}\cdot(\vec{t}_{g'} - p_g\vec{t}_{g'})]
 
     factorSystemNS[[i,1]]: 
  	The 3 dimensional k-vector of the i-th k-point.
     factorSystemNS[[i,2,m,n]]: 
  	The nonsymmorphic part of the factor system 
  	between the m-th symmetry operator and n-th one.
     	If the m-th or n-th symmetry is not an element of G_k on the i-th k-point,
  	this value becomes "Null".
   
2-7. factor_system_spin_import.txt: imported to "factorSystemSpin".
     The spin part of the factor system of G_k/T on each k-point.

     factorSystemSpin[[i,1]]: 
  	The 3 dimensional k-vector of the i-th k-point.
     factorSystemSpin[[i,2,m,n]]:
 	The spin part of the factor system
        between the m-th symmetry operator and n-th one.
	If the m-th or n-th symmetry is not an element of G_k on the i-th k-point,
	this value becomes "Null".
						      
======================================
