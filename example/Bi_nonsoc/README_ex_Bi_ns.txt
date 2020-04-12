# Copyright (c) 2020 Akishi Matsugatani, Seishiro Ono, Yusuke Nomura, Haruki Watanabe

======================================

Here, qeirreps/example/Bi_nonsoc is the directory 
which stores the input files for DFT calculation of bismuth
without spin-orbit coupling.

Bi.scf.in:
The input file of QE for self-consistent first-principle (scf) calculation

Bi.band.in:
The input file of QE for non-self-consistent first-principle (nscf) calculation
to obtain the band structure

Bi.bands.in:
The input file of QE to obtain the plot of the band structure

Bi.rep.in:
The input file of QE for non-self-consistent first-principle (nscf) calculation
to obtain the irreducible representation

dir-wfn:
The empty directory to store the input files for "qeirreps"
"qe2respack" will write the result files of it here

output:
The empty directory to store the result files of "qeirreps"
"qeirreps" will write the result files of it here

======================================

# To obtain band structure

Use ordinary functions of QE as follows

1. scf calculation by QE

   Set "outdir" in "Bi.scf.in" for your situation
   (i.e., outdir='./work/scf')
   
   Run "pw.x" in QE with the input file "Bi.scf.in"

2. nscf calculation by QE
   
   Set "outdir" in "Bi.band.in" for your situation
   (i.e., outdir='./work/band')

   Copy the output files of scf calculation 
   to the reading directory of "Bi.band.in"
   (i.e., cp -r ./work/scf/* ./work/band/)

   Run "pw.x" in QE with the input file "Bi.band.in"

3. band calculation by QE

   Set "outdir" in "Bi.bands.in" for your situation
   (i.e., outdir='./work/band')

   Run "band.x" in QE with the input file "Bi.bands.in"

   "Bi.band", "Bi.band.gnu", "Bi.band.rap" will be generated

   You can plot the band structure from "Bi.band.gnu" by using gnuplot

======================================

# To obtain character table

Use QE, qe2respack, and qeirreps as follows


1. scf calculation by QE

   Set "outdir" in "Bi.scf.in" for your situation
   (i.e., outdir='./work/scf')
   
   Run "pw.x" in QE with the input file "Bi.scf.in"

2. nscf calculation by QE
   
   Set "outdir" in "Bi.rep.in" for your situation
   (i.e., outdir='./work/rep')

   Copy the output files of scf calculation 
   to the reading directory of "Bi.band.in"
   (i.e., cp -r ./work/scf/* ./work/rep/)

   Run "pw.x" in QE with the input file "Bi.rep.in"

3. Preparation of input files for qeirreps by qe2respack

   Run qe2respack by typing as follows
   	$ PATH_OF_qe2respack/qe2respack OUTDIR/PREFIX.save
   (i.e. "~/respackDev-master/util/qe2respack/ ./work/rep/Bi.save")
   "PATH_OF_qe2respack" is the directory which has executable file of qe2respack
   "OUTDIR/PREFIX.save" is the directory produced by QE in step 2
   
   Some files (i.e. "dat.wfn"), will be generated in the directory "dir-wfn"
   qeirreps reads these files in the latter step


4. Calculation of the character tables by qeirreps

   Run qeirreps by typing as follows
	$ ./../../src/qeirreps.x . > qeirreps.log

   where "../../src/" is the location of the qeirreps executable file,
   "." is the current directory "Bi_nonsoc" which contains "dir-wfn" and "output",
   and "qeirreps.log" is the log file of standard output

   Some text files will be exported in "output", for example, "character_import.txt".
   They have various data to see irreducible representation.
   Standard output also show some data, for example, the lattice vectors, 
   the reciprocal lattice vectors, type of the symmetry operations, and so on.
   Check the document in reference directory "qeirreps/reference/README_ref.txt"
   for more information.

======================================
