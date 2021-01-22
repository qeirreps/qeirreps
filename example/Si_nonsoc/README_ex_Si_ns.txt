# Copyright (c) 2020 Akishi Matsugatani, Seishiro Ono, Yusuke Nomura, Haruki Watanabe

======================================

Here, qeirreps/example/Si_nonsoc is the directory which stores the input files for DFT calculation of silicon without spin-orbit coupling.

See also qeirreps/reference/Si_nonsoc to confirm the result.

======================================

# Contents

Si.scf.in:
The input file of QE for self-consistent first-principles (scf) calculation

Si.band.in:
The input file of QE for non-self-consistent first-principles (nscf) calculation
to obtain the band structure

Si.bands.in:
The input file of QE to obtain the plot of the band structure

Si.rep.in:
The input file of QE for non-self-consistent first-principles (nscf) calculation
to obtain the irreducible representation

output:
The empty directory to store the result files of "qeirreps"
"qeirreps" will write the result files here.

======================================

# To obtain band structure

Use original functions of QE as follows.

1. scf calculation by QE

   Set "outdir" in "Si.scf.in."
   (e.g., outdir='./work/scf')
   
   Run "pw.x" in QE with the input file "Si.scf.in."

2. nscf calculation by QE
   
   Set "outdir" in "Si.band.in."
   (e.g., outdir='./work/band')

   Copy the output files of scf calculation to the reading directory of "Si.band.in."
   (e.g., cp -r ./work/scf/* ./work/band/)

   Run "pw.x" in QE with the input file "Si.band.in."

3. band calculation by QE

   Set "outdir" in "Si.bands.in."
   (e.g., outdir='./work/band')

   Run "bands.x" in QE with the input file "Si.bands.in."

   "Si.band," "Si.band.gnu," "Si.band.rap" will be generated.

   One can plot the band structure from "Si.band.gnu" by using gnuplot.

======================================

# To obtain character table

Use QE, qe2respack, and qeirreps as follows.


1. scf calculation by QE

   Set "outdir" in "Si.scf.in."
   (e.g., outdir='./work/scf')

   Run "pw.x" in QE with the input file "Si.scf.in."

   REMARK:
   Norm-conserving calculations are necessary.
   Set the option "wf_collect = .TRUE."
   Use the pseudopotentials optimized for norm-conserving calculations. See README_psd in the directory "pseudo" for detailed information.

2. nscf calculation by QE
   
   Set "outdir" in "Si.rep.in."
   (e.g., outdir='./work/rep')

   Copy the output files of scf calculation to the reading directory of "Si.rep.in"
   (e.g., cp -r ./work/scf/* ./work/rep/)

   Run "pw.x" in QE with the input file "Si.rep.in."

3. Preparation of input files for qeirreps by qe2respack

   Run qe2respack by typing as

   	$ python PATH_OF_qe2respack/qe2respack.py OUTDIR/PREFIX.save

   (e.g. "python ~/respackDev-master/util/qe2respack/qe2respack.py ./work/rep/Si.save")

   "PATH_OF_qe2respack" is the directory which has the file qe2respack.py.
   "OUTDIR/PREFIX.save" is the directory produced by QE in step 2.
   
   Nine files including "dat.wfn" will be generated in the directory "dir-wfn." qeirreps reads these files in step 4.


4. Calculation of the character tables by qeirreps

   Run qeirreps by typing as follows.

	$ ../../src/qeirreps.x . > qeirreps.log

   "../../src/" is the location of the qeirreps executable file.
   "." is the current directory "Si_nonsoc" which contains "dir-wfn" and "output," and "qeirreps.log" is the log file of standard output.

   Some text files will be exported in "output," for example, "character.dat."
   Check the document and files in reference directory for more information.

======================================
