# Copyright (c) 2020 Akishi Matsugatani, Seishiro Ono, Yusuke Nomura, Haruki Watanabe

======================================

Here, qeirreps/example/PbPt3_nonsoc is the directory which stores the input files for DFT calculation of PbPt3 without spin-orbit coupling.

See also qeirreps/reference/PbPt3_nonsoc to confirm the result.

======================================

# Contents

PbPt3.scf.in:
The input file of QE for self-consistent first-principles (scf) calculation

PbPt3.band.in:
The input file of QE for non-self-consistent first-principles (nscf) calculation
to obtain the band structure

PbPt3.bands.in:
The input file of QE to obtain the plot of the band structure

PbPt3.rep.in:
The input file of QE for non-self-consistent first-principles (nscf) calculation
to obtain the irreducible representation

dir-wfn:
The empty directory to store the input files for "qeirreps"
"qe2respack" will write the result files here

output:
The empty directory to store the result files of "qeirreps"
"qeirreps" will write the result files here

======================================

# To obtain band structure

Use original functions of QE as follows.

1. scf calculation by QE

   Set "outdir" in "PbPt3.scf.in".
   (e.g., outdir='./work/scf')
   
   Run "pw.x" in QE with the input file "PbPt3.scf.in".

2. nscf calculation by QE
   
   Set "outdir" in "PbPt3.band.in".
   (e.g., outdir='./work/band')

   Copy the output files of scf calculation to the reading directory of "PbPt3.band.in".
   (i.e., cp -r ./work/scf/* ./work/band/)

   Run "pw.x" in QE with the input file "PbPt3.band.in".

3. band calculation by QE

   Set "outdir" in "PbPt3.bands.in".
   (i.e., outdir='./work/band')

   Run "band.x" in QE with the input file "PbPt3.bands.in".

   "PbPt3.band", "PbPt3.band.gnu", "PbPt3.band.rap" will be generated.

   You can plot the band structure from "PbPt3.band.gnu" by using gnuplot.

======================================

# To obtain character table

Use QE, qe2respack, and qeirreps as follows.


1. scf calculation by QE

   Set "outdir" in "PbPt3.scf.in".
   (e.g., outdir='./work/scf')

   Run "pw.x" in QE with the input file "PbPt3.scf.in".

   REMARK:
   Norm-conserving calculations are necessary.
   Set the option "wf_collect = .TRUE."
   Use the pseudo potentials optimized for norm-conserving calculations. See README_psd in the directory "pseudo" for detailed information.

2. nscf calculation by QE
   
   Set "outdir" in "PbPt3.rep.in".
   (e.g., outdir='./work/rep')

   Copy the output files of scf calculation to the reading directory of "PbPt3.band.in"
   (i.e., cp -r ./work/scf/* ./work/rep/)

   Run "pw.x" in QE with the input file "PbPt3.rep.in".

3. Preparation of input files for qeirreps by qe2respack

   Run qe2respack by typing as

   	$ PATH_OF_qe2respack/qe2respack OUTDIR/PREFIX.save

   (i.e. "~/respackDev-master/util/qe2respack/qe2respack ./work/rep/PbPt3.save")

   "PATH_OF_qe2respack" is the directory which has executable file of qe2respack.
   "OUTDIR/PREFIX.save" is the directory produced by QE in step 2.
   
   Eight files including "dat.wfn" will be generated in the directory "dir-wfn". qeirreps reads these files in step 4.


4. Calculation of the character tables by qeirreps

   Run qeirreps by typing as follows

	$ ../../src/qeirreps.x . > qeirreps.log

   "../../src/" is the location of the qeirreps executable file.
   "." is the current directory "PbPt3_nonsoc" which contains "dir-wfn" and "output", and "qeirreps.log" is the log file of standard output.

   Some text files will be exported in "output", for example, "character_import.txt".
   Check the document in reference directory "qeirreps/reference/README_ref.txt" for more information.

======================================
