# Copyright (c) 2020 Akishi Matsugatani, Seishiro Ono, Yusuke Nomura, Haruki Watanabe

======================================

Here, qeirreps/example/PbPt3_soc is the directory which stores the input files for DFT calculation of PbPt3 with spin-orbit coupling.

See also qeirreps/reference/PbPt3_soc to confirm the result.

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

output:
The empty directory to store the result files of "qeirreps"
"qeirreps" will write the result files here.

REMARK:
In the case of PbPt3, the processes of qeirreps are done with the "ctm" option.
The choice of k-points in PbPt3.rep.in, PbPt3.rep.out, and dir-wfn is specified for the calculation to export the input file of CheckTopologicalMat.
This input file is named Bilbao.txt and stored in output directory.

======================================

# To obtain band structure

Use original functions of QE as follows.

1. scf calculation by QE

   Set "outdir" in "PbPt3.scf.in."
   (e.g., outdir='./work/scf')
   
   Run "pw.x" in QE with the input file "PbPt3.scf.in."

2. nscf calculation by QE
   
   Set "outdir" in "PbPt3.band.in."
   (e.g., outdir='./work/band')

   Copy the output files of scf calculation to the reading directory of "PbPt3.band.in."
   (e.g., cp -r ./work/scf/* ./work/band/)

   Run "pw.x" in QE with the input file "PbPt3.band.in."

3. band calculation by QE

   Set "outdir" in "PbPt3.bands.in."
   (e.g., outdir='./work/band')

   Run "bands.x" in QE with the input file "PbPt3.bands.in."

   "PbPt3.band," "PbPt3.band.gnu," "PbPt3.band.rap" will be generated.

   One can plot the band structure from "PbPt3.band.gnu" by using gnuplot.

======================================

# To obtain character table

Use QE, qe2respack, and qeirreps as follows.


1. scf calculation by QE

   Set "outdir" in "PbPt3.scf.in."
   (e.g., outdir='./work/scf')

   Run "pw.x" in QE with the input file "PbPt3.scf.in."

   REMARK:
   Norm-conserving calculations are necessary.
   Set the option "wf_collect = .TRUE."
   Use the pseudopotentials optimized for norm-conserving calculations. See README_psd in the directory "pseudo" for detailed information.

2. nscf calculation by QE
   
   Set "outdir" in "PbPt3.rep.in."
   (e.g., outdir='./work/rep')

   Copy the output files of scf calculation to the reading directory of "PbPt3.rep.in"
   (e.g., cp -r ./work/scf/* ./work/rep/)

   Run "pw.x" in QE with the input file "PbPt3.rep.in."

3. Preparation of input files for qeirreps by qe2respack

   Run qe2respack by typing as

   	$ python PATH_OF_qe2respack/qe2respack.py OUTDIR/PREFIX.save

   (e.g. "python ~/respackDev-master/util/qe2respack/qe2respack.py ./work/rep/PbPt3.save")

   "PATH_OF_qe2respack" is the directory which has the file qe2respack.py.
   "OUTDIR/PREFIX.save" is the directory produced by QE in step 2.
   
   Nine files including "dat.wfn" will be generated in the directory "dir-wfn." qeirreps reads these files in step 4.


4. Calculation of the character tables by qeirreps

   Run qeirreps by typing as follows.

	$ ../../src/qeirreps.x . 68 ctm > qeirreps.log

   "../../src/" is the location of the qeirreps executable file.
   "." is the current directory "PbPt3_soc" which contains "dir-wfn" and "output," and "qeirreps.log" is the log file of standard output.
   "68" is the filling of the PbPt3. This value is produced in the standard output of scf calculation by QE as "number of electrons = 68.00."
   "ctm" is the option to export the input file of CheckTopologicalMat.

   Some text files will be exported in "output," for example, "character.dat."
   Check the document and files in reference directory for more information.
   Especially, the file named "Bilbao.txt" will be generated.
   This can be used as an input file of CheckTopologicalMat.

   REMARK:
   One can find a warning in the standard output as

	Warning: number of irreps is not an integer
	 at k=  0.000000000000000E+000  0.000000000000000E+000  0.000000000000000E+000 
	 and level=          28 :
	 (0.499999999999489,3.469446951953614E-017)

   This is the warning in the processes to count the irreducible representations in each energy level.
   This is the warning for the highest energy level in the unoccupied part. This indicates that the highest energy level has degeneracy, but only a part of degenerate bands are computed in the nscf calculation.
   There is no problem with the other energy levels, such as the occupied bands.


======================================
