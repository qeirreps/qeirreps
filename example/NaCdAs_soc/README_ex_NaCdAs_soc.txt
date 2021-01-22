# Copyright (c) 2020 Akishi Matsugatani, Seishiro Ono, Yusuke Nomura, Haruki Watanabe

======================================

Here, qeirreps/example/NaCdAs_soc is the directory which stores the input files for DFT calculation of NaCdAs with spin-orbit coupling.

See also qeirreps/reference/NaCdAs_soc to confirm the result.

======================================

# Contents

NaCdAs.scf.in:
The input file of QE for self-consistent first-principles (scf) calculation

NaCdAs.band.in:
The input file of QE for non-self-consistent first-principles (nscf) calculation
to obtain the band structure

NaCdAs.bands.in:
The input file of QE to obtain the plot of the band structure

NaCdAs.rep.in:
The input file of QE for non-self-consistent first-principles (nscf) calculation
to obtain the irreducible representation

output:
The empty directory to store the result files of "qeirreps"
"qeirreps" will write the result files here.

for_CTM:
The directory for the calculation with the "ctm" option to export the input file of CheckTopologicalMat.
This directory also has NaCdAs.rep.in and output.
In the case of NaCdAs, the choice of k-points is the same in both cases to export z4.dat and Bilbao.txt.

======================================

# To obtain band structure

Use original functions of QE as follows.

1. scf calculation by QE

   Set "outdir" in "NaCdAs.scf.in."
   (e.g., outdir='./work/scf')
   
   Run "pw.x" in QE with the input file "NaCdAs.scf.in."

2. nscf calculation by QE
   
   Set "outdir" in "NaCdAs.band.in."
   (e.g., outdir='./work/band')

   Copy the output files of scf calculation to the reading directory of "NaCdAs.band.in."
   (e.g., cp -r ./work/scf/* ./work/band/)

   Run "pw.x" in QE with the input file "NaCdAs.band.in."

3. band calculation by QE

   Set "outdir" in "NaCdAs.bands.in."
   (e.g., outdir='./work/band')

   Run "bands.x" in QE with the input file "NaCdAs.bands.in."

   "NaCdAs.band," "NaCdAs.band.gnu," "NaCdAs.band.rap" will be generated.

   One can plot the band structure from "NaCdAs.band.gnu" by using gnuplot.

======================================

# To obtain character table

Use QE, qe2respack, and qeirreps as follows.


1. scf calculation by QE

   Set "outdir" in "NaCdAs.scf.in."
   (e.g., outdir='./work/scf')

   Run "pw.x" in QE with the input file "NaCdAs.scf.in."

   REMARK:
   Norm-conserving calculations are necessary.
   Set the option "wf_collect = .TRUE."
   Use the pseudopotentials optimized for norm-conserving calculations. See README_psd in the directory "pseudo" for detailed information.

2. nscf calculation by QE
   
   Set "outdir" in "NaCdAs.rep.in."
   (e.g., outdir='./work/rep')

   Copy the output files of scf calculation to the reading directory of "NaCdAs.rep.in"
   (e.g., cp -r ./work/scf/* ./work/rep/)

   Run "pw.x" in QE with the input file "NaCdAs.rep.in."

3. Preparation of input files for qeirreps by qe2respack

   Run qe2respack by typing as

   	$ python PATH_OF_qe2respack/qe2respack.py OUTDIR/PREFIX.save

   (e.g. "python ~/respackDev-master/util/qe2respack/qe2respack.py ./work/rep/NaCdAs.save")

   "PATH_OF_qe2respack" is the directory which has the file qe2respack.py.
   "OUTDIR/PREFIX.save" is the directory produced by QE in step 2.
   
   Nine files including "dat.wfn" will be generated in the directory "dir-wfn." qeirreps reads these files in step 4.


4. Calculation of the character tables and Z4 index by qeirreps

   Run qeirreps by typing as follows.

	$ ../../src/qeirreps.x . 176 z4 > qeirreps.log

   "../../src/" is the location of the qeirreps executable file.
   "." is the current directory "NaCdAs_soc" which contains "dir-wfn" and "output," and "qeirreps.log" is the log file of standard output.
   "176" is the filling of the NaCdAs. This value is produced in the standard output of scf calculation by QE as "number of electrons = 176.00."
   "z4" is the option to obtain the Z4-index.

   Some text files will be exported in "output," for example, "character.dat."
   Check the document and files in reference directory for more information.

5. Exporting the input file for CheckTopologicalMat

   Move to the directory "for_CTM" and run nscf calculation as explained in step 2.
   Note that the result file of the scf calculation is needed for the nscf calculation, so that one should copy the file OUTDIR/PREFIX.save to the work directory.
   (e.g., cp -r ../work/scf/* ./work/rep/)
   The choice of k-points is the same in both cases to export z4.dat and Bilbao.txt.

   Run qe2respack as described in step 3.

   Run qeirreps by typing as follows.

	$ ../../../src/qeirreps.x . 176 ctm > qeirreps.log

   "ctm" is the option to export the input file of CheckTopologicalMat.

   The additional file named "Bilbao.txt" will be generated in "output."
   This can be used as an input file of CheckTopologicalMat.


======================================
