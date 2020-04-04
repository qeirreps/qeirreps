#!/bin/sh
#QSUB -queue i18cpu
#QSUB -node 1
#QSUB -mpi 24
#QSUB -omp 1
#QSUB -over false
#QSUB -place pack
#PBS -l walltime=00:30:00
#PBS -N Bi_ns
cd ${PBS_O_WORKDIR}
date
mkdir /work//scf/Bi_ns
mpijob ~/research/dft/qe-6.2/bin/pw.x < Bi.scf.in > Bi.scf.out
date
mkdir /work/band/Bi_ns
cp -r /work/scf/Bi_ns/* /work/band/Bi_ns
mpijob ~/research/dft/qe-6.2/bin/pw.x < Bi.band.in > Bi.band.out
date
mpijob ~/research/dft/qe-6.2/bin/bands.x < Bi.bands.in > Bi.bands.out
date
