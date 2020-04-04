#!/bin/sh
#QSUB -queue i18cpu
#QSUB -node 1
#QSUB -mpi 24
#QSUB -omp 1
#QSUB -over false
#QSUB -place pack
#PBS -l walltime=00:30:00
#PBS -N Si
cd ${PBS_O_WORKDIR}
date
mkdir /work/scf/Si
mpijob ~/research/dft/qe-6.2/bin/pw.x < Si.scf.in > Si.scf.out
date
mkdir /work/band/Si
cp -r /work/scf/Si/* /work/band/Si
mpijob ~/research/dft/qe-6.2/bin/pw.x < Si.band.in > Si.band.out
date
mpijob ~/research/dft/qe-6.2/bin/bands.x < Si.bands.in > Si.bands.out
date
