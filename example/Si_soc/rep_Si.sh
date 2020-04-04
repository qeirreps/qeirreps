#!/bin/sh
#QSUB -queue i18cpu
#QSUB -node 1
#QSUB -mpi 24
#QSUB -omp 1
#QSUB -over false
#QSUB -place pack
#PBS -l walltime=00:30:00
#PBS -N Si_rep
cd ${PBS_O_WORKDIR}
date
mkdir /work/rep/Si
cp -r /work/scf/Si/* /work/rep/Si
mpijob ~/research/dft/qe-6.2/bin/pw.x < Si.rep.in > Si.rep.out
date
