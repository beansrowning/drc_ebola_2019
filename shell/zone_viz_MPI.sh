#!/bin/bash
#$ -N ebola_health_zone_viz
#$ -cwd -V
#$ -M <email>
#$ -m e
#$ -l mem_free=16G,h_vmem=16G
#$ -q all.q
#$ -pe openmpi 32
#$ -e ../output/$JOB_NAME_error.log
#$ -o ../output/$JOB_NAME_out.log
#$ -R y

export R_LIBS="~/R/x86_64-pc-linux-gnu-library/3.3"
module load openmpi # Base MPI bin and headers
module load openmpi/1.8.1-ib-tcp-gcc-64 # Infiniband vocab
mpirun -np 32 R CMD BATCH ../R/zone_viz_MPI.R ../output/zone_viz_MPI.Rout
