#!/bin/bash
#$ -N ebola_simulation
#$ -cwd -V
#$ -M <email>
#$ -m e
#$ -l mem_free=32G,h_vmem=32G
#$ -q all.q
#$ -pe smp 48
#$ -e output/logs/$JOB_NAME_error.log
#$ -o output/logs/$JOB_NAME_out.log
#$ -R y
export R_LIBS="~R/x86_64-redhat-linux-gnu-library/3.6"
export N_CORES=48

module load R/3.6.1
R CMD BATCH R/simulation.R output/logs/ebola_simulation.Rout
