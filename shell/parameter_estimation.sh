#!/bin/bash
#$ -N ebola_parameter_estimation
#$ -cwd -V
#$ -M <email>
#$ -m e
#$ -l mem_free=32G,h_vmem=32G
#$ -q short.q
#$ -pe smp 32
#$ -e ../output/logs/$JOB_NAME_error.log
#$ -o ../output/logs/$JOB_NAME_out.log
#$ -R y
export R_LIBS="~R/x86_64-redhat-linux-gnu-library/3.6"
export N_CORES=32
R CMD BATCH ../R/parameter_estimation.R ../output/logs/parameter_estimation.Rout
