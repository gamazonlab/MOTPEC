#!/bin/bash
#SBATCH --mail-user=yuexu0312@qq.com
#SBATCH --mail-type=ALL
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --constraint=haswell|skylake #sandybridge|
#SBATCH --time=4:00:00
#SBATCH --mem-per-cpu=32G
#SBATCH --job-name=common_MAGMA_run
#SBATCH --array=1-49
##SBATCH --account=g_gamazon_lab

ml GCC/8.2.0  OpenMPI/3.1.4 Intel/2019.1.144  IntelMPI/2018.4.274 R/3.6.0
cd /gpfs52/data/g_gamazon_lab/zhoud2/cross_tissue_prediction/R_code/BMI/predixcan
Rscript run_predixcan.R ${SLURM_ARRAY_TASK_ID}
