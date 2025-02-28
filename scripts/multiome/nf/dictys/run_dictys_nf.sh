#!/bin/bash -l
#SBATCH --output=/scratch/users/%u/%j.out
#SBATCH --nodes=1
#SBATCH -n 1
#SBATCH --mem=10G
#SBATCH --time=48:00:00

#Load required packages
conda create --prefix ./nfEnv -y
source ~/.bashrc
source activate ./nfEnv
conda install -y -c bioconda nextflow=23.04.1
conda install -y -c conda-forge openjdk=17.0.9

#Make tmpdir
mkdir $PWD/singularity_tmp
mkdir $PWD/tmp_static

#Export some environmental variable to nextflow for singularity usage
export SINGULARITY_TMPDIR=$PWD/singularity_tmp
export dataDir=$PWD/data
export tmpDir=$PWD/tmp_static/
export subsetDir=$PWD/data/subsets

nextflow -C scripts/nf/dictys/dictys_nf.config run scripts/nf/dictys/dictys.nf
