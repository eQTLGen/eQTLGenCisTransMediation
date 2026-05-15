#!/bin/bash

#SBATCH --time=12:00:00
#SBATCH -N 1
#SBATCH --ntasks-per-node=1
#SBATCH --mem=6G
#SBATCH --mail-type=BEGIN
#SBATCH --mail-type=END
#SBATCH --mail-type=FAIL
#SBATCH --job-name="MRLink2"

# Here load needed system tools (Java 1.8 is required, one of singularity or anaconda - python 2.7 are needed,
# depending on the method for dependancy management)

module load Java/11.0.16

set -f

nextflow_path=[path to nextflow executable]

finemapped=[File with fine-mapped lead variants]
vcf_folder=[Folder with vcf.gz files]
expression_matrix=[Processed expression matrix]
covariate_matrix=[Covariate matrix]
covariates_to_test="GenPC1,GenPC2,GenPC3,GenPC4,ExpPC1,ExpPC2,ExpPC3,ExpPC4,ExpPC5,ExpPC6,ExpPC7,ExpPC8,ExpPC9,ExpPC10,ExpPC11,ExpPC12,ExpPC13,ExpPC14,ExpPC15,ExpPC16,ExpPC17,ExpPC18,ExpPC19,ExpPC20"

NXF_VER=23.04.4 ${nextflow_path}/nextflow run main.nf \
--finemapped ${finemapped} \
--vcf ${vcf_folder} \
--expression_matrix ${expression_matrix} \
--covariate_matrix ${covariate_matrix}\
--covariates ${covariates_to_test}\
--OutputDir results \
-profile slurm,singularity \
-resume
