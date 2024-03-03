#!/usr/bin/env bash
#SBATCH --nodes=1
#SBATCH --qos=medium
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH -J vcf
#SBATCH -o "%x"."%j".out
#SBATCH -e "%x"."%j".err
#SBATCH --time=2-00:00:00

module load build-env/2020
module load vcftools/0.1.16-foss-2018b-perl-5.28.0

cd /users/ronit.chakraborty/snp_project

vcftools --vcf All_20180418.vcf --remove-indels --recode --recode-INFO-all --out All_SNPs.vcf
