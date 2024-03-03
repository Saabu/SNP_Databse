#!/usr/bin/env bash
#SBATCH --nodes=1
#SBATCH --cpus-per-task=1
#SBATCH --ntasks-per-node=1
#SBATCH --mem=6G
#SBATCH --array=1-21
#SBATCH -J vcf_filter
#SBATCH -o vcf_filter_%a.out
#SBATCH -e vcf_filter_%a.err
#SBATCH --time=8:00:00


module load build-env/2020
module load vcftools/0.1.16-foss-2018b-perl-5.28.0

#wd="/scratch/ronit.chakraborty/snp_project/1000genome/1kGP"
od1="/scratch/ronit.chakraborty/snp_project/1000genome/filtered_vcf/no_indel"
od2="/scratch/ronit.chakraborty/snp_project/1000genome/filtered_vcf/final_vcf"

#cd ${wd}

#file=$(ls *.filtered.SNV_INDEL_SV_phased_panel.vcf.gz | sed -n ${SLURM_ARRAY_TASK_ID}p)

#echo ${file}
#vcftools --gzvcf ${file} --remove-indels --recode --recode-INFO-all --stdout | gzip -c > ${file%.filtered.SNV_INDEL_SV_phased_panel.vcf.gz}.out_snps.vcf.gz
#mv ${wd}/*.out_snps.vcf.gz ${od1}

cd ${od1}

file=$(ls *.out_snps.vcf.gz | sed -n ${SLURM_ARRAY_TASK_ID}p)

echo ${file}
vcftools --gzvcf ${file} --remove-filtered-all --recode --stdout | gzip -c > ${file%.out_snps.vcf.gz}.output_PASS_only.vcf.gz
mv ${od1}/*.output_PASS_only.vcf.gz ${od2}
