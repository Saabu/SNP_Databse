#!/usr/bin/env bash
#SBATCH --nodes=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH -J download
#SBATCH -o "%x"."%j".out
#SBATCH -e "%x"."%j".err
#SBATCH --time=8:00:00


for chr in {1..22}
do
wget -P /home/user/chakraborty/snp_project/files http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000G_2504_high_coverage/working/20220422_3202_phased_SNV_INDEL_SV/1kGP_high_coverage_Illumina.chr${chr}.filtered.SNV_INDEL_SV_phased_panel.vcf.gz
done
