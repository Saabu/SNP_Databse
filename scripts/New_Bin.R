library(VariantAnnotation)
library(GenomicFeatures)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
library(ggplot2)
library(biomaRt)
library(dplyr)
library(tidyr)
library(writexl)

setwd("/scratch/kovarik/Ronit/snp_project/final_vcf")

fi = "1kGP_high_coverage_Illumina.chr3.output_PASS_only.vcf"
hdr = scanVcfHeader(fi)
param = ScanVcfParam(info = NA, geno=NA)
vcf = readVcf(fi, "hg38", param = param)

print("VCF file loaded successfully")

txdb = TxDb.Hsapiens.UCSC.hg38.knownGene

intersect(seqlevels(vcf), seqlevels(txdb))

loc_all <- locateVariants(vcf, txdb, ThreeUTRVariants())

threeutr = unique(as.data.frame(loc_all@ranges))

gene_ids = na.omit(as.data.frame(unique(loc_all$GENEID)))

ensembl=useMart("ensembl")
datasets <- listDatasets(ensembl)
ensembl = useDataset("hsapiens_gene_ensembl", mart = ensembl)
entrzID = gene_ids$`unique(loc_all$GENEID)`
filters = listFilters(ensembl)
attributes = listAttributes(ensembl)

print("ensembl loaded successfully")

#Count transcript----
gene_list = na.omit(unique(loc_all$GENEID))

#gene_list = gene_list[1:10]

df_175 = data.frame()

print(length(gene_list))

for (i in gene_list) {
  result_id <- getBM(attributes = c("entrezgene_id", "hgnc_symbol", "entrezgene_description"), 
                     filters = "entrezgene_id", values = i, mart = ensembl)
  
  result_3p <- na.omit(getBM(attributes = c("3_utr_start", "3_utr_end"), 
                             filters = "entrezgene_id", values = i, mart = ensembl))
  
  if (nrow(result_3p) > 0) {
    result_3p$count_50bp <- 0
    result_3p$count_175_225bp <- 0
    result_3p$utr_length <- 0
    
    for (j in 1:nrow(result_3p)) {
      utr_start <- result_3p$`3_utr_start`[j]
      utr_end <- result_3p$`3_utr_end`[j]
      
      utr_length <- utr_end - utr_start
      if (utr_length >= 300 && utr_length <= 600) {
        
        # Add error handling using tryCatch
        tryCatch({
          if (!is.na(utr_start) && !is.na(utr_end)) {
            overlapping_rows <- threeutr %>%
              filter(!((end <= utr_start) | (start >= utr_end)))
            
            count_50bp <- sum(overlapping_rows$start >= utr_start & overlapping_rows$start <= utr_start + 50)
            count_175_225bp <- sum(overlapping_rows$start >= utr_start + 175 & overlapping_rows$start <= utr_start + 225)
            
            result_3p$count_50bp[j] <- count_50bp
            result_3p$count_175_225bp[j] <- count_175_225bp
            result_3p$utr_length[j] <- utr_length
            
          }
        }, error = function(e) {
          # Handle the error, e.g., print a message or log it
          cat("Error in loop iteration:", j, "for gene:", i, "-", e$message, "\n")
        })
      }
    }
    
    hgnc_symbol <- result_id$hgnc_symbol[!is.na(result_id$hgnc_symbol)][1]
    
    if (!is.na(hgnc_symbol)) {
      result_3p$hgnc_symbol <- hgnc_symbol
      result_3p <- result_3p %>%
        select(hgnc_symbol, utr_length, count_50bp, count_175_225bp)
      df_175 <- bind_rows(df_175, result_3p)
      df_175 <- df_175 %>% filter(count_50bp != 0 | count_175_225bp != 0)
    }
  }
}

print("df_175 generated successfully")

setwd("/lisc/user/chakraborty/snp_project/new_bin/50_175bp_bins/output_rds")
saveRDS(df_175, file="Bin_175_chr3.Rds")

print("RDS File 50-175 for chr 3 saved successfully")
