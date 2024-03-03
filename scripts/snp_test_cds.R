setwd("/lisc/scratch/kovarik/Ronit/snp_project/final_vcf")
library(VariantAnnotation)
library(GenomicFeatures)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
library(ggplot2)
library(biomaRt)
library(dplyr)
library(tidyr)

fi = "1kGP_high_coverage_Illumina.chr22.output_PASS_only.vcf"
hdr = scanVcfHeader(fi)
#info_list = hdr@header@listData$INFO@rownames
param = ScanVcfParam(info = NA, geno=NA)
vcf = readVcf(fi, "hg38", param = param)

txdb = TxDb.Hsapiens.UCSC.hg38.knownGene

intersect(seqlevels(vcf), seqlevels(txdb))

loc_all <- locateVariants(vcf, txdb, CodingVariants())

threeutr = unique(as.data.frame(loc_all@ranges))

gene_ids = na.omit(as.data.frame(unique(loc_all$GENEID)))

ensembl=useMart("ensembl")
datasets <- listDatasets(ensembl)
ensembl = useDataset("hsapiens_gene_ensembl", mart = ensembl)
entrzID = gene_ids$`unique(loc_all$GENEID)`
filters = listFilters(ensembl)
attributes = listAttributes(ensembl)

#biomartCacheClear()

#Normalised count by transcript----
gene_list = na.omit(unique(loc_all$GENEID))

#gene_list = gene_list[801:841]

df_all <- data.frame()

for (i in gene_list) {
  result_id <- getBM(attributes = c("entrezgene_id", "hgnc_symbol", "entrezgene_description"), 
                     filters = "entrezgene_id", values = i, mart = ensembl)
  
  result_cds <- na.omit(getBM(attributes = c("genomic_coding_start", "genomic_coding_end"), 
                             filters = "entrezgene_id", values = i, mart = ensembl))
  
  if (nrow(result_cds) > 0) {
    result_cds$count <- 0
    result_cds$utr_length <- 0
    result_cds$norm_count <- 0
    
    for (j in 1:nrow(result_cds)) {
      utr_start <- result_cds$`genomic_coding_start`[j]
      utr_end <- result_cds$`genomic_coding_end`[j]
      
      # Add error handling using tryCatch
      tryCatch({
        if (!is.na(utr_start) && !is.na(utr_end)) {
          overlapping_rows <- threeutr %>%
            filter(!((end <= utr_start) | (start >= utr_end)))
          
          count <- nrow(overlapping_rows)
          
          utr_length <- utr_end - utr_start
          if (utr_length > 0) {
            norm_count <- round(count / (utr_length / 1000))
            result_cds$utr_length[j] <- utr_length
            result_cds$count[j] <- count
            result_cds$norm_count[j] <- round(norm_count)
          }
        }
      }, error = function(e) {
        # Handle the error, e.g., print a message or log it
        cat("Error in loop iteration:", j, "for gene:", i, "-", e$message, "\n")
      })
      
    }
    
    hgnc_symbol <- result_id$hgnc_symbol[!is.na(result_id$hgnc_symbol)][1]
    
    if (!is.na(hgnc_symbol)) {
      result_cds$hgnc_symbol <- hgnc_symbol
      df_all <- bind_rows(df_all, result_cds)
    }
  }
}

df_all <- df_all %>% filter(count != 0)
df_all = subset(df_all, select = -count)
colnames(df_all) = c("genomic_coding_start", "genomic_coding_end", "cds_length", "norm_count", "hgnc_symbol")
#df_all = subset(df_all, select = -entrezgene_id)
#df_all = subset(df_all, select = -entrezgene_description)

#df_all$freq <- round((df_all$count / sum(df_all$count)) * 100, 2)

setwd("/scratch/kovarik/Ronit/snp_project/output_Rds_CDS")
saveRDS(df_all, file="1kGP_high_coverage_Illumina.chr22.output_PASS_only.vcf.Rds")






















#Old Script---------------------------------------------
vcf = readVcf("result_new.tsv", "hg38")
txdb = TxDb.Hsapiens.UCSC.hg38.knownGene

intersect(seqlevels(vcf), seqlevels(txdb))

vcf <- renameSeqlevels(vcf, paste0("chr", seqlevels(vcf)))
intersect(seqlevels(vcf), seqlevels(txdb))

loc_all <- locateVariants(vcf, txdb, CodingVariants())

threeutr = unique(as.data.frame(loc_all@ranges))

gene_ids = na.omit(as.data.frame(unique(loc_all$GENEID)))

ensembl=useMart("ensembl")
datasets <- listDatasets(ensembl)
ensembl = useDataset("hsapiens_gene_ensembl", mart = ensembl)
entrzID = gene_ids$`unique(loc_all$GENEID)`
filters = listFilters(ensembl)
attributes = listAttributes(ensembl)

gene_list = na.omit(unique(loc_all$GENEID))

df_all <- data.frame()

total_snps = nrow(threeutr)

filter_count = 0

for (i in gene_list) {
  
  result_id <- getBM(attributes = c("entrezgene_id", "hgnc_symbol", "entrezgene_description"), 
                     filters = "entrezgene_id", values = i, mart = ensembl)
  
  result_cds <- na.omit(getBM(attributes = c("genomic_coding_start", "genomic_coding_end"), filters = "entrezgene_id", values = i, mart = ensembl))
  
  
  mean_start = round(mean(result_cds$genomic_coding_start))
  mean_end = round(mean(result_cds$genomic_coding_end))
  
  result_id$cds_start = mean_start
  result_id$cds_end = mean_end

  count <- threeutr %>%
    filter(!((end <= result_id$cds_start) | (start >= result_id$cds_end))) %>%
    nrow()
  
  #freq = (count/total_snps)*100
  
  cds_length <- result_id$cds_end - result_id$cds_start
  norm_count <- round(count / (cds_length / 1000))
  
  filter_count = filter_count + norm_count  
}

for (j in gene_list) {
  
  result_id <- getBM(attributes = c("entrezgene_id", "hgnc_symbol", "entrezgene_description"), 
                     filters = "entrezgene_id", values = j, mart = ensembl)
  
  result_cds <- na.omit(getBM(attributes = c("genomic_coding_start", "genomic_coding_end"), filters = "entrezgene_id", values = j, mart = ensembl))
  
  
  mean_start = round(mean(result_cds$genomic_coding_start))
  mean_end = round(mean(result_cds$genomic_coding_end))
  
  result_id$cds_start = mean_start
  result_id$cds_end = mean_end
  
  count <- threeutr %>%
    filter(!((end <= result_id$cds_start) | (start >= result_id$cds_end))) %>%
    nrow()
  
  #freq = (count/total_snps)*100
  
  cds_length <- result_id$cds_end - result_id$cds_start
  norm_count <- count / (cds_length / 1000)
  norm_freq <- (norm_count / filter_count) * 100
  
  result_id$cds_length = round(cds_length)
  result_id$count = count
  #result_id$freq = round(freq, 2)
  result_id$norm_count = round(norm_count)
  result_id$norm_freq = round(norm_freq, 2)
  
  df_all = rbind(df_all, result_id)
} 

df_all = subset(df_all, count != 0)  
df_all = df_all[order(df_all$norm_count, decreasing = TRUE), ]

p <- ggplot(df_all, aes(x = reorder(hgnc_symbol, -norm_count), y = norm_count)) +
  ggtitle("SNPs distribution across genes in Y chromosome in CDS") +
  geom_bar(stat = "identity", fill = "steelblue") +
  labs(x = "Gene", y = "Normalised Count") +
  ylim(0, 355) +
  geom_text(aes(label = norm_freq), vjust = -0.5, color = "black", size = 4, angle = 45, hjust = -0.05) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggsave("cds_gene.pdf", plot = p, width = 10, height = 8, units = "in")


