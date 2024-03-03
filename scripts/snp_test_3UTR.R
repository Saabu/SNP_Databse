library(VariantAnnotation)
library(GenomicFeatures)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
library(ggplot2)
library(biomaRt)
library(dplyr)
library(tidyr)

#df_3utr = readRDS(file="data_3utr.Rds")

setwd("/scratch/kovarik/Ronit/snp_project/final_vcf")

fi = "1kGP_high_coverage_Illumina.chr22.output_PASS_only.vcf"
hdr = scanVcfHeader(fi)
#info_list = hdr@header@listData$INFO@rownames
param = ScanVcfParam(info = NA, geno=NA)
vcf = readVcf(fi, "hg38", param = param)

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

#Normalised count by transcript----
gene_list = na.omit(unique(loc_all$GENEID))

#gene_list = gene_list[801:841]

df_all <- data.frame()

for (i in gene_list) {
  result_id <- getBM(attributes = c("entrezgene_id", "hgnc_symbol", "entrezgene_description"), 
                     filters = "entrezgene_id", values = i, mart = ensembl)
  
  result_3p <- na.omit(getBM(attributes = c("3_utr_start", "3_utr_end"), 
                             filters = "entrezgene_id", values = i, mart = ensembl))
  
  if (nrow(result_3p) > 0) {
    result_3p$count <- 0
    result_3p$utr_length <- 0
    result_3p$norm_count <- 0
    
    for (j in 1:nrow(result_3p)) {
      utr_start <- result_3p$`3_utr_start`[j]
      utr_end <- result_3p$`3_utr_end`[j]
      
      # Add error handling using tryCatch
      tryCatch({
        if (!is.na(utr_start) && !is.na(utr_end)) {
          overlapping_rows <- threeutr %>%
            filter(!((end <= utr_start) | (start >= utr_end)))
          
          count <- nrow(overlapping_rows)
          
          utr_length <- utr_end - utr_start
          if (utr_length > 0) {
            norm_count <- round(count / (utr_length / 1000))
            result_3p$utr_length[j] <- utr_length
            result_3p$count[j] <- count
            result_3p$norm_count[j] <- round(norm_count)
          }
        }
      }, error = function(e) {
        # Handle the error, e.g., print a message or log it
        cat("Error in loop iteration:", j, "for gene:", i, "-", e$message, "\n")
      })
      
    }
    
    hgnc_symbol <- result_id$hgnc_symbol[!is.na(result_id$hgnc_symbol)][1]
    
    if (!is.na(hgnc_symbol)) {
      result_3p$hgnc_symbol <- hgnc_symbol
      df_all <- bind_rows(df_all, result_3p)
    }
  }
}

df_all <- df_all %>% filter(count != 0)
df_all = subset(df_all, select = -count)
#df_all = subset(df_all, select = -entrezgene_id)
#df_all = subset(df_all, select = -entrezgene_description)

#df_all$freq <- round((df_all$count / sum(df_all$count)) * 100, 2)

setwd("/scratch/kovarik/Ronit/snp_project/output_Rds_3UTR")
saveRDS(df_all, file="1kGP_high_coverage_Illumina.chr22.output_PASS_only.vcf.Rds")





#----------------------------------------------------------

non_na_indices <- !is.na(result_id$hgnc_symbol)
if (any(non_na_indices)) {
  hgnc_symbol <- result_id$hgnc_symbol[non_na_indices][1]
  result_3p$hgnc_symbol <- hgnc_symbol
  df_all <- bind_rows(df_all, result_3p)
}

#----------------------------------------------------------

#ylim_val = max(df_all$count) + 10

p <- ggplot(df_all, aes(x = reorder(hgnc_symbol, -count), y = count)) +
  ggtitle("SNPs distribution across genes in X chromosome in 3'UTR") +
  geom_bar(stat = "identity", fill = "steelblue") +
  labs(x = "Gene", y = "Normalised Count") +
  ylim(0, ylim_val) +
  geom_text(aes(label = freq), vjust = -0.5, color = "black", size = 4, angle = 45, hjust = -0.05) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

p


#Normalised Count by mean UTR length----

gene_list = na.omit(unique(loc_all$GENEID))

df_all <- data.frame()

total_snps = nrow(threeutr)

filter_count = 0

for (i in gene_list) {
  
  result_id <- getBM(attributes = c("entrezgene_id", "hgnc_symbol", "entrezgene_description"), 
                     filters = "entrezgene_id", values = i, mart = ensembl)
  
  result_3p <- na.omit(getBM(attributes = c("3_utr_start", "3_utr_end"), 
                             filters = "entrezgene_id", values = i, mart = ensembl))
  
  mean_start = round(mean(result_3p$`3_utr_start`))
  mean_end = round(mean(result_3p$`3_utr_end`))
  
  result_id$`3_utr_start` = mean_start
  result_id$`3_utr_end` = mean_end
  
  count <- threeutr %>%
    filter(!((end <= result_id$`3_utr_start`) | (start >= result_id$`3_utr_end`))) %>%
    nrow()
  
  #freq = (count/filter_snps)*100
  
  utr_length <- result_id$`3_utr_end` - result_id$`3_utr_start`
  norm_count <- round(count / (utr_length / 1000))
  
  filter_count = filter_count + norm_count
}

for (j in gene_list) {
  
  result_id <- getBM(attributes = c("entrezgene_id", "hgnc_symbol", "entrezgene_description"), 
                     filters = "entrezgene_id", values = j, mart = ensembl)
  
  result_3p <- na.omit(getBM(attributes = c("3_utr_start", "3_utr_end"), 
                             filters = "entrezgene_id", values = j, mart = ensembl))
  
  mean_start = round(mean(result_3p$`3_utr_start`))
  mean_end = round(mean(result_3p$`3_utr_end`))
  
  result_id$`3_utr_start` = mean_start
  result_id$`3_utr_end` = mean_end
  
  count <- threeutr %>%
    filter(!((end <= result_id$`3_utr_start`) | (start >= result_id$`3_utr_end`))) %>%
    nrow()
  
  #freq = (count/filter_snps)*100
  
  utr_length <- result_id$`3_utr_end` - result_id$`3_utr_start`
  norm_count <- count / (utr_length / 1000)
  norm_freq <- (norm_count / filter_count) * 100
  
  result_id$utr_length = round(utr_length)
  result_id$count = count
  #result_id$freq = round(freq, 2)
  result_id$norm_count = round(norm_count)
  result_id$norm_freq = round(norm_freq, 2)
  
  df_all = rbind(df_all, result_id)
}

df_all = subset(df_all, count != 0)  
df_all = df_all[order(df_all$norm_count, decreasing = TRUE), ]

p <- ggplot(df_all, aes(x = reorder(hgnc_symbol, -norm_count), y = norm_count)) +
  ggtitle("SNPs distribution across genes in Y chromosome in 3'UTR") +
  geom_bar(stat = "identity", fill = "steelblue") +
  labs(x = "Gene", y = "Normalised Count") +
  ylim(0, 425) +
  geom_text(aes(label = norm_freq), vjust = -0.5, color = "black", size = 4, angle = 45, hjust = -0.05) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggsave("3_utr_gene.pdf", plot = p, width = 10, height = 8, units = "in")
#---------------------------------------------------------------------------------------------------------  





result_3p$entrezgene_id <- result_id$entrezgene_id

result_3p$count = 0

result_3p$freq = 0

total_snps = nrow(threeutr)

for (j in 1:nrow(result_3p)) {
  count <- threeutr %>%
    filter(!((end <= result_3p$`3_utr_start`[j]) | (start >= result_3p$`3_utr_end`[j]))) %>%
    nrow()
  
  freq = (count/total_snps)*100
  
  utr_length <- result_3p$`3_utr_end`[j] - result_3p$`3_utr_start`[j]
  norm_count <- (count / utr_length) * 1000
  norm_freq <- (norm_count / total_snps) * 100
  
  result_3p$count[j] <- count
  result_3p$freq[j] <- round(freq, 2)
  
  result_3p$utr_length[j] <- utr_length
  result_3p$norm_count[j] <- round(norm_count)
  result_3p$norm_freq[j] <- round(norm_freq, 2)
}

result_final = result_3p %>%
  group_by(entrezgene_id) %>%
  summarize(`3_utr_start` = paste(`3_utr_start`, collapse = ", "), `3_utr_end` = paste(`3_utr_end`, collapse = ", "), utr_length = paste(utr_length, collapse = ", "), count = paste(count, collapse = ",  "), freq = paste(freq, collapse = ",  "), norm_count = paste(norm_count, collapse = ", "), norm_freq = paste(norm_freq, collapse = ", "))

final_table = merge(result_id, result_final, by="entrezgene_id", all = TRUE)

df_all <- rbind(df_all, final_table)

#---------------------------------------------------------------------------------------------------------

snps_df <- data.frame(region = names(table(loc_all$LOCATION)), count = table(loc_all$LOCATION))
snps_df$region <- factor(snps_df$region, levels = snps_df$region[order(snps_df$count.Freq, decreasing = TRUE)])

total_snps <- sum(snps_df$count.Freq)
snps_df$freq <- snps_df$count.Freq/total_snps

plot <- ggplot(snps_df, aes(x = reorder(region, -freq), y = freq, fill = region)) +
  geom_bar(stat = "identity") +
  labs(title = "SNP Frequencies by Genomic Region", x = "Genomic Region", y = "SNP Frequency") +
  geom_text(aes(label = count.Freq), vjust = -0.5, color = "black", size = 4, hjust = 0.5) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

plot

#ggsave("snps_plot_all.pdf", plot = plot, width = 10, height = 8, units = "in")

#My Code----

df_all <- data.frame()

for (i in gene_list) {
  
  result_id <- getBM(attributes = c("entrezgene_id", "hgnc_symbol", "entrezgene_description"), 
                     filters = "entrezgene_id", values = i, mart = ensembl)
  
  result_3p <- na.omit(getBM(attributes = c("3_utr_start", "3_utr_end"), 
                             filters = "entrezgene_id", values = i, mart = ensembl))
  
  if (nrow(result_3p) > 0) {
    result_3p$count <- 0
    
    for (j in 1:nrow(result_3p)) {
      count <- threeutr %>%
        filter(!((end <= result_3p$`3_utr_start`[j]) | (start >= result_3p$`3_utr_end`[j]))) %>%
        nrow()
      
      utr_length <- result_3p$`3_utr_end`[j] - result_3p$`3_utr_start`[j]
      if (utr_length > 0) {
        
        norm_count <- round(count / (utr_length / 1000))
        result_3p$utr_length[j] <- utr_length
        result_3p$count[j] <- count
        result_3p$norm_count[j] <- round(norm_count)
      }
      
    }
    result_id$count <- round(mean(result_3p$norm_count))
    df_all <- rbind(df_all, result_id)
  }
  
}

df_all = subset(df_all, count != 0)  
#df_all = df_all[order(df_all$count, decreasing = TRUE), ]

#df_all = df_all[!duplicated(df_all[, "entrezgene_id"]), ]

for (k in 1:nrow(df_all)) {
  freq = (df_all$count[k] / sum(df_all$count)) * 100
  df_all$freq[k] = round(freq, 2)
}
