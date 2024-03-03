library(readxl)
library(dplyr)
library(tidyr)
library(writexl)


setwd("/lisc/user/chakraborty/snp_project/3UTR/")

three_utr = read_xlsx("Database/final_df_3UTR.xlsx")
five_utr = read_xlsx("/lisc/user/chakraborty/snp_project/5UTR/final_df_5UTR.xlsx")
cds = read_xlsx("/lisc/user/chakraborty/snp_project/CDS/final_df_CDS.xlsx")
category_1 = read_xlsx("Database/Category_1.xlsx")
category_2 = read_xlsx("Database/Category_2.xlsx")
category_3 = read_xlsx("Database/Category_3.xlsx")

hk_genes = read.csv("Database/housekeeping_genes.txt", header = F)
colnames(hk_genes) = "hgnc_symbol"

cat_1_hk = merge(hk_genes, category_1, by = "hgnc_symbol")
cat_2_hk = merge(hk_genes, category_2, by = "hgnc_symbol")
cat_3_hk = merge(hk_genes, category_3, by = "hgnc_symbol")
three_utr_hk = merge(hk_genes, three_utr, by = "hgnc_symbol")
five_utr_hk = merge(hk_genes, five_utr, by = "hgnc_symbol")
cds_hk = merge(hk_genes, cds, by = "hgnc_symbol")

write.xlsx(cds_hk, "Housekeeping/cds_hk.xlsx")


cytokine_genes = read_xlsx("Database/Cytokine.xlsx")
cytokine_genes = drop_na(cytokine_genes)
cytokine_genes = cytokine_genes[, 1:2]

cytokine_genes <- cytokine_genes %>%
  mutate(Name = gsub("-", "", Name),  # Remove hyphens from Name column
         `Synonym(s)` = gsub("-", "", `Synonym(s)`)) %>%  # Remove hyphens from Synonym(s) column
  separate_rows(`Synonym(s)`, sep = ", ")  # Separate values in Synonym(s) column based on commas

cytokine_genes <- cytokine_genes %>%
  mutate(Name = gsub("α", "A", Name),
         Name = gsub("β", "B", Name),
         `Synonym(s)` = gsub("α", "A", `Synonym(s)`),
         `Synonym(s)` = gsub("β", "B", `Synonym(s)`))

name = as.data.frame(cytokine_genes$Name)
symb = as.data.frame(cytokine_genes$`Synonym(s)`)
colnames(name) = "name"
colnames(symb) = "name"

cytokine_new = rbind(name, symb)
colnames(cytokine_new) = "hgnc_symbol"
cytokine_new$hgnc_symbol = toupper(cytokine_new$hgnc_symbol)

cat_1_cyt = merge(cytokine_new, category_1, by = "hgnc_symbol")
cat_2_cyt = merge(cytokine_new, category_2, by = "hgnc_symbol")
cat_3_cyt = merge(cytokine_new, category_3, by = "hgnc_symbol")
three_utr_cyt = merge(cytokine_new, three_utr, by = "hgnc_symbol")
five_utr_cyt = merge(cytokine_new, five_utr, by = "hgnc_symbol")
cds_cyt = merge(cytokine_new, cds, by = "hgnc_symbol")

write.xlsx(cds_cyt, "Cytokine/cds_cyt.xlsx")



chemokine_genes = read_xlsx("Database/Chemokine.xlsx")
chemokine_genes = drop_na(chemokine_genes)
chemokine_genes = chemokine_genes[, 1:2]

chemokine_genes <- chemokine_genes %>%
  mutate(Name = gsub("-", "", Name),  # Remove hyphens from Name column
         `Synonym(s)` = gsub("-", "", `Synonym(s)`)) %>%  # Remove hyphens from Synonym(s) column
  separate_rows(`Synonym(s)`, sep = ", ")  # Separate values in Synonym(s) column based on commas

chemokine_genes <- chemokine_genes %>%
  mutate(Name = gsub("α", "A", Name),
         Name = gsub("β", "B", Name),
         `Synonym(s)` = gsub("α", "A", `Synonym(s)`),
         `Synonym(s)` = gsub("β", "B", `Synonym(s)`),
         `Synonym(s)` = gsub(" ", "", `Synonym(s)`))

chemokine_genes$`Synonym(s)` <- sub("\\+$", "", chemokine_genes$`Synonym(s)`)

name = as.data.frame(chemokine_genes$Name)
symb = as.data.frame(chemokine_genes$`Synonym(s)`)
colnames(name) = "name"
colnames(symb) = "name"

chemokine_new = rbind(name, symb)
colnames(chemokine_new) = "hgnc_symbol"
chemokine_new$hgnc_symbol = toupper(chemokine_new$hgnc_symbol)

cat_1_che = merge(chemokine_new, category_1, by = "hgnc_symbol")
cat_2_che = merge(chemokine_new, category_2, by = "hgnc_symbol")
cat_3_che = merge(chemokine_new, category_3, by = "hgnc_symbol")
three_utr_che = merge(chemokine_new, three_utr, by = "hgnc_symbol")
five_utr_che = merge(chemokine_new, five_utr, by = "hgnc_symbol")
cds_che = merge(chemokine_new, cds, by = "hgnc_symbol")

write.xlsx(cds_che, "Chemokine/cds_che.xlsx")
