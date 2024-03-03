library(readxl)
library(dplyr)
library(tidyr)
library(writexl)


setwd("/lisc/user/chakraborty/snp_project/new_bin/")

df_175 = read_xlsx("50_175bp_bins/50_175bp_bins.xlsx")
df_mid = read_xlsx("50_middle_bins/50_middle_bins.xlsx")

setwd("/lisc/user/chakraborty/snp_project/3UTR/")
hk_genes = read.csv("Database/housekeeping_genes.txt", header = F)
colnames(hk_genes) = "hgnc_symbol"

df_175_hk = merge(hk_genes, df_175, by = "hgnc_symbol")
df_mid_hk = merge(hk_genes, df_mid, by = "hgnc_symbol")

setwd("/lisc/user/chakraborty/snp_project/new_bin/")
write.xlsx(df_mid_hk, "Database/mid_hk.xlsx")

setwd("/lisc/user/chakraborty/snp_project/3UTR/")
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

df_175_cyt = merge(cytokine_new, df_175, by = "hgnc_symbol")
df_mid_cyt = merge(cytokine_new, df_mid, by = "hgnc_symbol")

setwd("/lisc/user/chakraborty/snp_project/new_bin/")
write.xlsx(df_mid_cyt, "Database/mid_cyt.xlsx")


setwd("/lisc/user/chakraborty/snp_project/3UTR/")
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

df_175_che = merge(chemokine_new, df_175, by = "hgnc_symbol")
df_mid_che = merge(chemokine_new, df_mid, by = "hgnc_symbol")

setwd("/lisc/user/chakraborty/snp_project/new_bin/")
write.xlsx(df_mid_che, "Database/mid_che.xlsx")
