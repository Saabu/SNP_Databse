library(openxlsx, lib.loc = "/lisc/app/R/4.3.1/lib64/R/library")
library(tidyr)

setwd("/lisc/user/chakraborty/snp_project/new_bin/50_middle_bins/output_rds/")

rds_files <- list.files(pattern = "\\.Rds$")

data_frames <- list()

for (file in rds_files) {
  data_frames[[file]] <- readRDS(file)
}

combined_df <- do.call(rbind, data_frames)
rownames(combined_df) = NULL

df_long <- pivot_longer(combined_df, cols = c(count_50bp, count_middle_50bp), names_to = "count_type", values_to = "count")

write.xlsx(combined_df, "/lisc/user/chakraborty/snp_project/new_bin/50_middle_bins")
