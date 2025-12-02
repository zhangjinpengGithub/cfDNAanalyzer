#!/usr/bin/env Rscript
args <- commandArgs(TRUE)
dir_path = args[1]

output_dir <- paste0(dir_path, "/filter_sample")
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}

PFE = read.csv(paste0(dir_path, "/PFE.csv"))
PFE_sample = PFE$sample

file_names <- list.files(dir_path, pattern = "\\.csv$", full.names = FALSE)

file_names <- file_names[file_names != "PFE.csv"]

for (i in file_names) {
  feature_csv = read.csv(paste0(dir_path, "/", i))
  filtered_data <- feature_csv[feature_csv$sample %in% PFE_sample, ]
  write.csv(filtered_data, file = paste0(output_dir, "/", i), quote = F, row.names = F)
  cat("Processed:", i, "\n")
}