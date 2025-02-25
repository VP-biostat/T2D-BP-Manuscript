library(arrow)
library(dplyr)

tissue <- "tissue id you want to extract"
parquet_folder <- "./GTEx_v8_allpair/GTEx_Analysis_v8_EUR_eQTL_all_associations/"
parquet_files <- list.files(path = parquet_folder, pattern = tissue)

parquet_data_frames <- list()

for (file in parquet_files) {
  parquet_data <- read_parquet(paste0(parquet_folder, file))
  parquet_data_frames[[file]] <- parquet_data
}

combined_data <- bind_rows(parquet_data_frames)