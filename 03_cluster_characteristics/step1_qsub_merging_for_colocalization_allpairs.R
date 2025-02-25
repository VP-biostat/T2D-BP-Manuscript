library(data.table)
library(arrow)

# Look at the input 
args = commandArgs(trailingOnly=TRUE)
if (length(args)==0) {
  warning("\nRequired arguments: 
       \n1. a numeric from 1 to 49\n", call.=FALSE)
} else {
  neq <- as.numeric(args[1])
}

parquet_folder <- "./GTEx_v8_allpair/GTEx_Analysis_v8_EUR_eQTL_all_associations/"

## Attempt of a colocalization analysis using merged datasets of both GWAS and eQTL, cis eqtl from GTEx v8 (n = 49 tissues) and SNPs from latest clustering

# List of tissues of interest: 
int_tissues <- unique(gsub(".v8.*", "", list.files(path = parquet_folder)))
#select the one of interest (1 per job submitted)
int_tissue <- int_tissues[neq]

#READING THE PARQUETS FILES 
cat("Reading eQTL all_pairs from GTEx (", int_tissue,")...\n")
parquet_files <- list.files(path = parquet_folder, pattern = int_tissue)

parquet_data_frames <- list()

for (file in parquet_files) {
  #excluding chr x 
  if (grepl("chrX", file)) {
    next
  }
  cat("\tReading eQTL all_pairs from GTEx -", file, "\n")
  parquet_data <- read_parquet(paste0(parquet_folder, file))
  parquet_data_frames[[file]] <- parquet_data
}

subeqtl <- rbindlist(parquet_data_frames)
rm(parquet_data)
rm(parquet_data_frames)
cat("\tAdding the necessary columns for eQTL file\n")
#add the chr # in a column , pos_38 and chrpos_38
subeqtl[, chr := as.numeric(gsub("chr(.+)_.*_.*_.*_b38", "\\1", variant_id))]
subeqtl[, Pos_b38 := as.numeric(gsub("chr.*_(.+)_.*_.*_b38", "\\1", variant_id))]
subeqtl[, chrpos_b38 := paste(chr, Pos_b38, sep = ":")]

# Reading the GWASs (which should be lift-overed to build 38!
#T2D
cat("reading T2D GWAS...\n")
t2d <- fread("./Mahajan.NatGenet2018b.T2D.European_with_b38.txt")
#remove the ones without liftover mapping
t2d <- t2d[liftOver_unMapped == F, ]
#add the chrpos_b38 column
t2d[, chrpos_b38 := paste(Chr, Pos_b38, sep = ":")]
#BP 
cat("reading BP GWASs...\n")
sbp <- fread("./UKB_SBP_summary_with_b38.txt")
sbp <- sbp[liftOver_unMapped == F, ]
sbp[, chrpos_b38 := paste(chromosome, Pos_b38, sep = ":")]
dbp <- fread("./UKB_DBP_summary_with_b38.txt")
dbp <- dbp[liftOver_unMapped == F, ]
dbp[, chrpos_b38 := paste(chromosome, Pos_b38, sep = ":")]
pp <- fread("./UKB_PP_summary_with_b38.txt")
pp <- pp[liftOver_unMapped == F, ]
pp[, chrpos_b38 := paste(chromosome, Pos_b38, sep = ":")]

# # Reading the fourier map created by lddetect 
# cat("reading Fourier map for EUR, created by lddetect to select the region of interest for each top SNP...\n")
# map <- fread("./colocalization/fourier_ls-all.bed")

#extract what we will use as df1: the GWASes from T2D and BP based on SNPs from the clustering
cat("reading Clustering...\n")
cluster <- fread("./clustering_result/imputeSCOPA/YOURCLUSTERINGRESULT_5clusters_table.txt") #replace by your output clustering name
subcluster <- cluster[, c("rsid", "GRCh37_chrpos", "outcome", "gcluster", "z_t2d", "z_SBP", "z_DBP", "z_PP")]
cat("\tAdding the chrpos_b38 column in clustering\n")
subcluster <- merge(subcluster, sbp[, c("rsid", "chrpos_b38")], by = "rsid")
ncl <- nrow(subcluster)

#assigning for each SNPs of the cluster a phenotype ID matching the name of locus in the eqtl
#  we assign it if the SNP is within 1 Mb of tss_b37
#  location adapted from the coloc FAQ: https://github.com/chr1swallace/coloc/blob/main/FAQ.md#if-i-understand-correctly-colocabf-can-be-run-with-correlated-variants-that-is-no-prerequisite-for-taking-through-ld-pruningclumping-is-required-am-i-correct-in-my-understanding-
for (j in 1:ncl) {
  snp <- as.character(subcluster[j, "chrpos_b38"])
  snp_chr <- as.numeric(gsub(":.*", "", snp))
  snp_pos <- as.numeric(gsub(".*:", "", snp))
  is_in_1mb <- which(abs(subeqtl$Pos_b38-snp_pos) <= 100000 & subeqtl$chr == snp_chr) #EDITED USED A 200KB WINDOW  (see above)
  # #identify good map row to take in our lddetect map, this will be the region of extraction for this top snp 
  # map_row <- map[chr == paste0("chr", snp_chr) & ((snp_pos-start) * (snp_pos-stop)<0)]
  # map_row <- map[j, ]
  # snp <- paste(map_row, collapse = "_")
  # snp_chr <- as.integer(gsub("chr", "", map_row$chr))
  # snp_pos <- as.numeric(map_row$start)
  # if (nrow(map_row) > 1) cat("Error in the mapping: 2 rows selected in map_row!\n")
  # is_in_1mb <- which((subeqtl$var_pos_b37 > map_row$start) & (subeqtl$var_pos_b37 < map_row$stop) & (subeqtl$chr == snp_chr))
  
  found <- length(is_in_1mb)
  if (found >= 1) {
    cat("\nfor", snp, "we found", found, "snps in eqtl, looking for them in GWAS")
    subeqtl2 <- subeqtl[is_in_1mb, ]
    
    temp_merge_eqtl_t2d <- merge(subeqtl2, t2d, by = "chrpos_b38")
    temp_merge_eqtl_t2d$cluster_snp <- snp
    cat("\n  we have", nrow(temp_merge_eqtl_t2d), "snps found both in eQTL and T2D")
    temp_merge_eqtl_sbp <- merge(subeqtl2, sbp, by = "chrpos_b38")
    temp_merge_eqtl_sbp$cluster_snp <- snp
    cat("\n  we have", nrow(temp_merge_eqtl_sbp), "snps found both in eQTL and SBP")
    temp_merge_eqtl_dbp <- merge(subeqtl2, dbp, by = "chrpos_b38")
    temp_merge_eqtl_dbp$cluster_snp <- snp
    cat("\n  we have", nrow(temp_merge_eqtl_dbp), "snps found both in eQTL and DBP")
    temp_merge_eqtl_pp <- merge(subeqtl2, pp, by = "chrpos_b38")
    temp_merge_eqtl_pp$cluster_snp <- snp
    cat("\n  we have", nrow(temp_merge_eqtl_pp), "snps found both in eQTL and PP")
    
    # #if we are the first try, we have to create the merge files, 
    # #if the 1st try was not conclusive (no association found), then we need to 
    # #create the merge files on the second or third attempts
    # if (j == 1 | !exists("merge_eqtl_t2d")) {
    #   merge_eqtl_t2d <- temp_merge_eqtl_t2d
    #   merge_eqtl_sbp <- temp_merge_eqtl_sbp
    #   merge_eqtl_dbp <- temp_merge_eqtl_dbp
    #   merge_eqtl_pp <- temp_merge_eqtl_pp
    # } else {
    #   merge_eqtl_t2d <- rbind(merge_eqtl_t2d, temp_merge_eqtl_t2d)
    #   merge_eqtl_sbp <- rbind(merge_eqtl_sbp, temp_merge_eqtl_sbp)
    #   merge_eqtl_dbp <- rbind(merge_eqtl_dbp, temp_merge_eqtl_dbp)
    #   merge_eqtl_pp <- rbind(merge_eqtl_pp, temp_merge_eqtl_pp)
    # }
    #EDIT: too many things to be merged, saving the temp merged gwas and eqtl 
    out_file <- paste0("./colocalization/merged_gwas_eqtl/t2d_bp_gwas_", int_tissue,"_eqtl_chr", snp_chr, "_", snp_pos, ".Rdata")
    cat("\n\nSaving the merged file for", snp, "in an Rdata with all t2d and bp gwas and", int_tissue, "in: \n\t", out_file,"\n")
    save(temp_merge_eqtl_t2d, temp_merge_eqtl_sbp, temp_merge_eqtl_dbp, temp_merge_eqtl_pp, file = out_file)
    
  }
  else {
    cat("\nno association found for", snp)
  }
}
cat("\n\nEnd of the process!")