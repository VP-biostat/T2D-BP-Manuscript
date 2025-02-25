library(data.table)
library(dplyr)
library(readxl)

out <- "YOUR_OUTPUT_NAME"

#READ DIRECTLY IF THE PREVIOUS COMMANDS HAVE BEEN DONE PREVIOUSLY 
signif_results <- fread(paste0('./colocalization/', out, 'signif_results_curated.csv'))

#READ the clustering - just to extract b38
subcluster <- fread("./clustering_result/imputeSCOPA/", out, "_5clusters_table.txt")

sorted_signif_results <- signif_results[PP.H4.abf > 0.8 & nsnps < 10000 & nsnps > 100 & PP.H3.abf < 0.5, ]
sorted_signif_results_with_gcluster <- merge(sorted_signif_results, subcluster, by.x = "top_snp", by.y = "mapping_chrpos", all.x = T)

prev_loci_scott <- read_excel("./db161253supplementarydata2.xlsx", sheet = 9, skip = 2)

prev_loci_suzuki <- read_excel("./41586_2024_7019_MOESM3_ESM.xlsx", sheet = 8, skip = 3)

for (g in unique(subcluster$gcluster)) {
  
  if (!dir.exists(paste0("./colocalization/g", g, "/"))) dir.create(paste0("./colocalization/g", g, "/"))
  
  #take only gcluster == g
  subdf <- sorted_signif_results_with_gcluster[gcluster == g, ]
  
  #removed the top_snp from clustering that were in the same 200kb region (taking 100kbp margin)
  subdf <- subdf[diff_pb == F, ]
  #removing the snps within the MHC region: from Leyden et al. https://www.cell.com/ajhg/fulltext/S0002-9297(21)00468-7
  # Variants within the MHC region (chr6: 25,000,000–35,000,000) were excluded from analyses
  subdf <- subdf[!(chr == 6 & as.numeric(gsub(".*:", "", GRCh37_chrpos)) > 25000000 & as.numeric(gsub(".*:", "", GRCh37_chrpos)) < 35000000), ]
  
  fwrite(subdf, file = paste0("./colocalization/g", g, "/", out, "_H4_0.8_H3_0.5_noMHC_indepSNV.tsv"),
         sep = "\t", quote = F)
  
  merge_wprevloci_scott <- merge(subdf, prev_loci_scott, by.x = "SYMBOL", by.y = "locus")
  fwrite(merge_wprevloci_scott, file = paste0("./colocalization/g", g, "/", out, "_merged_wprevloci_scott_et_al.tsv"),
          sep = "\t", quote = F)
  
  merge_wprevloci_suzuki <- merge(subdf, prev_loci_suzuki, by.x = "SYMBOL", by.y = "...1")
  fwrite(merge_wprevloci_suzuki, file = paste0("./colocalization/g", g, "/", out, "_merged_wprevloci_suzuki_et_al.tsv"),
          sep = "\t", quote = F)
  
}
