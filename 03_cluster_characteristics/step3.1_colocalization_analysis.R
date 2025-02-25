library(data.table)
library(pheatmap)
library(ggplot2)
library(dplyr)

out <- "YOUR_OUTPUT_NAME"

coloc_files <- list.files(paste0('./colocalization/', out), pattern = "result_table.csv", full.names = T)

# read all output from step2 and concatenate them all
rm(signif_results)
for (file in coloc_files) {
  coloc <- fread(file)
  if (!exists("signif_results")) {
    signif_results <- coloc
  } else {
    signif_results <- rbind(signif_results, coloc)
  }
}

sorted_signif_results <- signif_results[order(-PP.H4.abf), ]

sorted_signif_results[, position := as.numeric(gsub('chr.*_', '', top_snp))]
sorted_signif_results[, chr := factor(gsub('chr(.+)_.*', '\\1', top_snp))]

#removing the snps within the MHC region: from Leyden et al. https://www.cell.com/ajhg/fulltext/S0002-9297(21)00468-7
# Variants within the MHC region (chr6: 25,000,000–35,000,000) were excluded from analyses
sorted_signif_results <- sorted_signif_results[!(chr == 6 & as.numeric(position) > 25000000 & as.numeric(position) < 35000000), ]

#Filter : PPH4 > 0.8 ; nsnps < 10000 > 100 
sorted_signif_results <- sorted_signif_results[PP.H4.abf > 0.8 & nsnps < 10000 & nsnps > 100 & PP.H3.abf < 0.5, ]

#Make a count table
count_table_signif <- table(sorted_signif_results$gwas, sorted_signif_results$tissue)

#plot overall of the number of colocalization per eQTL/GWAS pairs
pheatmap::pheatmap(t(count_table_signif), 
                   color = colorRampPalette(c("white", "#b2182b"))(6), 
                   labels_col = c("DBP", "PP", "SBP", "T2D"), 
                   cluster_cols = F, cluster_rows = F, 
                   display_numbers = T, number_format = "%1.0f", 
                   filename = paste0("./colocalization/", out, "/heatmap_coloc_tissues_wTIGER.png"))

cluster <- fread("./clustering_result/imputeSCOPA/bp_t2d_TRUEnold0.2_wlocus_new_phenos2_hierarchical_clustering_5clusters_table_for_colocalization.txt")
#Adding the chrpos_b38 column in clustering
sbp <- fread("./gwas/warren_UKB_BP/UKB_SBP_summary_with_b38.txt")
sbp <- sbp[liftOver_unMapped == F, ]
sbp[, chrpos_b38 := paste(chromosome, Pos_b38, sep = ":")]
subcluster <- merge(cluster, sbp[, c("rsid", "chrpos_b38")], by = "rsid")
subcluster[, mapping_chrpos1 := paste0('chr',  gsub(":", "_", GRCh37_chrpos))] #reformat top b37 snp to match the eqtl sorted_signif format (for TIGER)
subcluster[, mapping_chrpos2 := paste0('chr',  gsub(":", "_", chrpos_b38.x))] #reformat top b38 snp to match the eqtl sorted_signif format (for the others)

for (cl in unique(subcluster$gcluster)) {
  sub_subcluster <- subcluster[gcluster == cl, ]
  sub_sorted_signig <- sorted_signif_results[top_snp %in% sub_subcluster$mapping_chrpos1 | 
                                               top_snp %in% sub_subcluster$mapping_chrpos2, ] #selecting only the top snp related to a cluster 
  
  # Calculate count for each tissue
  count_per_tissue <- sub_sorted_signig %>%
    count(tissue) %>%
    arrange(desc(n)) %>%
    pull(tissue)
  
  # Reorder the levels of tissue based on the count
  sub_sorted_signig$tissue <- factor(sub_sorted_signig$tissue, levels = count_per_tissue)
  
  #try the barplot
  ggplot(sub_sorted_signig, aes(x = tissue, fill = gwas)) +
    geom_bar(stat = "count") +
    labs(title = "Counts in Different Tissues", x = "Tissue", y = "Count") +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
  ggsave(filename = paste0("./colocalization/", out, "/barplot_coloc_tissues_cluster_", cl, ".png"), width = 10, height = 6)
  
  #Make a count table
  sub_count_table_signif <- table(sub_sorted_signig$gwas, sub_sorted_signig$tissue)
  
  #plot of the number of colocalization per eQTL/GWAS pairs among 1 cluster
  pheatmap::pheatmap(t(sub_count_table_signif), 
                     color = colorRampPalette(c("white", "#b2182b"))(6), 
                     labels_col = c("DBP", "PP", "SBP", "T2D"), 
                     cluster_cols = F, cluster_rows = F, 
                     display_numbers = T, number_format = "%1.0f", 
                     filename = paste0("./colocalization/", out, "/heatmap_coloc_tissues_cluster_", cl, ".png"))
}


##SAME for limited number of tissues
int_tissues <- c("Adipose_Subcutaneous", "Adipose_Visceral_Omentum", "Adrenal_Gland", "Artery_Aorta", "Heart_Left_Ventricle", "Liver", "Muscle_Skeletal", "Pancreas", "Whole_Blood", "Thyroid", "Artery_Tibial")

sorted_signif_results2 <- sorted_signif_results[tissue %in% int_tissues,]

#plot overall of the number of colocalization per eQTL/GWAS pairs
pheatmap::pheatmap(t(table(sorted_signif_results2$gwas, sorted_signif_results2$tissue)), 
                   color = colorRampPalette(c("white", "#b2182b"))(6), 
                   labels_col = c("DBP", "PP", "SBP", "T2D"), 
                   cluster_cols = F, cluster_rows = F, 
                   display_numbers = T, number_format = "%1.0f", 
                   filename = paste0("./colocalization/", out, "/heatmap_coloc_reduced_tissues.png"))

for (cl in unique(subcluster$gcluster)) {
  sub_subcluster <- subcluster[gcluster == cl, ]
  sub_sorted_signig <- sorted_signif_results2[top_snp %in% sub_subcluster$mapping_chrpos, ] #selecting only the top snp related to a cluster 
  
  #plot of the number of colocalization per eQTL/GWAS pairs among 1 cluster
  pheatmap::pheatmap(t(table(sub_sorted_signig$gwas, sub_sorted_signig$tissue)), 
                     color = colorRampPalette(c("white", "#b2182b"))(6), 
                     labels_col = c("DBP", "PP", "SBP", "T2D"), 
                     cluster_cols = F, cluster_rows = F, 
                     display_numbers = T, number_format = "%1.0f", 
                     filename = paste0("./colocalization/", out, "/heatmap_coloc_reduced_tissues_cluster_", cl, ".png"))
  
}

fwrite(sorted_signif_results, file = paste0("./colocalization/", out, "/signif_results_curated.tsv", sep = "\t", quote = F))
