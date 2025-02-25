library(pheatmap)
library(ggplot2)
library(RColorBrewer)

#set the output name
out <- "bp_t2d_imputeSCOPAimp"

#load your harmonised TABLE contening all the z-scores (the subsequent loaded table is named z_table2)
load("./z_table_harmonised.RData")
cat("there are ",sum(is.na(z_table2))," missing values")

 #QC1: checking column missingess
#number of missingness by GWAS - only step2
who_is_na <- colSums(is.na(z_table2[,9:ncol(z_table2)]))
print(who_is_na)
#remove FG (redundancy with RG), tot_chol (redundancy with LDL)
df <- z_table2[,-c("z_FG","z_total_cholesterol", "z_FGluadjBMI")]
df <- df[!duplicated(rsid)]

 #QC2: checking the same for rows 
#number of missingness by variant,
who_is_na2 <- rowSums(is.na(df[,9:ncol(df)]))
print(who_is_na2)
#only keep 20% or lower rate of missingness 
df <- df[!who_is_na2>0.2*(ncol(df)-8),]

 #QC3: for each variants, there is overlapping by including the proxy (in the previous table harmonisation step)
#this for loop takes each snp and check if the proxy is doing better (less NA)
for (snp in unique(df$rsid)) {
  rows <- which(df$rsid == snp | df$proxy == snp)
  
  if (length(rows) <= 1) next
  
  df_subset <- df[rows, ]
  who_is_na_temp <- rowSums(is.na(df_subset[,9:ncol(df_subset)]))
  who_is_min <- which.min(who_is_na_temp)
  rows_to_remove <- rows[-who_is_min]
  
  df <- df[-rows_to_remove, ]
}

#table with no annotations
df_no_labels <- cbind(df[,2], df[,6], df[,9:ncol(df)])
#table with numeric z-scores only
df_no_labels_no_na <- df_no_labels[,3:ncol(df_no_labels)]


#use imputeSCOPA to run multiple standard imputation
#this first line load imputeSCOPA, it requires:
# - ranger
# - Rcpp
# - tictoc 
# - you should have installed imputed scopa in your ~/ folder 
# see https://github.com/ImperialStatGen/imputeSCOPA#readme
source("~/imputeSCOPA/R/csImputeSCOPA.R") 
imp <- csImputeSCOPA(df_no_labels_no_na, maxiter = 1, num.trees = 10, verbose = T, seed = 7)
df_no_labels_no_na <- as.data.frame(imp)


#truncate first all column to max 5.45 (~= 5.45 sd) + mean imputation
for (j in names(df_no_labels_no_na)) {
  if (Inf %in% df_no_labels_no_na[,j]) {
    ind <- which(df_no_labels_no_na[, j] == Inf)
    df_no_labels_no_na[ind, j] <- max(df_no_labels_no_na[- ind, j])
  }
  
  # sd2 <- 5.45 #used to plot, but not to cluster
  sd2 <- 2*sd(as.matrix(df_no_labels_no_na[,j]), na.rm = T) #used to cluster
  for (i in 1:nrow(df_no_labels_no_na)) {
    if (df_no_labels_no_na[i,j] > sd2) {
      df_no_labels_no_na[i,j] <- sd2
    } else if (df_no_labels_no_na[i,j] < -sd2) {
      df_no_labels_no_na[i,j] <- -sd2
    }
  }
}

#reput in df_no_labels the imputation
df_no_labels[,5:ncol(df_no_labels)] <- df_no_labels_no_na


#choose the best number of cluster, k here
k = 5


 #clustering using pheatmap
cluster <- pheatmap(df_no_labels_no_na, file = paste0("./clustering_result/imputeSCOPA/",out,"_hierarchical_clustering_",k,"clusters_pheatmap.png"), 
         clustering_method = "ward.D2", 
         # scale = "column",
         breaks = seq(-5.45,5.45,length.out=8), color = colorRampPalette(rev(brewer.pal(n = 7,
         name = "RdBu")))(7), show_rownames = F, width = 12, height = 12, 
         cutree_rows = k)
		 
#obtain the trees in a table, need to put the number of cluster 
#determine the desired distance threshold
threshold <- max(cluster$tree_row$height) / 2 #typically half of the distance here
df_no_labels$cluster <- cutree(cluster$tree_row, h = threshold)
k <- length(unique(df_no_labels$cluster)) #iterative: you can redo line 87 based on this new k to generate another clustering
write.table(df_no_labels, paste0("./clustering_result/imputeSCOPA/",out,"_hierarchical_clustering_",k,"clusters_table.csv"))

#boxplot of phenos
for (cluster in unique(df_no_labels$cluster)) {
  
  temp_df <- df_no_labels_no_na[which(df_no_labels$cluster == cluster),]
  
  reformatted_df_no_labels_no_na <- data.frame()
  
  for (col in names(temp_df)) {
    
    temp_reformatted_df <- data.frame(z = temp_df[,col], 
                                      pheno = rep(col, nrow(temp_df)), 
                                      mean_value = rep(mean(temp_df[,col]),nrow(temp_df)))
    reformatted_df_no_labels_no_na <- rbind(reformatted_df_no_labels_no_na, temp_reformatted_df)
    
  }
  
  ggplot(reformatted_df_no_labels_no_na, aes(x = reorder(pheno, mean_value, mean), y = z, fill = mean_value)) +
    geom_boxplot() + 
    scale_fill_gradient2(low = "#2166AC", mid = "#F7F7F7", high = "#B2182B", midpoint = 0) +
    theme_minimal() + 
    labs(title = paste("Boxplot of cluster",cluster),y = "Z-score",x = "Phenotype") +
    theme(axis.text.x = element_text(angle = 90))
  ggsave(paste0("./clustering_result/imputeSCOPA/",out,"_hierarchical_clustering_", k, "cluster_",cluster,"_boxplot.png"), width = 10, height = 10)
  
}
