library(data.table)

#read the soft clustering
bnmf_path = "insert your path of the bNMF clustering"
bnmf_W = fread(bnmf_path)

#read the hard clustering
hieclu_path = "insert your path of the hierarchical clustering to compare"
hieclu = fread(hieclu_path)

#merge the clusteing 
comp = merge(bnmf_W, hieclu[, c("rsid", "gcluster")], by.x = "variant", by.y = "rsid")

#out name
out = "bp_t2d_imputeSCOPAimp"

### LOGISTIC REGRESSION ----------
library(logistf)

#test the stat logist assoc with one cluster from hieclu to the weights obtained
#from bNMF results by using logit link function 
#then compare the deviances of the two models by an anova test of enrichment (chisq test)

cls <- na.omit(unique(hieclu$gcluster))

get_sign <- function(x) {
  if (x > 0) {
    return(1)
  } else if (x < 0) {
    return(-1)
  } else {
    return(0)
  }
}

logit_list <- list()
if (exists("logit_df")) rm(logit_df)

for (cl in cls) {
  
  cat("Cluster", cl, "\n")
  
  subcomp <- comp[, current_cl := gcluster == cl]
  
  formula_constrained <- as.formula(paste("current_cl ~", paste(grep(pattern = "^X", names(subcomp), value = T), collapse = " + "), "- 1"))
    
  lm <- glm(formula_constrained, data = subcomp)
  
  print(summary(lm))
    
  logit_list[[cl]] <- summary(lm)$coefficients
  
  fwrite(logit_list[[cl]], file = paste0("./clustering_result/imputeSCOPA/comparison_bNMF_hieclu/", out, "_cl", cl, "_hieclu_bnmf_comparison_glm_table.csv"), sep = ";", quote = F)
  
  #write a table to do a heatmap
  temp_df <- -log(logit_list[[cl]][, 4])*unlist(lapply(logit_list[[cl]][, 1], get_sign))
  temp_bnmf <- names(temp_df)
  temp <- data.frame(hieclu = rep(cl, length(temp_df)), softclu = temp_bnmf, enrich = unlist(temp_df))
  if (!exists("logit_df")) {
    logit_df <- temp
  } else {
    logit_df <- rbind(logit_df, temp)
  }
  
}

p = ggplot(logit_df, aes(y = softclu, x = as.character(hieclu), fill = enrich)) +
  geom_tile(color = "white", linewidth = 1) +
  scale_fill_gradient2(low = "#2166AC", mid = "white", high = "#B2182B") +  # You can adjust the color gradient as needed
  theme_minimal() +
  xlab("T2D-BP hierarchical clusters (k = 5)") + ylab("T2D-BP bNMF clusters (k = 8)") +
  scale_x_discrete(labels = c("Inverse T2D-BP risk", "Metabolic Syndrome", "Higher adiposity", "Vascular dysfunction", "Reduced beta-cell function")) +  # Custom x-axis labels
  scale_y_discrete(labels = c("Cluster 1", "Cluster 2", "Cluster 3", "Cluster 4", "Cluster 5", "Cluster 6", "Cluster 7", "Cluster 8")) +  # Custom y-axis labels
  theme(text = element_text(size = 20),
        axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1),
        axis.title.x = element_text(size = 24),
        axis.title.y = element_text(size = 24), 
        legend.title = element_blank(), 
        legend.key.height = unit(2, "cm")) 
ggsave(p, filename = paste0("./comparison_bNMF_hieclu_", out, "_heatmap_of_comparison.png"), width = 8, height = 8, dpi = 600)

