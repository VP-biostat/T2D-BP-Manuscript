library(data.table)
library(ggplot2)
library(dplyr)

out <- "bp_t2d_imputeSCOPAimp"

cluster = fread(paste0("./clustering_result/imputeSCOPA/", out, "_hierarchical_clustering_5clusters_table_for_colocalization.txt"))
z_col = grep("z_", names(cluster), value = T)

if (exists("gg_heat")) rm(gg_heat)
for (z in z_col) {
  
  n_cl= length(unique(cluster$gcluster))
  
  formula = as.formula(paste0(z, " ~ c1+c2+c3+c4+c5-1"))
  
  lm = summary(glm(data = cluster, formula))$coefficients
  
  temp = data.table(Z = rep(z, n_cl), 
                    Cluster = c(1,2,3,4,5), 
                    Beta = lm[,1], 
                    SE = lm[,2],
                    P = lm[,4])
  
  if (exists("gg_heat")) {
    gg_heat = rbind(gg_heat, temp)
  } else {
    gg_heat = temp
  }
    
}

# Plot heatmap
desired_order <- c("z_heart_rate", 
                   "z_IGF_1", 
                   "z_small_vessel_stroke", 
                   "z_PAI_1", 
                   "z_ostradiol_m", 
                   "z_oestradiol_f",
                   "z_serum_level_ACE_protein", 
                   "z_TNF_related_apoptosis_inducing_ligand_receptor_2", 
                   "z_leptin", 
                   "z_renin", 
                   "z_proinsulin", 
                   "z_menopause",
                   "z_IL10",
                   "z_IL6", 
                   "z_TNF_alpha",
                   "z_testosterone_f",
                   "z_HOMA_B", 
                   "z_HOMA_IR", 
                   "z_fasting_insulin",
                   "z_2hrGluadjBMI", 
                   "z_modified_stumvoll_ISI",
                   "z_HbA1C",
                   "z_TG",
                   "z_RG", 
                   "z_heart_failure",
                   "z_CAD", 
                   "z_stroke", 
                   "z_atrial_fibrillation", 
                   "z_aspartate_aminotransferase", 
                   "z_alanine_aminotransferase",
                   "z_WHRadjBMI",
                   "z_WBC", 
                   "z_heel_bmd", 
                   "z_CRP", 
                   "z_BMI", 
                   "z_t2d", 
                   "z_SBP", 
                   "z_PP", 
                   "z_DBP", 
                   "z_LDL", 
                   "z_HDL", 
                   "z_insulin_fold_change",
                   "z_adiponectin", 
                   "z_menarche", 
                   "z_birth_weight", 
                   "z_SHBG_m", 
                   "z_testosterone_m", 
                   "z_SHBG_f", 
                   "z_height")
#change the name for cosmetic labels
desired_order = gsub("z_", "", desired_order)
desired_order = gsub("aspartate_aminotransferase", "AST", desired_order)
desired_order = gsub("alanine_aminotransferase", "ALT", desired_order)
desired_order = gsub("TNF_related_apoptosis_inducing_ligand_receptor_2", "TRAIL-R2", desired_order)
desired_order = gsub("_", " ", desired_order)
desired_order = gsub("t2d", "T2D", desired_order)

#apply same cosmetic changes in the dt 
gg_heat$Z = gsub("z_", "", gg_heat$Z)
gg_heat$Z = gsub("aspartate_aminotransferase", "AST", gg_heat$Z)
gg_heat$Z = gsub("alanine_aminotransferase", "ALT", gg_heat$Z)
gg_heat$Z = gsub("TNF_related_apoptosis_inducing_ligand_receptor_2", "TRAIL-R2", gg_heat$Z)
gg_heat$Z = gsub("_", " ", gg_heat$Z)
gg_heat$Z = gsub("t2d", "T2D", gg_heat$Z)

#put them as factor to order them like the original hierarchical clustering
gg_heat$Z <- factor(gg_heat$Z, levels = desired_order)

fwrite(gg_heat, file = paste0("./clustering_result/imputeSCOPA/", out, "_table_per_cluster.txt"), sep = "\t", quote = F)

gg_heat[P < 5e-8, P := 5e-8] #truncate to avoid blurred pictures

#EDIT: remove some pale values that does not add much 
gg_heat <- gg_heat[!Z %in% c("TNF alpha", "IL6", "IL10", "menopause", "PAI 1", "ostradiol m", "oestradiol f"),]

p = ggplot(gg_heat, aes(y = Z, x = Cluster, fill = -log(P)*sign(Beta))) +
  geom_tile(color = "white", linewidth = 1) +
  scale_fill_gradient2(low = "#2166AC", mid = "white", high = "#B2182B") +  # You can adjust the color gradient as needed
  theme_minimal() +
  xlab("") + ylab("") +
  theme(text = element_text(size = 20),
        axis.title.x = element_text(size = 24),
        axis.title.y = element_text(size = 24), 
        legend.title = element_blank(), 
        legend.key.height = unit(2, "cm")) 
ggsave(p, filename = paste0("./clustering_result/imputeSCOPA/", out, "_heatmap_subset_per_cluster.png"), width = 10, height = 12, dpi = 600)

