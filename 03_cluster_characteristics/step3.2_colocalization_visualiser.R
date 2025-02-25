library(data.table)
library(ggplot2)
library(dplyr)
library(ggrepel)

out <- "YOUR_OUTPUT_NAMEYOUR_OUTPUT_NAME"

#function to visualise the colocalisation per cluster/eQTL pairs using a cool Manhattan plot 
manhattan_colocalization <- function(tissues = c("Adipose_Subcutaneous", "Thyroid", "Artery_Tibial"), 
                                     coloc_table = NA, 
                                     text_threshold = 0.8) {
  ###A function to plot colocalization analysis
  ### the 'tissues' parameter needs a list of characters corresponding to tissue names in 'coloc_table'
  ### the 'coloc_table' needs to have the following columns: tissue, PP.H3.abf, chr, position, H3+H4, SYMBOL
  ### the text_threshold is a numeric between 0.8 and 1 to show the associated gene to one signal
  
  coloc_table <- as.data.table(coloc_table)
  coloc_table <- coloc_table[, .(`tissue`, `PP.H4.abf`, `chr`, `position`, `H3+H4`, `SYMBOL`)]
  sub_coloc_table <- coloc_table[`tissue` %in% tissues,] #select only the limited number of tissues
  sub_coloc_table <- sub_coloc_table %>% arrange(chr, position)
  sub_coloc_table[, cumsum := .I]
  
  custom_palette <- c("#1f78b4", "#33a02c", "#e31a1c", "#ff7f00", "#6a3d9a",
                      "#a6cee3", "#b2df8a", "#fb9a99", "#fdbf6f", "#cab2d6",
                      "#1f4424", "#7fc97f", "#b15928", "#fee08b", "#377eb8",
                      "#984ea3", "#4daf4a", "#377eb8", "#ff7f00", "#984ea3",
                      "#a65628", "#f781bf")
  
  # Filter data to keep only non-duplicated labels for the strongest signals
  filtered_data <- sub_coloc_table %>%
    filter(`PP.H4.abf` > text_threshold) %>%
    group_by(SYMBOL) %>%
    filter(row_number() == 1) %>%
    ungroup()
  
  p <- ggplot(data = sub_coloc_table, mapping = aes(x = cumsum, y = `PP.H4.abf`, shape = tissue, col = factor(chr))) +
        geom_point(size = 2) +
        scale_color_manual(values = custom_palette[1:length(unique(sub_coloc_table$chr))]) +
        ylim(0,1.10) + 
        labs(x = "Genomic Position",
             y = "Posterior Probability of colocalization (PP.H4)",
             color = "Chromosome", 
             shape = "Tissue") +
        geom_hline(yintercept = 0.8, linetype = "dotted", color = "gray") +
        theme_minimal() +
        theme(axis.text.x = element_blank()) + 
        geom_text_repel(data = filtered_data,
                          aes(x = cumsum, y = `PP.H4.abf`, label = SYMBOL),
                          box.padding = 0.5,
                          point.padding = 0.1,
                          force = 2,
                          segment.size = 0.2,
                          segment.color = "grey", 
                          max.overlaps = 40)
  
  return(p)
  
}

#Apply the Manhattan Plot visualiser on 3 tissues and each cluster
for (cl in unique(subcluster$gcluster)) {
  for (tiss in c("Adipose_Subcutaneous", "Thyroid", "Artery_Tibial")) {
    
    cat("Trying to build Manhattan plot for cluster", cl, "and tissue", tiss, "\n")
    sub_subcluster <- subcluster[gcluster == cl, ]
    
    #MERGE THE INFO ON THE CLUSTER AND THE COLOCALIZATION TABLE
    sub_sorted_signig <- merge(signif_results, sub_subcluster, by.x = "top_snp", by.y = "mapping_chrpos")
    cat(nrow(sub_sorted_signig), "signals to plot \n")
    
    #RUN THE PLOT FUNCTION 
    p <- manhattan_colocalization(coloc_table = sub_sorted_signig, tissues = tiss, text_threshold = 0.85)
    ggsave(p, filename = paste0("./colocalization/", out, "/manhattan_plot_colocalisation_tissue_", tiss, "_cluster_", cl, ".png"), width = 20, height = 7)
    
  }
}

