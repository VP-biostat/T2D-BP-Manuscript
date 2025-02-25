library(data.table)
library(ggplot2)
library(ggrepel)

path = "./pgs_after_clustering/PHEWAS/" #PheWAS files + partitioned PGS extracted from UKB
res_files = list.files(path = path, pattern = "\\d+_to_\\d+") #PheWAS files + partitioned PGS extracted from UKB

truncate_name <- function(name, n_words = 3) {
  words <- unlist(strsplit(name, "\\s+"))  # Split the name into words
  if (length(words) > n_words) {
    paste(words[1:n_words], collapse = " ")  # Join the first n_words with space
  } else {
    name  # Return the original name if it has less than n_words
  }
}

results = data.frame()
for (file in res_files) {
  res = fread(paste0(path, file))
  if (ncol(res) == 13) results = rbind(results, res)
  else cat("\n",file, "has an issue")
}

diseases = fread("./MEGAPHEWAS_icd10_list.tsv")
df = merge(results, diseases[, c("icd10", "name")], by.x = "Phenotype", by.y = "icd10")
# df$truncated_name <- sapply(df$name, truncate_name)
df$truncated_name <- gsub("\\,.*", "", df$name)
df$truncated_name <- gsub("without.*", "", df$truncated_name)

df$PGS = gsub("1.z.T2D", "Inverse T2D-BP risk", df$PGS)
df$PGS = gsub("2.z.T2D", "Metabolic Syndrome", df$PGS)
df$PGS = gsub("3.z.T2D", "Higher adiposity", df$PGS)
df$PGS = gsub("4.z.T2D", "Vascular dysfunction", df$PGS)
df$PGS = gsub("5.z.T2D", "Reduced beta-cell function", df$PGS)
df$PGS = factor(df$PGS, levels = c("Inverse T2D-BP risk", "Metabolic Syndrome", "Higher adiposity", "Vascular dysfunction", "Reduced beta-cell function"))

color_legend = data.frame(Cat = c("A", "B", "C", "D", "E", "F", "G", "H", "I", 
                                  "J", "K", "L", "M", "N", "O", "P", "Q", "R", 
                                  "S", "T", "U", "V", "W", "X", "Y", "Z"),
                          Cat_name = c("Infectious and parasitic diseases", 
                                        "Infectious and parasitic diseases", 
                                        "Cancer", 
                                        "Neoplasms, blood, and blood-forming organs", 
                                        "Endocrine, nutritional, or metabolic", 
                                        "Mental and behavioral disorders", 
                                        "Nervous system", 
                                        "Eyes, ears, nose, and throat", 
                                        "Circulatory system",
                                        "Respiratory system",
                                        "Digestive system",
                                        "Skin",
                                        "Musculoskeletal system",
                                        "Genitourinary system",
                                        "Pregnancy and childbirth",
                                        "Perinatal conditions",
                                        "Congenital and chromosomal abnormalities",
                                        "Abnormal clinical and lab findings",
                                        "Injury, poisoning, and other external causes",
                                        "Injury, poisoning, and other external causes",
                                        "Used for emergency designation",
                                        "External causes of morbidity",
                                        "External causes of morbidity",
                                        "External causes of morbidity",
                                        "External causes of morbidity",
                                        "Factors influencing health status and contact with health services"))

df = merge(df[, Cat := substr(Phenotype, 1, 1)], color_legend, by = "Cat")
df$Shape <- ifelse(df$Effect > 1, "Positive", "Negative")

crop_df = df[-log10(P_value) > 40, P_value := 1e-40]
crop_df = crop_df[!Cat %in% c("A", "B", "C", "D", "R", "S", "T", "U", "V", "W", "X", "Y", "Z"),]

p_thres = 0.05/nrow(crop_df)

p <- ggplot(crop_df, aes(x = Phenotype, y = -log10(P_value), color = Cat_name, shape = Shape)) +
  geom_point(size = 1.2) +
  scale_shape_manual(values = c("Positive" = 2, "Negative" = 6)) + # "+" for positive, "-" for negative
  scale_colour_manual(values = c("#E69F00", "#56B4E9", "#009E73", "#355E3B", 
                                 "#0072B2", "#D55E00", "#CC79A7", "#F0A3FF", 
                                 "#40E0D0", "#8B4513")) +
  geom_text_repel(aes(label = ifelse(P_value < p_thres*1e-5, truncated_name, "")),
                  max.overlaps = 6, # Adjust this to limit the number of overlapping labels
                  size = 4, # Adjust text size
                  box.padding = 0.35, # Padding around text boxes
                  point.padding = 0.3, # Padding around points
                  segment.color = "grey50", # Color of the line connecting text to points
                  segment.size = 0.5) + # Thickness of the line
  labs(x = "Phenotype", y = expression(-log[10](p)), color = "Disease category", shape = "Effect direction") +
  theme_minimal() +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    legend.position = "none",
    axis.title.x = element_text(size = 20),
    axis.title.y = element_text(size = 20), 
    axis.text.x.bottom = element_blank()
  ) +
  facet_wrap(~ PGS, ncol = 1) +
  geom_hline(yintercept = -log10(p_thres), color = "gray", linetype = "dashed")

ggsave(p, filename = "./pgs_after_clustering/PHEWAS/legend_new_phewas_plot_for_fig4c.png", dpi = 600, height = 16, width = 12)

p1 <- ggplot(crop_df, aes(x = Phenotype, y = -log10(P_value), color = Cat_name, shape = Shape)) +
  geom_point(alpha = 0.8, size = 1.2) +
  scale_shape_manual(values = c("Positive" = 2, "Negative" = 6)) + # "+" for positive, "-" for negative
  scale_colour_manual(values = c("#E69F00", "#56B4E9", "#009E73", "#355E3B", 
                                 "#0072B2", "#D55E00", "#CC79A7", "#F0A3FF", 
                                 "#40E0D0", "#8B4513")) +
  geom_text_repel(aes(label = ifelse(P_value < p_thres, truncated_name, "")),
                  max.overlaps = 8, # Adjust this to limit the number of overlapping labels
                  size = 4, # Adjust text size
                  box.padding = 0.35, # Padding around text boxes
                  point.padding = 0.3, # Padding around points
                  segment.color = "grey50", # Color of the line connecting text to points
                  segment.size = 0.5) + # Thickness of the line
  labs(x = "Phenotype", y = expression(-log[10](p)), color = "Disease category", shape = "Effect direction") +
  theme_minimal() +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    # axis.text.x.top = element_text(size = 24),
    axis.title.x = element_text(size = 20),
    axis.title.y = element_text(size = 20), 
    # strip.text.x.top = element_text(size = 10), # Increase facet label font size here
    axis.text.x.bottom = element_blank(),
  ) +
  facet_wrap(~ PGS, ncol = 1) +
  geom_hline(yintercept = -log10(p_thres), color = "gray", linetype = "dashed")


ggsave(p1, filename = "./pgs_after_clustering/PHEWAS/phewas_plot_for_fig5c.png", dpi = 600, height = 32, width = 12)
