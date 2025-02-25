library(data.table)
library(readxl)
library(ggplot2)
library(dplyr)
library(gridExtra)

st = data.table(read_excel("BP_T2D_Supplementary_Table_240816.xlsx", sheet = 8, skip = 1))
st$`eQTL tissue` = gsub("Not Sun Exposed", "Sun-", st$`eQTL tissue`)
st$`eQTL tissue` = gsub("Sun Exposed", "Sun+", st$`eQTL tissue`)
st$`eQTL tissue` = gsub("Esophagus Gastroesophageal", "Esophagus GE", st$`eQTL tissue`)

# Assuming 'st' is your data frame and it is already loaded
# Aggregating the data by eQTL tissue and Cluster
agg_data <- st %>%
  group_by(`eQTL tissue`, `Cluster`, `GWAS trait`) %>%
  summarise(Intersection_Size = n()) %>%
  ungroup()

argg_data <- agg_data %>%
  group_by(`eQTL tissue`) %>%
  summarise(Tot_coloc_per_tissue = sum(Intersection_Size)) %>%
  filter(Tot_coloc_per_tissue >= 100) %>%
  ungroup()

clu_order = c("Inverse T2D-BP risk", "Metabolic Syndrome", "Higher adiposity", "Vascular dysfunction", "Reduced beta-cell function")


# Creating the bar plot for Intersection size by eQTL tissue
p1 <- ggplot(agg_data[agg_data$`eQTL tissue` %in% argg_data$`eQTL tissue`, ], aes(x = reorder(`eQTL tissue`, -Intersection_Size), y = Intersection_Size, fill = factor(Cluster, levels = clu_order))) +
  geom_bar(stat = "identity") +
  theme_minimal() +
  # geom_text(aes(label = paste0(round(Percentage, 1), "%")), 
  #           position = position_stack(vjust = 0.5), size = 5) +
  labs(x = "", y = "No of colocalisations", fill = "Cluster") +
  theme(axis.text.x = element_text(angle = 90),
        text = element_text(size = 20),
        axis.title.x = element_text(size = 24),
        axis.title.y = element_text(size = 24),
        legend.position = "none")
ggsave(p1, file = './colocalization/fig3a_tot_count_per_tissue.png', width = 10, height = 8, dpi = 600)


# Do the same plot, but using GWAS trait as colour 
p2 <- ggplot(agg_data[agg_data$`eQTL tissue` %in% argg_data$`eQTL tissue`, ], aes(x = factor(`eQTL tissue`, levels = tissue_order), y = Intersection_Size, fill = `GWAS trait`)) +
  geom_bar(stat = "identity") +
  theme_minimal() +
  labs(x = "eQTL Tissue", y = "Number per GWAS of colocalisations", fill = "GWAS trait") +
  theme(axis.text.x = element_blank(), 
        text = element_text(size = 20),
        axis.title.x = element_blank(),
        axis.title.y = element_text(size = 20), 
        legend.position = "none")
ggsave(p2, file = './colocalization/fig3a_tot_count_per_gwas.png', width = 10, height = 8, dpi = 600)


subst = st[`eQTL tissue` %in% argg_data$`eQTL tissue`,]

#find the SNPs that colocalise multiple times 
snp_duplicated = unique(st[duplicated(st$`SNV from clustering`) & duplicated(st$`GWAS trait`), ]$`SNV from clustering`)
#remove them in this st: 
singlest = st[!`SNV from clustering` %in% snp_duplicated,]
# Aggregating the data by eQTL tissue and Cluster
agg_singlest <- singlest %>%
  group_by(`eQTL tissue`, `Cluster`, `GWAS trait`) %>%
  summarise(Intersection_Size = n()) %>%
  ungroup()

# Extract the order of eQTL tissue from p1
tissue_order <- levels(reorder(agg_data[agg_data$`eQTL tissue` %in% argg_data$`eQTL tissue`, ]$`eQTL tissue`, -agg_data[agg_data$`eQTL tissue` %in% argg_data$`eQTL tissue`, ]$Intersection_Size))

# Manual addition of the non specific tissue to match with p1 
agg_singlest <- rbind(agg_singlest, c("Skin Sun- Suprapubic", "Inverse T2D-BP risk", "SBP", 0))
agg_singlest <- rbind(agg_singlest, c("Adipose Visceral Omentum", "Metabolic Syndrome", "T2D", 0))
agg_singlest <- rbind(agg_singlest, c("Colon Sigmoid", "Higher adiposity", "T2D", 0))
agg_singlest <- rbind(agg_singlest, c("Esophagus GE Junction", "Vascular dysfunction", "SBP", 0))
agg_singlest <- rbind(agg_singlest, c("Pancreas", "Reduced beta-cell function", "T2D", 0))
agg_singlest <- rbind(agg_singlest, c("Brain Caudate basal ganglia", "Vascular dysfunction", "SBP", 0))

# Ensure the same order in p2
p3 <- ggplot(agg_singlest[agg_singlest$`eQTL tissue` %in% argg_data$`eQTL tissue`, ], aes(x = factor(`eQTL tissue`, levels = tissue_order), y = as.numeric(Intersection_Size), fill = factor(Cluster, levels = clu_order))) +
  geom_bar(stat = "identity") +
  theme_minimal() +
  labs(x = "eQTL Tissue", y = "No single-tissue colocalisations", fill = "Cluster") +
  scale_y_reverse() +
  theme(axis.text.x = element_text(angle = 90), 
        
        text = element_text(size = 20),
        axis.title.x = element_text(size = 20),
        axis.title.y = element_text(size = 20), 
        legend.position = "none")
ggsave(p3, file = './colocalization/fig3a_singlecoloc_count_per_tissue.png', width = 10, height = 8, dpi = 600)

# Use grid.arrange to align the plots
grid.arrange(p1, p3, ncol = 1)




# Aggregating the data by eQTL tissue and Cluster
argg_singlest <- agg_singlest %>%
  group_by(`eQTL tissue`) %>%
  summarise(Tot_singlecoloc_per_tissue = sum(Intersection_Size)) %>%
  ungroup()

argg = merge(argg_data, argg_singlest, by = "eQTL tissue")
argg$perc = argg$Tot_singlecoloc_per_tissue/argg$Tot_coloc_per_tissue

agg_data = merge(agg_data, argg, by = "eQTL tissue", all.x = T)

agg_data$Tistis = paste0(agg_data$`eQTL tissue`, " ", round(agg_data$perc*100, 1), "%")

# Creating the bar plot for Intersection size by eQTL tissue
p4 <- ggplot(agg_data[agg_data$`eQTL tissue` %in% argg_data$`eQTL tissue`, ], aes(y = reorder(Tistis, Intersection_Size), x = Intersection_Size, fill = Cluster)) +
  # geom_text(aes(label = paste0(round(perc, 1), "%")), 
  #           position = position_stack(vjust = 0.5), size = 5) +
  geom_bar(stat = "identity") +
  theme_minimal() +
  labs(x = "eQTL Tissue", y = "No of colocalisations", fill = "Cluster") +
  theme(text = element_text(size = 18),
        axis.title.x = element_text(size = 24),
        axis.title.y = element_text(size = 24)) 
# legend.position = "none")
ggsave(p4, file = './colocalization/fig3a_tot_count_per_tissue.png', width = 10, height = 8, dpi = 600)
