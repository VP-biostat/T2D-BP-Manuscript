library(readxl)
library(plotly)
library(dplyr)
library(tidyr)
library(tibble)

# Read the Supplementary Table 4
df = read_xlsx("YOUR_PATH_TO_SUPPLEMENTARY_DATA_3", sheet = 4, skip = 2)

## STARTING WITH SUZUKI ET AL 2024
count_df = table(df$`Attributed cluster`, df$Cluster...8)

# Define the order of source and target nodes
source_nodes <- rownames(count_df)
# [1] "Higher adiposity"           "Inverse T2D-BP risk"        "Metabolic Syndrome"        
# [4] "Reduced beta-cell function" "Vascular dysfunction"
target_nodes <- colnames(count_df)
# [1] "Beta cell -PI"          "Beta cell +PI"          "Body fat"               "Lipodystrophy"         
# [5] "Liver/lipid metabolism" "Metabolic syndrome"     "Obesity"                "Residual glycaemic" 

# Define nodes and links
nodes <- list(
  label = c(source_nodes, target_nodes),
  color = c("#51d3a6", "#f98b83", "#a9ab00", "#f471ff", "#00b4fc",
            "#e41a1c", "#3c82bc", "#4daf4a", "#9c50a8", "#ff8300", "#fefe34", "#aa5829", "#f585be"),
  pad = 15,
  thickness = 20,
  line = list(color = "black", width = 0.5)
)

# Create numeric indices
# (0-based indexing for Plotly: sources 0–4, targets 5–12)
source_indices <- setNames(seq_along(source_nodes) - 1, source_nodes)
target_indices <- setNames(seq_along(target_nodes) - 1 + length(source_nodes), target_nodes)

# Convert count_df to long format
links_df <- count_df %>%
  as.data.frame() %>%
  filter(Freq > 0) %>%
  mutate(
    source_id = source_indices[Var1],
    target_id = target_indices[Var2]
  )

# Build the links list for Plotly Sankey
links <- list(
  source = links_df$source_id,
  target = links_df$target_id,
  value = links_df$Freq
)

# Create Sankey diagram
fig <- plot_ly(
  type = "sankey",
  orientation = "h",
  node = nodes,
  link = links
)

# fig <- fig %>% layout(title = "Comparison of T2D-BP hierarchical clusters with\nSuzuki et al. (2024) hierarchical clusters")
fig
htmlwidgets::saveWidget(fig, "Suzuki_hieclu_Sankey.html")
# Save to PNG (approx. 300 dpi)
webshot::webshot(
  "Suzuki_hieclu_Sankey.html",         # input HTML
  file = "Suzuki_hieclu_Sankey.png",   # output PNG
  vwidth = 700,                # virtual browser width
  vheight = 700,               # virtual browser height
  zoom = 1                      # higher zoom ≈ higher DPI
)


## MOVING TO SMITH ET AL 2024
count_df = table(df$`Attributed cluster`, df$Cluster...11)

# Define the order of source and target nodes
source_nodes <- rownames(count_df)
# [1] "Higher adiposity"           "Inverse T2D-BP risk"        "Metabolic Syndrome"        
# [4] "Reduced beta-cell function" "Vascular dysfunction"
target_nodes <- colnames(count_df)
# [1] "ALP Neg"         "Beta Cell 1"     "Beta Cell 2"     "Cholesterol"     "Hyper Insulin"  
# [6] "Lipodystrophy 1" "Lipodystrophy 2" "Liver-Lipid"     "Obesity"         "Proinsulin"     
# [11] "SHBG-LpA"   

# Define nodes and links
nodes <- list(
  label = c(source_nodes, target_nodes),
  color = c("#51d3a6", "#f98b83", "#a9ab00", "#f471ff", "#00b4fc",
            "#9e0142", "#d63e4f", "#f56e44", "#fdae61", "#fee08b", "#ffffc0", "#f0fe9e", "#abdda4", "#66c2a5", "#328cc4", "#5e4fa2"),
  pad = 15,
  thickness = 20,
  line = list(color = "black", width = 0.5)
)

# Create numeric indices
# (0-based indexing for Plotly: sources 0–4, targets 5–12)
source_indices <- setNames(seq_along(source_nodes) - 1, source_nodes)
target_indices <- setNames(seq_along(target_nodes) - 1 + length(source_nodes), target_nodes)

# Convert count_df to long format
links_df <- count_df %>%
  as.data.frame() %>%
  filter(Freq > 0) %>%
  mutate(
    source_id = source_indices[Var1],
    target_id = target_indices[Var2]
  )

# Build the links list for Plotly Sankey
links <- list(
  source = links_df$source_id,
  target = links_df$target_id,
  value = links_df$Freq
)

# Create Sankey diagram
fig <- plot_ly(
  type = "sankey",
  orientation = "h",
  node = nodes,
  link = links
)

# fig <- fig %>% layout(title = "Comparison of T2D-BP hierarchical clusters with\nSmith et al. (2024) 'soft' clusters")
fig
htmlwidgets::saveWidget(fig, "Smith_hieclu_Sankey.html")
# Save to PNG (approx. 300 dpi)
webshot::webshot(
  "Smith_hieclu_Sankey.html",         # input HTML
  file = "Smith_hieclu_Sankey.png",   # output PNG
  vwidth = 700,                # virtual browser width
  vheight = 700,               # virtual browser height
  zoom = 1                      # higher zoom ≈ higher DPI
)


## MOVING TO VAURA ET AL 2024
count_df = table(df$`Attributed cluster`, df$Cluster...14)

# Define the order of source and target nodes
source_nodes <- rownames(count_df)
# [1] "Higher adiposity"           "Inverse T2D-BP risk"        "Metabolic Syndrome"        
# [4] "Reduced beta-cell function" "Vascular dysfunction"
target_nodes <- colnames(count_df)
# [1] "Dyslipidemia N Loci = 14"  "Hypolipidemia N Loci = 27" "Obesity N Loci = 28"      
# [4] "Short Stature N Loci = 25"

# Define nodes and links
nodes <- list(
  label = c(source_nodes, target_nodes),
  color = c("#51d3a6", "#f98b83", "#a9ab00", "#f471ff", "#00b4fc",
            "#8dd3c7", "#ffffb3", "#bebada", "#fb8072"),
  pad = 15,
  thickness = 20,
  line = list(color = "black", width = 0.5)
)

# Create numeric indices
# (0-based indexing for Plotly: sources 0–4, targets 5–12)
source_indices <- setNames(seq_along(source_nodes) - 1, source_nodes)
target_indices <- setNames(seq_along(target_nodes) - 1 + length(source_nodes), target_nodes)

# Convert count_df to long format
links_df <- count_df %>%
  as.data.frame() %>%
  filter(Freq > 0) %>%
  mutate(
    source_id = source_indices[Var1],
    target_id = target_indices[Var2]
  )

# Build the links list for Plotly Sankey
links <- list(
  source = links_df$source_id,
  target = links_df$target_id,
  value = links_df$Freq
)

# Create Sankey diagram
fig <- plot_ly(
  type = "sankey",
  orientation = "h",
  node = nodes,
  link = links
)

# fig <- fig %>% layout(title = "Comparison of T2D-BP hierarchical clusters with\nVaura et al. (2024) 'soft' clusters")
fig
htmlwidgets::saveWidget(fig, "Vaura_softclu_Sankey.html")
# Save to PNG (approx. 300 dpi)
webshot::webshot(
  "Vaura_softclu_Sankey.html",         # input HTML
  file = "Vaura_softclu_Sankey.png",   # output PNG
  vwidth = 700,                # virtual browser width
  vheight = 700,               # virtual browser height
  zoom = 1                      # higher zoom ≈ higher DPI
)
