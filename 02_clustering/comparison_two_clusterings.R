library(data.table)
library(ggplot2)

#insert plot of clustering 1 (referred to as "old")
old <- fread("the path of your first clustering to compare")
#insert plot of clustering 2 (referred to as "latest")
latest <- fread("the path of your second clustering")

#merge them 
m <- merge(latest, old, by = "rsid")

#check their clustering
t <- table(m$gcluster.x, m$gcluster.y)
#format for ggplot 
t2 <- reshape2::melt(t)
names(t2) <- c("Clustering 2", "Clustering 1", "N")
t2$`Clustering 2` <- factor(t2$`Clustering 2`)
t2$`Clustering 1`<- factor(t2$`Clustering 1`)
t2$Proportion <- NA
for (old_latest_cluster in unique(t2$`Clustering 1`)) {
  tot_variant_cluster <- sum(t2[which(t2$`Clustering 1` == old_latest_cluster), "N"])
  t2[which(t2$`Clustering 1` == old_latest_cluster), "Proportion"] <- t2[which(t2$`Clustering 1` == old_latest_cluster), "N"] / tot_variant_cluster
}

#do a non clustering heatmap
p <- ggplot(data = t2, aes(x = `Clustering 2`, y = `Clustering 1`, fill = Proportion)) +
       geom_tile() + 
       geom_text(aes(label = N)) +
       scale_fill_gradient(low = "white", high = "#B2182B") +
       labs(fill = "Proportion of Adjusted\nclustering variants found\nin the Unadjusted Clustering", 
            x = "Adjusted Clustering", y = "Unadjusted Clustering")
ggsave(p, filename = "./clustering_result/comparison_two_clusterings.png", width = 8, height = 6)
