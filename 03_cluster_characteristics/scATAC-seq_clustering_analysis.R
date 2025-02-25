library(data.table)

### DATA GATHERING AND MERGING -------

# cCRE_by_cell_type from http://catlas.org/catlas_downloads/humantissues/cCRE_by_cell_type/
# =================
#   
#   This directory stores cCRE by cell type, containing three files:
#   
#   1. matrix.mtx.gz: matrix in matrix market format, each line represents a cCRE-cell type pair. For example, "2 22" means cCRE No. 2 is accessibile in cell type No. 22.
#   2. celltypes.txt.gz: column names (cell types) of the matrix.
#   3. cCREs.bed.gz: row names (cCREs) of the matrix.
# 
# (NOTE: the original files were wrong and the correct files were uploaded on Nov 26, 2021.)
# README: http://catlas.org/catlas_downloads/humantissues/README.md
ccre <- fread("./scATAC-seq_zhang/cCREs.bed.gz")
cell_type <- fread("./scATAC-seq_zhang/celltypes.txt.gz", sep = "\t", header = F)
matrix <- fread("./scATAC-seq_zhang/matrix.tsv.gz", sep = " ", skip = 3)

#create a matrix of the clustered SNPs per cell type
# ###CHANGE IF NEEDED
out <- "YOUR_OUTPUT_NAME"
cluster <- fread(paste0("./clustering_result/imputeSCOPA/", out, "_5clusters_table.txt"))

df <- cluster[, c("rsid", "GRCh37_chrpos", "gcluster")]

#MIGHT TAKE TIME, this is why this step from here is not necessary everytime
library(biomaRt)
mart <- useMart("ENSEMBL_MART_SNP", dataset = "hsapiens_snp")
snp_info <- getBM(attributes = c('refsnp_id', 'chr_name', 'chrom_start', 
                                 'consequence_type_tv'),
                  filters = "snp_filter",
                  values = df$rsid,
                  mart = mart, useCache = F)

snp_info_b38 <- data.table(snp_info[!duplicated(snp_info$refsnp_id), ])
snp_info_b38[, chrpos := paste0(chr_name, ":", chrom_start)]
names(snp_info_b38) <- c("rsid", "chr", "position", "GRCh38_chrpos")

df <- merge(df, snp_info_b38, by = 'rsid')
fwrite(df, file = paste0("./clustering_result/imputeSCOPA/", out, "_snps_b37_b38.txt"), sep = "\t", quote = F)

#now the big table to count enrichment 
df <- cbind(df, data.table(matrix(data = 0, nrow = nrow(df), ncol = nrow(cell_type))))
setnames(df, paste0("V", seq(1:nrow(cell_type))), cell_type$V1)

#and do a for loop (easier, but can be computed in a function for more efficiency)
for (i in 1:nrow(df)) {
  if (i %% 1000 == 0) cat('Have done the', i, 'first SNP - cCRE associations\n')
  current_snp <- data.table(chr = as.numeric(gsub(":.*", "", df[i , ]$GRCh38_chrpos)), 
                            position = as.numeric(gsub(".*:", "", df[i, ]$GRCh38_chrpos)))
  
  current_ccre_ind <- which(ccre$V1 == paste0("chr", current_snp$chr) & ccre$V2 <= current_snp$position & ccre$V3 >= current_snp$position)
  if (length(current_ccre_ind) == 1) {
    current_matrix <- matrix[V1 == current_ccre_ind, ]
    current_cell_type <- cell_type$V1[current_matrix$V2]
    df[i, (current_cell_type) := lapply(.SD, function(x) x + 1), .SDcols = current_cell_type]
  } else if (length(current_ccre_ind) > 1) cat(i, "indice has more than one match?!")
}
fwrite(df, file = paste0("./scATAC-seq_zhang/scATAQ-seq_zhang_", out, "_clustering_enrichment.csv"), sep =";", quote = F)


### LOGISTIC REGRESSION ----------
library(logistf)

#the idea is to use logit link function 
#to put constained logit using EXON and 3'UTR and 5'UTR as covariates
#then to performed the logit unconstrained with 1 cluster and 1 tissue cCRE count
#then compare the deviances of the two models by an anova test of enrichment (chisq test)

cts <- cell_type$V1
cls <- na.omit(unique(df$gcluster))

if (exists("logit_df")) rm(logit_df)

for (cl in cls) {
  
  cat("Cluster", cl, "\n")
  
  rsid_in_cluster <- df[gcluster == cl, ]$rsid
  
  subdf <- df[rsid %in% rsid_in_cluster | associated_index_snp %in% rsid_in_cluster, ]
  
  for (ct in cts) {
    
    cat("\tCell type", ct, "\n")
    
    formula_constrained <- as.formula("is_index_snp ~ EXON + `3_prime_UTR` + `5_prime_UTR`")
    formula_unconstrained <- as.formula(paste0("is_index_snp ~ EXON + `3_prime_UTR` + `5_prime_UTR` + `", ct, "`"))
    
    lcon <- logistf(data = subdf, formula_constrained) #logistf model with delta = 0
    luncon <- logistf(data = subdf, formula_unconstrained) #logistf model with unconstrained delta
    
    lfr <- anova(lcon, luncon) #comparison of the deviances models integrated for logistf models, will give a chi-squares statistics + associated p-value
    
    temp <- data.table(cluster = cl, 
                       cell_type = ct, 
                       logistf_coef = luncon$coefficients[5], 
                       deviance_diff_chisq = lfr$chisq, 
                       deviance_diff_pval = lfr$pval)
    
    if (exists("logit_df")) {
      logit_df <- rbind(logit_df, temp)
    } else {
      logit_df <- temp
    }
    
  }
  
}
fwrite(logit_df, file = paste0("./scATAC-seq_zhang/", out, "_logistf_table.csv"), sep = ";", quote = F)


### VISUALISATION -------

library(ggplot2)
library(dplyr)

#visualisation for logistf pval:
gg_heat <- logit_df
pval_bonferroni <- 0.05/length(cts)
pval_int <- 0.05/15
interesting_cell_type <- unique(gg_heat[deviance_diff_pval < pval_bonferroni, ]$cell_type)
interesting_cell_type2 <- unique(gg_heat[deviance_diff_pval < pval_int, ]$cell_type)
# gg_heat[deviance_diff_pval < 5e-8, deviance_diff_pval := 5e-8] #truncate to avoid blurred pictures
gg_heat <- gg_heat[cell_type %in% interesting_cell_type2, ] #take only the one of interest


p = ggplot(gg_heat, aes(y = cell_type, x = cluster, fill = -log(deviance_diff_pval))) +
  geom_tile(color = "white") +
  scale_fill_gradient(low = "white", high = "#B2182B") +  # You can adjust the color gradient as needed
  theme_minimal() +
  xlab("") + ylab("") +
  # geom_text(data = filter(gg_heat, deviance_diff_pval < pval_int),
  #           aes(label = "*"), color = "white", size = 6) +
  geom_text(data = filter(gg_heat, deviance_diff_pval < pval_bonferroni),
            aes(label = "***"), color = "white", size = 6) +
  theme(text = element_text(size = 20),
        axis.title.x = element_text(size = 24),
        axis.title.y = element_text(size = 24), 
        legend.title = element_blank(), 
        legend.key.height = unit(2, "cm")) 
ggsave(p, filename = paste0("./scATAC-seq_zhang/", out, "_logistf_", pval_int, "_heatmap_per_cluster.png"), width = 8, height = 12)
