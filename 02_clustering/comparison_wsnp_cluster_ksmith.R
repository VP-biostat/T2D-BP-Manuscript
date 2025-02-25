library(readxl)
library(LDlinkR)
library(data.table)

token <- 'your LDlink token'

#use the manuscript supplementary tables
cluster <- data.table(read_excel("./BP_T2D_Supplementary_Table_240816.xlsx", sheet = 6, skip =1))
cluster[, KSmith_snp := "none"]
cluster[, KSmith_r2 := 0]
cluster[, KSmith_cluster := "unassigned"]

#use k smith et al. supplementary tables
ksmith <- data.table(read_excel("./41591_2024_2865_MOESM2_ESM.xlsx", sheet = 5, skip = 3, range = "A4:O654"))
ksmith[, Chromosome := as.numeric(gsub("\\_.*", "", VAR_ID_hg19))] #fill the chromosome col

for (i in 1:nrow(cluster)) {
  if (exists("mat")) rm(mat)
  
  snp = cluster[i, ]$rsid
  chr = as.numeric(gsub("\\:.*", "", cluster[i, ]$`chr:pos (b37)`))
  
  perc_loading = round(i / nrow(cluster)*50)
  perc_remaining = 50-perc_loading
  loading_bar = paste0(strrep("=", perc_loading), strrep("-", perc_remaining))
  cat("\n|",loading_bar , "|\nSNP: ", snp, "\tin chromosome: ", chr, "\tin cluster: ", cluster[i, ]$`Attributed cluster`)
  
  tempsnps = c(snp, ksmith[Chromosome == chr, ]$rsID)
  
  #try to query LDmatrix from LDlink
  mat = try(LDmatrix(tempsnps, pop = "CEU", r2d = "r2", token = token))
  
  # but there might be several options:
  #EDIT VP 23/11/22: 
  # 0. LDlink server was not available while making the LDproxy request
  while (inherits(mat, "try-error")) {
    #just try it again until there is an answer from LDlink, 
    #here we are trying to fix the error " Bad Gateway (HTTP 502) " --> no object created
    mat = try(LDmatrix(tempsnps, pop = "CEU", r2d = "r2", token = token))
  }
  while (is.null(mat)) {
    #just try it again until there is an answer from LDlink, 
    #LDlink will return NULL if there is no answer from the server " Server is down "
    mat = try(LDmatrix(tempsnps, pop = "CEU", r2d = "r2", token = token))
  }
  
  rind_osnp = which(mat$RS_number == snp)
  cmat = data.table(mat[rind_osnp, ])
  cmat = cmat[, !..snp]
  rmat = data.table(rsid = names(cmat)[-1], `r2` = as.vector(t(cmat))[-1])
  rmat[is.na(r2), r2 := 0]
  
  max = which.max(rmat$r2)
  
  cat("\nClosest SNP r2 in ksmith: ", rmat[max, ]$rsid, "\twith r2: ", rmat[max, ]$r2)
  
  if (rmat[max, ]$r2 > 0.6) {
    smi_snp = rmat[max, ]$rsid
    smi_clus = ksmith[rsID == smi_snp, !c("VAR_ID_hg19", "rsID", "locus", "Chromosome")]
    smi_rclus = data.table(cluster = names(smi_clus), weight = as.vector(t(smi_clus)))
    max_clus = which.max(smi_rclus$weight)
    
    if (smi_rclus[max_clus, ]$weight > 0.75) {
      
      cat("\nThis SNP is attributed to cluster (using weight>0.75):", smi_rclus[max_clus, ]$cluster)
      
      cluster[rsid == snp, `:=`(ksmith_snp = smi_snp,
                                ksmith_r2 = as.numeric(rmat[max, ]$r2), 
                                ksmith_cluster = smi_rclus[max_clus, ]$cluster)]
      
    } else {
      cat("\nUnclustered SNP...")
    }
    
  }
}

fwrite(cluster, file = "./ST5_wksmith_cluster_r2_0.6.tsv", sep = "\t", quote = F)