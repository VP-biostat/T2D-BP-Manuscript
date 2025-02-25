library(readxl)
library(LDlinkR)
library(data.table)

token <- 'your LDlink token'

#use the manuscript supplementary tables
cluster <- data.table(read_excel("./BP_T2D_Supplementary_Table_240816.xlsx", sheet = 6, skip =1))
cluster[, Suzuki_snp := "none"]
cluster[, Suzuki_r2 := 0]
cluster[, Suzuki_cluster := "unassigned"]

#load suzuki et al. supplementary tables
suzuki <- data.table(read_excel("./41586_2024_7019_MOESM3_ESM.xlsx", sheet = 6, skip = 2))
suzuki <- suzuki[!is.na(`Cluster assignment`),] #remove NA
suzuki[, `:=`(Chromosome = nafill(Chromosome, type = "locf"))] #fill the chromosome col

for (i in 1:nrow(cluster)) {
  if (exists("mat")) rm(mat)
  
  snp = cluster[i, ]$rsid
  chr = as.numeric(gsub("\\:.*", "", cluster[i, ]$`chr:pos (b37)`))
  
  perc_loading = round(i / nrow(cluster)*10)
  perc_remaining = 10-perc_loading
  loading_bar = paste0(strrep("=", perc_loading), strrep("-", perc_remaining))
  cat("\n|",loading_bar , "|\nSNP: ", snp, "\tin chromosome: ", chr, "\tin cluster: ", cluster[i, ]$`Attributed cluster`)
  
  tempsnps = c(snp, suzuki[Chromosome == chr, ]$`Index SNV`)
  
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
  
  cat("\nClosest SNP r2 in Suzuki: ", rmat[max, ]$rsid, "\twith r2: ", rmat[max, ]$r2)
  
  if (rmat[max, ]$r2 > 0.6) {
    suz_snp = rmat[max, ]$rsid
    suz_clu = suzuki[`Index SNV` == suz_snp, ]$`Cluster assignment`
    
    cluster[rsid == snp, `:=`(Suzuki_snp = suz_snp,
                              Suzuki_r2 = as.numeric(rmat[max, ]$r2), 
                              Suzuki_cluster = suz_clu)]
    
    cat("\nAttributing this Suzuki snp: ", suz_snp, "\tfrom this cluster: ", suz_clu)
  }
}

fwrite(cluster, file = "./ST5_wSuzuki_cluster_r2_0.6.tsv", sep = "\t", quote = F)