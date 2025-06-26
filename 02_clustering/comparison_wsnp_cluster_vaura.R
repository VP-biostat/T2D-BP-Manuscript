library(readxl)
library(LDlinkR)
library(data.table)

token <- 'your LDlink token'

#use the manuscript supplementary tables 
cluster <- data.table(read_excel("./BP_T2D_Supplementary_Table_240816.xlsx", sheet = 6, skip =1))
cluster[, Vaura_snp := "none"]
cluster[, Vaura_r2 := 0]
cluster[, Vaura_cluster := "unassigned"]

#use vqurq et al. supplementary table
vaura <- data.table(read_excel("./vaura_et_al_bp_clustering.xlsx", sheet = 1))
vaura <- vaura[!is.na(Group),] #remove NA
vaura[, `:=`(Chromosome = nafill(Chr, type = "locf"))] #fill the chromosome col
vaura[, rsid := gsub(" \\(.*", "", `ID (Locus)`)] #format the rsid col
vaura[, rsid := gsub("(.+)\\:(.+)\\_.*_.*", "chr\\1\\:\\2", rsid)] #format the rsid col

for (i in 1:nrow(cluster)) {
  if (exists("mat")) rm(mat)
  
  snp = cluster[i, ]$rsid
  chr = as.numeric(gsub("\\:.*", "", cluster[i, ]$`chr:pos (b37)`))
  
  perc_loading = round(i / nrow(cluster)*10)
  perc_remaining = 10-perc_loading
  loading_bar = paste0(strrep("=", perc_loading), strrep("-", perc_remaining))
  cat("\n|",loading_bar , "|\nSNP: ", snp, "\tin chromosome: ", chr, "\tin cluster: ", cluster[i, ]$`Attributed cluster`)
  
  tempsnps = c(snp, vaura[Chromosome == chr, ]$rsid)
  
  if (length(tempsnps) < 2) {
    cat("\nSkipping here, no SNP found in", chr)
    next
  }
  
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
  
  cat("\nClosest SNP r2 in vaura: ", rmat[max, ]$rsid, "\twith r2: ", rmat[max, ]$r2)
  
  if (rmat[max, ]$r2 > 0.6) {
    vau_snp = rmat[max, ]$rsid
    #select the variant(s) in vaura et al
    vaura_temp = vaura[rsid == vau_snp, ]
    if (nrow(vaura_temp) == 0) {
      tempsnps_collapsed = paste(tempsnps, collapse = ',')
      cat("\n\tPROBLEM: the best SNP is in a list but unable to find it, use this ", vau_snp, "to look at this list", tempsnps_collapsed, "to find the actual cluster")
      
      cluster[rsid == snp, `:=`(vaura_snp = vau_snp,
                                vaura_r2 = as.numeric(rmat[max, ]$r2), 
                                vaura_cluster = tempsnps_collapsed)]
    } else {
      max_clus = which.max(vaura_temp$`bNMF Weight`)
      
      #if there is more than 1 variant, select the max cluster
      if (vaura_temp[max_clus, ]$`bNMF Weight` > 0.75) {
        
        cat("\nThis SNP is attributed to cluster (using weight>0.75):", vaura_temp[max_clus, ]$Group)
        
        cluster[rsid == snp, `:=`(vaura_snp = vau_snp,
                                  vaura_r2 = as.numeric(rmat[max, ]$r2), 
                                  vaura_cluster = vaura_temp[max_clus, ]$Group)]
        
      } else {
        cat("\nUnclustered SNP...")
      }
      
    }
    
  }
}

fwrite(cluster, file = "./ST5_wVaura_cluster_r2_0.6.tsv", sep = "\t", quote = F)
