library(data.table)

df <- fread("./pgs_after_clustering/bp_t2d_clustering_new_phenos2_5clusters_ukb_imp_6PCs_unweighted_df.tsv") #Pheno file + partitioned PGS extracted from UKB
df[, comorbidity := (`T2D.all` == 1 & hypertension == T)]
df[, `1_over0.9` := (`1.z.T2D` > quantile(`1.z.T2D`, 0.9))]
df[, `1test_over0.9` := (`1.z.T2D.test` > quantile(`1.z.T2D.test`, 0.9))]
df[, `1testnodoubledosage_over0.9` := (`1.z.T2D.test.nodoubledosage` > quantile(`1.z.T2D.test.nodoubledosage`, 0.9))]
df[, `1nodoubledosage_over0.9` := (`1.z.T2D.nodoubledosage` > quantile(`1.z.T2D.nodoubledosage`, 0.9))]
df[, `2_over0.9` := (`2.z.T2D` > quantile(`2.z.T2D`, 0.9))]
df[, `3_over0.9` := (`3.z.T2D` > quantile(`3.z.T2D`, 0.9))]
df[, `4_over0.9` := (`4.z.T2D` > quantile(`4.z.T2D`, 0.9))]
df[, `5_over0.9` := (`5.z.T2D` > quantile(`5.z.T2D`, 0.9))]

prs = c(grep(pattern = ".z.", names(df), value = T), 
        "2.z.T2DAND3.z.T2D",
        "2.z.T2DAND4.z.T2D",
        "2.z.T2DAND5.z.T2D",
        "2.z.T2DAND3.z.T2DAND5.z.T2D", 
        "overall")
phenotype = c("T2D.all", "hypertension")

povt2d = nrow(df[T2D.all == 1,])/nrow(df)
povhyp = nrow(df[hypertension == T,])/nrow(df)

povcom = nrow(df[comorbidity == T,])/nrow(df)

if (exists("ptable")) rm(ptable)

for (pr in prs) {
  
  #probability and relative risk 
  if (pr == "overall") {
    subdf = df
  } else if (grepl("AND", pr)) {
    prbis = gsub("T2D", "", pr)
    prcl = str_extract_all(prbis, "\\d+")
    prcl = unlist(prcl)

    subdf = df
    for (cl in prcl) {
      subdf = subdf[get(paste0(cl, "_over0.9")) == T, ]
    }
    
  } else {
    subdf = df[get(paste0(str_extract(pr, "\\d+"), "_over0.9")) == T, ]
  }
    
  
  pt2d = nrow(subdf[T2D.all == 1,])/nrow(subdf)
  phyp = nrow(subdf[hypertension == T,])/nrow(subdf)
  
  pcom = nrow(subdf[comorbidity == T,])/nrow(subdf)
  
  RRt2d = pt2d/povt2d
  RRhyp = phyp/povhyp
  RRcom = pcom/povcom
  
  temp = data.table(pgs = pr, 
                    n = nrow(subdf), 
                    pt2d = pt2d, 
                    phyp = phyp, 
                    pcom = pcom, 
                    RRt2d = RRt2d,
                    RRhyp = RRhyp,
                    RRcom = RRcom)
  
  if (!exists("ptable")) {
    ptable <- temp
  } else {
    ptable <- rbind(ptable, temp)
  }
    
}

fwrite(ptable, file = "./pgs_after_clustering/probability_of_comorbidity_per_cluster_0.9_quantile.csv", sep = ";", quote = F)
