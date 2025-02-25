library(data.table)

#changing baseline ldscores to add SNP as chrpos

files <- list.files(path = "~/ldsc/1000G_Phase3_ldscores/LDscore/", pattern = "l2.ldscore.gz")

for (file in files) {
  
  df <- fread(paste0("~/ldsc/1000G_Phase3_ldscores/LDscore/",file))
  
  #store rsid in rsid column
  # df$rsid <- df$SNP
  #create a chrpos column that will be used, to note it is build 37
  df$SNP <- paste0(df$CHR,":",df$BP)
  # names(df)[4] <- "L2"
  
  fwrite(df, file = paste0("./ldsc/baselineLD/", file), sep = "\t", quote = F)
}


#changing BP GWAS SNP as chrpos
sbp <- fread("./data/ICBP_SBP_02082017.txt.gz")
dbp <- fread("./data/ICBP_DBP_02082017.txt.gz")
pp <- fread("./ICBP_PP_02082017.txt.gz")

sbp$SNP <- gsub(":SNP", "", sbp$markername)
dbp$SNP <- gsub(":SNP", "", dbp$markername)
pp$SNP <- gsub(":SNP", "", pp$markername)

fwrite(sbp, file = "./ldsc/sumstats/ICBP_SBP_02082017.txt.gz", sep = "\t", quote = F)
fwrite(dbp, file = "./ldsc/sumstats/ICBP_DBP_02082017.txt.gz", sep = "\t", quote = F)
fwrite(pp, file = "./ldsc/sumstats/ICBP_PP_02082017.txt.gz", sep = "\t", quote = F)
