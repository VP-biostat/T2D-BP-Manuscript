library(data.table)

#ref gwas is the less powerful 
pp <- fread("./ldsc/sumstats/icbp.pp.eur.sumstats.gz")

#the other sumstats 
sbp <- fread("./ldsc/sumstats/icbp.sbp.eur.sumstats.gz")
dbp <- fread("./ldsc/sumstats/icbp.dbp.eur.sumstats.gz")
t2d <- fread("./ldsc/sumstats/mahajan.t2d.eur.sumstats.gz")

#replace sbp, t2d, dbp by the good one, line to edit have a #*
merged <- merge(pp, sbp, by = "SNP") #*
allele_pb <- !(merged$A1.x==merged$A1.y & merged$A2.x==merged$A2.y)

if (sum(allele_pb) == 0 ) {
  cat("No merge-alleles required!")
  sub_merged <- merged[, c("SNP","A1.y","A2.y","Z.y","N.y")]
  names(sub_merged) <- c("SNP", "A1", "A2", "Z", "N")
  
  fwrite(sub_merged, file = "./ldsc/sumstats/icbp.sbp.eur.merge-alleles.sumstats.gz", sep = "\t", quote = F)
  
} else if (sum(allele_pb) == nrow(merged)) {
  cat("systematic inversion of A1 and A2, swapping all of them in second gwas!")
  sub_merged <- merged[, c("SNP","A2.y","A1.y","Z.y","N.y")]
  names(sub_merged) <- c("SNP", "A1", "A2", "Z", "N")
  sub_merged$Z <- -sub_merged$Z 
  
  #then check again the problems
  merged <- merge(pp, sub_merged, by = "SNP") #*
  allele_pb <- !(merged$A1.x==merged$A1.y & merged$A2.x==merged$A2.y)
}

if (sum(allele_pb) > 0) {
  cat("still some issues, there are", sum(allele_pb), "SNPs with multialleles, skipping them")
  
  sub_merged <- merged[!allele_pb, c("SNP","A1.y","A2.y","Z.y","N.y")]
  names(sub_merged) <- c("SNP", "A1", "A2", "Z", "N")
  
  fwrite(sub_merged, file = "./ldsc/sumstats/mahajan.t2d.eur.merge-alleles.sumstats.gz", sep = "\t", quote = F)
}
