library(data.table)
library(comorbidPRS)

df <- fread("./UKB_phenotypes_and_T2D_BP_PGS.tsv")

#list the phenotype of interest in the df 
phens <- c("T2D.all","dbp","pp","sbp","hypertension")
#same for prs dist (identfied by the pattern "Zscore"
prs <- grep("Zscore", names(df), value = T)
#expand the associations
assoc_table <- expand.grid(prs, phens)
#same for covar
covs <- c("sex","array","bmi","age","PC1","PC2","PC3","PC4","PC5","PC6")

#run comorbidPRS
results <- multiassoc(df, assoc_table, covar_col = covs)

write.table(results, file = "./bp_t2d_comorbisPGS_results.tsv", sep = "\t", quote = F, row.names = F, col.names = T)
