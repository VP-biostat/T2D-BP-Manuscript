library(data.table)
library(comorbidPGS)

# Putting the step (-J 1-60) TO PARALLELISE the R script
args = commandArgs(trailingOnly=TRUE)
if (length(args)==0) {
  stop("\nRequired arguments:
       \n1. a numeric from 1 to 60\n", call.=FALSE)
} else {
  job <- as.numeric(args[1]) #looping on gwas
}
step = 783 #to adjust with your data: based on your cluster and the total analysis to do 
min_ind = (job-1)*step+1
max_ind = job*step

# Taking all phen, pgs, covar
diseases = fread("./MEGAPHEWAS_icd10_list.tsv") #The icd10 available within the UK Biobank and their description 
icd_col = diseases$icd10 
pgs_col = c("1.z.T2D", "2.z.T2D", "3.z.T2D", "4.z.T2D", "5.z.T2D")
assoc_table = expand.grid(pgs_col, icd_col)

# Subtaking the phen, pgs and covar for current job
sub_assoc_table = assoc_table[min_ind:max_ind, ]
covar_col = c("age", "sex", "array", "PC1", "PC2", "PC3", "PC4", "PC5", "PC6")

cur_pgs_col = as.character(unique(sub_assoc_table$Var1))
cur_phen_col = as.character(unique(sub_assoc_table$Var2))
all_col = c(covar_col, cur_pgs_col, cur_phen_col)

df = fread("./MEGAPHEWAS_bp_t2d_phenofile.tsv.gz", select = all_col) #Pheno file + partitioned PGS extracted from UKB

# Find columns with only one unique value
cols_with_one_value <- sapply(df[, ..cur_phen_col], function(x) length(unique(x)) == 1)
# Get the names of those columns
sub_phen_col <- cur_phen_col[!cols_with_one_value]

sub_cur_assoc_table = sub_assoc_table[sub_assoc_table$Var2 %in% sub_phen_col, ]

res = multiassoc(df, sub_cur_assoc_table, scale = T, covar_col, parallel = T, num_cores = 6)

fwrite(res, file = paste0("./pgs_after_clustering/PHEWAS/pgs_association_PHEWAS_", min_ind, "_to_", max_ind, ".tsv"), sep = "\t", quote = F, na = "")
