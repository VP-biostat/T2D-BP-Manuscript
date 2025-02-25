library(data.table)
library(coloc)
library(tibble)

# Allele alignment
allele_alignment <- function(eqtl_gwas, gwas_t2d) {
  if (gwas_t2d) {
    beta_col <- "Beta"
    ea_col <- "EA"
    nea_col <- "NEA"
  } else {
    beta_col <- "beta"
    ea_col <- "allele_effect"
    nea_col <- "allele_ref"
  }
  
  eqtl_EA <- gsub(".*_.*_.*_(.+)_b38", "\\1", eqtl_gwas$variant_id)
  eqtl_NEA <- gsub(".*_.*_(.+)_.*_b38", "\\1", eqtl_gwas$variant_id)
  
  allele_pb <- !((eqtl_EA == eqtl_gwas[, get(ea_col)]) & (eqtl_NEA == eqtl_gwas[, get(nea_col)]))
  eqtl_gwas[allele_pb, c(beta_col, ea_col, nea_col) := .(-get(beta_col), get(nea_col), get(ea_col))]
  still_allele_pb <- !((eqtl_EA == eqtl_gwas[, get(ea_col)]) & (eqtl_NEA == eqtl_gwas[, get(nea_col)]))
  
  return(eqtl_gwas[!still_allele_pb, ])
}

# Runing coloc.abf per gene function
run_coloc_per_gene <- function(gene_eqtl_gwas, gwas_t2d) {
  #check the values in the df
  #remove the impossible values here:
  gene_eqtl_gwas <- gene_eqtl_gwas[slope != Inf,] #slope with inf
  gene_eqtl_gwas <- gene_eqtl_gwas[slope != -Inf,] #slope with -inf
  gene_eqtl_gwas <- gene_eqtl_gwas[pval_nominal < 1,] #slope with inf, and pval < 1
  
  #if we have less than 30 signals, skipping
  if (nrow(gene_eqtl_gwas) < 30) {
    
    return(NULL)
    
  } else {
    
    # FIRST we are dividing the merge table into df1 and df2
    df1 <- gene_eqtl_gwas[, c("chrpos_b38", "slope", "slope_se", "maf", "Pos_b38.x", "pval_nominal")]
    df1$maf <- as.numeric(df1$maf)
    df1$pval_nominal <- as.numeric(df1$pval_nominal)
    df1 <- df1[!duplicated(df1$chrpos_b38),]
    ldf1 <- list(
      type = "quant", 
      beta = df1$slope, 
      varbeta = df1$slope_se^2,
      pvalues = df1$pval_nominal,
      N = 838, 
      MAF = df1$maf, 
      snp = df1$chrpos_b38, 
      sdY = 1 #external information
    )
    
    # check_dataset(df1)
    
    # Make the gwas dataset2 for coloc
    if (gwas_t2d) {
      df2 <- gene_eqtl_gwas[, c("chrpos_b38", "Beta", "SE", "Pos_b38.y", "Neff", "EAF", "Pvalue")]
      #modify EAF into MAF
      df2[EAF > 0.5, EAF := 1 - EAF]
    } else {
      df2 <- gene_eqtl_gwas[, c("chrpos_b38", "beta", "SE", "Pos_b38.y", "N", "MAF", "p")]
    }
    names(df2) <- c("snp", "beta", "SE", "position", "N", "MAF", "pvalues")
    df2$MAF <- as.numeric(df2$MAF)
    df2$pvalues <- as.numeric(df2$pvalues)
    df2 <- df2[!duplicated(df2$snp),]
    
    # If MAF is 0 or 1, replace it with very low numbers
    df2$MAF[which(df2$MAF == 0)] <- 1e-12
    df2$MAF[which(df2$MAF == 1)] <- 0.9999999
    
    if (gwas_t2d) {
      ldf2 <- list(
        type = "cc", 
        beta = df2$beta, 
        varbeta = df2$SE^2,
        pvalues = df2$pvalues, 
        N = 425802+33447, 
        s = 33447/(425802+33447),
        MAF = df2$MAF, 
        snp = df2$snp
      )
    } else {
      ldf2 <- list(
        type = "quant", 
        beta = df2$beta, 
        varbeta = df2$SE^2,
        pvalues = df2$pvalues, 
        N = 459247, 
        MAF = df2$MAF, 
        snp = df2$snp, 
        sdY = 1
      )
    }
    
    # check_dataset(df2)
    
    res <- coloc.abf(ldf1, ldf2)
    return(res$summary)
    
  }
}

# Look at the input that would be the gwas/eqtl pairs we are looking into, 
#   EXAMPLE: if we have 9 eQTLs and 4 GWAS, we need to run 1-36 job (9*4)
args = commandArgs(trailingOnly=TRUE)
if (length(args)==0) {
  warning("\nRequired arguments:
       \n1. a numeric from 1 to 4\n", call.=FALSE)
} else {
  job <- as.numeric(args[1]) #looping on gwas and eqtl
}

out <- "YOUR_OUTPUT_NAME"
dir.create(paste0("./colocalization/", out))

if (exists("result_table")) rm(result_table)

## Fourth attempt of a colocalization analysis using coloc package, cis eqtl from GTEx v8 (n = 9 tissues) allpairs, and T2D and BP GWASes

# this one is based on the only script qsub_merging_for_colocalization.sh, ran multiple time (49) for each eQTL
# we aim to target, for each SNP in the clustering, the Loci around and available in eQTL, for each loci, 
# we add the top SNPs associated with what we found in GWAS (T2D and BP) and store it in a big table including
# both GWAS, eQTL data and assigned to SNPs of the cluster 

# List of gwases of interest: 
gwases <- c("t2d", "sbp", "dbp", "pp")

# List of tissues of interest: 
parquet_folder <- "/rds/general/project/eph-prokopenko-lab-silk/live/data/GTEx_v8_allpair/GTEx_Analysis_v8_EUR_eQTL_all_associations/"
int_tissues <- unique(gsub(".v8.*", "", list.files(path = parquet_folder)))

# Select the gwas to do in this job, and the tissue 
all_assoc <- expand.grid(gwases, int_tissues)
# Select the eQTL and GWAS pair
int_tissue <- as.character(all_assoc[job, 2])
gwas <- as.character(all_assoc[job, 1])

# Load all the files related to this gwas + tissues
merged_files <- list.files(path = "./colocalization/merged_gwas_eqtl/", pattern = int_tissue, full.names = T)
n_merged_files <- length(merged_files)

# Window in COLOC? https://github.com/chr1swallace/coloc/blob/main/FAQ.md#if-i-understand-correctly-colocabf-can-be-run-with-correlated-variants-that-is-no-prerequisite-for-taking-through-ld-pruningclumping-is-required-am-i-correct-in-my-understanding-
# Coloc is designed to address whether two traits share causal variant(s) in a genomic region. 
# It leaves the definition of "region" up to the user. You need to break the genome into smaller regions, 
# within which it is reasonable to assume there is at most one (coloc.abf) or a small number (coloc.signals)
# of causal variants per trait. One way to do this is to use the boundaries defined by recombination hotspots, 
# proxied by this map created by lddetect.
# 
# How big should a region be? Big enough that all variants in LD with a lead SNP are included; 
# small enough that only one or a small number of causal variants might exist in it. 
# I have found using densely genotyped studies and the lddetect boundaries above that 
# regions typically contain 1,000-10,000 SNPs.

for (k in 1:n_merged_files) {
  
  merged_file <- merged_files[k]
  #load the 4 tables with GWAS-1eQTL
  load(merged_file)
  
  #find in the file name the top_SNP name
  top_snp <- gsub(".*chr(.+)\\.Rdata", "chr\\1", merged_file)
  
  cat("-----\n")
  cat("Doing coloc for", top_snp, "(", k, "/", n_merged_files, ") in", gwas, "GWAS and", int_tissue, "eQTL\n\n")
  
  #select only the merged value using 
  merged_dt <- get(paste0("temp_merge_eqtl_", gwas))
  
  #remove potential NA in the allpairs files (slope + slope_se) 
  merged_dt <- merged_dt[!is.na(slope), ]
  merged_dt <- merged_dt[!is.na(slope_se), ]
  
  #ALLELE alignment: check if eqtl and gwas snps have the same allele
  merged_dt <- allele_alignment(merged_dt, (gwas == "t2d"))
  
  #then split by gene
  #the naming has also transcript (.##), that we are getting rid of
  #EDIT: removed this line to keep the transcript info
  # merged_dt$phenotype_id <- gsub("\\..*", "", merged_dt$phenotype_id)
  #and we split by gene
  genes_eqtl_gwas <- split(merged_dt, by = "phenotype_id")
  
  #apply the coloc per gene
  results <- as.data.frame(do.call(rbind, lapply(genes_eqtl_gwas, run_coloc_per_gene, (gwas == "t2d"))))
  #change the rownames (it is by default the genes, but might overlap later)
  results <- rownames_to_column(results, var = "phenotype_id")
  #add the gwas, tissue and snp targeted in this analysis 
  n_res <- nrow(results)
  
  if(n_res > 0) {
    results$top_snp <- rep(top_snp, n_res)
    results$gwas <- rep(gwas, n_res)
    results$tissue <- rep(int_tissue, n_res)
    
    if (!exists("result_table")) {
      result_table <- results
    } else {
      result_table <- rbind(result_table, results)
    }
  }
  
  # There is colocalisation if the posterior probability of the model sharing a single causal variant was larger than 80%
  if (nrow(result_table) > 0) {
    is_coloc <- result_table$PP.H4.abf>0.8
    
    if (sum(is_coloc)>0) {
  
      # subset(res$results,SNP.PP.H4>0.01)
      found_message <- paste(" - Colocalization found between", top_snp, ",", gwas, ",", int_tissue, "---\n")
      found_file <- paste0("./colocalization/", out, "/", "association_found.txt")
      cat(found_message)
      if (!file.exists(found_file)) file.create(found_file)
      cat(found_message, file = paste0("./colocalization/", out, "/", "association_found.txt"), append = T)
    }
  }
  
  
}



fwrite(result_table, file = paste0("./colocalization/", out, "/", gwas, "_", int_tissue, "_coloc_cluster_snp_result_table.csv"), sep = ";", quote = F, row.names = F)