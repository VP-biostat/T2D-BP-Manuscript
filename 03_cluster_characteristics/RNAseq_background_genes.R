library(data.table)

path = "./RNAseq/YOUR_FILE" #replace with RNAseq data (from GTEx in this study)
tissue = gsub(".*v8_(.+).gct.gz", "\\1", path)
RNAseq = fread(path)

ind_col = grep(names(RNAseq), pattern = "GTEX", value = T)

RNAseq[, med_TPM := apply(.SD, 1, median), .SDcols = ind_col]

subdf = RNAseq[med_TPM > 120, ]

write(gsub("\\..*", "", subdf$Name), file = paste0("./RNAseq/", length(gsub("\\..*", "", subdf$Name)), "_background_genes_in_", tissue, ".txt"))

