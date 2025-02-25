library(data.table)
library(readxl)

st = data.table(read_excel("BP_T2D_Supplementary_Table_240816.xlsx", sheet = 8, skip = 1))

#select tissue
tissue = "Adipose Subcutaneous"
# tissue = "Artery Tibial"
# tissue = "Nerve Tibial"
# tissue = "Skin Sun Exposed Lower leg"
# tissue = "Thyroid"

#select cell type
cell = "Adipocyte"
# cell = "Vasc Sm Muscle 1"
# cell = "Schwann General"
# cell = "Keratinocyte 1"
# cell = "Follicular"

ccre = fread("./scATAC-seq_zhang/scATAQ-seq_zhang_bp_t2d_nold0.6_thres250_new_phenos2_clustering_enrichment_withEnsemblAnnot.csv")

thyr = st[`eQTL tissue` == tissue, ]
thyr2 = ccre[get(cell) == 1, ]

thyr3 = merge(thyr, thyr2, by.x = "SNV from clustering", by.y = "rsid")

fwrite(thyr3, paste0("./", tissue, "_tissue_fig2a_", cell, "_celltype_fig2b.csv"))
