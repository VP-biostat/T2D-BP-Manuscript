## Multiomic characterisation of T2D-BP clusters

### Changes in gene expression associated with the clustered SNVs

Protocol for colocalization analysis:

1. Extract from the GTEx multitissue eQTLs (n = 48) all pairs. The script `GTEx_all_pairs_parquet_reader.R` was used to extract and read them per tissue

2. lift-over from build 37 to build 38 the four GWAS (T2D, SBP, DBP and PP) using [_liftOver_](https://genome-store.ucsc.edu/) setup in your environment. The chain file `hg19ToHg38.over.chain` should be stored in the same directory as _liftOver_. 

3. Extract and merge GWAS and eQTLs information around each clustered SNV using `step1_qsub_merging_for_colocalization_allpairs.R`. It is advised to parallelised this step.

4. Use `step2_colocalization_coloc_by_gene.R` to perform coloc.abg per clustered SNV region. It is advised to parallelised this step (info in comments).

5. Several representation of the colocalization results are provided with:
   - `step3.1_colocalization_analysis.R` which gives a summary table per clustered SNV per GWAS and per tissue-specific eQTL and the heatmap from SF6
   - `step3.2_colocalization_visualiser.R` provides a cool visualisation like a Manhattan plot per cluster per tissue-specific eQTL
   - `step3.3_fig3a.R` gives the plots (+additional ones) from Fig. 3a
   
6. Compare with previous T2D loci (from [Suzuki et al. 2024](https://doi.org/10.1038/s41586-024-07019-6) and [K. Smith et al. 2024](https://doi.org/10.1038/s41591-024-02865-3) using `step4_colocalization_allg_colocalization.R`

For further investigation we used:
   - `top5_tissues_in_fig3a_and_b.R` identifies for a given tissues and associated cell type the common clustered SNVs.
   - `RNAseq_background_genes.R` extracts from RNAseq data from GTEx the top genes per tissues (with TPM>120) to use as background genes in _metascape_.

### Changes in regulatory elements associated with the clustered SNVs

scATAC-seq can be downloaded from [CATLAS](http://catlas.org/humanenhancer/) and analysed using `scATAC-seq_clustering_analysis.R`
