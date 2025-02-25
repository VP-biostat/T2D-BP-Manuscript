## Hierarchical clustering correlation between T2D and BP

1. Load and harmonise the GWASes included in the clustering 
   - All GWAS should have the same rsID or chr:pos format
   - The first GWAS (referred to as the reference GWAS) only have Beta > 0 (flip the alleles when there is a variants with a negative effect)
   - The subsequent GWASes should have their effect alleles aligned with the reference GWAS
   - Compute for each variant in each GWAS their Z-score = Beta / SE
   - In the end of this step, have an *harmonised table with one variant per row and its associated Z-score for the given GWAS traits in the columns*
   
2. Adapt and run hierarchical clustering R script:
   - Using [_imputeSCOPA_](https://github.com/ImperialStatGen/imputeSCOPA) imputation (recommended for high-dimensional analysis): `hierarchical_clustering_imputeSCOPA_imputation.R`
   - Using [_mice_](https://github.com/amices/mice) imputation: `hierarchical_clustering_MICE_imputation.R`

3. Convert the actual hierarchical clustering into a nice and workable plot with `hierarchical_clustering_visualizer.R`

4. Compare this clustering with other methods such as [_bNMF clustering_](https://github.com/gwas-partitioning/bnmf-clustering/) with `comparison_bnmf_hieclu.R`

5. Compare the clusters with latest T2D clustering from the litterature: 
   - [T2D hierarchical clustering published](https://doi.org/10.1038/s41586-024-07019-6) with this R script ``
   - [T2D bNMF clustering](https://doi.org/10.1038/s41591-024-02865-3) with this R script ``
