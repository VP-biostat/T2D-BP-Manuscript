## T2D-BP comorbidity using partitioned PGSs

This section relies on the UK Biobank individual-level datasets which are not publicly available. Unweighted PGS were done per cluster, giving the partitioned PGS. Several R scripts to build the plots are provided.

### Comorbidity Relative Risk 

Obtain the relative risk of T2D-BP comorbidity per being in the top 10% distribution of partitioned PGS: `pgs_ptable_after_clustering.R`


### Comorbidity survival curve 

Perform survival analysis with Cox proportional hazards model per being in the top 33% distribution of partitioned PGS: `pgs_after_cluster_survcurv_wcovar.R` 


### PheWAS per partitioned PGS 

1. Run PheWAS across all ICD10 codes available in the UKB per partitioned PGS using [_comorbidPGS_](https://github.com/VP-biostat/comorbidPGS) with: `pgs_after_clustering_phewas.R`
It is advised to parallelised this script, specially if run in the UK Biobank

2. Visualise PheWAS per partitioned PGS with: `pgs_after_clustering_phewas_visualizer.R`
