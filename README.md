# T2D-BP Manuscript

## Overview
This repository contains the data and analysis scripts supporting our manuscript on the genetic relationship between type 2 diabetes (T2D) and blood pressure (BP). The study uses genome-wide association studies (GWAS), partitioned polygenic scores (PGS), colocalization analyses, and enrichment in cis-regulatory elements (cCREs) to identify genetic drivers underlying T2D-BP comorbidities.

## Repository Structure
The repository is organised into four main folders, each corresponding to a key section of the results:

### 1. **Genetic overlap between T2D and BP** (`01_T2D_BP_shared_genetics/`)
   - Evaluation of the genetic correlation between T2D and Systolic Blood Pressure (SBP)/Diastolic Blood Pressure (DBP)/Pulse Pressure (PP) using [_ldsc_](https://github.com/bulik/ldsc)
   - Build classic PGS for T2D, SBP, DBP, and PP using robust and reduced SNV lists to ensure we are limitating their pleiotropic effect on the disease of interest
   - Use [_comorbidPGS_](https://github.com/VP-biostat/comorbidPGS) to assess whether the genetic predisposition towards one condition can predict the risk of the other

### 2. **Clusters of pathogenetic processes** (`02_clustering/`)
   - Define five pathogenetic clusters using hierarchical clustering based on 39 GWAS: 
     - *Inverse T2D-BP Risk*
     - *Metabolic Syndrome*
     - *Higher Adiposity*
     - *Vascular Dysfunction*
     - *Reduced Beta-Cell Function*
   - Replicate our results using [_MRClust_](https://github.com/cnfoley/mrclust)
   - Replicate our results using [_bNMF clustering_](https://github.com/gwas-partitioning/bnmf-clustering/)
   - Validate the clusters comparing with the latest [T2D hierarchical clustering published](https://doi.org/10.1038/s41586-024-07019-6) and [T2D bNMF clustering](https://doi.org/10.1038/s41591-024-02865-3) respectively

### 3. **Multiomic characterisation of T2D-BP clusters** (`03_cluster_characteristics/`)
   - Evaluate the changes in gene expression associated with the clustered SNVs using eQTLs and [_coloc_](https://github.com/chr1swallace/coloc)
   - Examine enrichment of cluster-specific SNVs in cCREs across relevant tissues
   - Provides functional annotation insights.

### 4. **T2D-BP comorbidity using partitioned PGSs** (`04_partitioned_PGS/`)
   - Build partitioned PGS based on the five clusters
   - Investigate the associations between partitioned PGS and clinical outcomes and complications within the UK Biobank using [_comorbidPGS_](https://github.com/VP-biostat/comorbidPGS)
   - Make cumulative hazard trends for T2D-BP complications on individuals stratified by partitioned PGS using Cox proportional hazards models from [_survival_](https://CRAN.R-project.org/package=survival)

## Study Design and Analysis Workflow
![Study design and Analysis Workflow](https://raw.githubusercontent.com/VP-biostat/T2D-BP-Manuscript/main/figures/Figure1.png)

## Contact
For inquiries, please contact the corresponding author(s) as indicated in the manuscript.
