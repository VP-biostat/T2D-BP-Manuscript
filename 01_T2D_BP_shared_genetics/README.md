## Genetic correlation between T2D and BP

1. Load your environment with 
   - [_ldsc_](https://github.com/bulik/ldsc) v1.0.1 along with the baselineLD LDscores
   - [_mungesumstat_](https://github.com/Al-Murphy/MungeSumstats) 
   - [_plink_](https://www.cog-genomics.org/plink2/) v1.9

2. Format T2D and BP GWASs to have `chr:pos` format using the R script `ldsc_chrpos_formatting.R`

3. Run mungesumstat on the 4 GWAS using the unix script `munge_sumstat_ldsc.sh` *BUT do not use --merge-alleles (it does not work)*

4. Run an R script that merges alleles by hand: `ldsc_merge-alleles_instead_of_mungesumstat.R`

5. Perform genetic correlation with _ldsc_ using a command such as: 
```
[ldsc path]/ldsc/ldsc.py \
--out [Output name].rg \
--rg [T2D sumstats path],[SBP sumstats path],[DBP sumstats path],[PP sumstats path] \
--w-ld-chr [ldsc path]/ldsc/baselineLD/LDscore. \
--ref-ld-chr [ldsc path]/ldsc/baselineLD/LDscore.
```
