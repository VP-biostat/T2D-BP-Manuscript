#!/bin/bash
#PBS -l walltime=72:00:00
#PBS -l select=1:ncpus=4:mem=64gb

module load anaconda3/personal
source activate ldsc

#for icbp.pp.eur:
~/ldsc/munge_sumstats.py \
--out ./ldsc/sumstats/icbp.pp.eur \
--snp SNP \
--a1 allele1 \
--a2 allele2 \
--N-col totalsamplesize \
--sumstats ./ICBP_PP_02082017.txt.gz \
--ignore markername 

#for icbp.sbp.eur
~/ldsc/munge_sumstats.py \
--out ./ldsc/sumstats/icbp.sbp.eur \
--snp SNP \
--a1 allele1 \
--a2 allele2 \
--N-col totalsamplesize \
--sumstats ./ICBP_SBP_02082017.txt.gz \
--ignore markername

#for icbp.dbp.eur:
~/ldsc/munge_sumstats.py \
--out ./ldsc/sumstats/icbp.dbp.eur \
--snp SNP \
--a1 allele1 \
--a2 allele2 \
--N-col totalsamplesize \
--sumstats ./ICBP_DBP_02082017.txt.gz \
--ignore markername

#for mahajan.t2d.eur:
~/ldsc/munge_sumstats.py \
--out  ./ldsc/sumstats/mahajan.t2d.eur \
--N 1114458.0 \
--sumstats ./Mahajan.NatGenet2018b.T2D.European.txt \
--a1 NEA \
--a2 EA 