Adapted from https://github.com/LieberInstitute/spatial_hpc/tree/main/code/enrichment_analysis/project_all &
https://github.com/LieberInstitute/spatialdACC/tree/main/code/17_LDSC

##############################################
step 1: compute cell type specificity score, create "bedfiles" 
##############################################

`01_prep_marker_bed.R`

##############################################
# step 3: run LDSC
##############################################

1. jobs for generating annotation files
perl ldsc_anno_jobs.pl

perl simpleArray.pl -j ldsc_anno_jobs.txt -n 20 -o ldsc_anno_jobs.pbs
change "module load conda_R/4.0" to "source activate ldsc" in above job submission file

sbatch ldsc_anno_jobs.pbs

2. jobs for generating ldscore files
perl ldsc_score_jobs.pl

```
python /dcs04/lieber/shared/statsgen/LDSC/base/scripts/ldsc.py \
    # Estimate l2
	 --l2 --thin-annot \
	 
	 # Specify the window size to be used for estimating LD Scores in units of centiMorgans (cM)
	 --ld-wind-cm 1 \
	 
	 # Prefix for Plink .bed/.bim/.fam file
	 --bfile /dcs04/lieber/shared/statsgen/LDSC/base/referencefiles/1000G_EUR_Phase3_plink/1000G.EUR.QC.${chr} \
	 
	 # Filename prefix for annotation file for partitioned LD Score estimation
	 --anno ${OUT_DIR}/${cluster}/chr.${chr}.annot.gz \
	 
	 # output
	 --out ${OUT_DIR}/${cluster}/chr.${chr} \
	 
	 # only print LD Scores for the SNPs listed
	 --print-snps /dcs04/lieber/shared/statsgen/LDSC/base/referencefiles/hapmap3_snps/hm.${chr}.snp

```



perl simpleArray.pl -j ldsc_score_jobs.txt -n 20 -o ldsc_score_jobs.pbs
change "module load conda_R/4.0" to "source activate ldsc" in above job submission file

sbatch ldsc_score_jobs.pbs

3. jobs for generating h2 files
perl ldsc_h2_jobs.pl

```
	python /dcs04/lieber/shared/statsgen/LDSC/base/scripts/ldsc.py \
	
	    # Filename for a .sumstats[.gz] file for one-phenotype LD Score regression
		--h2 /dcs04/lieber/shared/statsgen/LDSC/base/gwas_brain/${gwas}.gz \
		
		#Filename prefix for file with LD Scores with sum r^2 taken over SNPs included
		--w-ld-chr /dcs04/lieber/shared/statsgen/LDSC/base/referencefiles/1000G_Phase3_weights_hm3_no_MHC/weights.hm3_noMHC. \
		
		# Name of a file that has a list of file name prefixes for cell-type-specific analysis
		--ref-ld-chr ${OUT_DIR}/${cluster}/chr.,/dcs04/lieber/shared/statsgen/LDSC/base/baseline2/baselineLD. \
		
		# partitioned LD Scores were generates using an annot matrix with overlapping categories
		--overlap-annot \
		
		# frqfile files split over chromosome
		--frqfile-chr /dcs04/lieber/shared/statsgen/LDSC/base/referencefiles/1000G_Phase3_frq/1000G.EUR.QC. \
		
		# Output filename prefix
		--out ${OUT_DIR}/${cluster}/${gwas}.out \
		
		# when categories are overlapping, print coefficients as well as heritabilities
		--print-coefficients
		
```

perl simpleArray.pl -j ldsc_h2_jobs.txt -n 20 -o ldsc_h2_jobs.pbs




4. results collection
perl ldsc_results.pl
Rscript ldsc_results2.R