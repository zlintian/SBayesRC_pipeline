This is a pipeline how I format any GWAS summary-level data into cojo format, and use SBayesRC method to generate predictors. Both the scripts to use the [R version SBayesRC](https://github.com/zhilizheng/SBayesRC) and [GCTB version SBayesRC](https://cnsgenomics.com/software/gctb/#SBayesRCTutorial) are included. 

Resource data is available for local users with the following path setting, and public for [downloading](https://cnsgenomics.com/software/gctb/#Download) too. 

```{bash, eval = F}
cd $workingpath

## this is where you placed file cojo_format_v7.R
exedir="./"

## this is the path to LD files and annotation file. 
LD_PATH1="/scratch/project_mnt/S0007/uqzzhen4/project/UKB/LD/ukb20k_7M_4cM/"
LD_PATH2="/scratch/project_mnt/S0007/uqzzhen4/project/UKB/LD/ukb20k_hm3_4cM/"
annot="/scratch/project_mnt/S0007/uqzzhen4/project/UKB/annot/annot_baseline2.2.txt" # annottion file

## you can put your GWAS file into a folder and name it with the trait
trait=MDD_01
gwas_file=PGC_UKB_depression_genome-wide.txt
ma_file=${trait}/${gwas_file}
```


# format a GWAS summary level data to cojo format

for example, if the raw gwas file has a header like this:

> MarkerName A1 A2 Freq LogOR StdErrLogOR P

you can put each element of the header name into the script per flag:

```{bash, eval = F}
cmd1="Rscript  ${exedir}/cojo_format_v7.R  \
  --file  ${trait}/${gwas_file}  \
  --out  ${trait}/${gwas_file}.ma   \
  --SNP  MarkerName  \
  --A1 A1  \
  --A2  A2 \
  --freq  Freq   \
  --pvalue P  \
  --beta  LogOR  \
  --se  StdErrLogOR   \
  --samplesize  500199   "

job_name="format_"${trait}
formatqsub=`qsubshcom  "$cmd1" 1 50G  $job_name  2:00:00  " "     `
```

This R script has several extra functions than formatting:
> 1. If dbSNP ID is missing, and chr and bp are provided. We will denote the SNP column name as "missing", and match to the dbSNP IDs using the snp.info file in LD reference.
> 2. If the alleles are in lower case, it converts them to uppercase.
> 3. Convert odds ratio to effect size if it's named OR.
> 4. Make a allele frequency comparison plot between data vs LD reference if it is provided
> 5. If allele frequency is missing, it fills it with the AF in LD reference.
> 6. If there is not a column of per-SNP-sample-size in raw data, the script can fill it with a unique number found in publication.
> 7. If the data has two columns for the sample size as Ncase and Ncontrol, the script adds them up to be N. Put them in as "Ncase,Ncontrol" behind --samplesize.
> 8. Allele frequency is also sometimes separately included for cases and controls. If that’s the case, we will do a calculation with the number of cases and controls.
> 9. If the data does not have sample size but "NCHROBS", it will be divided by 2 to be sample size.
> 10. If SE is missing, we will estimate it from BETA and P.
> 11. If beta is missing but z is provided, we will calculate effect size with z score, allele frequency and sample size. 

 
This is an exmple of comparison of allele frequency in gwas summary data vs. reference data:

<img src="GIANT_HEIGHT_YENGO_2022_GWAS_SUMMARY_STATS_EUR_AF_plot.png" width="50%" height="50%" />


# R version

The 3 key steps of QC, imputation and SBayesRC are adopted derectly from [Zhili's method](https://github.com/zhilizheng/SBayesRC).

## Tidy 

optional step, tidy summary data

```{bash, eval = F}
## "log2file=TRUE" means the messages will be redirected to a log file 
job_name="tidy_"${trait} # your job name, customize
tidyqsub=`qsubshcom "Rscript -e \"SBayesRC::tidy(mafile='${ma_file}.ma', LDdir='$LD_PATH1', output='${ma_file}_tidy.ma', log2file=TRUE) \"" 1 50G $job_name 10:00:00 "  -wait=$formatqsub  " `
```

Best practice: read the log to check issues in your GWAS summary data.  

## choose LD matrix 

We don't want to impute more than 30% SNPs in the LD matrix. If you GWAS summary stat has less than 70% of the SNPs in the 7.3M LD reference, we will switch to the HapMap3 reference. 

```{bash, eval = F}
## choose LD matrix based on number of SNPs
if [ $(wc -l  ${trait}/${gwas_file}  | awk '{print $1}'  ) -gt  5149563 ]; then   LD_PATH=$LD_PATH1 ; else  LD_PATH=$LD_PATH2 ; fi
```
For convinience, we moved this step into the beginning by checking rows in original data in the bash pipeline file. 

## Impute 

optional step if your summary data doesn't cover the SNP panel

```{bash, eval = F}
job_name="imputation_"${trait}  # customize
imputesub=`qsubshcom "Rscript -e \"SBayesRC::impute(mafile='${ma_file}_tidy.ma', LDdir='$LD_PATH', output='${ma_file}_imp.ma', log2file=TRUE) \"" 4 150G $job_name 12:00:00 " -wait=$tidyqsub  "   `
```

## SBayesRC: main function for SBayesRC

```{bash, eval = F}
job_name="sbr_eig_"${trait}
sbrcsub=`qsubshcom "Rscript -e \"SBayesRC::sbayesrc(mafile='${ma_file}_imp.ma', LDdir='$LD_PATH', outPrefix='${ma_file}_sbrc', annot='$annot', log2file=TRUE) \"" 10 150G $job_name 25:00:00 " -wait=$imputesub " `
```

## plot SBayesRC effect size

At last we compare the marginal effect size with the effect size from SBayesRC with a simple plot. 

```{bash, eval = F}
plotcmd="Rscript  ${exedir}/effect_size_plot_updated.R    $trait   ${gwas_file}.imp.ma  ${gwas_file}.imp_sbrc.txt  "
jobname="effect_plot_"${trait}
plotsub=`qsubshcom "$plotcmd"  1  50G  $jobname  1:00:00  " -wait=$sbrcsub  " `
```

As an example:

<img src="Anorexia_01_pgcAN2.2019-07.modified.vcf.tsv_sbrc.txt_compare_marginal_effect_vs_SBayesRC_20231103_10_18.png" width="50%" height="50%" />

# GCTB version

The same method is implemented into GCTB. Although the output files appear in different format, and the content of output are different, the core computation is using the same SBayesRC method in GCTB version and R version. 

```{bash, eval = F}
ldm=/QRISdata/Q3895/ldm/eigen/ukbEUR_Imputed/
```
## QC and Impute: GCTB does the two things in one step.

```{bash, eval = F}
job_name="imputation_"${trait}  
imputesub=`qsubshcom "gctb --ldm-eigen $ldm --gwas-summary ${ma_file}.ma --impute-summary --out ${ma_file}_imp.ma --thread 4"  4 150G $job_name 12:00:00 " -wait=$impqsub  "   `
```

## SBayesRC
```{bash, eval = F}
job_name="sbrc_gctb_"${trait}
sbrcsub=`qsubshcom "gctb  --sbayes RC  --ldm-eigen  ${ldm} --gwas-summary ${ma_file}_imp.ma.imputed.ma  --annot $annot --out  ${ma_file}_sbrc  --thread 4"  4 150G $job_name 72:00:00 " -wait=$imputesub  "   `
```

# Useful links:

GCTB:  
> https://cnsgenomics.com/software/gctb/#SBayesRCTutorial

SBayesRC:  
> https://github.com/zhilizheng/SBayesRC  

qsubshcom:  
> https://github.com/zhilizheng/qsubshcom  




