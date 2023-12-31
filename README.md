# SBayesRC_pipeline
This is a pipeline how I use SBayesRC method to generate predictor. The 3 key steps of QC, imputation and SBayesRC are adopted derectly from [Zhili's method](https://github.com/zhilizheng/SBayesRC).

## define path and input/out 

```{bash, eval = F}
cd $workingpath

## this is where you placed file cojo_format_v3.R
exedir="./"

## this is where you placed LD files and annotation file. They are available for downloading from GCTB website:
## https://cnsgenomics.com/software/gctb/#Download

LD_PATH1="/scratch/project_mnt/S0007/uqzzhen4/project/UKB/LD/ukb20k_7M_4cM/"
LD_PATH2="/scratch/project_mnt/S0007/uqzzhen4/project/UKB/LD/ukb20k_hm3_4cM/"
annot="/scratch/project_mnt/S0007/uqzzhen4/project/UKB/annot/annot_baseline2.2.txt" # annottion file

## you can put your GWAS file into a folder and name it with the trait
trait=MDD_01
gwas_file=PGC_UKB_depression_genome-wide.txt

ma_file=${trait}/${gwas_file}

```


## format to cojo

for example, if the raw gwas file has a header like this:

> MarkerName A1 A2 Freq LogOR StdErrLogOR P

you can put each element of the header name into the script per flag:

```{bash, eval = F}
cmd1="Rscript  ${exedir}/cojo_format_v6.R  \
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
> 1. Convert odds ratio to effect size if it's named OR.
> 2. Make a allele frequency comparison plot between data vs LD reference
> 3. If allele frequency is missing in your data, it fills it with the AF in LD reference.
> 4. If there is not a column of per-SNP-sample-size in your data, you can fill it with a unique number.
> 5. If your data has two columns for the sample size as Ncase and Ncontrol, the script adds them up to be N. Put them in as "Ncase,Ncontrol" behind --samplesize.
> 6. If the data does not have sample size but "NCHROBS", it will be divided by 2 to be sample size.
> 7. If SE is missing, we will estimate it from BETA and P.
> 8. If SNP ID is in chr:pos or chr:pos:A1:A2 format, we will replace it by denoting the chr column and bp column, and the SNP column name as "missing".
> 9. If sample size and allele frequency are give separately for cases and controls, we can input them together and separate with ",". (make sure the AF of case/control are put in the same order as their sample size)

 
As an exmple:

<img src="GIANT_HEIGHT_YENGO_2022_GWAS_SUMMARY_STATS_EUR_AF_plot.png" width="50%" height="50%" />





## Tidy: optional step, tidy summary data

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

## Impute: optional step if your summary data doesn't cover the SNP panel

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


## Useful links:

GCTB:  
> https://cnsgenomics.com/software/gctb/#Download  

SBayesRC:  
> https://github.com/zhilizheng/SBayesRC  

qsubshcom:  
> https://github.com/zhilizheng/qsubshcom  




