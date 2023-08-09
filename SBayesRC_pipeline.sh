#!/bin/bash

##############################################################################
## this is the part to modify for each trait
## you will need to put your GWAS file into a folder and name it with the trait

trait="CNT_02"
gwas_file="chronotype_raw_BOLT.output_HRC.only_plus.metrics_maf0.001_hwep1em12_info0.3.txt"

## tell the script how the GWAS summary statistic named each column of essential information:
SNP="SNP"
A1="ALLELE1"
A2="ALLELE0"
AF="A1FREQ"
PVAL="P_BOLT_LMM"
BETA="BETA"
SE="SE"
Nsize=449734

# special issues:
# put in "missing" for AF if it's not available in sumstat
# put in column name of Nsize if it's available as per snp sample size
# put in both column name separated with "," if there are two columns of sample size of cases and controls. 
# put in the general sample size as a number if you find can only find it from publication. 

####################################################################################
#### define path and input/out 
#### this is the part to modify when you change your environment.

cd /scratch/project/genetic_data_analysis/uqtlin5/Predictor_Generation

## this is where you placed file cojo_format_v3.R
exedir="./"

## this is where you placed LD files and annotation file. They are available for downloading from GCTB website:
## https://cnsgenomics.com/software/gctb/#Download
LD_PATH1="/scratch/project_mnt/S0007/uqzzhen4/project/UKB/LD/ukb20k_7M_4cM/"
LD_PATH2="/scratch/project_mnt/S0007/uqzzhen4/project/UKB/LD/ukb20k_hm3_4cM/"
annot="/scratch/project_mnt/S0007/uqzzhen4/project/UKB/annot/annot_baseline2.2.txt" # annottion file

#####################################################################################
## Below this line does not need modification most of time. 

ma_file=${trait}/${gwas_file}

## format to cojo

cmd1="Rscript  ${exedir}/cojo_format_v3.R  \
  --file  ${ma_file}  \
  --out  ${ma_file}.ma   \
  --SNP  $SNP   \
  --A1  $A1    \
  --A2   $A2    \
  --freq  $AF   \
  --pvalue  $PVAL  \
  --beta  $BETA  \
  --se  $SE   \
  --samplesize  $Nsize   "

job_name="format_"${trait}
formatqsub=`qsubshcom  "$cmd1" 1 50G  $job_name  2:00:00  " "     `



## check row numbers
cmd2="if [ $(wc -l  ${ma_file}.ma  | awk '{print $1}'  )  -ne   $(wc -l   ${ma_file}  | awk '{print $1}' ) ]； then   echo \"formatted file could be truncated or filtered with allele frequency. Double check!\"  ; fi  "
job_name="check1_"${trait}
checkqsub=`qsubshcom "$cmd2"   1 1G  $job_name  1:00:00  " -wait=$formatqsub " ` 


## Tidy: optional step, tidy summary data
job_name="tidy_"${trait} 
tidyqsub=`qsubshcom "Rscript -e \"SBayesRC::tidy(mafile='${ma_file}.ma', LDdir='$LD_PATH1', output='${ma_file}_tidy.ma', log2file=TRUE) \"" 1 50G $job_name 10:00:00 "  -wait=$checkqsub  " `
## Best practice: read the log to check issues in your GWAS summary data.  


## choose LD matrix based on number of SNPs
job_name="ldpick_"${trait}
ldpick=`qsubshcom "if [ $(wc -l  ${trait}/${gwas_file}_tidy.ma  | awk '{print $1}'  ) -gt  5149563 ]; then      echo "yes"; else         echo "no" ; fi"   1 1G  $job_name  1:00:00  " -wait=$tidyqsub "  `


## Impute: optional step if your summary data doesn't cover the SNP panel
job_name="imputation_"${trait}  
imputesub=`qsubshcom "Rscript -e \"SBayesRC::impute(mafile='${ma_file}_tidy.ma', LDdir='$LD_PATH', output='${ma_file}_imp.ma', log2file=TRUE) \"" 4 150G $job_name 12:00:00 " -wait=$ldpick  "   `


## SBayesRC: main function for SBayesRC
job_name="sbr_eig_"${trait}
sbrcsub=`qsubshcom "Rscript -e \"SBayesRC::sbayesrc(mafile='${ma_file}_imp.ma', LDdir='$LD_PATH', outPrefix='${ma_file}_sbrc', annot='$annot', log2file=TRUE) \"" 10 150G $job_name 25:00:00 " -wait=$imputesub " `


## plot SBayesRC effect size
plotcmd=" Rscript  ${exedir}/effect_size_plot.R  $trait   $gwas_file "
jobname="effect_plot_"${trait}
plotsub=`qsubshcom "$plotcmd"  1  50G  $jobname  1:00:00  " -wait=$sbrcsub  " `


