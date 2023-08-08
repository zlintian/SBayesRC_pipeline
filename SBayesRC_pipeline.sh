#!/bin/bash

####################################################################################
#### define path and input/out 

cd /scratch/user/uqtlin5/Predictor_Generation/

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


#####################################################################################33

## format to cojo

## for example, if the raw gwas file has a header like this:
## MarkerName A1 A2 Freq LogOR StdErrLogOR P
## you can put each element of the header name into the script per flag:


cmd1="Rscript  ${exedir}/cojo_format_v3.R  \
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


## check row numbers
cmd2="if [ $(wc -l  ${trait}/${gwas_file}.ma  | awk '{print $1}'  )  -ne   $(wc -l   ${trait}/${gwas_file}  | awk '{print $1}' ) ]； then   echo "formatted file could be truncated or filtered with allele frequency. Double check!"  ; fi  "
job_name="check1_"${trait}
checkqsub=`qsubshcom "$cmd2"   1 1G  $job_name  1:00:00  " -wait=$formatqsub " ` 


## Tidy: optional step, tidy summary data
## "log2file=TRUE" means the messages will be redirected to a log file 
ma_file=${trait}/${gwas_file}
job_name="tidy_"${trait} # your job name, customize
tidyqsub=`qsubshcom "Rscript -e \"SBayesRC::tidy(mafile='${ma_file}.ma', LDdir='$LD_PATH1', output='${ma_file}_tidy.ma', log2file=TRUE) \"" 1 50G $job_name 10:00:00 "  -wait=$checkqsub  " `
## Best practice: read the log to check issues in your GWAS summary data.  


## choose LD matrix based on number of SNPs
job_name="ldpick_"${trait}
ldpick=`qsubshcom "if [ $(wc -l  ${trait}/${gwas_file}_tidy.ma  | awk '{print $1}'  ) -gt  5149563 ]; then      echo "yes"; else         echo "no" ; fi"   1 1G  $job_name  1:00:00  " -wait=$tidyqsub "  `


## Impute: optional step if your summary data doesn't cover the SNP panel
job_name="imputation_"${trait}  # customize
imputesub=`qsubshcom "Rscript -e \"SBayesRC::impute(mafile='${ma_file}_tidy.ma', LDdir='$LD_PATH', output='${ma_file}_imp.ma', log2file=TRUE) \"" 4 150G $job_name 12:00:00 " -wait=$ldpick  "   `


## SBayesRC: main function for SBayesRC
job_name="sbr_eig_"${trait}
sbrcsub=`qsubshcom "Rscript -e \"SBayesRC::sbayesrc(mafile='${ma_file}_imp.ma', LDdir='$LD_PATH', outPrefix='${ma_file}_sbrc', annot='$annot', log2file=TRUE) \"" 10 150G $job_name 25:00:00 " -wait=$imputesub " `


## plot SBayesRC effect size
plotcmd=" Rscript  ${exedir}/effect_size_plot.R  $trait   $gwas_file "
jobname="effect_plot_"${trait}
plotsub=`qsubshcom "$plotcmd"  1  50G  $jobname  1:00:00  " -wait=$sbrcsub  " `


