#!/usr/bin/Rscript

# load required libraries
library(optparse)
library(data.table)
library(ggplot2)

freq.file="/scratch/project_mnt/S0007/uqzzhen4/project/UKB/LD/ukb20k_7M_4cM/snp.info"

# Define the command line argument parser
option_list <- list(
  make_option(c("-f", "--file"),         type="character",  help="path to summary statistics"),
  make_option(c("-i", "--SNP"),          type="character", default= "SNP"  , help="column name of SNP ID used in the sum stat file
"),
  make_option(c("-A", "--A1"),           type="character", default= "A1", help="column name of A1 used in the sum stat file"),
  make_option(c("-B", "--A2"),           type="character", default= "A2", help="column name of A2 used in the sum stat file"),
  make_option(c("-q", "--freq"),         type="character", default= "freq", help="column name of A1 allele frequency used in the s
um stat file"),
  make_option(c("-b", "--beta"),         type="character", default= "beta", help="column name of effect size used in the sum stat 
file"),
  make_option(c("-s", "--se"),           type="character", default= "se", help="column name of standard error used in the sum stat
 file"),
  make_option(c("-p", "--pvalue"),       type="character", default= "pvalue", help="column name of p.value used in the sum stat fi
le"),
  make_option(c("-n", "--samplesize"),  type="character", default= "N", help="column name of sample size used in the sum stat file
, or general sample size if it's missing in the file"),
  make_option(c("-o", "--out"),       type="character", default= "formatted_with_cojo_format.txt", help="file name of ouput")
  )


args <- parse_args(OptionParser(option_list=option_list))



# Access the values of the command line arguments

file.name   <- args$file
colname.SNP <- args$SNP 
colname.A1  <- args$A1  
colname.A2  <- args$A2  
colname.freq<- args$freq
colname.b   <- args$beta   
colname.se  <- args$se  
colname.p   <- args$pvalue   
colname.N   <- args$samplesize   
output   <- args$out




# print on screen for error check


cat("file name is ", file.name, "\n")
cat("column name used for SNP ID in this file is ", colname.SNP, "\n")
cat("column name used for A1 in this file is ", colname.A1, "\n")
cat("column name used for A2 in this file is ", colname.A2, "\n")
cat("column name used for freq in this file is ", colname.freq, "\n")
cat("column name used for effect size in this file is ", colname.b, "\n")
cat("column name used for SE in this file is ", colname.se, "\n")
cat("column name used for p value in this file is ", colname.p, "\n")
cat("column name used for sample size or the overal sample size in this file is ", colname.N, "\n")



# make a new table

gwas = data.frame(fread(file.name))
n.snp = nrow(gwas)

formatted.gwas = data.frame(
 matrix(NA, nrow = n.snp, ncol = 8)
)

colnames(formatted.gwas) = c("SNP", "A1", "A2", "freq", "b", "se", "p", "N")

formatted.gwas$SNP = gwas[,colname.SNP]
formatted.gwas$A1= toupper(gwas[,colname.A1])
formatted.gwas$A2 = toupper(gwas[,colname.A2])
formatted.gwas$se = gwas[,colname.se]
formatted.gwas$p = gwas[,colname.p]

## effect size is put in as it is, but if it's OR, we log it

if(colname.b == "OR"){
	formatted.gwas$b = log(gwas[,colname.b])
	}else{
	formatted.gwas$b = gwas[,colname.b]
}


## make a plot to check AF 

ref.freq = data.frame(fread(freq.file))
ref.freq$gwasA1 = formatted.gwas[match(ref.freq$ID, formatted.gwas$SNP),"A1"]
ref.freq$sign = sign(( as.numeric(ref.freq$gwasA1  == ref.freq$A1))-0.5)
ref.freq$gwasA1freq= abs(as.numeric(ref.freq$gwasA1  != ref.freq$A1) - ref.freq$A1Freq)


## allele frequency is sometimes missing

if(colname.freq == "missing"){
 # freq.file = "/afm01/UQ/Q3046/Reference//HRC.r1-1.GRCh37.wgs.mac5.sites.tab"
  cat("AF was not supplied in the file, so we will borrow allel frequency from data ", freq.file, "\n")
  
  formatted.gwas$freq = ref.freq[match( formatted.gwas$SNP, ref.freq$ID),"gwasA1freq"]
  
}else{
    
  formatted.gwas$freq = gwas[,colname.freq]
  
  ref.freq$freq.in.gwas = formatted.gwas[match(ref.freq$ID, formatted.gwas$SNP) , "freq"]
  freq.plot = ggplot(data = ref.freq, aes(x = gwasA1freq, y =freq.in.gwas )) + geom_point(size = 0.2) + xlab("AF in reference data") + ylab("AF in GWAS summary statistic")
  ggsave(paste0(file.name, "_AF_plot.png"), freq.plot, height = 8, width = 8)
  
  }




## sample size is sometimes missing

sample.size.names=unlist(strsplit(colname.N, ","))



if(length(sample.size.names)==1){
	if(is.na(as.numeric(colname.N) )== T) {
  
    		formatted.gwas$N = gwas[, colname.N]

	}else{
   	formatted.gwas$N = as.numeric(colname.N)
 	}
}else{
	formatted.gwas$N= gwas[, sample.size.names[1]] + gwas[,sample.size.names[2]]
}


formatted.gwas = formatted.gwas[which(is.na(formatted.gwas$freq) == F) ,]


write.table(formatted.gwas, file=output, quote = F, sep ="\t", row.names = F)


