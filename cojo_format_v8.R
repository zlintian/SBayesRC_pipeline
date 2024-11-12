#!/usr/bin/Rscript

# load required libraries
library(optparse)
library(data.table)
library(ggplot2)


#####################################################################################


# Define the command line argument parser
option_list <- list(
  make_option(c("-f", "--file"),         type="character",  help="path to summary statistics"),
  make_option(c("-i", "--SNP"),          type="character", default= "SNP"  , help="column name of SNP ID used in the sum stat file"),
  make_option(c("-A", "--A1"),           type="character", default= "A1", help="column name of A1 used in the sum stat file"),
  make_option(c("-B", "--A2"),           type="character", default= "A2", help="column name of A2 used in the sum stat file"),
  make_option(c("-q", "--freq"),         type="character", default= "freq", help="column name of A1 allele frequency used in the sum stat file"),
  make_option(c("-b", "--beta"),         type="character", default= "beta", help="column name of effect size used in the sum stat file"),
  make_option(c("-s", "--se"),           type="character", default= "se", help="column name of standard error used in the sum stat file"),
  make_option(c("-p", "--pvalue"),       type="character", default= "pvalue", help="column name of p.value used in the sum stat file"),
  make_option(c("-n", "--samplesize"),  type="character", default= "N", help="column name of sample size used in the sum stat file, or general sample size if it's missing in the file"),
  make_option(c("-o", "--out"),       type="character", default= "formatted_with_cojo_format.txt", help="file name of ouput"),
  make_option(c("-c", "--chr"),       type="character", default= "Chrom", help="chromosome of the SNP"),
  make_option(c("-l", "--pos"),       type="character", default= "POS", help="position of the SNP"),
  make_option(c("-z", "--z"),       type="character", default= "missing", help="z score"),
  make_option(c("-r", "--ref"),       type="character", default= "/scratch/project_mnt/S0007/uqzzhen4/project/UKB/LD/ukb20k_7M_4cM/snp.info", help="snp.info file of reference data")
)


args <- parse_args(OptionParser(option_list=option_list))

#####################################################################################


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
colname.chr <- args$chr
colname.pos <- args$pos
colname.z   <- args$z
output      <- args$out
freq.file   <- args$ref


#####################################################################################


# print on screen for error check


cat("file name is ", file.name, "\n")
cat("reference file is ", freq.file, "\n")
cat("column name used for SNP ID in this file is ", colname.SNP, "\n")
cat("column name used for A1 in this file is ", colname.A1, "\n")
cat("column name used for A2 in this file is ", colname.A2, "\n")
cat("column name used for freq in this file is ", colname.freq, "\n")
cat("column name used for effect size in this file is ", colname.b, "\n")
cat("column name used for SE in this file is ", colname.se, "\n")
cat("column name used for p value in this file is ", colname.p, "\n")
cat("column name used for sample size or the overal sample size in this file is ", colname.N, "\n")


#####################################################################################

# read in the data

gwas = data.frame(fread(file.name, check.names = FALSE))
ref.freq = data.frame(fread(freq.file))




#####################################################################################
# make a new table

n.snp = nrow(gwas)

formatted.gwas = data.frame(
  matrix(NA, nrow = n.snp, ncol = 8)
)

colnames(formatted.gwas) = c("SNP", "A1", "A2", "freq", "b", "se", "p", "N")




######################################################################################

# remap SNP ID

if(   (colname.SNP != "missing")  ){
  
  formatted.gwas$SNP = gwas[,colname.SNP]
  
}else{
  
  
  gwas$A = pmin(toupper(gwas[,colname.A1]), toupper(gwas[,colname.A2]))
  gwas$B = pmax(toupper(gwas[,colname.A1]), toupper(gwas[,colname.A2]))
  
  
  gwas$chrbpAB = paste0(gwas[,colname.chr], "_", gwas[,colname.pos], "_", gwas$A, "_", gwas$B)
  
  
  ref.freq$A = with(ref.freq, pmin(A1, A2))
  ref.freq$B = with(ref.freq, pmax(A1, A2))
  ref.freq$chrbpAB = paste0(ref.freq$Chrom, "_", ref.freq$PhysPos, "_", ref.freq$A, "_", ref.freq$B)
  
  gwas$matched_RSID = ref.freq[match(gwas$chrbpAB , ref.freq$chrbpAB), "ID"] 
  
  gwas[ !(gwas$chrbpAB %in%ref.freq$chrbpAB), "matched_RSID"]  =  gwas[ !(gwas$chrbpAB %in%ref.freq$chrbpAB), "chrbpAB" ]  
  
  formatted.gwas$SNP = gwas$matched_RSID
  
}


#####################################################################################

# forward the columns that are usually available

formatted.gwas$A1= toupper(gwas[,colname.A1])
formatted.gwas$A2 = toupper(gwas[,colname.A2])
formatted.gwas$p = as.numeric(gwas[,colname.p])



#####################################################################################



## sample size is sometimes missing, and sometimes as two columns with case and control separately
## if it is NCHROBS, we will divide it by 2. 
## if it's case and control separate, we will add them up.
## if it's missing, we will find it from the publication. 

sample.size.names=unlist(strsplit(colname.N, ","))


if(length(sample.size.names)==1){
  if(is.na(as.numeric(colname.N) )== T) {
    if(colname.N == "NCHROBS"){
      formatted.gwas$N = 0.5* (gwas[, colname.N])
    }else{
      formatted.gwas$N = gwas[, colname.N]
    }
  }else{
    formatted.gwas$N = as.numeric(colname.N)
  }
}else{


	if(is.na(as.numeric(sample.size.names[1]) )== T) {
  		formatted.gwas$N= gwas[, sample.size.names[1]] + gwas[,sample.size.names[2]]
	}else{ formatted.gwas$N=  as.numeric(sample.size.names[1]) + as.numeric(sample.size.names[2])   }


}


#####################################################################################

## Allele Frequency


allele.frequency.names=unlist(strsplit(colname.freq, ","))



ref.freq$gwasA1 = formatted.gwas[match(ref.freq$ID, formatted.gwas$SNP),"A1"]
ref.freq$sign = sign(( as.numeric(ref.freq$gwasA1  == ref.freq$A1))-0.5)
ref.freq$gwasA1freq= abs(as.numeric(ref.freq$gwasA1  != ref.freq$A1) - ref.freq$A1Freq)



if(length(allele.frequency.names) == 2) {
  ## Allele frequency is sometimes separately included for cases and controls. 
  freq.case = allele.frequency.names[1]
  freq.control = allele.frequency.names[2]
  N.case = sample.size.names[1]
  N.control = sample.size.names[2]
  
  formatted.gwas$freq = ( (gwas[,freq.case] * gwas[,N.case]) + (gwas[,freq.control] * gwas[,N.control]) ) /(gwas[,N.case] + gwas[,N.control])
  
  ref.freq$freq.in.gwas = formatted.gwas[match(ref.freq$ID, formatted.gwas$SNP) , "freq"]
  freq.plot = ggplot(data = ref.freq, aes(x = gwasA1freq, y =freq.in.gwas )) + geom_point(size = 0.2) + xlab("AF in reference data") + ylab("AF in GWAS summary statistic")
  ggsave(paste0(file.name, "_AF_plot.png"), freq.plot, height = 8, width = 8)
  
  
}else{
  
  
  ## allele frequency is sometimes missing
  
  if(colname.freq == "missing"){
    # freq.file = "/afm01/UQ/Q3046/Reference//HRC.r1-1.GRCh37.wgs.mac5.sites.tab"
    cat("AF was not supplied in the file, so we will borrow allel frequency from data ", freq.file, "\n")
    formatted.gwas$freq = ref.freq[match( formatted.gwas$SNP, ref.freq$ID),"gwasA1freq"]
    
  }else{
    
    # if freq is available, we will compare it with reference data and make a plot
    
    formatted.gwas$freq = as.numeric( gwas[,colname.freq] )
    
    ref.freq$freq.in.gwas = formatted.gwas[match(ref.freq$ID, formatted.gwas$SNP) , "freq"]
    freq.plot = ggplot(data = ref.freq, aes(x = gwasA1freq, y =freq.in.gwas )) + geom_point(size = 0.2) + xlab("AF in reference data") + ylab("AF in GWAS summary statistic")
    ggsave(paste0(file.name, "_AF_plot.png"), freq.plot, height = 8, width = 8)
    
  }
  
  
}


#####################################################################################

## effect size is put in as it is, but if it's OR, we log it
## If only z score is supplied but not beta, we will calculate


if(colname.b == "missing" & colname.z != "missing"){
  formatted.gwas$z = gwas[,colname.z]
  formatted.gwas$b = ( formatted.gwas$z )/( ((2 *(formatted.gwas$freq) * (1 - formatted.gwas$freq) * ( formatted.gwas$N + (formatted.gwas$z)^2)))^0.5 )
}  else if(colname.b == "OR" | colname.b == "odds_ratio"){
  formatted.gwas$b = log(gwas[,colname.b])
}else{
  formatted.gwas$b = gwas[,colname.b]
}




#####################################################################################

## if SE is missing, we will estimate it from BETA and P


if(colname.se == "missing"){
  
  cat("SE was not supplied in the file, so we will calculate it from known information")
  
  formatted.gwas$z = sign(formatted.gwas$b) * abs( qnorm(formatted.gwas$p/2) )
#  formatted.gwas$se = formatted.gwas$b/formatted.gwas$z
  formatted.gwas$se =1/sqrt(2 * formatted.gwas$freq *(1 - formatted.gwas$freq) *(formatted.gwas$N +  (formatted.gwas$z ^ 2)  )  )
 
  
}else{
  
  formatted.gwas$se = gwas[,colname.se]
  
  
}


#####################################################################################


formatted.gwas = formatted.gwas[which(is.na(formatted.gwas$freq) == F) ,]


write.table(formatted.gwas[,1:8], file=output, quote = F, sep ="\t", row.names = F)




