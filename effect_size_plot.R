args=commandArgs(trailingOnly = TRUE)

trait=args[1]
gwas.file = args[2]


library(data.table)
library(ggplot2)

gwas = data.frame(fread(paste0(trait, "/", gwas.file, ".ma")))
predictor=data.frame(fread(paste0(trait, "/", gwas.file, "_sbrc.txt")))


predictor$effect.in.gwas = gwas[match(predictor$SNP, gwas$SNP),"b"]
predictor$A1.in.gwas =  gwas[match(predictor$SNP, gwas$SNP),"A1"]

predictor$marginal.effect = predictor[,"effect.in.gwas"] *(sign((as.numeric(predictor$A1 == predictor$A1.in.gwas) - 0.5)))

effect.plot = ggplot(data = predictor, aes(x = marginal.effect, y = BETA)) + geom_point(size = 0.4) +  
  geom_abline(intercept=0, slope=1, color="blue")    +
  labs(title=paste0(nrow(predictor[which(predictor$BETA != 0),])  ," SNPs != 0 in predictor" ),  x="GWAS marginal effect", y = "SBayesR Effect size") 


ggsave(effect.plot, filename = paste0(trait, "/", trait,  "_compare_marginal_effect_vs_SBayesRC.png"), height = 8, width = 8)
