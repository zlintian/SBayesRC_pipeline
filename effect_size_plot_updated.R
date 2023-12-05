args=commandArgs(trailingOnly = TRUE)

trait=args[1]
gwas.file = args[2]
predictor.file=args[3]

library(data.table)
library(ggplot2)

gwas = data.frame(fread(paste0(trait, "/", gwas.file)))
predictor=data.frame(fread(paste0(trait, "/", predictor.file)))


predictor$effect.in.gwas = as.numeric(gwas[match(predictor$SNP, gwas$SNP),"b"])
predictor$A1.in.gwas =  gwas[match(predictor$SNP, gwas$SNP),"A1"]

predictor$marginal.effect = predictor[,"effect.in.gwas"] *(sign((as.numeric(predictor$A1 == predictor$A1.in.gwas) - 0.5)))

predictor$gwas.se =  gwas[match(predictor$SNP, gwas$SNP),"se"]
predictor$threshold = abs(predictor$effect.in.gwas) + 3*(predictor$gwas.se)
predictor$outlier = (abs(predictor$BETA) > abs(predictor$threshold))

st=format(Sys.time(), "%Y%m%d_%H_%M")

write.table(predictor[which(predictor$outlier == TRUE),], paste0("Figures/", trait,"_"  , predictor.file ,  "_effect_size_outliers.txt"  ), quote = F, sep ="\t"  )

effect.plot = ggplot(data = predictor, aes(x = marginal.effect, y = BETA,  color= outlier)) + geom_point(size = 0.4) +  
  geom_abline(intercept=0, slope=1, color="blue")    +
  labs(title=predictor.file,  x="GWAS marginal effect", y = "SBayesR Effect size") 


ggsave(effect.plot, filename = paste0("Figures/", trait, "_"  , predictor.file ,  "_compare_marginal_effect_vs_SBayesRC_", st, ".png"), height = 8, width = 8)




#pip.plot = ggplot(data = predictor[which(predictor$PIP > 0.1),], aes(x = BETA, y = PIP )) + geom_point()

#ggsave(pip.plot, filename = paste0(trait, "/", gwas.file ,  "_PIP_vs_BETA.png"), height = 6, width = 8)


