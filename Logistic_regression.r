#######################################################################################################################################
## Genome-wide association study of Pancreas cancer risk                                                                             ##
## Evangelina López de Maturana                                                                                                      ##
## Logistic regression models to assess the association between SNPs and pancreas cancer risk                                                                                                                                  ##
## 
#######################################################################################################################################

rm(list=ls())
###############################################################################
require(data.table)
library(parallel)
###############################################################################
#caco variable refers to case control status
#data.df is the dataframe containing the adjustment variables
#geno is the dataframe with the individuals in the rows and the SNPs in the columns
#geno and data.df dataframes are ordered by the individual id

  mylm <- function(i){
    
    GenoMod  <- glm(as.factor(caco) ~  age + gender_Phenotype + as.factor(Country.collapse)+ PC1+PC2+PC3+PC4+PC5+ geno[,i],data = data.df, family=binomial(logit))
    
    c(
      colnames(geno)[i], 
      exp(summary(GenoMod)$coefficients['geno[, i]',"Estimate"]),
      summary(GenoMod)$coefficients['geno[, i]',"Pr(>|z|)"],
      summary(GenoMod)$coefficients['geno[, i]',"Estimate"],
      summary(GenoMod)$coefficients['geno[, i]',"Std. Error"],
      dim(data.df)[1]
    )
  }
  
  Out <- mclapply(1:dim(geno)[2] , mylm, mc.cores=16)  
  Result <- data.frame(do.call(rbind,Out) ,stringsAsFactors = FALSE)
  colnames(Result) <- c('SNP','OR','pvalue','Estimate', 'SE','N')
