#######################################################################################################################################
## Genome-wide association study of Pancreas cancer risk                                                                             ##
## Evangelina López de Maturana                                                                                                      ##
## Logistic regression models to assess the association between SNPs and pancreas cancer risk                                                                                                                                  ##
## 
#######################################################################################################################################

rm(list=ls())
###############################################################################
## Set script timer
initTime <- proc.time()
require(data.table)
library(parallel)

###############################################################################
###Logistic regression to assess the SNP association with PC risk
  mylm <- function(i){
    
    GenoMod  <- glm(as.factor(caco) ~  age + gender_Phenotype + as.factor(Country.collapse)+ PC1+PC2+PC3+PC4+PC5+ geno2[,i],data = data.df, family=binomial(logit))
    
    c(
      colnames(geno)[i], 
      exp(summary(GenoMod)$coefficients['geno',"Estimate"]),
      summary(GenoMod)$coefficients['geno[, i]',"Pr(>|z|)"],
      summary(GenoMod)$coefficients['geno[, i]',"Estimate"],
      summary(GenoMod)$coefficients['geno[, i]',"Std. Error"],
      dim(data.df)[1]
    )
  }
  
  Out <- mclapply(1:dim(geno2.pangen.isblac.epicuro.t.round.MAF0.05)[2] , mylm, mc.cores=16)  
  Result <- do.call(rbind,Out) # make list a matrix
  Result <- data.frame (Result,stringsAsFactors = FALSE)
  colnames(Result) <- c('SNP','OR_PANGEN_ISBLAC_EPICURO','pvalue_PANGEN_ISBLAC_EPICURO','Estimate_PANGEN_ISBLAC_EPICURO','SE_PANGEN_ISBLAC_EPICURO','N_PANGEN_ISBLAC_EPICURO')

  #Logistic regression to assess the OR of the SNPs # PANGEN+ISBLAC
  mylm <- function(i){
    
    GenoMod  <- glm(as.factor(caco) ~  age + gender_Phenotype + as.factor(Country.collapse)+ PC1+PC2+PC3++PC4+PC5+  geno2.pangen.isblac.t.round.MAF0.05[,i],data = final.data.pangen1.isblac1, family=binomial(logit))
    #rownames(summary(GenoMod)$coefficients)[6] = colnames(geno2.t.round)[i]
    #return the vector
    
    c(
      colnames(geno2.pangen.isblac.t.round.MAF0.05)[i], 
      exp(summary(GenoMod)$coefficients['geno2.pangen.isblac.t.round.MAF0.05[, i]',"Estimate"]),
      summary(GenoMod)$coefficients['geno2.pangen.isblac.t.round.MAF0.05[, i]',"Pr(>|z|)"],
      summary(GenoMod)$coefficients['geno2.pangen.isblac.t.round.MAF0.05[, i]',"Estimate"],
      summary(GenoMod)$coefficients['geno2.pangen.isblac.t.round.MAF0.05[, i]',"Std. Error"],
      dim(final.data.pangen1.isblac1)[1]
    )
  }
  
  Out <- mclapply(1:dim(geno2.pangen.isblac.t.round.MAF0.05)[2] , mylm, mc.cores=16)  
  Result_PANGEN_ISBLAC <- do.call(rbind,Out) # make list a matrix
  Result_PANGEN_ISBLAC <- data.frame (Result_PANGEN_ISBLAC,stringsAsFactors = FALSE)
  colnames(Result_PANGEN_ISBLAC) <- c('SNP','OR_PANGEN_ISBLAC','pvalue_PANGEN_ISBLAC','Estimate_PANGEN_ISBLAC','SE_PANGEN_ISBLAC','N_PANGEN_ISBLAC')
  
  #Logistic regression to assess the OR of the SNPs # only PANGEN
  mylm <- function(i){
    
    GenoMod  <- glm(as.factor(caco) ~  age + gender_Phenotype + as.factor(Country.collapse)+ PC1+PC2+PC3+ +PC4+PC5+ geno2.pangen.t.round.MAF0.05[,i],data = final.data.pangen1, family=binomial(logit))

    c(
      colnames(geno2.pangen.t.round.MAF0.05)[i], 
      exp(summary(GenoMod)$coefficients['geno2.pangen.t.round.MAF0.05[, i]',"Estimate"]),
      summary(GenoMod)$coefficients['geno2.pangen.t.round.MAF0.05[, i]',"Pr(>|z|)"],
      summary(GenoMod)$coefficients['geno2.pangen.t.round.MAF0.05[, i]',"Estimate"],
      summary(GenoMod)$coefficients['geno2.pangen.t.round.MAF0.05[, i]',"Std. Error"],
      dim(final.data.pangen1)[1]
    )
  }
  
  Out <- mclapply(1:dim(geno2.pangen.t.round.MAF0.05)[2] , mylm, mc.cores=16)  
  Result_PANGEN <- do.call(rbind,Out) # make list a matrix
  Result_PANGEN <- data.frame (Result_PANGEN,stringsAsFactors = FALSE)
  colnames(Result_PANGEN) <- c('SNP','OR_PANGEN','pvalue_PANGEN','Estimate_PANGEN','SE_PANGEN','N_PANGEN')

  Result <- merge(Result,Result_PANGEN_ISBLAC,by='SNP',sort=FALSE,all.x=TRUE,all.y=TRUE) #I choose all the SNPs in both datasets Meeting with Núria 27/10/2016
  Result <- merge(Result,Result_PANGEN,by='SNP',sort=FALSE,all.x=TRUE,all.y=TRUE)  

  dim(Result)
  Result_MAF <- merge(Result,rs_MAF_info_geno.pangen [,c('SNP','MAF_PANGEN')],by='SNP',sort=FALSE,all.x=TRUE)
  dim(Result_MAF)
  Result_MAF <- merge(Result_MAF,rs_MAF_info_geno.isblac[,c('SNP','MAF_ISBLAC')],by='SNP',sort=FALSE,all.x=TRUE)
  dim(Result_MAF)
   Result_MAF <- merge(Result_MAF,rs_MAF_info_geno.pangen.isblac[,c('SNP','MAF_PANGEN_ISBLAC')],by='SNP',sort=FALSE,all.x=TRUE)
  dim(Result_MAF)
  Result_MAF <- merge(Result_MAF,rs_MAF_info_geno.pangen.isblac.epicuro[,c('SNP','MAF_PANGEN_ISBLAC_EPIC')],by='SNP',sort=FALSE,all.x=TRUE)
  dim(Result_MAF)
  Result_MAF <- merge(Result_MAF,rs_MAF_info_geno.ESP [,c('SNP','MAF_ESP')],by='SNP',sort=FALSE,all.x=TRUE)
  dim(Result_MAF)
  Result_MAF <- merge(Result_MAF,rs_MAF_info_geno.ESP.pangen[,c('SNP','MAF_ESP.pangen')],by='SNP',sort=FALSE,all.x=TRUE)
  dim(Result_MAF)
  Result_MAF <- merge(Result_MAF,rs_MAF_info_geno.ESP.isblac[,c('SNP','MAF_ESP.isblac')],by='SNP',sort=FALSE,all.x=TRUE)
  dim(Result_MAF)
  Result_MAF <- merge(Result_MAF,rs_MAF_info_geno.UK.IRE.NIRE[,c('SNP','MAF_UK.IRE.NIRE')],by='SNP',sort=FALSE,all.x=TRUE)
  dim(Result_MAF)
  Result_MAF <- merge(Result_MAF,rs_MAF_info_geno.GERMANY[,c('SNP','MAF_GER')],by='SNP',sort=FALSE,all.x=TRUE)
  dim(Result_MAF)
  Result_MAF <- merge(Result_MAF,rs_MAF_info_geno.ITA[,c('SNP','MAF_ITA')],by='SNP',sort=FALSE,all.x=TRUE)
  dim(Result_MAF)
  Result_MAF <- merge(Result_MAF,rs_MAF_info_geno.SWE[,c('SNP','MAF_SWE')],by='SNP',sort=FALSE,all.x=TRUE)
  
  Result_MAF <- merge(Result_MAF,rs_MAF_info_geno.epicuro[,c('SNP','MAF_EPIC','pos_alleles')],by='SNP',sort=FALSE,all.x=TRUE)
  
  Result_MAF <- merge(Result_MAF,rs_MAF_info_geno.pangen [,c('SNP','chr','position', 'Major_allele', 'Minor_allele')],by='SNP',sort=FALSE,all.x=TRUE)
  
  
  print('dim(Result_MAF)')
  print(dim(Result_MAF))
  
  Result_MAF <- merge(Result_MAF,SNPs_EPICURO_chr[,c('pos_alleles','info','type','array','flip_del')],by='pos_alleles',sort=FALSE,all.x=TRUE)
  print('dim(Result_MAF)')
  print(dim(Result_MAF))

  Result_MAF1 <- rbind(Result_MAF, Result_MAF1)
}

Result_MAF1$pvalue=as.numeric(Result_MAF1$pvalue_PANGEN_ISBLAC)
Result_MAF2 <- Result_MAF1[order(Result_MAF1$pvalue_PANGEN_ISBLAC_EPICURO,decreasing=FALSE),]

save(Result_MAF2,file='Results.RData') 
###############################################################################

###############################################################################
## Display time elapsed
proc.time() - initTime
###############################################################################