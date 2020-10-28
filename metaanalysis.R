##################################################################################################################################################
##################################################################################################################################################
#######################Metaanalyses with PanGenEU, PanScanI-II, PanScanIII, Panc4################################################################
rm(list=ls())

library(metafor)

#res is the dataframe including the summary statistics of the different studies  with the fields c("study", "studyname", "CHR_BP", "CHR", "SNP", "BP", "logOR", "OR", "SE")
snplist<-unique(res$CHR_BP)

list.of.fit<-lapply(snplist, function(x) rma.uni(yi=logOR, sei=SE, data=res[res$CHR_BP==x, ], method="REML", slab=res[res$CHR_BP==x, "studyname"]))

names(list.of.fit)<-snplist
list.of.results<-lapply(list.of.fit, function(est) unlist(est[c("b", "se", "pval", "ci.lb", "ci.ub", "I2")]))
results.meta.analysis.REML<-data.frame( do.call(rbind, list.of.results))
results.meta.analysis.REML$OR<-exp(results.meta.analysis.REML$b)

pdf(file.path(main.dir, "results", "forest_plot.pdf"))

lapply(rownames(results.meta.analysis.REML) , function(snp) {
  forest(list.of.fit[[snp]], main=snp)
})

dev.off()