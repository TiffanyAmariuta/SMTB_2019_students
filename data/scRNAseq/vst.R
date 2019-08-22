data = readRDS("~/R/SMTB_2019_students/data/scRNAseq/sparse_data.rds")

library("Rcpp")
library("Matrix")

setwd("~/R/SMTB_2019_students/data/scRNAseq/")
sourceCpp("papa_mojet_v_c.cpp")

vst = function(data, loess.spam=0.3, power=2) {
  expressors = t(data)@p[-1]
  
  real_means = asS4(as.vector(rowSums(data)/expressors, mode="numeric"))
  
  varience = SparseRowVar2(data, real_means, power)

  to_loess = data.frame(mean=real_means, var=varience, exp=expressors)
  loess_output = loess(data=to_loess,
                       formula=log10(var) ~ log10(mean),
                       spam=loess.spam)
  
  expected_varience = 10 ^ loess_output$fitted
  
  vst_scores = SparseRowVarStd2(data, 
                                real_means, 
                                expected_varience, 
                                ncol(data)^(power/2),
                                power)
  
  return(vst_scores)
}

vst_scores = vst(data, power=2)

library("Seurat")

vst_output = order(as.numeric(vst_scores), decreasing=TRUE)
seur_output = order(seur_out$vst.variance.standardized, decreasing = TRUE)

vst_best = vst_output[1:1000]
seur_best = seur_output[1:1000]



in_common = intersect(vst_best, seur_best)
length(in_common)
setwd("~/R/SMTB_2019_students/data/scRNAseq/vst/")
saveRDS(vst_output, "2nzm.rds")


anti_jinny = function(expression) {
  expression = sort(expression, decreasing = TRUE)
  under_curve = 0
  i = 1
  
  for (num in expression) {
    under_curve = under_curve + i*num
    i = i+1
  }
  
  return((sum(expression)*length(expression))/under_curve)
}