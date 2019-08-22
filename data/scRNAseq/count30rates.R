count_rate30 = function(devs) {
  return(sqrt((sum(devs[1:30]^2))/(sum(devs^2))))
}

setwd("~/R/SMTB_2019_students/data/scRNAseq/rates/")
fnames = list.files(pattern = ".rds")

prev_version = readRDS("30rates")
prev_version = c()
out = c(length=length(fnames))

i = 1
for (file in fnames) {
  if (sum(names(prev_version) == file) > 0) {
    out[i] = prev_version[file]
    i = i+1
    next
  }
  range=readRDS(file)
  sign_genes = 1000
  sign = data[range[1:sign_genes],]
  
  pca = prcomp(t(sign), center = T, scale. = T)
  
  out[i] = count_rate30(pca$sdev)
  i = i + 1
}

names(out) = fnames
saveRDS(out, "30rates")

