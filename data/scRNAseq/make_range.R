a = readRDS("~/R/SMTB_2019_students/data/scRNAseq/shrinked_table.rds")
our_av_dev = function(expression) {
  deviation = expression - median(expression)
  return(mean(deviation^2)^(1/2))
}
av_devs = c()
step = 300
for (i in seq(1, nrow(a), by=step)) {
  last = min(c(step + i - 1, nrow(a)))
  av_devs = c(av_devs, apply(a[i:last,], 1, our_av_dev))
}

install.packages("Seurat")

length(av_devs)
nrow(a)
dev_rate = order(av_devs, decreasing=TRUE)
saveRDS(dev_rate, "~/R/SMTB_2019_students/data/scRNAseq/2power_median_rate.rds")
dev_range = rownames(a)[dev_rate]
grep("PAX5", dev_range)
dev_range[1:5]
