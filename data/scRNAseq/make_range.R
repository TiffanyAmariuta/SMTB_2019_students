a = readRDS("~/R/SMTB_2019_students/data/scRNAseq/shrinked_table.rds")
powers = c(2, 4, 8, 16, 32)

for (power in powers) {
  method_name = sprintf("nzm%d", power)
  
  our_av_dev = function(expression) {
    deviation = expression - (sum(expression))/(sum(expression > 0))
    return(mean(deviation^power)^(1/power))
  }
  av_devs = c()
  step = 300
  for (i in seq(1, nrow(a), by=step)) {
    last = min(c(step + i - 1, nrow(a)))
    av_devs = c(av_devs, apply(a[i:last,], 1, our_av_dev))
  }
  
  dev_rate = order(av_devs, decreasing=TRUE)
  saveRDS(dev_rate, sprintf("~/R/SMTB_2019_students/data/scRNAseq/rates/%s.rds", method_name))
}

dev_rate = readRDS("jinny.rds")
lengths = c()

for (gene_index in dev_rate) {
  lengths = c(lengths, sum(a[gene_index,] > 0))
}

plot(1:nrow(a), lengths)
dev_range[1:5]
dev_rate[1]
exp = a[dev_rate[6],]

count_nz = function(x) {
  return(sum(x > 0))
}

nz = c()
step = 150
for (i in seq(1, nrow(a), by=step)) {
  last = min(c(step + i - 1, nrow(a)))
  nz = c(nz, apply(a[i:last,], 1, count_nz))
}
unicell = nz == 1
one_cell_genes_exp = a[unicell,]
one_cell_genes_exp = t(one_cell_genes_exp)
sum(rowSums(one_cell_genes_exp) > 0)
  