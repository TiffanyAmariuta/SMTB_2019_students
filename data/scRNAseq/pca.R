a = readRDS("~/R/SMTB_2019_students/data/scRNAseq/shrinked_table.rds")
range = readRDS("~/R/SMTB_2019_students/data/scRNAseq/8power_median_rate.rds")

significant = a[range[1:1000],]
rm(a)

pca_res = prcomp(t(significant), center = T, scale. = T)

df = as.data.frame(pca_res$x)

clusters = kmeans(t(significant), 5)
colors = c("blue", "red", "green", "black", "magenta")
plot(df$PC1, df$PC2, col = colors[clusters$cluster])
  