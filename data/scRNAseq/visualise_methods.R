marker_genes = c("PDGFRA","ISLR","CD2","CD3D","CD79A","RALGPS2","CD14","C1QA")
colors = c("grey", "grey", "orange", "orange", "blue", "blue", "pink", "pink")
marker_gene_indices = match(marker_genes, rownames(a))

setwd("~/R/SMTB_2019_students/data/scRNAseq/rates/")
fnames = list.files(pattern = ".rds")
out_matrix = matrix(nrow = length(marker_genes), ncol = length(fnames))
colnames(out_matrix) = fnames
rownames(out_matrix) = marker_genes
for (name in fnames) {
  range = readRDS(name)
  indices = match(marker_gene_indices, range)
  out_matrix[,name] = indices
}

colnames(out_matrix) = gsub(".rds", "", fnames)

library("DT")
install.packages("DT")

df = as.data.frame(out_matrix)
out = datatable(df) %>% formatStyle(
  "2mean",
  target="row",
  backgroundColor=styleEqual(df[,"2mean"], colors)
)

saveWidget(out, "marker_genes_rate.html")
