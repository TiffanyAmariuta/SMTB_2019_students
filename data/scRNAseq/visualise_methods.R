marker_genes = c("PDGFRA","CD2", "CD79A",  "CD14")
colors = c("rgb(150, 150, 150)", "orange", "rgb(110, 117, 216)", "pink")

marker_gene_indices = match(marker_genes, rownames(data))

setwd("~/R/SMTB_2019_students/data/scRNAseq/rates/")
fnames = list.files(pattern = ".rds")
fnames
fnames = fnames[c(2, 3, 7, 12, 14)]
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

rates30 = readRDS("30rates")[fnames]
rates30 = round(rates30*100, digits = 2)
rates30 = sprintf("<b>%g%%</b>", rates30)

rates30 = t(as.matrix(rates30))
rownames(rates30) = c("<b>30 score</b>")

df = as.data.frame(rbind(out_matrix, rates30))
saveRDS(df, "~/R/SMTB_2019_students/data/scRNAseq/simp_comparison.rds")
colors = c(colors, "white")

desc = data.frame(marker_gene = c("PDGFRA(fib)","CD2(Tc)", "CD79A(Bc)",  "CD14(monoc.)", "30 score"))

out = datatable(desc, 
                rownames=FALSE,
                width=85*ncol(desc),
                options=list(searching=FALSE,
                             ordering=FALSE,
                             paging=FALSE), 
                class = "cell-border stripe", 
                caption = htmltools::tags$caption(
                          style = 'caption-side: top; text-align: left;',
                          htmltools::em("Comparison of variance measures for choosing variable genes")), 
                escape = FALSE) %>% 
  formatStyle(
              "marker_gene",
              target="row",
              backgroundColor=styleEqual(desc[,"marker_gene"], colors))

saveWidget(out, "~/Downloads/desc.html")


