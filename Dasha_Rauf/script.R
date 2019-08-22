library(umap)
#A summary of our tasks thusfar. 
"We have learned how to load gene expression data into R. Then we have learned how 
to perform quality control on the genes, removing those that are not expressed (e.g.)
zero total expression over the cells. We also prioritize nuclear mRNA and genes and therefore 
remove all mitochondrial genes (44 genes). Now that we have reduced the dimensionality of 
our data in terms of the number of genes, we can continue to perform dimensionality 
reduction. Next we learned how to compute the variance of gene expression of each gene across all cells. 
The idea is that the most variable genes will help us distinguish cell types. For example, 
if a gene has low variance, meaning that it is similarly expressed in all cells, it will 
not serve to help us distinguish cell types. Then, we learned to select the most variable genes- 
in pratice, my colleagues and I tend to select the top 1000 most variable genes. Feel free to modify 
this approach as you wish. Then, we learned to perform principal components analysis (PCA) to find
the major axes of variation and structure within our data. This allows us to visualize our data 
in two dimensions, for example the PC1 (the most important dimension) by PC2 (the next most important 
dimension). We can also use a cool algorithm called UMAP to plot more than 2 PCs in 2D space. The math
behind it is complicated, but it can compress the information of many PCs into a 2D representation.

Our next tasks, are to figure out how many candidate cell types are represented among the RA cells 
and among the OA cells, separately. Once you have separated the gene expression matrix into two matrices, 
1 for RA cells, and 1 for OA cells, you must separately look for variable genes, because the most variable genes
may differ between RA and OA cells. Therefore to do this step you must start with the matrix of gene expression
where the only dimensionality reduction that has been done is on unexpressed and mitochondrial genes. Then, 
you may run PCA and plot the first 2 PCs or rather run UMAP and plot the data based on the top 30, for example, 
PCs (as shown below.)

While we may visually hypothesize that there are X cell types for RA cells and Y cell types for OA cells, we need a 
way to figure out which cells belong to which groups. For this we will learn how to perform clustering. After clustering, 
we will learn how to hypothesize what cell type is represented by the cells from a given cluster, using 
differential gene expression analysis. After this task is done, each pair of students will work on separate tasks, regarding 
DNA (GWAS), protein (ChIP-seq), and enhancer (IMPACT) data."


################################################################################

#our progress analyzing scRNAseq data 
rna <- readRDS("/home/tmpam10/SMTB_2019_students/data/scRNAseq/celseq_synovium_log2_5265cells_paper.rds") #log2tpm already! 

#remove genes with 0 expression 
rs <- rowSums(rna)
rna <- rna[-which(rs == 0),]

#identify mitochondrial genes
mito.genes_1 <- grep(pattern = "^MT-", x = rownames(rna))
mito.genes_2 <- grep(pattern = "^MTRNR", x = rownames(rna))

#plot percentage of mitochdondrial reads per cell 
mito_read_percentage <- colSums(rna[c(mito.genes_1, mito.genes_2),]) / colSums(rna)
hist(mito_read_percentage) #as long as %mito reads / cell is < 5%, we consider the cells to be of good quality

#remove mitochondrial genes; we care about nuclear genes and mRNA
rna <- rna[-c(mito.genes_1, mito.genes_2),] 

#Is there batch effect? (individual, plate, disease)
meta <- read.table("/home/tmpam10/SMTB_2019_students/data/scRNAseq/celseq_synovium_meta_5265cells_paper_v2.txt.gz", 
                   header = T, sep = " ")
rv <- sapply(1:nrow(rna), function(x) var(rna[x,]))
#library(metaMA)
#rv <- rowVars(rna)
var_thresh <- sort(rv,decreasing = T)[1000]
rna_mostvariablegenes <- rna[which(rv >= var_thresh),]
pca <- prcomp(t(rna_mostvariablegenes), center = T, scale. = T)
PCs <- pca$x
plot(PCs[,1],PCs[,2])

#color by disease

diseases <- unique(meta$disease)
mycolors <- rainbow((length(diseases)))
plot(PCs[,1],PCs[,2], col = mycolors[match(meta$disease, diseases)])
#color by person
people <- unique(meta$sample)
mycolors <- rainbow((length(people)))
plot(PCs[,1],PCs[,2], col = mycolors[match(meta$sample, people)])
#color by lane
lane <- unique(meta$lane)
mycolors <- rainbow((length(lane)))
plot(PCs[,1],PCs[,2], col = mycolors[match(meta$lane, lane)])
#different way to plot groups (compresses N PCs into 2D space)
library(umap)
prc_umap <- umap(PCs[,1:30])
plot(prc_umap$layout[,1],prc_umap$layout[,2]) #a lot more information than the 1st 2 PCs. 

###### OA versus RA ######
#1. split the matrix into two matrices: OA and RA
all(meta$cell_name == colnames(rna)) #make sure cells are in the same order in both matrices
rna_oa <- rna[,which(meta$disease == "OA")]
rna_ra <- rna[,which(meta$disease == "RA")]

#2. select different variable genes per OA or RA 
type <- "ra" #change between oa and ra 
myrna <- get(paste0("rna_",type))
rv <- sapply(1:nrow(myrna), function(x) var(myrna[x,]))
rv <- rowVars(myrna)
var_thresh <- sort(rv,decreasing = T)[1500]
myrna_mostvariablegenes <- myrna[which(rv >= var_thresh),]

#3. plot the PCA and UMAP, count number of distinct clusters. 
pca <- prcomp(t(myrna_mostvariablegenes), center = T, scale. = T)
PCs <- pca$x
plot(PCs[,1],PCs[,2],main = type)
library(umap)
prc_umap <- umap(PCs[,1:30])
plot(prc_umap$layout[,1],prc_umap$layout[,2], main = type) 
#How many clusters for OA and RA?
prc_umap_ra <- prc_umap

#

markergenes <- c("PDGFRA", "ISLR", "CD2", "CD3D", "CD79A", "RALGPS2", "CD14", "C1QA")
markercelltypes <- c(rep("Fibroblast",2), rep("T cell", 2), rep("B cell", 2), rep("Monocyte"))
markergene_matrix <- cbind(markergenes, markercelltypes)
normgenes <- markergenes[markergenes %in% rownames(myrna_mostvariablegenes)]
markercelltypes <- markercelltypes[markergenes %in% rownames(myrna_mostvariablegenes)]
matrica <- matrix( nrow=length(normgenes), ncol=length(unique(kms)))
for (j in 1:length(unique(kms))) {
  cells <- which(kms ==j)
  for (i in 1:length(normgenes)) {
    need_wc_x <- myrna_mostvariablegenes[rownames(myrna_mostvariablegenes) == normgenes[i], cells]
    need_wc_y <- myrna_mostvariablegenes[rownames(myrna_mostvariablegenes) == normgenes[i], -cells]
    #WCTEST <- wilcox.test(x = as.numeric(need_wc_x), y = as.numeric(need_wc_y), alternative = "two.sided", mu = 0, paired = FALSE, exact = T, correct = TRUE, conf.int = T, conf.level = 0.95)
    WCTEST <- wilcox.test(x = as.numeric(need_wc_x), y = as.numeric(need_wc_y), alternative = "greater")
    matrica [i, j] <- WCTEST$p.value
  }
}
rownames(matrica) <- normgenes
colnames(matrica) <- paste0("cluster",1:length(unique(kms)))
apply(matrica, c(1,2), function(i) i<0.05/5000)

library(tiger)


red <- rgb(red = 1, green = 0, blue = 0, alpha = 0.5)
blue <- rgb(red = 0, green = 0, blue = 1, alpha = 0.5)
green <- rgb(red= 0, green = 1, blue = 0, alpha = 0.5)
col1 <- rgb(red=0.5, green = 0.5, blue = 0.5, alpha = 0.5)
col2 <- rgb(red=0.2, green = 0.2, blue = 0.2, alpha = 0.5)


listofcols <- rainbow(length(normgenes)) #c(red, green, blue, col1, col2)
pdf("~/gene_expression_plot.pdf", width = 10, height = 15)
par(mfrow = c(3,2))
plot(prc_umap_oa$layout[,1],prc_umap_oa$layout[,2], col = "gray", xlab = "UMAP 1", ylab = "UMAP 2", cex.axis = 1.5, cex.lab = 1.5, pch = 19)
mtext(text = paste0("All ",toupper(type)," Cells"), side = 3, cex = 2)
for (i in 1:length(normgenes)) {
  
  #plot(prc_umap_oa$layout[,1],prc_umap_oa$layout[,2], col = "white", xlab = "UMAP 1", ylab = "UMAP 2", cex.axis = 1.5, cex.lab = 1.5, pch = 19)

  cols <- color.factor(color = listofcols[i], myrna_mostvariablegenes[normgenes[i], ], max(myrna_mostvariablegenes[normgenes[i], ]))
  points(prc_umap_oa$layout[cols != "#FFFFFF",1],prc_umap_oa$layout[cols != "#FFFFFF",2], col = cols, xlab = "UMAP 1", ylab = "UMAP 2", cex.axis = 1.5, cex.lab = 1.5, pch = 19)
  mtext(text = paste0(normgenes[i], " ",markercelltypes[i]), side = 3, cex = 2)
  #legend_image <- as.raster(matrix(sort(cols,decreasing = T), ncol=1))
  #rasterImage(legend_image, xleft = -10, ybottom = 8, xright = -8, ytop = 10)
}
dev.off()




cols <- color.factor(listofcols[1], myrna_mostvariablegenes[normgenes[1], ], max(myrna_mostvariablegenes[normgenes[1], ]))
plot(prc_umap_oa$layout[,1],prc_umap_oa$layout[,2], col =cols)

par(new=TRUE)
cols <- color.factor(listofcols[3], myrna_mostvariablegenes[normgenes[3], ], max(myrna_mostvariablegenes[normgenes[3], ]))
plot(prc_umap_oa$layout[,1],prc_umap_oa$layout[,2], col =cols)

pdf("~/gene_expression_ggplot.pdf", width = 10, height = 15)
par(mfrow = c(3,2))

a <- as.data.frame(prc_umap_ra$layout)
a$V3 <- myrna_mostvariablegenes[normgenes[1],]
colnames(a) <- c("UMAP 1", "UMAP 2", "gene expression")

allcells <- ggplot(a, aes(x = `UMAP 1`, y = `UMAP 2`, colour = `gene expression`))+
  geom_point()+theme_dark()+
  scale_color_gradient(low = "#FFFFFF", high = "#FFFFFF")+
  ggtitle("All RA cells")+
  theme(plot.title = element_text(face = "bold", hjust = 0.5, size = 16))+
  annotate(geom = "text", x = -7, y = 0.5, label = "T-cells", color = "#FFFFFF", size = 6)+
  annotate(geom = "text", x = 7, y = 10, label = "Fibroblasts", color = "#FFFFFF", size = 6)+
  annotate(geom = "text", x = 8, y = -1, label = "Monocytes", color = "#FFFFFF", size = 6)+
  annotate(geom = "text", x = 2, y = -5, label = "B-cells", color = "#FFFFFF", size = 6)


  
PDGFRA <- ggplot(a, aes(x = `UMAP 1`, y = `UMAP 2`, colour = `gene expression`))+
  geom_point()+theme_dark()+
  scale_color_gradient(low = "#E6E6FA", high = "#6A5ACD")+
  ggtitle("Distribution of PDGFRA expression")+
  theme(plot.title = element_text(face = "bold", hjust = 0.5, size = 16))+
  annotate(geom = "text", x = -7, y = 0.5, label = "T-cells", color = "#FFFFFF", size = 6)+
  annotate(geom = "text", x = 7, y = 10, label = "Fibroblasts", color = "#FFFFFF", size = 6)+
  annotate(geom = "text", x = 8, y = -1, label = "Monocytes", color = "#FFFFFF", size = 6)+
  annotate(geom = "text", x = 2, y = -5, label = "B-cells", color = "#FFFFFF", size = 6)

a$`gene expression` <- myrna_mostvariablegenes[normgenes[2],]

ISLR <- ggplot(a, aes(x = `UMAP 1`, y = `UMAP 2`, colour = `gene expression`))+
  geom_point()+theme_dark()+
  scale_color_gradient(low = "#8BC8AA", high = "#007348")+
  ggtitle("Distribution of ISLR expression")+
  theme(plot.title = element_text(face = "bold", hjust = 0.5, size = 16))+
  annotate(geom = "text", x = -7, y = 0.5, label = "T-cells", color = "#FFFFFF", size = 6)+
  annotate(geom = "text", x = 7, y = 10, label = "Fibroblasts", color = "#FFFFFF", size = 6)+
  annotate(geom = "text", x = 8, y = -1, label = "Monocytes", color = "#FFFFFF", size = 6)+
  annotate(geom = "text", x = 2, y = -5, label = "B-cells", color = "#FFFFFF", size = 6)

 a$`gene expression` <- myrna_mostvariablegenes[normgenes[3],]
 CD2 <- ggplot(a, aes(x = `UMAP 1`, y = `UMAP 2`, colour = `gene expression`))+
   geom_point()+theme_dark()+
   scale_color_gradient(low = "#FAD3E0", high = "#7D1A3B")+
   ggtitle("Distribution of CD2 expression")+
   theme(plot.title = element_text(face = "bold", hjust = 0.5, size = 16))+
   annotate(geom = "text", x = -7, y = 0.5, label = "T-cells", color = "#FFFFFF", size = 6)+
   annotate(geom = "text", x = 7, y = 10, label = "Fibroblasts", color = "#FFFFFF", size = 6)+
   annotate(geom = "text", x = 8, y = -1, label = "Monocytes", color = "#FFFFFF", size = 6)+
   annotate(geom = "text", x = 2, y = -5, label = "B-cells", color = "#FFFFFF", size = 6)

 a$`gene expression` <- myrna_mostvariablegenes[normgenes[4],]

 CD79A <- ggplot(a, aes(x = `UMAP 1`, y = `UMAP 2`, colour = `gene expression`))+
  geom_point()+theme_dark()+
  scale_color_gradient(low = "#EAE6CA", high = "#9E9764")+
  ggtitle("Distribution of CD79A expression")+
  theme(plot.title = element_text(face = "bold", hjust = 0.5, size = 16))+
  annotate(geom = "text", x = -7, y = 0.5, label = "T-cells", color = "#FFFFFF", size = 6)+
  annotate(geom = "text", x = 7, y = 10, label = "Fibroblasts", color = "#FFFFFF", size = 6)+
  annotate(geom = "text", x = 8, y = -1, label = "Monocytes", color = "#FFFFFF", size = 6)+
  annotate(geom = "text", x = 2, y = -5, label = "B-cells", color = "#FFFFFF", size = 6)

a$`gene expression` <- myrna_mostvariablegenes[normgenes[5],]

 CD14 <- ggplot(a, aes(x = `UMAP 1`, y = `UMAP 2`, colour = `gene expression`))+
  geom_point()+theme_dark()+
  scale_color_gradient(low = "#FEDEBE", high = "#FD5602")+
  ggtitle("Distribution of CD14 expression")+
  theme(plot.title = element_text(face = "bold", hjust = 0.5, size = 16))+
  annotate(geom = "text", x = -7, y = 0.5, label = "T-cells", color = "#FFFFFF", size = 6)+
  annotate(geom = "text", x = 7, y = 10, label = "Fibroblasts", color = "#FFFFFF", size = 6)+
  annotate(geom = "text", x = 8, y = -1, label = "Monocytes", color = "#FFFFFF", size = 6)+
  annotate(geom = "text", x = 2, y = -5, label = "B-cells", color = "#FFFFFF", size = 6)
 ggarrange(allcells, CD14, CD2, CD79A, ISLR, PDGFRA, ncol = 2, nrow = 3)
 dev.off()
 
 
 #OA
 #2. select different variable genes per OA or RA 
 type <- "oa" #change between oa and ra 
 myrna_oa <- get(paste0("rna_",type))
 rv_oa <- sapply(1:nrow(myrna_oa), function(x) var(myrna_oa[x,]))
 var_thresh_oa <- sort(rv_oa,decreasing = T)[1500]
 myrna_mostvariablegenes_oa <- myrna_oa[which(rv_oa >= var_thresh_oa),]
 
 #3. plot the PCA and UMAP, count number of distinct clusters. 
 pca_oa <- prcomp(t(myrna_mostvariablegenes_oa), center = T, scale. = T)
 PCs_oa <- pca_oa$x
 library(umap)
 
 prc_umap_oa <- umap(PCs_oa[,1:30])
 plot(prc_umap_oa$layout[,1],prc_umap_oa$layout[,2], main = type) 
 #How many clusters for OA and RA?

 
 #k-means (k is a number of clusters we think that we have)
 km_ra <- kmeans(t(myrna_mostvariablegenes_oa), 5, iter.max = 100, algorithm = c("Lloyd"))
 length(km_ra$cluster)
 unique(km_ra$cluster)
 km_ra <- km_ra$cluster
 cell1_oa <- which(km_ra ==1)
 
 markergenes <- c("PDGFRA", "ISLR", "CD2", "CD3D", "CD79A", "RALGPS2", "CD14", "C1QA")
 markercelltypes <- c(rep("Fibroblast",2), rep("T cell", 2), rep("B cell", 2), rep("Monocyte"))
 markergene_matrix <- cbind(markergenes, markercelltypes)
 normgenes <- markergenes[markergenes %in% rownames(myrna_mostvariablegenes_oa)]
 markercelltypes <- markercelltypes[markergenes %in% rownames(myrna_mostvariablegenes_oa)]
 matrica <- matrix( nrow=length(normgenes), ncol=length(unique(km_ra)))
 for (j in 1:length(unique(km_ra))) {
   cells <- which(km_ra ==j)
   for (i in 1:length(normgenes)) {
     need_wc_x <- myrna_mostvariablegenes_oa[rownames(myrna_mostvariablegenes_oa) == normgenes[i], cells]
     need_wc_y <- myrna_mostvariablegenes_oa[rownames(myrna_mostvariablegenes_oa) == normgenes[i], -cells]
     #WCTEST <- wilcox.test(x = as.numeric(need_wc_x), y = as.numeric(need_wc_y), alternative = "two.sided", mu = 0, paired = FALSE, exact = T, correct = TRUE, conf.int = T, conf.level = 0.95)
     WCTEST <- wilcox.test(x = as.numeric(need_wc_x), y = as.numeric(need_wc_y), alternative = "greater")
     matrica [i, j] <- WCTEST$p.value
   }
 }
 rownames(matrica) <- normgenes
 colnames(matrica) <- paste0("cluster",1:length(unique(kms_oa)))
 apply(matrica, c(1,2), function(i) i<0.05/5000)
 
 #cluster 5 for OA is T cells 
 
 
 library(tiger)
 library(ggplot2)
 
 pdf("~/gene_expression_ggplot_oa.pdf", width = 10, height = 15)
 par(mfrow = c(3,2))
 
 a_oa <- as.data.frame(prc_umap_oa$layout)

 a_oa$V3 <- myrna_mostvariablegenes_oa[normgenes[1],]
 
 colnames(a_oa) <- c("UMAP 1", "UMAP 2", "gene expression")
 
 allcells_oa <- ggplot(a_oa, aes(x = `UMAP 1`, y = `UMAP 2`, colour = `gene expression`))+
   geom_point()+theme_dark()+
   scale_color_gradient(low = "#FFFFFF", high = "#FFFFFF")+
   ggtitle("All RA cells")+
   theme(plot.title = element_text(face = "bold", hjust = 0.5, size = 16))+
   annotate(geom = "text", x = -7, y = 0, label = "T-cells", color = "#FFFFFF", size = 6)+
   annotate(geom = "text", x = 10, y = 1, label = "Fibroblasts", color = "#FFFFFF", size = 6)+
   annotate(geom = "text", x = 6, y = 5, label = "Monocytes", color = "#FFFFFF", size = 6)+
   annotate(geom = "text", x = -8, y = 12.5, label = "B-cells", color = "#FFFFFF", size = 6)
 
 
 PDGFRA_oa <- ggplot(a_oa, aes(x = `UMAP 1`, y = `UMAP 2`, colour = `gene expression`))+
   geom_point()+theme_dark()+
   scale_color_gradient(low = "#E6E6FA", high = "#6A5ACD")+
   ggtitle("Distribution of PDGFRA expression")+
   theme(plot.title = element_text(face = "bold", hjust = 0.5, size = 16))+
   annotate(geom = "text", x = -7, y = 0, label = "T-cells", color = "#FFFFFF", size = 6)+
   annotate(geom = "text", x = 10, y = 1, label = "Fibroblasts", color = "#FFFFFF", size = 6)+
   annotate(geom = "text", x = 6, y = 5, label = "Monocytes", color = "#FFFFFF", size = 6)+
   annotate(geom = "text", x = -8, y = 12.5, label = "B-cells", color = "#FFFFFF", size = 6)
 
 a_oa$`gene expression` <- myrna_mostvariablegenes_oa[normgenes[2],]
 
 ISLR_oa <- ggplot(a_oa, aes(x = `UMAP 1`, y = `UMAP 2`, colour = `gene expression`))+
   geom_point()+theme_dark()+
   scale_color_gradient(low = "#8BC8AA", high = "#007348")+
   ggtitle("Distribution of ISLR expression")+
   theme(plot.title = element_text(face = "bold", hjust = 0.5, size = 16))+
   annotate(geom = "text", x = -7, y = 0, label = "T-cells", color = "#FFFFFF", size = 6)+
   annotate(geom = "text", x = 10, y = 1, label = "Fibroblasts", color = "#FFFFFF", size = 6)+
   annotate(geom = "text", x = 6, y = 5, label = "Monocytes", color = "#FFFFFF", size = 6)+
   annotate(geom = "text", x = -8, y = 12.5, label = "B-cells", color = "#FFFFFF", size = 6)
 
 
 a_oa$`gene expression` <- myrna_mostvariablegenes_oa[normgenes[3],]
 CD2_oa <- ggplot(a_oa, aes(x = `UMAP 1`, y = `UMAP 2`, colour = `gene expression`))+
   geom_point()+theme_dark()+
   scale_color_gradient(low = "#FAD3E0", high = "#7D1A3B")+
   ggtitle("Distribution of CD2 expression")+
   theme(plot.title = element_text(face = "bold", hjust = 0.5, size = 16))+
   annotate(geom = "text", x = -7, y = 0, label = "T-cells", color = "#FFFFFF", size = 6)+
   annotate(geom = "text", x = 10, y = 1, label = "Fibroblasts", color = "#FFFFFF", size = 6)+
   annotate(geom = "text", x = 6, y = 5, label = "Monocytes", color = "#FFFFFF", size = 6)+
   annotate(geom = "text", x = -8, y = 12.5, label = "B-cells", color = "#FFFFFF", size = 6)
 
 a_oa$`gene expression` <- myrna_mostvariablegenes_oa[normgenes[4],]
 
 CD79A_oa <- ggplot(a_oa, aes(x = `UMAP 1`, y = `UMAP 2`, colour = `gene expression`))+
   geom_point()+theme_dark()+
   scale_color_gradient(low = "#EAE6CA", high = "#9E9764")+
   ggtitle("Distribution of CD79A expression")+
   theme(plot.title = element_text(face = "bold", hjust = 0.5, size = 16))+
   annotate(geom = "text", x = -7, y = 0, label = "T-cells", color = "#FFFFFF", size = 6)+
   annotate(geom = "text", x = 10, y = 1, label = "Fibroblasts", color = "#FFFFFF", size = 6)+
   annotate(geom = "text", x = 6, y = 5, label = "Monocytes", color = "#FFFFFF", size = 6)+
   annotate(geom = "text", x = -8, y = 12.5, label = "B-cells", color = "#FFFFFF", size = 6)
 
 a_oa$`gene expression` <- myrna_mostvariablegenes_oa[normgenes[5],]
 
 CD14_oa <- ggplot(a_oa, aes(x = `UMAP 1`, y = `UMAP 2`, colour = `gene expression`))+
   geom_point()+theme_dark()+
   scale_color_gradient(low = "#FEDEBE", high = "#FD5602")+
   ggtitle("Distribution of CD14 expression")+
   theme(plot.title = element_text(face = "bold", hjust = 0.5, size = 16))+
   annotate(geom = "text", x = -7, y = 0, label = "T-cells", color = "#FFFFFF", size = 6)+
   annotate(geom = "text", x = 10, y = 1, label = "Fibroblasts", color = "#FFFFFF", size = 6)+
   annotate(geom = "text", x = 6, y = 5, label = "Monocytes", color = "#FFFFFF", size = 6)+
   annotate(geom = "text", x = -8, y = 12.5, label = "B-cells", color = "#FFFFFF", size = 6)
 
 library(ggpubr)
 ggarrange(allcells_oa, CD14_oa, CD2_oa, CD79A_oa, ISLR_oa, PDGFRA_oa, ncol = 2, nrow = 3)
 dev.off()

 #Plsamoblasts 
  
 
 #GWAS
 
 snitch <- read.table("/home/tmpam10/SMTB_2019_students/data/GWAS/PASS_Rheumatoid_Arthritis.sumstats.LDscores.gz", stringsAsFactors = F, header = T, sep = " ") 
 genesdata <- read.table("/home/tmpam10/SMTB_2019_students/data/allgenes.txt", stringsAsFactors = F, header = F, sep = "\t")

 
 #t-cell RA
 myrna_TRA <- which(km_ra == 2)
 myrna_TRAdata <- rna_ra[ , myrna_TRA]

 
 #t-cell OA
 myrna_TOA <- which(km_oa == 5)
 myrna_TOAdata <- rna_oa[ , myrna_TOA]

 #DE analysis
 #matrica_DE_RAOA <- sapply(1:nrow(myrna_TOAdata), function(x) wilcox.test(x = as.numeric(myrna_TRAdata[x, ]), y = as.numeric(myrna_TOAdata[x, ]), alternative = "greater")$p.val)
 # matrica_DE_RAOA <- matrix(nrow=length(rownames(myrna_TOAdata)), ncol=1)
 #   for (i in 1:length(rownames(myrna_TOAdata))) {
 #     need_wc_x_DE <- myrna_TRAdata[i, ]
 #     need_wc_y_DE <- myrna_TOAdata[i, ]
 #     #WCTEST <- wilcox.test(x = as.numeric(need_wc_x), y = as.numeric(need_wc_y), alternative = "two.sided", mu = 0, paired = FALSE, exact = T, correct = TRUE, conf.int = T, conf.level = 0.95)
 #     WCTEST_DE <- wilcox.test(x = as.numeric(need_wc_x_DE), y = as.numeric(need_wc_y_DE), alternative = "greater")
 #     matrica_DE_RAOA [i, 1] <- WCTEST_DE$p.value
 #   }

matrica_DE_RAOA <- sapply(1:nrow(rna), function(x) wilcox.test(x = as.numeric(rna_ra[x, ]), y = as.numeric(rna_oa[x, ]), alternative = "greater")$p.val) 
write.table(matrica_DE_RAOA, "matrica_DE_RAOA_allcells.txt", row.names = F, col.names = F, quote = F)

 
w_sig <- which(matrica_DE_RAOA < 0.05/nrow(rna)) 
m <- match(rownames(rna)[w_sig], genesdata$V4)
m_remove <- which(is.na(m))
m <- m[-m_remove]
perri <- genesdata[m,]
write.table(perri, "perri.txt", row.names = F, col.names = F, quote = F)
 #1396 DE genes
#perri <- read.table("/home/tmpam10/perri.txt", stringsAsFactors = F, header = F, sep = " ") 
#perri <- perri[-which(perri$V1 == "chrY"),]


# rownames(matrica_DE_RAOA) <- rownames(myrna_TRAdata) 
# colnames(matrica_DE_RAOA) <- 1
# TRUEMAtrix <-apply(matrica_DE_RAOA, c(1,2), function(i) i<0.05/nrow(myrna_TOAdata))
# 
# truetrue <- TRUEMAtrix[TRUEMAtrix == 1,]
# m <- match(names(truetrue), genesdata$V4)
# perri <- genesdata[m,]


#GWAS


#perri <- na.omit(perri)

#bezNA <- perri$V1[!(is.na(perri$V1))]

# snit <- list()
# 
# snit_r <- list()
# snit_govno <- list()
# snit_1 <- list()
# snit_2 <- list()
# 
#  for (i in c(1:(length(unique(snitch$CHR))-1), "X")) {
#    ind <- i
#    snit[[ind]] <- snitch[snitch$CHR == ind, "BP"]
#    snit_r[[ind]] <- snitch[snitch$CHR== ind, "SNP"]
#    snit_1[[ind]] <- snitch[snitch$CHR == ind, "Z"]
#    snit_2[[ind]] <- snitch[snitch$CHR== ind, "LDscore"]
#    snit_govno[[ind]] <- cbind(snit[[ind]], snit_r[[ind]], snit_1[[ind]], snit_2[[ind]])
#  } 

pir <- list() #start
pir_r <- list() #end
pir_str <- list() #strand
ka <- list() #start and end
for (i in c(1:(length(unique(perri$V1))-1), "X")) {
  ind <- paste0("chr", i, "")
  pir[[ind]] <- perri[perri$V1 == ind, "V2"]
  pir_r[[ind]] <- perri[perri$V1== ind, "V3"]
  pir_str[[ind]] <- perri[perri$V1== ind, "V6"]
  ka[[ind]] <- cbind(pir[[ind]], pir_r[[ind]], pir_str[[ind]])
}

gene_SNPs_chisquared <- c()
count <- 1
for(j in c(1:(length(unique(perri$V1))-1), "X")) {
  print(j)
  ind <- paste0("chr", j, "")
  table <- ka[[ind]] #genes starts and ends
  if(nrow(table) > 0){
    snitch_small <- snitch[snitch$CHR == j,]
    for(i in 1:nrow(table)) {
      #name <- paste0(ind, "_", table[i, 1], "_", table[i, 2])
      if (table[i, 3] == "+") { 
        window_plus <- which(snitch_small$BP > (as.numeric(table[i, 1])-2000) & snitch_small$BP < (as.numeric(table[i, 1])))
        strvar <- window_plus
      }
      else {
        window_minus <- which(snitch_small$BP > (as.numeric(table[i, 2])) & snitch_small$BP < (as.numeric(table[i, 2])+2000))
        strvar <- window_minus
      }
      chisquared <- max(((snitch_small$Z[strvar])^2)/snitch_small$LDscore[strvar])
      gene_SNPs_chisquared[count] <- chisquared
      count <- count + 1
    }
  }
}

wInf <- which(gene_SNPs_chisquared == "-Inf")
if (length(wInf)> 0){
  gene_SNPs_chisquared <- gene_SNPs_chisquared[-wInf]
}
mean(gene_SNPs_chisquared, na.rm = T)
write.table(gene_SNPs_chisquared, "/home/tmpam10/SMTB_2019_students/gene_SNPs_chisquared_prom.txt", quote = F, row.names = F, col.names = F)

######
# genesdata_nonDE <- genesdata
# DEgenes <- perri$V4
# for (i in 1:length(DEgenes)){
#   w <- which(genesdata$V4 == DEgenes[i])
#   genesdata_nonDE <- genesdata_nonDE[-w,]
# }
# genesdata_nonDE <- genesdata_nonDE[-which(nchar(genesdata_nonDE$V1) > 5),]
# genesdata_nonDE <- genesdata_nonDE[-which(genesdata_nonDE$V1 == "chrY"),]
# 
# nonDEgenes <- genesdata_nonDE[sample(1:nrow(genesdata_nonDE), size = nrow(perri), replace = F),]
# perri_notsig_ds <- nonDEgenes
# # w_not_sig <- which(matrica_DE_RAOA > 0.9) #9977
# # m <- match(rownames(myrna_TOAdata)[w_not_sig], genesdata$V4)
# # m_remove <- which(is.na(m))
# # m <- m[-m_remove]
# # perri_notsig <- genesdata[m,]
# # write.table(perri_notsig, "/home/tmpam10/SMTB_2019_students/perri_notsig.txt", row.names = F, col.names = F, quote = F)
# # perri_notsig <- perri_notsig[-which(perri_notsig$V1 == "chrY"),]
# # #perri_notsig <- perri_notsig[-which(nchar(perri_notsig$V1) > 5),] 
# # #perri_notsig_ds <- perri_notsig
# # perri_notsig_ds <- perri_notsig[sample(1:nrow(perri_notsig), size = nrow(perri), replace = F),]
# pir <- list() #start
# pir_r <- list() #end
# ka <- list() #start and end
# for (i in c(1:(length(unique(perri_notsig_ds$V1))-1), "X")) {
#   ind <- paste0("chr", i, "")
#   pir[[ind]] <- perri_notsig_ds[perri_notsig_ds$V1 == ind, "V2"]
#   pir_r[[ind]] <- perri_notsig_ds[perri_notsig_ds$V1== ind, "V3"]
#   pir_str[[ind]] <- perri_notsig_ds[perri_notsig_ds$V1== ind, "V6"]
#   ka[[ind]] <- cbind(pir[[ind]], pir_r[[ind]], pir_str[[ind]])
# }

mean_overtrials <- c()
for (tr in 1:100){
  print(paste0("trials_",tr))
  nonDEgenes <- genesdata_nonDE[sample(1:nrow(genesdata_nonDE), size = nrow(perri), replace = F),]
  perri_notsig_ds <- nonDEgenes
  
  pir <- list() #start
  pir_r <- list() #end
  ka <- list() #start and end
  for (i in c(1:(length(unique(perri_notsig_ds$V1))-1), "X")) {
    ind <- paste0("chr", i, "")
    pir[[ind]] <- perri_notsig_ds[perri_notsig_ds$V1 == ind, "V2"]
    pir_r[[ind]] <- perri_notsig_ds[perri_notsig_ds$V1== ind, "V3"]
    pir_str[[ind]] <- perri_notsig_ds[perri_notsig_ds$V1== ind, "V6"]
    ka[[ind]] <- cbind(pir[[ind]], pir_r[[ind]], pir_str[[ind]])
  }
  
  gene_SNPs_chisquared_negative <- c()
  count <- 1
  for(j in c(1:(length(unique(perri_notsig_ds$V1))-1), "X")) {
    print(j)
    ind <- paste0("chr", j, "")
    table <- ka[[ind]] #genes starts and ends
    if(nrow(table) > 0){
      snitch_small <- snitch[snitch$CHR == j,]
      for(i in 1:nrow(table)) {
        #name <- paste0(ind, "_", table[i, 1], "_", table[i, 2])
        if (table[i, 3] == "+") { 
          window_plus <- which(snitch_small$BP > (as.numeric(table[i, 1])-2000) & snitch_small$BP < (as.numeric(table[i, 1])))
          strvar <- window_plus
        }
        else {
          window_minus <- which(snitch_small$BP > (as.numeric(table[i, 2]))  & snitch_small$BP < (as.numeric(table[i, 2])+2000) )
          strvar <- window_minus
        }
        chisquared <- max(((snitch_small$Z[strvar])^2)/snitch_small$LDscore[strvar])
        gene_SNPs_chisquared_negative[count] <- chisquared
        count <- count + 1
      }
    }
  }
  
  wInf <- which(gene_SNPs_chisquared_negative == "-Inf")
  if (length(wInf)> 0){
    gene_SNPs_chisquared_negative <- gene_SNPs_chisquared_negative[-wInf]
  }
  
  mean_overtrials[tr] <- mean(gene_SNPs_chisquared_negative, na.rm = T)
}



#write.table(gene_SNPs_chisquared_negative, "/home/tmpam10/SMTB_2019_students/gene_SNPs_chisquared_negative_prom.txt", quote = F, row.names = F, col.names = F)

obs <- mean(gene_SNPs_chisquared, na.rm = T)
pdf("permutation_promoters_small.pdf", height = 4, width = 4)
hist(mean_overtrials, col = "lavender", main = "Disease associations in RA specific genes", ylab = "Number of Genes", xlab = "max SNP chi-squared statistic / gene")
abline(v = obs, col = "darkred", lwd = 2)
p <- length(which(mean_overtrials >= obs)) / (length(mean_overtrials)+1)
text(0.065, 15, "p = 0.37")
dev.off()
#gene_SNPs_chisquared_negative <- gene_SNPs_chisquared_negative[-which(gene_SNPs_chisquared_negative == "-Inf")]

#t.test(gene_SNPs_chisquared, gene_SNPs_chisquared_negative, alternative = "greater")

#boxplot(gene_SNPs_chisquared, gene_SNPs_chisquared_negative)
