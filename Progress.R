
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
library(metaMA)
rv <- rowVars(rna)
var_thresh <- sort(rv,decreasing = T)[1000]
rna_mostvariablegenes <- rna[which(rv >= var_thresh),]
pca <- prcomp(t(rna_mostvariablegenes), center = T, scale. = T)
PCs <- pca$x
plot(PCs[,1],PCs[,2])

#color by disease
diseases <- unique(meta$disease)
mycolors <- hue_pal()(length(diseases))
plot(PCs[,1],PCs[,2], col = mycolors[match(meta$disease, diseases)])
#color by person
people <- unique(meta$sample)
mycolors <- hue_pal()(length(people))
plot(PCs[,1],PCs[,2], col = mycolors[match(meta$sample, people)])
#color by lane
lane <- unique(meta$lane)
mycolors <- hue_pal()(length(lane))
plot(PCs[,1],PCs[,2], col = mycolors[match(meta$lane, lane)])
#different way to plot groups (compresses N PCs into 2D space)
prc_umap <- umap(PCs[,1:30])
plot(prc_umap$layout[,1],prc_umap$layout[,2]) #a lot more information than the 1st 2 PCs. 

###### OA versus RA ######
#1. split the matrix into two matrices: OA and RA
all(meta$cell_name == colnames(rna)) #make sure cells are in the same order in both matrices
rna_oa <- rna[,which(meta$disease == "OA")]
rna_ra <- rna[,which(meta$disease == "RA")]

#2. select different variable genes per OA or RA 
type <- "oa" #change between oa and ra 
myrna <- get(paste0("rna_",type))
rv <- rowVars(myrna)
var_thresh <- sort(rv,decreasing = T)[1000]
myrna_mostvariablegenes <- myrna[which(rv >= var_thresh),]

#3. plot the PCA and UMAP, count number of distinct clusters. 
pca <- prcomp(t(myrna_mostvariablegenes), center = T, scale. = T)
PCs <- pca$x
plot(PCs[,1],PCs[,2],main = type)
prc_umap <- umap(PCs[,1:30])
plot(prc_umap$layout[,1],prc_umap$layout[,2], main = type) 





