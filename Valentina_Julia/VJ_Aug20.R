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
pdf("histmito_read_percentage.pdf", height = 4, width =  4)
hist(mito_read_percentage, main = "Cell quality assessment", col = "darkseagreen3", ylab = "Cell count", xlab = "mitochondrial read %") #as long as %mito reads / cell is < 5%, we consider the cells to be of good quality
text(x = 0.03, y = 800, cex = 0.6, "Result: All cells < 3% mito reads \n Conclusion: High quality", adj = 1)
dev.off()
#remove mitochondrial genes; we care about nuclear genes and mRNA
rna <- rna[-c(mito.genes_1, mito.genes_2),] 


total_mean_genexp <- rowSums(rna) / ncol(rna)
pdf("hist_avggeneexp.pdf", height = 4, width =  4)
hist(total_mean_genexp, main = "Total Gene Expression", col = "mistyrose3", ylab = "Cell count", xlab = "Mean log2 UMI count") #as long as %mito reads / cell is < 5%, we consider the cells to be of good quality
text(x = 0.03, y = 800, cex = 0.6, "Result: All cells < 3% mito reads \n Conclusion: High quality", adj = 1)
dev.off()

#Is there batch effect? (individual, plate, disease)
meta <- read.table("/home/tmpam10/SMTB_2019_students/data/scRNAseq/celseq_synovium_meta_5265cells_paper_v2.txt.gz", 
                   header = T, sep = " ")


#install.packages(library(metaMA))
# rv <- sapply (1:nrow(rna), function(x) var(rna[x,]))
# var_thresh <- sort(rv,decreasing = T)[1500]
# rna_mostvariablegenes <- rna[which(rv >= var_thresh),]
# pca <- prcomp(t(rna_mostvariablegenes), center = T, scale. = T)
# PCs <- pca$x
# plot(PCs[,1],PCs[,2])

#color by disease
# diseases <- unique(meta$disease)
# library(scales)
# rainbow(n = length(diseases))
# mycolors <- hue_pal()(length(diseases))
# plot(PCs[,1],PCs[,2], col = mycolors[match(meta$disease, diseases)])
# #color by person
# people <- unique(meta$sample)
# mycolors <- hue_pal()(length(people))
# pdf("personcolcells.pdf", width = 4, height = 4)
# plot(PCs[,1],PCs[,2], col = mycolors[match(meta$sample, people)])
# dev.off()
# #color by lane
# lane <- unique(meta$lane)
# 
# mycolors <- hue_pal()(length(lane))
# plot(PCs[,1],PCs[,2], col = mycolors[match(meta$lane, lane)])
# #different way to plot groups (compresses N PCs into 2D space)
# #install.packages("umap")
# library(umap)
# prc_umap <- umap(PCs[,1:30])
# plot(prc_umap$layout[,1],prc_umap$layout[,2]) #a lot more information than the 1st 2 PCs. 

###### OA versus RA ######
#1. split the matrix into two matrices: OA and RA
all(meta$cell_name == colnames(rna)) #make sure cells are in the same order in both matrices
rna_oa <- rna[,which(meta$disease == "OA")]
rna_ra <- rna[,which(meta$disease == "RA")]

#2. select different variable genes per OA or RA 
types <- c("ra","oa")
for (i in 1:length(types)){
  myrna <- get(paste0("rna_",types[i]))
  rv <- sapply(1:nrow(myrna), function(x) var(myrna[x,]))
  var_thresh <- sort(rv,decreasing = T)[1500]
  myrna_mostvariablegenes <- myrna[which(rv >= var_thresh),]
  
  #3. plot the PCA and UMAP, count number of distinct clusters. 
  pca <- prcomp(t(myrna_mostvariablegenes), center = T, scale. = T)
  PCs <- pca$x
  #plot(PCs[,1],PCs[,2],main = types[i])
  prc_umap <- umap(PCs[,1:30])
  assign(paste0("prc_umap_",types[i]),prc_umap)
  plot(prc_umap$layout[,1],prc_umap$layout[,2], main = types[i]) 
}

########################################################

types <- c("ra","oa")
for (i in 1:length(types)){
  
myrna <- get(paste0("rna_",types[i]))
rv <- sapply(1:nrow(myrna), function(x) var(myrna[x,]))
var_thresh <- sort(rv,decreasing = T)[1500]
myrna_mostvariablegenes <- myrna[which(rv >= var_thresh),]

x <- myrna_mostvariablegenes
cen <- 5
km <- kmeans(x = t(x), centers = cen, nstart = cen)
cvet <- hue_pal()(cen) #rainbow(n=cen)
ckm <- cvet[km$cluster]
#pdf(paste0(types[i],"_CLUSTJV.pdf"), width = 4, height = 4)
prc_umap <- get(paste0("prc_umap_",types[i]))
plot(prc_umap$layout[,1],prc_umap$layout[,2], main = types[i], col = ckm, xlab = "Dim 1 UMAP", ylab = "Dim 2 UMAP") 
#dev.off()
#i = 1
mygenes <- rownames(x)
markergenes <- c("PDGFRA","ISLR","CD2","CD3D","CD79A","RALGPS2","CD14","C1QA")
m <- match(markergenes, mygenes) #some of these were never in the list of most variable genes, so we are never going to find them. 
markergenes <- markergenes[!is.na(m)]
markercelltypes <- c("Fibroblast","Fibroblast","T","T","B","B","Monocyte","Monocyte")
markercelltypes <- markercelltypes[!is.na(m)]

pvalue_mat <- matrix(0,5,length(markergenes))
for (j in 1:length(unique(km$cluster))){ #iterate over clusters
  w <- which(km$cluster == j)
  for (k in 1:length(markergenes)){ #iterate over marker genes
    w_gene <- which(mygenes == markergenes[k])
    boxplot(x[w_gene,w], x[w_gene,-w])
    pvalue_mat[j,k] <- wilcox.test(x[w_gene,w], x[w_gene,-w], alternative = "greater")$p.val
  }
}
pvalue_mat <- round(pvalue_mat, digits = 2)
pvalue_mat <- ifelse(pvalue_mat < 0.05/(length(pvalue_mat)), 0, 1) #0 = significant, 1 = not significant
rownames(pvalue_mat) <- paste0("Cluster ",1:length(unique(km$cluster)))
colnames(pvalue_mat) <- paste0(markergenes, ":", markercelltypes)

n_cl <- km$cluster #cluster assignment
n_cl <- as.data.frame(n_cl)
n_cl$n <- rownames(n_cl)

bcell_clusters <- as.numeric(which(pvalue_mat[,which(markercelltypes == "B")] == 0))

#pdf(paste0(types[i],"_CLUSTJV.pdf"), width = 4, height = 4)
aa <- get(paste0("prc_umap_",types[i]))
plot(aa$layout[,1],aa$layout[,2], main = types[i], col = ckm, xlab = "Dim 1 UMAP", ylab = "Dim 2 UMAP") 
mycelltypes <- c()
for (cluster in 1:nrow(pvalue_mat)){
  ww <- which(pvalue_mat[cluster,] < 0.05/(nrow(pvalue_mat)*ncol(pvalue_mat)))
  mytext <- unique(markercelltypes[ww])
  mycelltypes[cluster] <- mytext
  meanx <- mean(aa$layout[km$cluster == cluster,1])
  meany <- mean(aa$layout[km$cluster == cluster,2])
  text(x = meanx, y = meany, mytext)
}
#dev.off()

#head(names_5k)
#dim(names_5k)
#length(names_5k)
##################################################################
for (j in 1:length(bcell_clusters)){
  b <- n_cl$n[n_cl$n_cl == bcell_clusters[j]]
  assign(paste0("b_",j),b)
  k_genes <- myrna_mostvariablegenes[,match(b, colnames(myrna_mostvariablegenes))]
  assign(paste0("k_genes_",j),k_genes)
}

thresh <- 0.05/nrow(k_genes)
p_value_1k <- sapply(1:nrow(k_genes), function(x)
  wilcox.test(k_genes_1[x,], k_genes_2[x,], alternative = "greater")$p.val)
names_1k <- rownames(myrna_mostvariablegenes[which(p_value_1k < thresh),])

p_value_2k <- sapply(1:nrow(k_genes), function(x)
  wilcox.test(k_genes_2[x,], k_genes_1[x,], alternative = "greater")$p.val)
names_2k <- rownames(myrna_mostvariablegenes[which(p_value_2k < thresh),])

#save(names_5k, "/home/tmpam10/Documents/name_2k.rda")
#ames_5k <- readRDS("/home/tmpam10/Documents/name_2k.rda")

write.table(names_1k, paste0("/home/tmpam10/Documents/name_1k_",types[i],".txt"), quote = FALSE)
write.table(names_2k, paste0("/home/tmpam10/Documents/name_2k_",types[i],".txt"), quote = FALSE)

IGgenes <- grep("^IGH",rownames(myrna_mostvariablegenes), value = T)
IGgene_count <- c()
for (j in 1:length(bcell_clusters)){
  m <- match(IGgenes, get(paste0("names_",j,"k")))
  IGgene_count[j] <- length(which(is.na(m)))
}
if (max(IGgene_count) == 0){
  plasmacluster <- NULL
}else{
  plasmacluster <- bcell_clusters[which(IGgene_count == max(IGgene_count))]
}

pdf(paste0(types[i],"_CLUSTJV_plasmablasts.pdf"), width = 4, height = 4)
aa <- get(paste0("prc_umap_",types[i]))
plot(aa$layout[,1],aa$layout[,2], main = toupper(types[i]), col = ckm, xlab = "Dim 1 UMAP", ylab = "Dim 2 UMAP") 
#mycelltypes[plasmacluster] <- "Plasmablasts"
mycelltypes[4] <- "B"
mycelltypes[3] <- "Plasmablasts"
for (cluster in 1:nrow(pvalue_mat)){
  meanx <- mean(aa$layout[km$cluster == cluster,1])
  meany <- mean(aa$layout[km$cluster == cluster,2])
  if (cluster == 3){
    text(x = meanx, y = meany, mycelltypes[cluster])
  }else{
    text(x = meanx, y = meany, mycelltypes[cluster], adj = 1)
  }
  
}
dev.off()

}
#skip PAX5 
#rerun OA
#for ra and oa separately, find out proportions of cell types. 

mycelltypes 
for (k in 1:length(mycelltypes)){
  print(length(which(km$cluster == k)) / length(km$cluster))
}







##############################################################################################33
kmo <- kmeans(x = t(x), centers = cen, nstart = cen)

mygenes <- rownames(x)
markergenes <- c("PDGFRA","ISLR","CD2","CD3D","CD79A","RALGPS2","CD14","C1QA")
m <- match(markergenes, mygenes) #some of these were never in the list of most variable genes, so we are never going to find them. 
markergenes <- markergenes[!is.na(m)]
markercelltypes <- c("Fibroblast","Fibroblast","T","T","B","B","Monocyte","Monocyte")
markercelltypes <- markercelltypes[!is.na(m)]
myrna <- get(paste0("rna_",type))
rv <- sapply(1:nrow(myrna), function(x) var(myrna[x,]))
var_thresh <- sort(rv,decreasing = T)[1500]
myrna_mostvariablegenes <- myrna[which(rv >= var_thresh),]


x <- myrna_mostvariablegenes
cen <- 5
kmo <- kmeans(x = t(x), centers = cen, nstart = cen)

tabkmoOA <- table(kmo$cluster)/length(kmo$cluster)
for (cluster in 1:nrow(pvalue_mat)){
  meanx <- mean(prc_umap$layout[km$cluster == cluster,1])
  meany <- mean(prc_umap$layout[km$cluster == cluster,2])
  text(x = meanx, y = meany, mycelltypes[cluster])
}
######################################################################################

i =1
myrna <- get(paste0("rna_",types[i]))
rv <- sapply(1:nrow(myrna), function(x) var(myrna[x,]))
var_thresh <- sort(rv,decreasing = T)[1500]
myrna_mostvariablegenes <- myrna[which(rv >= var_thresh),]



x <- myrna_mostvariablegenes
cen <- 5
km <- kmeans(x = t(x), centers = cen, nstart = cen)
tab1 <- table(km$cluster)/length(km$cluster)

rv <- sapply(1:nrow(myrna), function(x) var(myrna[x,]))
var_thresh <- sort(rv,decreasing = T)[1500]
myrna_mostvariablegenes <- myrna[which(rv >= var_thresh),]
x <- myrna_mostvariablegenes
cen <- 5
km <- kmeans(x = t(x), centers = cen, nstart = cen)
tab2 <- table(km$cluster)/length(km$cluster)
write.table(data.frame(t(tab2))[,3], "/home/tmpam10/Documents/tab2.txt", append = FALSE, quote = FALSE, row.names = mycelltypes, sep =" ")
b <- read.table("/home/tmpam10/Documents/tab2.txt")


write.table(data.frame(t(tab1))[,3], "/home/tmpam10/Documents/tab1.txt", append = FALSE, quote = FALSE, row.names = mycelltypes, sep =" ")
a <- read.table("/home/tmpam10/Documents/tab1.txt")


library(ggplot2)
mycelltypes[1] <- "Plasmablast"
#df <- data.frame(cbind(tab1,mycelltypes))
df <- data.frame(cell_proportion=as.numeric(tab1),cell_type=mycelltypes)

#df[1,2] <- "Plasmablast"
#colnames(df) <- c("cell_proportion", "cell_type")
#df1 <- data.frame(cbind(tab2,mycelltypes))
df1 <- data.frame(cell_proportion=as.numeric(tab2),cell_type=mycelltypes)

#df1[1,2] <- "Plasmablast"
#colnames(df1) <- c("cell_proportion", "cell_type")

df2 <- rbind(df, df1)
df2$desease <- c(rep("OA", 5), rep("RA", 5))
#df2[,1] <- sapply(as.numeric(levels(df2[,1])[df2[,1]]), function(x) round(x, 2))
#df0 <- as.data.frame(df2$mycelltypes) 
#df0$RA <- df2$tab1
#df0$OA <- df2$tab2
#colnames(df0) <- c("Cell types", "RA", "OA")
#df[,1] <- as.numeric(df[,1])
ggplot(df, aes(x=mycelltypes, y=tab1)) + geom_bar(stat="identity")
desease <- c("RA", "OA")
ggplot(df2, aes(x=cell_type, y=cell_proportion, fill=desease)) + geom_bar(stat="identity",  position=position_dodge()) + theme(text = element_text(size=20), axis.text.x = element_text(angle = 90))

                                                                                                                               