# k_5 <- n_cl$n[n_cl$n_cl == 5] #list of cells in cluster 5
# 
# cn <-colnames(myrna_mostvariablegenes)
# w_g <- which(mygenes == markergenes[5])
# 
# k_5_genes <- myrna_mostvariablegenes[,match(k_5, colnames(myrna_mostvariablegenes))]
# k_genes <- myrna_mostvariablegenes[,-match(k_5, colnames(myrna_mostvariablegenes))]
# 
# p_value_5k <- sapply(1:nrow(k_genes), function(x)
#   wilcox.test(k_5_genes[x,], k_genes[x,], alternative = "greater")$p.val)
#mulitple hypothesis correction
#we performed a statistical test 1000 times (i.e. number of genes)
#therefore we must make our threshold for significance stricter

# thresh <- 0.05 / nrow(k_genes)
# which(p_value_5k < thresh) #differentially expressed genes.
# 
# names_5k <- rownames(myrna_mostvariablegenes[which(p_value_5k < thresh),])
# write.table(names_5k, paste0("my_DEgenes_",type,".txt", row.names = F, col.names = F, quote = F))
# 
# pdf(paste0(type,"_CLUSTJV.pdf"), width = 4, height = 4)
# aa <- get(paste0("prc_umap_",type))
# plot(aa$layout[,1],aa$layout[,2], main = type, col = ckm, xlab = "Dim 1 UMAP", ylab = "Dim 2 UMAP") 
# mycelltypes <- c()
# for (cluster in 1:nrow(pvalue_mat)){
#   ww <- which(pvalue_mat[cluster,] < 0.05/(nrow(pvalue_mat)*ncol(pvalue_mat)))
#   mytext <- unique(markercelltypes[ww])
#   mycelltypes[cluster] <- mytext
#   meanx <- mean(aa$layout[km$cluster == cluster,1])
#   meany <- mean(aa$layout[km$cluster == cluster,2])
#   text(x = meanx, y = meany, mytext)
# }
# dev.off()
# 
# head(names_5k)
# dim(names_5k)
# length(names_5k)

##################################################################
# k_3 <- n_cl$n[n_cl$n_cl == 3]
# k_2 <- n_cl$n[n_cl$n_cl == 2]
# 
# k_3_genes <- myrna_mostvariablegenes[,match(k_3, colnames(myrna_mostvariablegenes))]
# k_2_genes <- myrna_mostvariablegenes[,match(k_2, colnames(myrna_mostvariablegenes))]
# 
# p_value_3k <- sapply(1:nrow(k_3_genes), function(x)
#   wilcox.test(k_3_genes[x,], k_genes[x,], alternative = "greater")$p.val)
# names_3k <- rownames(myrna_mostvariablegenes[which(p_value_3k < thresh),])
# 
# kgenes <- myrna_mostvariablegenes[,-match(k_2, colnames(myrna_mostvariablegenes))]
# p_value_2k <- sapply(1:nrow(k_2_genes), function(x)
#   wilcox.test(k_2_genes[x,], kgenes[x,], alternative = "greater")$p.val)
# names_2k <- rownames(myrna_mostvariablegenes[which(p_value_2k < thresh),])
# 
# save(names_2k, "/home/tmpam10/Documents/name_2k.rda")
# names_2k <- readRDS("/home/tmpam10/Documents/name_2k.rda")
# 
# write.table(names_3k, "/home/tmpam10/Documents/name_3k.txt", quote = FALSE)
# write.table(names_2k, "/home/tmpam10/Documents/name_2k.txt", quote = FALSE)
# common_genes <- intersect(names_2k,names_3k)
# variable_genes_names_2k <- names_2k[which(is.na(match(names_2k, common_genes)))]
# variable_genes_names_3k <- names_3k[which(is.na(match(names_3k, common_genes)))]

######################################################################################
# OA_2 <- variable_genes_names_2k
# OA_3 <- variable_genes_names_3k
# RA_1 <- names_DEgenes
# RA_3 <- read.table("/home/tmpam10/Documents/name_3k.txt")
# 
# allgenes <- read.table("/home/tmpam10/SMTB_2019_students/data/allgenes.txt", header = F, stringsAsFactors = F, sep = "\t")
# matOA1 <- match(OA_1$V2, allgenes$V4)
# matOA5 <-match(OA_5$V2, allgenes$V4)
# matRA1 <-match(RA_1$x, allgenes$V4)
# matRA3 <-match(RA_3$x, allgenes$V4)
# 
# tab_mat_OA1 <- allgenes[matOA1, c("V2", "V3", "V4")]

#########################################################################################

type <- "ra" #change between oa and ra 
myrna <- get(paste0("rna_",type))
rv <- sapply(1:nrow(myrna), function(x) var(myrna[x,]))
var_thresh <- sort(rv,decreasing = T)[1500]
myrna_mostvariablegenes <- myrna[which(rv >= var_thresh),]
#pca <- prcomp(t(myrna_mostvariablegenes), center = T, scale. = T)
#PCs <- pca$x
#prc_umap <- umap(PCs[,1:30])
#plot(prc_umap$layout[,1],prc_umap$layout[,2], main = type) 
#assign(paste0("prc_umap",type),prc_umap)

x <- myrna_mostvariablegenes
cen <- 5
km <- kmeans(x = t(x), centers = cen, nstart = cen)
cvet <- rainbow(n=cen)
ckm <- cvet[km$cluster]
#pdf(paste0(type,"_CLUSTJV.pdf"), width = 4, height = 4)
aa <- get(paste0("prc_umap_",type))
plot(aa$layout[,1],aa$layout[,2], main = type, col = ckm, xlab = "Dim 1 UMAP", ylab = "Dim 2 UMAP") 
#dev.off()
#i = 1
mygenes <- rownames(x)
markergenes <- c("PDGFRA","ISLR","CD2","CD3D","CD79A","RALGPS2","CD14","C1QA")
m <- match(markergenes, mygenes) #some of these were never in the list of most variable genes, so we are never going to find them. 
markergenes <- markergenes[!is.na(m)]
markercelltypes <- c("Fibroblast","Fibroblast","T","T","B","B","Monocyte","Monocyte")
markercelltypes <- markercelltypes[!is.na(m)]

pvalue_mat <- matrix(0,5,length(markergenes))
for (i in 1:length(unique(km$cluster))){ #iterate over clusters
  w <- which(km$cluster == i)
  for (j in 1:length(markergenes)){ #iterate over marker genes
    w_gene <- which(mygenes == markergenes[j])
    boxplot(x[w_gene,w], x[w_gene,-w])
    pvalue_mat[i,j] <- wilcox.test(x[w_gene,w], x[w_gene,-w], alternative = "greater")$p.val
  }
}
pvalue_mat <- round(pvalue_mat, digits = 2)
pvalue_mat <- ifelse(pvalue_mat < 0.05/(length(pvalue_mat)), 0, 1) #0 = significant, 1 = not significant
rownames(pvalue_mat) <- paste0("Cluster ",1:length(unique(km$cluster)))
colnames(pvalue_mat) <- paste0(markergenes, ":", markercelltypes)

#list of cells in cluster 5

#cn <-colnames(myrna_mostvariablegenes)
#w_g <- which(mygenes == markergenes[5])

#k_5_genes <- myrna_mostvariablegenes[,match(k_5, colnames(myrna_mostvariablegenes))]
#k_genes <- myrna_mostvariablegenes[,-match(k_1, colnames(myrna_mostvariablegenes))]

#p_value_5k <- sapply(1:nrow(k_genes), function(x)
# wilcox.test(k_5_genes[x,], k_genes[x,], alternative = "greater")$p.val)
#mulitple hypothesis correction
#we performed a statistical test 1000 times (i.e. number of genes)
#therefore we must make our threshold for significance stricter

#thresh <- 0.05 / nrow(k_genes)
#which(p_value_5k < thresh) #differentially expressed genes.

#names_5k <- rownames(myrna_mostvariablegenes[which(p_value_5k < thresh),])
#write.table(names_5k, paste0("my_DEgenes_",type,".txt", row.names = F, col.names = F, quote = F))

type <- "oa" #change between oa and ra 
myrna <- get(paste0("rna_",type))
rv <- sapply(1:nrow(myrna), function(x) var(myrna[x,]))
var_thresh <- sort(rv,decreasing = T)[1500]
myrna_mostvariablegenes <- myrna[which(rv >= var_thresh),]


x <- myrna_mostvariablegenes
cen <- 5
kmo <- kmeans(x = t(x), centers = cen, nstart = cen)
cvet <- rainbow(n=cen)
ckm <- cvet[km$cluster]
#pdf(paste0(type,"_CLUSTJV.pdf"), width = 4, height = 4)
aa <- get(paste0("prc_umap_",type))
plot(aa$layout[,1],aa$layout[,2], main = type, col = ckm, xlab = "Dim 1 UMAP", ylab = "Dim 2 UMAP") 
#dev.off()
#i = 1
mygenes <- rownames(x)
markergenes <- c("PDGFRA","ISLR","CD2","CD3D","CD79A","RALGPS2","CD14","C1QA")
m <- match(markergenes, mygenes) #some of these were never in the list of most variable genes, so we are never going to find them. 
markergenes <- markergenes[!is.na(m)]
markercelltypes <- c("Fibroblast","Fibroblast","T","T","B","B","Monocyte","Monocyte")
markercelltypes <- markercelltypes[!is.na(m)]

pvalue_mat <- matrix(0,5,length(markergenes))
for (i in 1:length(unique(km$cluster))){ #iterate over clusters
  w <- which(km$cluster == i)
  for (j in 1:length(markergenes)){ #iterate over marker genes
    w_gene <- which(mygenes == markergenes[j])
    boxplot(x[w_gene,w], x[w_gene,-w])
    pvalue_mat[i,j] <- wilcox.test(x[w_gene,w], x[w_gene,-w], alternative = "greater")$p.val
  }
}
pvalue_mat <- round(pvalue_mat, digits = 2)
pvalue_mat <- ifelse(pvalue_mat < 0.05/(length(pvalue_mat)), 0, 1) #0 = significant, 1 = not significant
rownames(pvalue_mat) <- paste0("Cluster ",1:length(unique(km$cluster)))
colnames(pvalue_mat) <- paste0(markergenes, ":", markercelltypes)

bcell_clusters <- as.numeric(which(pvalue_mat[,which(markercelltypes == "B")] == 0))

n_cl <- km$cluster #cluster assignment
n_cl <- as.data.frame(n_cl)
n_cl$n <- rownames(n_cl)

#list of cells in cluster 5





#p_value_5k <- sapply(1:nrow(k_genes), function(x)
# wilcox.test(k_5_genes[x,], k_genes[x,], alternative = "greater")$p.val)
#mulitple hypothesis correction
#we performed a statistical test 1000 times (i.e. number of genes)
#therefore we must make our threshold for significance stricter




#write.table(names_5k, paste0("my_DEgenes_",type,".txt", row.names = F, col.names = F, quote = F))

#pdf(paste0(type,"_CLUSTJV.pdf"), width = 4, height = 4)
aa <- get(paste0("prc_umap_",type))
plot(aa$layout[,1],aa$layout[,2], main = type, col = ckm, xlab = "Dim 1 UMAP", ylab = "Dim 2 UMAP") 
mycelltypes <- c()
for (cluster in 1:nrow(pvalue_mat)){
  ww <- which(pvalue_mat[cluster,] < 0.05/(nrow(pvalue_mat)*ncol(pvalue_mat)))
  mytext <- unique(markercelltypes[ww])
  mycelltypes[cluster] <- mytext
  meanx <- mean(aa$layout[km$cluster == cluster,1])
  meany <- mean(aa$layout[km$cluster == cluster,2])
  text(x = meanx, y = meany, mytext)
}

##################################################################
for (i in 1:length(bcell_clusters)){
  b <- n_cl$n[n_cl$n_cl == bcell_clusters[i]]
  assign(paste0("b_",i),b)
  k_genes <- myrna_mostvariablegenes[,match(b, colnames(myrna_mostvariablegenes))]
  assign(paste0("k_genes_",i),k_genes)
}

thresh <- 0.05/nrow(k_genes)
p_value_1k <- sapply(1:nrow(k_genes), function(x)
  wilcox.test(k_genes_1[x,], k_genes_2[x,], alternative = "greater")$p.val)
genes_pvalues_1 <- cbind(rownames(myrna_mostvariablegenes),p_value_1k)
names_1k <- rownames(myrna_mostvariablegenes[which(p_value_1k < thresh),])

p_value_2k <- sapply(1:nrow(k_genes), function(x)
  wilcox.test(k_genes_2[x,], k_genes_1[x,], alternative = "greater")$p.val)
genes_pvalues_2 <- cbind(rownames(myrna_mostvariablegenes),p_value_2k)
names_2k <- rownames(myrna_mostvariablegenes[which(p_value_2k < thresh),])


#IGgenes <- grep("^IGH",rownames(myrna_mostvariablegenes), value = T)
#boxplot(-log(as.numeric(genes_pvalues_1[IGgenes,2]),10), -log(as.numeric(genes_pvalues_2[IGgenes,2]),10))
#wilcox.test(-log(as.numeric(genes_pvalues_1[IGgenes,2]),10), -log(as.numeric(genes_pvalues_2[IGgenes,2]),10))

write.table(names_1k, paste0("/home/tmpam10/Documents/name_1k_",type,".txt"), quote = FALSE)
write.table(names_2k, paste0("/home/tmpam10/Documents/name_2k_",type,".txt"), quote = FALSE)

plasmacluster <- c()
IGgenes <- grep("^IGH",rownames(myrna_mostvariablegenes), value = T)
for (i in 1:length(bcell_clusters)){
  m <- match(IGgenes, get(paste0("names_",i,"k")))
  print(m)
  if(length(which(is.na(m))) > 0){
    plasmacluster <- bcell_clusters[i]
  }
}

pdf(paste0(type,"_CLUSTJV_plasmablasts.pdf"), width = 4, height = 4)
aa <- get(paste0("prc_umap_",type))
plot(aa$layout[,1],aa$layout[,2], main = type, col = ckm, xlab = "Dim 1 UMAP", ylab = "Dim 2 UMAP") 
mycelltypes[plasmacluster] <- "Plasmablasts"
for (cluster in 1:nrow(pvalue_mat)){
  meanx <- mean(aa$layout[km$cluster == cluster,1])
  meany <- mean(aa$layout[km$cluster == cluster,2])
  text(x = meanx, y = meany, mycelltypes[cluster])
}
dev.off()

##################################################################################################

rna_oa <-rna[,which(meta$disease == "OA")]

lenOA <- length(colnames(rna_oa))

type <- "oa" #change between oa and ra 
kmo <- kmeans(x = t(x), centers = cen, nstart = cen)

mygenes <- rownames(kmo)
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

kmo1 <- which(kmo$cluster == 1)
length(kmo1)
pros1 <- length(kmo1)/lenOA*100

kmo2 <- which(kmo$cluster == 2)
length(kmo2)
pros2 <- length(kmo2)/lenOA*100

kmo3 <- which(kmo$cluster == 3)
length(kmo3)
pros3 <- length(kmo3)/lenOA*100

kmo4 <- which(kmo$cluster == 4)
length(kmo4)
pros4 <- length(kmo4)/lenOA*100

kmo5 <- which(kmo$cluster == 5)
length(kmo5)
pros5 <- length(kmo5)/lenOA*100
######################################################################################
rna_ra <- rna[,which(meta$disease == "RA")]
lenRA <- length(colnames(rna_ra))

type <- "ra" #change between oa and ra 
myrna <- get(paste0("rna_",type))
rv <- sapply(1:nrow(myrna), function(x) var(myrna[x,]))
var_thresh <- sort(rv,decreasing = T)[1500]
myrna_mostvariablegenes <- myrna[which(rv >= var_thresh),]



x <- myrna_mostvariablegenes
cen <- 5
km <- kmeans(x = t(x), centers = cen, nstart = cen)
km1 <- which(km$cluster == 1)
length(km1)
pros1r <- length(km1)/lenRA*100



km2 <- which(km$cluster == 2)
length(km2)
pros2r <- length(km2)/lenRA*100

km3 <- which(km$cluster == 3)
length(km3)

pros3r <- length(km3)/lenRA*100

km4 <- which(km$cluster == 4)
length(km4)
pros4r <- length(km4)/lenRA*100

km5 <- which(km$cluster == 5)
length(km5)
pros5r <- length(km5)/lenRA*100