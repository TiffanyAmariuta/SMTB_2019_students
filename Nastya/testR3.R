#simulation of gene expression data 

ngenes <- 5
ncells <- 5 
#expression matrix
myvalues <- sample(x = 1:100, size = 1000, replace = T)
hist(myvalues)
myvalues <- sample(x = 1:100, size = 1000, replace = T,
                   prob = 1/(1:100))
hist(myvalues)

myvalues <- sample(x = 1:100, size = ngenes*ncells, 
                   replace = T,
                   prob = 1/(1:100))
fakegeneexp <- matrix(myvalues,nrow = ngenes, 
                      ncol = ncells)

#goal: want total gene expression from each cell to be equal 
rs <- rowSums(t(fakegeneexp)) #sum per cell (column)
#cell sums * ? = 1e6
multiplication_factor <- 1e6/rs
fakegeneexp_cpm <- sapply(1:ncol(fakegeneexp), function(x) fakegeneexp[,x]*multiplication_factor[x])
colSums(fakegeneexp_cpm) #all 1 million
# log normalization to handle outlier gene expression
fakegeneexp_cpm_log2 <- log(x = fakegeneexp_cpm + 1, base = 2)

pdf("~/myfile.pdf")
hist(myvalues)
dev.off()


#load data
realdata <- readRDS("/home/tmpam10/SMTB_2019_students/data/scRNAseq/
                    celseq_synovium_log2_5265cells_paper.rds")

#find genes with zero expression across cells
rs <- rowSums(realdata) #sum of each gene's expression over cells
w <- which(rs == 0) #which equal 0?
realdata_update <- realdata[-w,]
#mitochrondrial gene expression is a sign of poor cell viability (dying cell)
#locate mito genes, determine if any cells are of poor quality 
gene_names <- rownames(realdata_update)
head(gene_names)
#which are mitochrondrial genes?
g1 <- grep(pattern = "^MT-", x = gene_names) # ^ indicates that the word must start with the following characters
g2 <- grep(pattern = "^MTRNR", x = gene_names)
g <- c(g1,g2) #contains locations of our mitochondrial genes 
#are there any cells where mito gene exp is > nuclear gene exp?
i <- 1 #for example, take cell 1 
sum(realdata_update[g,i]) / sum(realdata_update[,i]) 
#sum mitochondrial genes / sum total gene exp

#for all cells
cs <- colSums(realdata_update)
total_mito_exp <- sapply(1:ncol(realdata_update)
                         , function(x) 
                           sum(realdata_update[g,x]))
mito_percentages <- total_mito_exp/cs
hist(mito_percentages)
mean(mito_percentages)
max(mito_percentages)
#now, we know that our cells are of good quality. 

#remove mito genes
realdata_update <- realdata_update[-g,]

dim(realdata_update) #currently, 29,121 genes, 5,265

#next task- how many genes are expressed in each cell?
i <- 1
length(which(realdata_update[,i] > 0)) #num genes expressed in cell 1
s <- sapply(1:ncol(realdata_update), function(i) length(which(realdata_update[,i] 
                                                              > 0)))
hist(s)
min(s)
max(s)

#filter out lowly variable genes 
#why? because highly variable genes will help us distinguish cell types.
i <- 1
sd(realdata_update[i,])
sds_overgenes <- sapply(1:nrow(realdata_update), 
                        function(x) sd(realdata_update[x,]))
hist(sds_overgenes)

sd_threshold <- sort(sds_overgenes, decreasing = T)[1000]
keep_these_genes <- which(sds_overgenes > sd_threshold)
realdata_update_vargenes <- realdata_update[keep_these_genes,]
#there are many ways to do this step, 
#you can preselect a threshold for sd of gene expression 
#you can preselect a threshold for mean and sd of gene expression, etc 
#later, you can update this step and see how gene selection changes

#more dimensionality reduction. 
#out strategy- principal components analysis (PCA)

pca_output <- prcomp(x = t(realdata_update_vargenes), 
                     center = T, scale. = T)
dim(pca_output$x) 
dim(t(realdata_update_vargenes))
pca_output$x[1:5,1:5]
PCs <- pca_output$x
plot(x = PCs[,1],y = PCs[,2], 
     xlab = "PC 1", ylab = "PC 2")
#question- is the separation of cells that we observe, 
#due to biology or technical artifact? 

#load meta data (info about each cell)
metadata <- read.table("/home/tmpam10/SMTB_2019_students/data/scRNAseq/celseq_synovium_meta_5265cells_paper_v2.txt.gz", 
                       header = T, sep = " ")
head(metadata)
metadata_small <- metadata[,c(1,3,4,5)]
#check that cells are ordered in the same way between matrices
all(metadata_small[,1] == rownames(PCs))
#goal: color PCA plot by disease
w_oa <- which(metadata_small$disease == "OA")
mycolors <- rep(0,nrow(PCs))
mycolors[w_oa] <- rgb(red = 0.2, green = 0.2, blue = 0.9, alpha = 0.5)
mycolors[-w_oa] <- rgb(red = 0.9, green = 0.2, blue = 0.2, alpha = 0.5)
plot(x = PCs[w_oa,1],y = PCs[w_oa,2], col = mycolors[w_oa],
     xlab = "PC 1", ylab = "PC 2", pch = 19)
plot(x = PCs[-w_oa,1],y = PCs[-w_oa,2], col = mycolors[-w_oa],
     xlab = "PC 1", ylab = "PC 2", pch = 19)
points(x = PCs[w_oa,1],y = PCs[w_oa,2], col = mycolors[w_oa],
       xlab = "PC 1", ylab = "PC 2", pch = 19)









#goal: color PCA plot by sample
> w_oa <- which(metadata_small$disease == "OA")
> mycolors <- rep(0,nrow(PCs))
> mycolors[w_oa] <- rgb(red = 0.2, green = 0.2, blue = 0.9, alpha = 0.5)
> mycolors[-w_oa] <- rgb(red = 0.9, green = 0.2, blue = 0.2, alpha = 0.5)
> plot(x = PCs[w_oa,1],y = PCs[w_oa,2], col = mycolors[w_oa],
       +      xlab = "PC 1", ylab = "PC 2", pch = 19)
> plot(x = PCs[-w_oa,1],y = PCs[-w_oa,2], col = mycolors[-w_oa],
       +      xlab = "PC 1", ylab = "PC 2", pch = 19)
> points(x = PCs[w_oa,1],y = PCs[w_oa,2], col = mycolors[w_oa],
         +        xlab = "PC 1", ylab = "PC 2", pch = 19)



www=unique(metadata_small$sample)
install.packages("scales")
library(scales)

mycolors <- hue_pal()(length(www)) #a vector of colors
plot(PCs[which(metadata_small$sample == www[1]),1],PCs[which(metadata_small$sample == www[1]),2], col = mycolors[1])
points(PCs[which(metadata_small$sample == www[2]),1],PCs[which(metadata_small$sample == www[2]),2], col = mycolors[2])

for (i in 1:length(www)){
  if (i == 1){plot(PCs[which(metadata_small$sample == www[i]),1],PCs[which(metadata_small$sample == www[i]),2], col = mycolors[i])
  }else{points(PCs[which(metadata_small$sample == www[i]),1],PCs[which(metadata_small$sample == www[i]),2], col = mycolors[i])
}
  
}

#try a loop



aaa <- unique(metadata_small$sample)
install.packages("scales")
library(scales)

mycolors <- hue_pal()(length(www)) #a vector of colors
plot(PCs[which(metadata_small$sample == www[1]),1],PCs[which(metadata_small$sample == www[1]),2], col = mycolors[1])
points(PCs[which(metadata_small$sample == www[2]),1],PCs[which(metadata_small$sample == www[2]),2], col = mycolors[2])

for (i in 1:length(www)){
  if (i == 1){plot(PCs[which(metadata_small$sample == www[i]),1],PCs[which(metadata_small$sample == www[i]),2], col = mycolors[i])
  }else{points(PCs[which(metadata_small$sample == www[i]),1],PCs[which(metadata_small$sample == www[i]),2], col = mycolors[i])}}




aaa <- unique(metadata_small$plate)

for (i in 1:length(aaa)){
  if (i == 1){plot(PCs[which(metadata_small$plate == aaa[i]),1],PCs[which(metadata_small$plate== aaa[i]),2], col = mycolors[i])
  }else{points(PCs[which(metadata_small$plate == aaa[i]),1],PCs[which(metadata_small$plate == aaa[i]),2], col = mycolors[i])
  }
  
}


