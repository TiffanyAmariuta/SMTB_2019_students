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





realdata <- readRDS("/home/tmpam10/SMTB_2019_students/data/scRNAseq/celseq_synovium_log2_5265cells_paper.rds")
dim(realdata) #view dimensions, genes by cells
realdata[1:5,1:5] #view some elements

length(which(realdata == 0)) / (nrow(realdata)*ncol(realdata))
#result 92% sparse.
length(which(rowSums(realdata) > 0))




