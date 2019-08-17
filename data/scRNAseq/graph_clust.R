a = readRDS("~/R/SMTB_2019_students/data/scRNAseq/shrinked_table.rds")

range = readRDS("~/R/SMTB_2019_students/data/scRNAseq/16power_median_rate.rds")
grep("PAX5", rownames(a)[range])
length(range)


sign_genes = 1000
sign = a[range[1:sign_genes],]

pca = prcomp(t(sign), center = T, scale. = T)
med4_pca = prcomp(t(med4_sign), center = T, scale. = T)
nmed4_pca = prcomp(t(nmed4_sign), center = T, scale. = T)
med8_pca = prcomp(t(med8_sign), center = T, scale. = T)
med16_pca = prcomp(t(med16_sign), center = T, scale. = T)

rm(a)
rm(sd_sign)
rm(med8_sign)
count_rate30 = function(devs) {
  return(sqrt((sum(devs[1:30]^2))/(sum(devs^2))))
}

count_rate30(pca$sdev)
count_rate30(sd_pca$sdev)
sd_pca = as.data.frame(sd_pca$x)
med4_pca = as.data.frame(med4_pca$x)
med8_pca = as.data.frame(med8_pca$x)
count_rate30(med8_pca$sdev)
sign_pc = 30
sd_pca_cut = sd_pca[,1:sign_pc]
med4_pca_cut = med4_pca[,1:sign_pc]
med8_pca_cut = med8_pca[,1:sign_pc]

sd_dist = as.matrix(dist(sd_pca_cut))
med4_dist = as.matrix(dist(med4_pca_cut))
med8_dist = as.matrix(dist(med8_pca_cut))

joint_rate = 0.2
cells = 5265

joint_rates = c(0.001, 0.002, 0.003, 0.004, seq(0.005, 0.7, 0.1), 0.75, 0.8, 0.9)
clusts = readRDS("~/R/SMTB_2019_students/data/scRNAseq/first13_joint_rates.rds")
clusts[14] = 1
for (joint_rate_i in joint_rates[14:14]) {
  sd_graph = graph.adjacency(make.kNNG(sd_dist, k=cells*joint_rate_i))
  clusts = c(clusts, length(cluster_infomap(sd_graph)))
}
saveRDS(clusts, "~/R/SMTB_2019_students/data/scRNAseq/first13_joint_rates.rds")
plot(joint_rates[1:13], clusts)
clusts
library("loe")
sd_graph = graph.adjacency(make.kNNG(sd_dist, k=cells*joint_rate))
med4_graph = graph.adjacency(make.kNNG(med4_dist, k=cells*joint_rate))
med8_graph = graph.adjacency(make.kNNG(med8_dist, k=cells*joint_rate))

sd_communities = cluster_infomap(sd_graph)
length(sd_communities)
med4_communities = cluster_infomap(med4_graph)
med8_communities = cluster_infomap(med8_graph)
    
sd_communities[1]
  