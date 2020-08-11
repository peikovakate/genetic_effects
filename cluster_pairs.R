library(tidyverse)
library(pheatmap)
library(RColorBrewer)
source("pipeline/visualisation_utils.R")

effects <- read_tsv("pipeline_data/effects_microarr.tsv")
effects <- read_tsv("../data/rnaseq/rnaseq_effects.tsv")
nrow(effects)

betas <- select(effects, starts_with('beta.')) %>% rename_all(function(x){sub("beta.", "", x)})
betas <- as.matrix(betas)
# SOMETHING IS WRONG
# betas <- betas*sign(betas[max.col(abs(betas))])

# betas <- t(apply(betas, 1, function(x){
#   abs_x <- abs(x)
#   x*sign(x[which.max(abs_x)])
# }))

# betas <- abs(betas)
# betas <- m$result$PosteriorMean
tissue_names <- colnames(betas)
    
# h = hclust(dist(t(betas)))
# plot(h)
# Pairwise correlation between tissues
m = "pearson"
cols.cor <- cor(betas, method = m)
hclust.col <- hclust(as.dist(1-cols.cor))
plot(hclust.col)
hmap(
  as.data.frame(cols.cor),
  fontSize = 7,
  colorPalette = 8,
  dendroLineSize = 0.2,
  scaleName = expression(paste("Pearson's r"))
  # main=paste("prior bitas", m)
)

pheatmap(cols.cor, cutree_cols = 7, colorRampPalette(brewer.pal(n = 7, name ="OrRd"))(100),
         filename="final_figures/rnaseq/pearson_raw_pheatmap.png", width=11, height=9)
# ggsave("final_figures/rnaseq/pearson_raw.png", scale=1, width = 7, height = 6)

# Pairwise correlation between rows eqtls
rows.cor <- cor(t(betas), method = "pearson")
hclust.row <- hclust(as.dist(1-rows.cor))
# plot(hclust.row)

# 
# dists <- dist(betas)
# tree <- hclust(cols.cor, method = "complete")
N_clusters = 10
cluster_cut <- cutree(hclust.row, N_clusters)

# sample <- sample_n(betas, 300)

lineplot_cluster <- function(cluster_ind, clusters, add_info){
  dir.create(paste0("final_figures/", add_info), recursive = T, showWarnings = F)
  c <- betas[clusters == cluster_ind,]
  size = dim(c)[1]
  print(dim(c))
  # tissue_names <- names(c[1, ])
  c <- aperm(as.array(c))
  
  png(sprintf("final_figures/%s/clt%i.png", add_info, cluster_ind), width=800, height = 800)
  par(mar=c(20,4,4,2)+0.1)
  # c <- aperm(sign(c) * c)
  matplot(
    c,
    type = "l",
    lwd = 0.1,
    lty = 1,
    main = paste("cluster", i, "size:", size, add_info),
    axes = F,
    ylab = "Prior betas",
    cex.axis=1.3, cex.lab=1.3
  )
  axis(2)
  axis(side=1,at=1:nrow(c),labels=tissue_names, las=2, cex.axis=1.3)
  dev.off()
}

for (i in 1:N_clusters) {
  lineplot_cluster(i, cluster_cut, "microarr/hierarchical_raw")
}

N_clusters = 8
clusters <- kmeans(betas, N_clusters, iter.max = 30, nstart = 5)
clusters$size


for (i in 1:N_clusters) {
  lineplot_cluster(i, clusters$cluster, "rnaseq/poster_k_means")
}

# which(clusters$size > 5)
# clusters$clusters
# new_betas <- betas[which(clusters$cluster %in% which(clusters$size > 5)), ]
# clusters <- kmeans(new_betas, 6, iter.max = 30, nstart = 5)

# pheatmap(betas, kmeans_k = 5, clustering_method="complete")
# 
# pca <- prcomp(betas)
# summary(pca)
# pca$x %>% 
#   as.data.frame %>%
#   ggplot(aes(x=PC3, y=PC4, color=tissue_names[max.col(betas)])) + geom_point(size=1) +
#   theme(legend.position="top")
# 
# t_betas <- t(betas)
# pca_tissues <- prcomp(t_betas)
# summary(pca_tissues)
# pca_tissues$x %>% 
#   as.data.frame %>%
#   ggplot(aes(x=PC1, y=PC2, color=rownames(t_betas))) + geom_point(size=3) +
#   geom_text(aes(label=rownames(t_betas)), position=position_jitter(width=1,height=1))+
#   theme(legend.position="top") +
#   labs(color="cell types")
