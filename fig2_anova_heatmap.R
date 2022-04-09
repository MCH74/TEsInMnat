library(pheatmap)
library(RColorBrewer)
library(tidyverse)
library(viridis)

countstable <- as.data.frame(read.table("z_score_anova.txt",sep="\t",header=TRUE))
countstable$TE<-NULL
colnames(countstable) <- c("Worker","Worker","Worker","Worker","T0","T0","T0","T0","T1","T1","T1","T1","T2","T2","T2","T3","T3","T3","T4","T4","T4","T4","King","King","King")

hmcol <- brewer.pal(11,"RdYlBu")
pdf("heatmap_anova_z_score.pdf",width=6, height=4)
pheatmap(as.matrix(countstable),clustering_distance_rows = dist(countstable,method="manhattan"),treeheight_row = 15,Rowv=FALSE,col=hmcol,mar=c(10,2),cexCol=0.5,legend_labels="Z-score",cluster_cols = F,show_rownames=F)
dev.off()

#for annotation of cluster
map<-pheatmap(as.matrix(countstable),clustering_distance_rows = dist(countstable,method="manhattan"),treeheight_row = 0,Rowv=FALSE,col=hmcol,mar=c(10,2),cexCol=0.5,legend_labels="Z-score",cluster_cols = F,show_rownames=F)

#import countstable again but do not remove TE column to retain the IDs for annotation
countstable <- as.data.frame(read.table("z_score_anova.txt",sep="\t",header=TRUE))
write.table(countstable[map$tree_row$order,],"anova_cluster_manhattan.txt")
