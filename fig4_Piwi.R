library(pheatmap)
library(viridis)

df<- as.data.frame(read.table("piwi.txt",sep="\t",header = T))
rownames(df)<-df$gene_id
df$gene_id<-NULL

hmcol <- viridis(10)
pdf("heatmap_piwi.pdf",width=6, height=4)
pheatmap(as.matrix(df),treeheight_row = 0,Rowv=FALSE,col=hmcol,mar=c(10,2),cexCol=0.5,legend_labels="score",cluster_cols = F,show_rownames=T)
dev.off()