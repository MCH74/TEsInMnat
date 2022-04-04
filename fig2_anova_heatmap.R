library(ComplexHeatmap)
library(RColorBrewer)
#heatmap
countstable <- as.data.frame(read.table("z_score_anova.txt",sep="\t",header=TRUE,row.names=1))
col_anno <- as.data.frame(read.table("col_anno.txt",sep="\t",header=TRUE))
colnames(countstable) <- c("Worker","Worker","Worker","Worker","Q0","Q0","Q0","Q0","Q1","Q1","Q1","Q1","Q2","Q2","Q2","Q3","Q3","Q3","Q4","Q4","Q4","Q4","King","King","King")

colAnn <- HeatmapAnnotation(df = col_anno,
                            col = list("Caste" = c("King" = "#F8766D", "Q0" = "#C39900", "Q1" = "#53B400", "Q2" = "#06C196", "Q3" = "#19BDED", "Q4" = "#A387FF", "Worker" = "#FB55D4")),
                            which = 'col')

pdf(file="anova_heatmap_fig2.pdf")
Heatmap(as.matrix(countstable),
        col=brewer.pal(n=11,name="RdYlBu"),
        top_annotation = colAnn,
        show_row_names = F,
        cluster_columns = T)
dev.off()


