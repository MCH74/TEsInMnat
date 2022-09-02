library(ggplot2)
library(dplyr)
library(ggokabeito)
library(viridis)

#bar chart TE class; castes
df <-as.data.frame(read.table("TE_families_caste.txt", sep = "\t", header =TRUE))
svg(filename="TE_families_caste.svg",width = 9, height = 4)
ggplot(df, aes(x=factor(Caste,levels=c("FW","Q0m","Q3m","Q9m","Q31m","Q20y","K20y")),y=Number, fill = factor(Element,levels = c("LINE","Other","DNA","LTR","RC")))) +
  geom_bar(position="dodge", stat="identity") +
  geom_text(label=df$Number, position=position_dodge(width=0.9), vjust=-0.25)+
  scale_fill_viridis()+
  xlab("Caste")+
  theme_classic()+
  theme(axis.text.y= element_blank(),
        axis.ticks.y=element_blank())
dev.off()

#Bar chart TE Class; Heatmap cluster
df <-as.data.frame(read.table("Heatmap_cluster_annotation_manhattan.txt", sep = "\t", header =TRUE))
svg(filename="TE_families_cluster_manhattan.svg",width = 12, height = 4)
ggplot(df, aes(x=factor(Caste,levels=c("FW","Q0m","Q3m","Q9m","Q31m","Q20y","K20y")),y=Number, fill = factor(Element,levels = c("LINE","Other","DNA","LTR","RC")))) +
  geom_bar(position="dodge", stat="identity") +
  geom_text(label=df$Number, position=position_dodge(width=0.9), vjust=-0.25)+
  scale_fill_viridis()+
  xlab("Caste")+
  theme_classic()+
  theme(axis.text.y= element_blank(),
        axis.ticks.y=element_blank())
dev.off()
+
  theme(axis.line = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        plot.title = element_text(hjust = 0.5, color = "#666666"))

#Heatmap cluster distance plot
df <-as.data.frame(read.table("heatmap_cluster_gene_ID_update.txt", sep = "\t", header =TRUE))
df <- na.omit(df)
pdf("anova_distance_manhattan_ito_log2_0.01.pdf",height = 2,width = 7)
ggplot(df[abs(df$Distance)<=1000000,],aes(factor(Caste,levels=c("FW","Q0m","Q3m","Q9m","Q31m","Q20y","K20y")),log2(as.numeric(abs(Distance+0.01)))))+
  geom_boxplot(aes(color=factor(Caste,levels=c("FW","Q0m","Q3m","Q9m","Q31m","Q20y","K20y"))))+
  scale_color_okabe_ito()+
  theme_minimal()+
  xlab("Caste")+
  ylab(expression(log2("Distance to gene [kb]")))+
  theme(legend.position = "none")
dev.off()

