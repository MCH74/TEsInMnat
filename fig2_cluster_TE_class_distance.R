#bar chart TE class; castes
df <-as.data.frame(read.table("TE_families_caste.txt", sep = "\t", header =TRUE))
svg(filename="TE_families_caste.svg",width = 9, height = 4)
ggplot(df, aes(x=factor(Caste,levels = c("Worker","Q0","Q1","Q2","Q3","Q4","King")),y=Number, fill = factor(Element,levels = c("LINE","Other","DNA","LTR","RC")))) +
  geom_bar(position="dodge", stat="identity") +
  geom_text(label=df$Number, position=position_dodge(width=0.9), vjust=-0.25)+
  scale_fill_brewer(palette="RdYlGn",name="Element")+
  xlab("Caste")+
  theme_classic()+
  theme(axis.text.y= element_blank(),
        axis.ticks.y=element_blank())
dev.off()

#Bar chart TE Class; Heatmap cluster
df <-as.data.frame(read.table("Heatmap_cluster_annotation_manhattan.txt", sep = "\t", header =TRUE))
svg(filename="TE_families_cluster_manhattan.svg",width = 12, height = 4)
ggplot(df, aes(x=factor(Caste,levels = c("Worker","Q0","Q1","Q2","Q3","Q4","King")),y=Number, fill = factor(Element,levels = c("LINE","Other","DNA","LTR","RC")))) +
  geom_bar(position="dodge", stat="identity") +
  geom_text(label=df$Number, position=position_dodge(width=0.9), vjust=-0.25)+
  scale_fill_brewer(palette="RdYlGn",name="Element")+
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
library(ggplot2)
library(dplyr)
library(viridis)
df <-as.data.frame(read.table("heatmap_cluster_gene_ID_update.txt", sep = "\t", header =TRUE))
df <- na.omit(df)
pdf("anova_distance_manhattan_1000kb_update.pdf",height = 2,width = 7)
ggplot(df,aes(factor(Caste,levels=c("Worker","Q0","Q1","Q2","Q3","Q4","King")), as.numeric(Distance)))+
  geom_point(aes(color=Distance,fill=Distance),color="black",position="jitter",shape=21,size=3)+
  scale_fill_viridis(option=1)+
  ylim(c(-1000000,1000000))+
  #scale_y_continuous(breaks=c(-1000000,-500000,0,500000,1000000),
  #                   labels=c("-1000000"="-1000 kb","-500000"="-500 kb","0","500000"="500 kb","1000000" = "1000 kb"))+
  theme_minimal()+
  xlab("Distance")+
  ylab("Cluster")+
  theme(legend.position = "none")
dev.off()

