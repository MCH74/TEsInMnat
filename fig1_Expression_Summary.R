#Summary of total expression and expressed TEs
library(ggplot2)
library(gridExtra)
sum_fig <- as.data.frame(read.table("expression_summary_combined.txt", sep = "\t", header =TRUE))
sum_fig2<- sum_fig[sum_fig$type == "Counts",]
sum_fig3<- sum_fig[sum_fig$type != "Counts",]

#%repeat expression
per_fig <- as.data.frame(read.table("percentage_TE_expression.txt", sep = "\t", header =TRUE))

svg(filename="percentage_TE_expression.svg",width = 5, height = 4)
ggplot(per_fig,aes(x=factor(caste,levels = c("Worker","Q0","Q1","Q2","Q3","Q4","King")),y=perc_TE*100, color=factor(caste,levels = c("Worker","Q0","Q1","Q2","Q3","Q4","King"))))+
  geom_violin()+
  geom_point()+
  ylab("%TE expression")+
  xlab("")+
  expand_limits(y = 0)+
  guides(color=guide_legend(title="Caste"))+
  theme_minimal()
dev.off()

svg(filename="summary_expression_combined.svg",width = 8, height = 6)
ggplot(sum_fig2,aes(x=sample,y=value, group = caste))+
  geom_bar(aes(fill = region,color=region), stat = "identity")+
  geom_text(aes(label=percent))+
  ylab("sum of counts")+
  xlab("")+
  theme_minimal()
dev.off()

p1<-ggplot(sum_fig2,aes(x=factor(caste,levels = c("Worker","Q0","Q1","Q2","Q3","Q4","King")),y=value,group=type, color=type))+
  geom_point(color="red")+
  ylab("sum of counts")+
  xlab("")+
  theme_minimal()+
  facet_wrap(~type)
p2<-ggplot(sum_fig3,aes(x=factor(caste,levels = c("Worker","Q0","Q1","Q2","Q3","Q4","King")),y=value,group=type, color=type))+
  geom_point(color="blue")+
  ylab("expressed TEs")+
  xlab("caste")+
  theme_minimal()+
  facet_wrap(~type)

svg(filename="summary_expression_10kb_intergenic.svg",width = 8, height = 4)
grid.arrange(p1,p2,ncol=1,nrow=2)
dev.off()