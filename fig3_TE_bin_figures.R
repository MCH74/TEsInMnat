library(gridExtra)
library(ggplot2)
library(cowplot)
library(lemon)
library(viridis)
library(ggokabeito)
library(ggalt)

exp <- as.data.frame(read.table("genome_bin_intron_merge_20.txt", sep = "\t", header =TRUE))
exp$gene_expression <- as.numeric(as.character(exp$gene_expression))
exp$sum_TE_expression <- as.numeric(as.character(exp$sum_TE_expression))
exp_label <- within(exp, exp$genome_bin_label <- ifelse(genome_bin == "overlap", "0", 
                    ifelse(genome_bin == "1 to 1000", "1",
                    ifelse(genome_bin == "-1 to -1000", "-1",
                    ifelse(genome_bin == "-1001 to -2000", "-2",
                    ifelse(genome_bin == "1001 to 2000", "2",
                    ifelse(genome_bin == "-2001 to -3000", "-3",
                    ifelse(genome_bin == "2001 to 3000", "3",
                    ifelse(genome_bin == "-3001 to -4000", "-4",
                    ifelse(genome_bin == "3001 to 4000", "4",
                    ifelse(genome_bin == "-4001 to -5000", "-5",
                    ifelse(genome_bin == "4001 to 5000", "5",
                    ifelse(genome_bin == "-5001 to -6000", "-6",
                    ifelse(genome_bin == "5001 to 6000", "6",
                    ifelse(genome_bin == "-6001 to -7000", "-7",
                    ifelse(genome_bin == "6001 to 7000", "7",
                    ifelse(genome_bin == "-7001 to -8000", "-8",
                    ifelse(genome_bin == "7001 to 8000", "8",
                    ifelse(genome_bin == "-8001 to -9000", "-9",
                    ifelse(genome_bin == "8001 to 9000", "9",
                    ifelse(genome_bin == "-9001 to -10000", "-10",
                    ifelse(genome_bin == "9001 to 10000", "10",
                    ifelse(genome_bin == "10001 to 11000", "11",
                    ifelse(genome_bin == "-10001 to -11000", "-11",
                    ifelse(genome_bin == "-11001 to -12000", "-12",
                    ifelse(genome_bin == "11001 to 12000", "12",
                    ifelse(genome_bin == "-12001 to -13000", "-13",
                    ifelse(genome_bin == "12001 to 13000", "13",
                    ifelse(genome_bin == "-13001 to -14000", "-14",
                    ifelse(genome_bin == "13001 to 14000", "14",
                    ifelse(genome_bin == "-14001 to -15000", "-15",
                    ifelse(genome_bin == "14001 to 15000", "15",
                    ifelse(genome_bin == "-15001 to -16000", "-16",
                    ifelse(genome_bin == "15001 to 16000", "16",
                    ifelse(genome_bin == "-16001 to -17000", "-17",
                    ifelse(genome_bin == "16001 to 17000", "17",
                    ifelse(genome_bin == "-17001 to -18000", "-18",
                    ifelse(genome_bin == "17001 to 18000", "18",
                    ifelse(genome_bin == "-18001 to -19000", "-19",
                    ifelse(genome_bin == "18001 to 19000", "19",
                    ifelse(genome_bin == "-19001 to -20000", "-20",
                    ifelse(genome_bin == "19001 to 20000", "20",
                    NA))))))))))))))))))))))))))))))))))))))))))
                                                        
levels(exp_label$exp$genome_bin_label) <- c("-1","-11","-2","-12","-13","-14","-15","-16","-17","-18",
                            "-19","-20","-3","-4","-5","-6","-7","-8","-9","-10",
                            "1","11","2","12","13","14","15","16","17","18",
                            "19","20","3","4","5","6","7","8","9","10","0")
exp_label$exp$genome_bin_label <- factor(exp_label$exp$genome_bin_label,c("-20","-19","-18","-17","-16","-15","-14","-13","-12","-11",
                                          "-10","-9","-8","-7","-6","-5","-4","-3","-2","-1","0",
                                          "1","2","3","4","5","6","7","8","9","10",
                                          "11","12","13","14","15","16","17","18","19","20"))

#linear model
explm <- exp_label[exp_label$exp$genome_bin_label %in% c("20","-20","-10","10","-1","1", "0"), ] 
explm <- na.omit(explm)
explm <- explm[explm$exp$genome_bin_label == "10",] #select bin here
lm_TE_gene <- lm(log(sum_TE_expression+1,2) ~ gene_expression*caste, data = explm)
summary(lm_TE_gene)

#linear model figure
pdf("hex_introns_linear_model_bin_selection_test.pdf")
ggplot(explm,aes(gene_expression,log(sum_TE_expression+1,2)))+
  geom_hex()+
  scale_fill_viridis(trans = "log", breaks = c(1,20,400,8000,160000))+
  geom_smooth(method = "lm", se=FALSE,aes(color = factor(Caste,levels=c("FW","Q0m","Q3m","Q9m","Q31m","Q20y","K20y"))))+
  scale_color_okabe_ito()+
  theme_minimal()+
  ylab("log2median of sum of TE counts")+
  xlab("median of gene log2counts")+
  facet_wrap(~genome_bin)+
  labs(color = "Caste")
dev.off()

#total TE counts per bin figure
exp <- na.omit(exp)
exp_sum <- as.data.frame(aggregate(exp$sum_TE_expression, by=list(genome_bin=exp$genome_bin,caste=exp$caste), FUN=sum))
exp_sum_label <- read.table("TE_bins_line_plot.txt",sep = "\t", header = T)
pdf(file="Wide_20kb_sum_bin_test.pdf",width = 8, height = 3)
ggplot(exp_sum_label,aes(x=genome_bin_no,y=log(x,2),color = factor(Caste,levels=c("FW","Q0m","Q3m","Q9m","Q31m","Q20y","K20y")),group=factor(caste,levels = c("Worker","T0","T1","T2","T3","T4","King"))))+
  geom_line()+
  scale_color_okabe_ito()+
  theme_minimal()+
  ylab("log2sum of TE counts")+
  xlab("distance to gene [1 kb bin]")+
  scale_x_continuous(breaks=seq(-20,20,by=2))+
  coord_cartesian(ylim = c(0, 20), xlim = c(-19,19))+
  labs(color = "Caste")
dev.off()

