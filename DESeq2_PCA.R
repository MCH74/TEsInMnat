library(DESeq2)
BiocManager::install("DESeq2")
#have .count files in a "counts" folder
directory <- "counts"
sampleFiles <- grep("counts",list.files(directory),value=TRUE)
sampleCondition <- c("Worker","Worker","Worker","T0","T0","T0","T0","T1","T1","T1","T1","T2","T4","T4","T4","T4","King","King","King","Worker","T3","T3","T3","T2","T2")
worker <- c(1,1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0)
T0 <- c(0,0,0,1,1,1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0)
T1 <- c(0,0,0,0,0,0,0,1,1,1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0)
T2 <- c(0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,1,1)
T3 <- c(0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,1,1,0,0)
T4 <- c(0,0,0,0,0,0,0,0,0,0,0,0,1,1,1,1,0,0,0,0,0,0,0,0,0)
King <- c(0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,1,1,0,0,0,0,0,0)
sampleTable <- data.frame(sampleName = sampleFiles,
                          fileName = sampleFiles,
                          condition = sampleCondition,
                          Worker = worker,
                          T0 = T0,
                          T1 = T1,
                          T2 = T2,
                          T3 = T3,
                          T4 = T4,
                          King = King)
ddsHTSeq <- DESeqDataSetFromHTSeqCount(sampleTable = sampleTable,
                                       directory = directory,
                                       design= ~ condition)
dds <- DESeq(ddsHTSeq)
keep <- rowSums(counts(ddsHTSeq)) >= 1
dds <- dds[keep,]
countstable <- as.data.frame(counts(dds, normalize = TRUE))
colnames(countstable) <- c("Worker1","Worker2","Worker3","T01","T02","T03","T04","T11","T12","T13","T14","T22","T41","T42","T43","T44","King1","King2","King3","Worker4","T31","T32","T33","T23","T24")
#create normalised countstable; order of samples not sorted
write.table(countstable, file = "normalised_counts.txt", sep="\t")

#PCA
vsd <- vst(dds, blind=FALSE)
svg(filename="TE_PCA.svg",width = 12, height = 4)
plotPCA(vsd,intgroup=c("condition"))+
  theme_minimal()
dev.off() 

#PCA from txt file
install.packages("stats")
BiocManager::install("prcomp")
library(prcomp)
library(ggplot2)

hm_df <- as.data.frame(read.table("10kb_TE_intergenic_region_counts.txt",header=T,sep="\t"))
rownames(hm_df) <- hm_df$TE
hm_df$TE<-NULL
prot.pca <- prcomp(data.matrix(data.frame(t(na.omit(hm_df)))), scale. = TRUE)
score.df <- as.data.frame(prot.pca$x)

#Skree plot

Contributions <- round(summary(prot.pca)$importance[2,] * 100, 2)

pdf(file="skree_plot_pcs_TE.pdf", width=7, height=6)

barplot(Contributions[1:20], las=3, main="Screeplot CCLE PCA", xlab="component", ylab="Contribution(%)", col="lightblue")

PC1contr <- Contributions[1]

PC2contr <- Contributions[2]

dev.off()

pdf(file="Intergenic_Repeat_Expression_PCA.pdf", width=6, height=5)
ggplot(data = score.df, aes(x = PC1, y = PC2, label = rownames(score.df))) +
  geom_point(aes(color=Caste,size = 4))+
  #scale_colour_manual(values = c("0" = viridis(20)[2],"15" = viridis(20)[7],"45" =viridis(20)[13],"120" =viridis(20)[18]))+
  xlab(sprintf("PC1 (%s%%)", PC1contr)) + ylab(sprintf("PC2 (%s%%)", PC2contr)) + theme_minimal()+
  labs(color = "Caste")
dev.off()

score.df$Caste <- c("Worker","Worker","Worker","Worker","T0","T0","T0","T0","T1","T1","T1","T1","T2","T2","T2","T3","T3","T3","T4","T4","T4","T4","King","King","King")
