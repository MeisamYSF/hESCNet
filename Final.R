setwd("Project/HumanNaive/")
library(GEOquery)
library(sva)
library(ggplot2)
library(DESeq2)
library(plyr)
library(limma)
library(pamr)
library(reshape2)

Arr.Hanna <- read.delim("Data/Array/Hanna2010.txt")
Arr.Gafni <- read.delim("Data/Array/Gafni2013.txt")
Arr.Qin <- read.delim("Data/Array/Qin2016.txt")
Arr.Vala <- read.delim("Data/Array/Valamehr2014.txt")
Arr.Theu <- read.csv("Data/Array/Theunissen.csv")
Seq.Guo <- read.csv("Data/RNAseq/Guo-count.csv")
# Seq.Tak <- read.csv("Data/RNAseq/Takashima-count.csv")
Seq.Warr <- read.csv("Data/RNAseq/Warrier-count.csv")

colnames(Arr.Vala)[18] <- "Gene"      # rep("N",5),"P",rep("N",4),rep("P",7)
colnames(Arr.Qin)[13] <- "Gene"       # rep("P",2),rep("N",3),"P","P","N","P","P","N","N"
colnames(Arr.Hanna)[19] <- "Gene"     # rep("P",12),rep("N",6)
colnames(Arr.Gafni)[13] <- "Gene"     # rep("N",4),"P","N","P",rep("N",4),"P"
colnames(Arr.Theu)[1] <- "Gene"       # rep("P",2),rep("N",5)
# colnames(Seq.Tak)[1] <- "Gene"        # rep("P",3),rep("N",3)
colnames(Seq.Warr)[1] <- "Gene"       # rep(c(rep("P",3),rep("N",9)),3)
                                      # rep("N",9)  Guo
Groups <- c(rep("N",5),"P",rep("N",4),rep("P",7),rep("P",2),rep("N",3),"P","P","N","P","P","N","N",rep("P",12),rep("N",6),rep("N",4),"P","N","P",rep("N",4),"P",rep("P",2),rep("N",5),rep(c(rep("P",3),rep("N",9)),3),rep("N",9),rep("N",3),rep("P",3))
Batch <- c(rep("Vala",17),rep("Qin",12),rep("Hanna",18),rep("Gafni",12),rep("Theu",7),rep("NCM",6),rep("NHSM",3),rep("RT",3),rep("NHSM",3),rep("NCM",3),rep("NHSM",3),rep("RT",6),rep("NCM",3),rep("NHSM",3),rep("RT",3),rep("Guo",9),rep("Tak",6))

Groups <- c(rep("P",2),rep("N",5),rep(c(rep("P",3),rep("N",9)),3),rep("N",9),rep("N",3),rep("P",3))
Batch <- c(rep("Theu",7),rep("War1",12),rep("War2",12),rep("War3",12),rep("Guo",9),rep("Tak",6))

Arr.Theu <- Arr.Theu[,-9:-10]
Seq.Guo <- Seq.Guo[,-1]
Seq.Guo <- Seq.Guo[,1:16]
# Seq.Tak <- Seq.Tak[,-2:-5]
# Seq.Tak <- Seq.Tak[,1:7]

gr <- c(rep("N",12),rep("P",3))
colData <- data.frame(condition=factor(gr), libType="single-end")
cds <- DESeqDataSetFromMatrix(Seq.Guo[,-1], colData, design=~condition)
cds <- DESeq(cds)
Seq.Guo[,-1] <- log2(1+counts(cds, normalized=T))

# Seq.Guo[,-1] <- log2(Seq.Guo[,-1] +1)
# Seq.Tak[,-1] <- log2(Seq.Tak[,-1] +1)
Seq.Warr[,-1] <- log2(Seq.Warr[,-1] +1)

Arr.Gafni <- Arr.Gafni[!Arr.Gafni$Gene=="", ]
Arr.Hanna <- Arr.Hanna[!Arr.Hanna$Gene=="", ]
Arr.Qin <- Arr.Qin[!Arr.Qin$Gene=="", ]
Arr.Vala <- Arr.Vala[!Arr.Vala$Gene=="", ]

Arr.Gafni <- aggregate(Arr.Gafni, list(Arr.Gafni$Gene), mean)
Arr.Hanna <- aggregate(Arr.Hanna, list(Arr.Hanna$Gene), mean)
Arr.Qin <- aggregate(Arr.Qin, list(Arr.Qin$Gene), mean)
Arr.Vala <- aggregate(Arr.Vala, list(Arr.Vala$Gene), mean)
Arr.Theu <- aggregate(Arr.Theu, list(Arr.Theu$Gene), mean)
Seq.Guo <- aggregate(Seq.Guo, list(Seq.Guo$Gene), mean)
# Seq.Tak <- aggregate(Seq.Tak, list(Seq.Tak$Gene), mean)
Seq.Warr <- aggregate(Seq.Warr, list(Seq.Warr$Gene), mean)

Arr.Gafni <- Arr.Gafni[,-14]
Arr.Hanna <- Arr.Hanna[,-20]
Arr.Qin <- Arr.Qin[,-14]
Arr.Theu <- Arr.Theu[,-2]
Arr.Vala <- Arr.Vala[,-19]
Seq.Guo <- Seq.Guo[,-2]
Seq.Warr <- Seq.Warr[,-2]
#Arr.Gafni[,13] <- make.names(Arr.Gafni[,13], unique = TRUE)
#Arr.Qin[,13] <- make.names(Arr.Qin[,13], unique = TRUE)
#Arr.Hanna[,19] <- make.names(Arr.Hanna[,19], unique = T)
#Arr.Vala[,18] <- make.names(Arr.Vala[,18], unique = T)

Arr.Hanna[,-1] <- log2(Arr.Hanna[,-1] +1)
Arr.Theu[,-1] <- log2(Arr.Theu[,-1] +1)

Trns <- merge(Arr.Vala, Arr.Qin, by="Group.1")
Trns <- merge(Trns, Arr.Hanna, by="Group.1")
Trns <- merge(Trns, Arr.Gafni, by="Group.1")
Trns <- merge(Trns, Arr.Theu, by="Group.1")
# Trns <- merge(Trns, Seq.Tak, by="Group.1")
Trns <- merge(Trns, Seq.Warr, by="Group.1")
Trns <- merge(Trns, Seq.Guo, by="Group.1")

Trns <- Arr.Theu
Trns <- merge(Trns, Seq.Warr, by="Group.1")
Trns <- merge(Trns, Seq.Guo, by="Group.1")

rownames(Trns) <- Trns$Group.1
Trns[,-1] <- normalizeQuantiles(as.matrix(Trns[,-1]))
scale.Tr <- t(scale(t(Trns[,-1]), scale = F))
pc <- prcomp(scale.Tr)
pcr <- data.frame(pc$rotation)
pcr$Sample <- Groups
pcr$Batch <- Batch
ggplot(pcr, aes(PC1,PC2, color=Batch, shape=Sample)) + geom_point(size=3) + theme_bw()

Transcription <- Trns[,-1]
Transcription <- normalizeQuantiles(as.matrix(Transcription))
condition <- data.frame("Sample"=colnames(Transcription),"con"=Groups,"batch"=Batch)
mod <- model.matrix(~con, data=condition)
#  batch <- model.matrix(~batch, data=condition)
combat_Trans <- ComBat(dat=Transcription, batch=Batch, mod=NULL, par.prior=T, prior.plots =F)
combat_Trans[combat_Trans < 0] <- 0
scale.combat <- t(scale(t(combat_Trans), scale = F))
pc <- prcomp(scale.combat)
pcr <- data.frame(pc$rotation)
pcr$Sample <- Groups
pcr$Batch <- Batch
ggplot(pcr, aes(PC1,PC2, color=Batch, shape=Sample)) + geom_point(size=3) + theme_bw()

gset <- data.frame(combat_Trans)
rownames(gset) <- Trns[,1]
design <- model.matrix(~-1 + factor(Groups))
colnames(design) <- c("N","P")
contrast <- makeContrasts(N-P, levels = design)
fit <- lmFit(gset, design)
fit2 <- contrasts.fit(fit, contrast)
fit2 <- eBayes(fit2,0.01)
tT <- topTable(fit2, adjust="fdr", sort.by="B", number=Inf)
tT <- subset(tT, select=c("adj.P.Val","logFC"))
tT$Gene <- rownames(tT)

ggplot(tT,aes(logFC,-log10(adj.P.Val)))+geom_point()+theme_bw()

write.table(tT, "Naive-Primed-Expression.txt", row.names=F, sep="\t", quote = F)
x <- data.frame(Gene=Trns$Group.1,combat_Trans)
write.table(x, "Trans.txt", row.names=F,sep = "\t", quote = F)
tT.Up <- subset(tT, logFC> 1 & adj.P.Val < 0.05)
write.table(as.character(strsplit2(unique(tT.Up$Gene),"///")),
            "UpRegulated.txt", row.names = F, sep = "\t", quote = F)
tT.Down <- subset(tT, logFC< -1 & adj.P.Val < 0.05)
write.table(as.character(strsplit2(unique(tT.Down$Gene),"///")),
            "DownRegulated.txt", row.names = F, sep = "\t", quote = F)

gsea <- combat_Trans
colnames(gsea) <- Groups
gsea <- cbind(gsea[,colnames(gsea)=="N"],gsea[,colnames(gsea)=="P"])
gsea <- data.frame(Gene=Trns$Group.1, gsea)
gsea[,-1] <- round(gsea[,-1],3)
write.table(gsea, "Naive_Trim_modNULL.txt", row.names =F, col.names = T, sep="\t", quote = F)

#ALL
cat("N", paste0("N.",1:73), "P", paste0("P.", 1:42),sep ="\n")
#TRIM
cat("N", paste0("N.",1:43), "P", paste0("P.", 1:13),sep ="\n")

zzz <- data.frame(Gene=gsea$Gene,N=rowMeans(gsea[,3:45]),P=rowMeans(gsea[,46:59]))
subset(zzz, abs(zzz$N-zzz$P)>2)
ggg <- read.delim("recon-genes.tsv")
colnames(ggg)[2] <- "Gene"
mmm <- merge(ggg, zzz, by="Gene")
mmm <- mmm[,-2:-17]




gsea <- Seq.Warr
colnames(gsea) <- rep(c(rep("P",3),rep("N",9)),3)
View(gsea)
colnames(gsea) <- c("Group.1",rep(c(rep("P",3),rep("N",9)),3))
gsea <- cbind(gsea[,colnames(gsea)=="N"],gsea[,colnames(gsea)=="P"])
gsea <- data.frame(Gene=Seq.Warr$Group.1, gsea)
gsea[,-1] <- round(gsea[,-1],3)
write.table(gsea, "gseaSEQwarr.txt", row.names =F, col.names = T, sep="\t", quote = F)
cat("N", paste0("N.",1:26), "P", paste0("P.", 1:8),sep ="\n")


z <- scale.combat[,c(1:5,7:10,20:22,25,28:29,42:51,53,55:58,62:66,70:78,82:90,94:114)]
batche <- c(rep("Vala",9),rep("Qin",6),rep("Hanna",6),rep("Gafni",9),rep("Theu",5),rep("War1",9),rep("War2",9),rep("War3",9),rep("Guo",9),rep("Tak",3))
pc <- prcomp(z)
pcr <- data.frame(pc$rotation)
cl<-kmeans(pcr,2)
out <- cbind(pcr, clusterNum = cl$cluster)
ggplot(out, aes(PC1,PC2, color=clusterNum)) + geom_point(size=3) +theme_bw()

pcr$batch <- batche
ggplot(pcr, aes(PC1,PC2, color=batche)) + geom_point(size=3)
