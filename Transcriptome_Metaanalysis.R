### The script was used to combine transcriptome data from multiple studies describing human naive pluripotent cells, and to analyze the data to be used as input for metabolic network reconstruction
### Meisam Yousefi et al. (2019) Cell and Bioscience [PMID:31485322]


# Setting directory and loading required packages

setwd("HumanNaive")
library(sva)
library(ggplot2)
library(limma)

# Reading expression data tables

Arr.Hanna <- read.delim("Data/Array/Hanna2010.txt")
Arr.Gafni <- read.delim("Data/Array/Gafni2013.txt")
Arr.Qin <- read.delim("Data/Array/Qin2016.txt")
Arr.Vala <- read.delim("Data/Array/Valamehr2014.txt")
Arr.Theu <- read.csv("Data/Array/Theunissen.csv")
Seq.Guo <- read.csv("Data/RNAseq/Guo-count.csv")
Seq.Tak <- read.csv("Data/RNAseq/Takashima-count.csv")
Seq.Warr <- read.csv("Data/RNAseq/Warrier-count.csv")

# Trimming/cleaning data tables for subsequent analysis

colnames(Arr.Vala)[18] <- "Gene"
colnames(Arr.Qin)[13] <- "Gene"
colnames(Arr.Hanna)[19] <- "Gene"
colnames(Arr.Gafni)[13] <- "Gene"
colnames(Arr.Theu)[1] <- "Gene"
colnames(Seq.Tak)[1] <- "Gene"
colnames(Seq.Warr)[1] <- "Gene"
Arr.Theu <- Arr.Theu[,-9:-10]
Seq.Guo <- Seq.Guo[,2:11]
Seq.Tak <- Seq.Tak[,c(1,6:11)]

TranscriptomeList <- list(Arr.Gafni,Arr.Hanna,Arr.Qin,Arr.Vala,Arr.Theu,Seq.Guo,Seq.Tak,Seq.Warr)
TranscriptomeList <- lapply(TranscriptomeList, function(x) {
  x <- x[!x$Gene=="", ]
})
names(TranscriptomeList) <- c("Arr.Gafni","Arr.Hanna","Arr.Qin","Arr.Vala","Arr.Theu","Seq.Guo","Seq.Tak","Seq.Warr")
list2env(TranscriptomeList, .GlobalEnv)

TranscriptomeList <- list(Arr.Gafni,Arr.Hanna,Arr.Qin,Arr.Vala,Arr.Theu,Seq.Guo,Seq.Tak,Seq.Warr)
TranscriptomeList <- lapply(TranscriptomeList, function(x) {
  x <- aggregate(x, list(x$Gene), mean)
})
names(TranscriptomeList) <- c("Arr.Gafni","Arr.Hanna","Arr.Qin","Arr.Vala","Arr.Theu","Seq.Guo","Seq.Tak","Seq.Warr")
list2env(TranscriptomeList, .GlobalEnv)

Arr.Gafni <- Arr.Gafni[,-14]
Arr.Hanna <- Arr.Hanna[,-20]
Arr.Qin <- Arr.Qin[,-14]
Arr.Theu <- Arr.Theu[,-2]
Arr.Vala <- Arr.Vala[,-19]
Seq.Guo <- Seq.Guo[,-2]
Seq.Warr <- Seq.Warr[,-2]
Seq.Tak <- Seq.Tak[,-2]

# Combining the expression tables from all studies

Trns <- merge(Arr.Vala, Arr.Qin, by="Group.1")
Trns <- merge(Trns, Arr.Hanna, by="Group.1")
Trns <- merge(Trns, Arr.Gafni, by="Group.1")
Trns <- merge(Trns, Arr.Theu, by="Group.1")
Trns <- merge(Trns, Seq.Warr, by="Group.1")
Trns <- merge(Trns, Seq.Guo, by="Group.1")
Trns <- merge(Trns, Seq.Tak, by="Group.1")

# Normalizing the gene expression values across all study samples

rownames(Trns) <- Trns$Group.1
Trns[,-1] <- normalizeQuantiles(as.matrix(Trns[,-1]))

# Principal components analysis without batch effects removal

Groups <- c(rep("N",5),"P",rep("N",4),rep("P",7),rep("P",2),rep("N",3),"P","P","N","P","P","N","N",rep("P",12),rep("N",6),rep("N",4),"P","N","P",rep("N",4),"P",rep("P",2),rep("N",5),rep(c(rep("P",3),rep("N",9)),3),rep("N",9),rep("N",3),rep("P",3))
Batch <- c(rep("Vala",17),rep("Qin",12),rep("Hanna",18),rep("Gafni",12),rep("Theu",7),rep("NCM",6),rep("NHSM",3),rep("RT",3),rep("NHSM",3),rep("NCM",3),rep("NHSM",3),rep("RT",6),rep("NCM",3),rep("NHSM",3),rep("RT",3),rep("Guo",9),rep("Tak",6))

scale.Tr <- t(scale(t(Trns[,-1]), scale = F))
pc <- prcomp(scale.Tr)
pcr <- data.frame(pc$rotation)
pcr$Sample <- Groups
pcr$Batch <- Batch
ggplot(pcr, aes(PC1,PC2, color=Batch, shape=Sample)) + geom_point(size=3) + theme_bw()

# Principal components analysis after batch effects removal using ComBat

Transcription <- Trns[,-1]
Transcription <- normalizeQuantiles(as.matrix(Transcription))
condition <- data.frame("Sample"=colnames(Transcription),"con"=Groups,"batch"=Batch)
mod <- model.matrix(~con, data=condition)
# batch <- model.matrix(~batch, data=condition) #(Commented for unsupervised batch effects removal)
combat_Trans <- ComBat(dat=Transcription, batch=Batch, mod=NULL, par.prior=T, prior.plots =F)
combat_Trans[combat_Trans < 0] <- 0
scale.combat <- t(scale(t(combat_Trans), scale = F))
pc <- prcomp(scale.combat)
pcr <- data.frame(pc$rotation)
pcr$Sample <- Groups
pcr$Batch <- Batch
ggplot(pcr, aes(PC1,PC2, color=Batch, shape=Sample)) + geom_point(size=3) + theme_bw()

# Principal components analysis after k-mer clustering on naive samples data

Naive.Tr <- scale.combat[,c(1:5,7:10,20:22,25,28:29,42:51,53,55:58,62:66,70:78,82:90,94:114)]
batche <- c(rep("Vala",9),rep("Qin",6),rep("Hanna",6),rep("Gafni",9),rep("Theu",5),rep("War1",9),rep("War2",9),rep("War3",9),rep("Guo",9),rep("Tak",3))
pc <- prcomp(Naive.Tr)
pcr <- data.frame(pc$rotation)
cl<-kmeans(pcr,2)
out <- cbind(pcr, clusterNum = cl$cluster)
ggplot(out, aes(PC1,PC2, color=clusterNum)) + geom_point(size=3) +theme_bw()

# Differentially expressed genes (DEGs) analysis using normalized data after batch effects removal 

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
write.table(tT, "Naive-Primed-Expression.txt", row.names=F, sep="\t", quote = F)
x <- data.frame(Gene=Trns$Group.1,combat_Trans)
write.table(x, "Trans.txt", row.names=F,sep = "\t", quote = F)
tT.Up <- subset(tT, logFC> 1 & adj.P.Val < 0.05)
write.table(as.character(strsplit2(unique(tT.Up$Gene),"///")),
            "UpRegulated.txt", row.names = F, sep = "\t", quote = F)
tT.Down <- subset(tT, logFC< -1 & adj.P.Val < 0.05)
write.table(as.character(strsplit2(unique(tT.Down$Gene),"///")),
            "DownRegulated.txt", row.names = F, sep = "\t", quote = F)