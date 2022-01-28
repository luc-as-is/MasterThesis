library (tidyr)
library(stringr)
library(dplyr)

library(DESeq2)
library(ggplot2)

setwd("~/lukeball/HostPref_R")

# LEGS
  # data prep
    # combine data tables ====
    # name columns with individual ID
fem_old <- read.table("HTSeq output/merged_htseq_legs_F1_MPCP_female_old_genecounts.txt",row.names = 1)
colnames(fem_old) <- c("142","215","259","327","359")

male_old <-read.table("HTSeq output/merged_htseq_legs_F1_MPCP_genecounts.txt",row.names = 1)
colnames(male_old) <- c("181","282","294","300","334","348")

fem_young <-read.table("HTSeq output/merged_htseq_legs_F1_MPCP_female_young_genecounts.txt",row.names = 1)
colnames(fem_young) <- c("265","314","322","60")

male_young <-read.table("HTSeq output/merged_htseq_legs_F1_MPCP_young_genecounts.txt",row.names = 1)
colnames(male_young) <- c("173","239","285","304","317")

legs <- cbind(fem_old,male_old,fem_young,male_young)

write.table(legs, "genecount PCA/PCA_legs.txt")

    # create coldata table ====
sample <- factor(colnames(legs))
sex <- factor(c(rep("F",5),rep("M",6),rep("F",4),rep("M",5)))
age <- factor(c(rep("old",11),rep("young",9)))
brood <- factor(c(rep("007",3),rep("011",2),rep("007",1),rep("011",8),rep("007",3),rep("011",3)))

coldata <- data.frame(row.names = colnames(legs),sample,sex,age,brood)

coldata$group <- as.factor(paste(coldata$age,coldata$brood,sep = "."))

write.table(coldata,"genecount PCA/coldata_legs.txt")

  # PCA ====

legs <- read.table("2. genecount PCA/PCA_legs.txt")
coldata <- read.table("2. genecount PCA/coldata_legs.txt")
  coldata$sex <- as.factor(coldata$sex)
  coldata$group<- as.factor(coldata$group)
  coldata$sample <- as.factor(coldata$sample)
  
obp <- read.table("obp.txt",sep = "")
  obp$ID <- substr(x = obp$V9,regexpr("=", obp$V9)+1, regexpr(";", obp$V9)-1)
  obp$gene <- substr(obp$V12, 1, regexpr(";", obp$V12)-1)
  obp <- obp[,13:14]
or <- read.table("or.txt",sep = "")
  or$ID <- substr(x = or$V9,regexpr("=", or$V9)+1, regexpr(";", or$V9)-1)
  or$gene <- substr(or$V11, 1, regexpr(";", or$V11)-1)
  or<-or[,12:13]
csp <- read.delim("csp.txt",header=FALSE)
  csp$ID <- substr(x = csp$V9,regexpr("=", csp$V9)+1, regexpr(";", csp$V9)-1)
  csp$gene <- substr(csp$V9, regexpr("CSP", csp$V9), regexpr(";s", csp$V9)-1)
  csp=csp[,10:11]

chemogenes <- rbind(or,obp)
chemogenes <- rbind(chemogenes,csp)
rownames(legs)[match(chemogenes$ID, rownames(legs))] <- chemogenes$gene
which(grepl(paste(chemogenes$ID,collapse = "|"),row.names(list)))
  
  # helpful guide: RNA-seq workflow: gene-level exploratory analysis and differential expression
    # https://www.bioconductor.org/packages/devel/workflows/vignettes/rnaseqGene/inst/doc/rnaseqGene.html#exploratory-analysis-and-visualization

dds <- DESeqDataSetFromMatrix(countData = legs, colData = coldata, design = ~ sex + group)

# transform data using vsd (use rlog for smaller datasets i.e. < 30)
vst <- vst(dds,blind = FALSE)

# create pca data
pcaData <- plotPCA(vst, intgroup = c("group", "sex"),returnData=TRUE)

# calculate the variance        
percentVar <- round(100 * attr(pcaData, "percentVar"))

color = c("#237a01","#8ce34d","#4a0475","#ab79c9")
grouptext <- c("mature & brood 07","mature & brood 11","young & brood 07", "young & brood 11")

# ggplot
ggplot(pcaData, aes(x = PC1, y = PC2,fill=group.1,shape=sex),color="black") +
  geom_point(size = 3,alpha=0.6) +
  scale_shape_manual(values=c(24,21),name="Sex") + 
  scale_fill_manual(values=color,name = "Age & Brood", labels = grouptext) +
  xlab(paste0("PC1: ", percentVar[1], "% variance")) +
  ylab(paste0("PC2: ", percentVar[2], "% variance")) +
  guides(fill = guide_legend(override.aes = list(color=color,shape=22,size=4))) +
  #  scale_y_continuous(minor_breaks = NULL, limits=c(-20,30)) +
  #  scale_x_continuous(minor_breaks = NULL, limits=c(-30,30)) +
  #  ggtitle("RNAseq gene expression") +
  coord_fixed()

    # PC 3 and 4====

  # modify 'plotPCA' function from the DESeq2 package by running PCA 3&4.R code
pcaData <- PCA_more(vst, intgroup = c("group", "sex"),returnData=TRUE)

# calculate the variance        
percentVar <- round(100 * attr(pcaData, "percentVar"))

color = c("#2e7512","#8ce34d","#7315ad","#ab79c9")
grouptext <- c("mature & brood 07","mature & brood 11","young & brood 07", "young & brood 11")
# ggplot

ggplot(pcaData, aes(x = PC3, y = PC4,fill=group.1,shape=sex),color="black") +
  geom_point(size = 3,alpha=0.78) +
  scale_shape_manual(values=c(24,21),name="Sex") + 
  scale_fill_manual(values=color,name = "Age & Brood", labels = grouptext) +
  xlab(paste0("PC3: ", percentVar[1], "% variance")) +
  ylab(paste0("PC4: ", percentVar[2], "% variance")) +
  guides(fill = guide_legend(override.aes = list(color=color,shape=22,size=4))) +
  #  scale_y_continuous(minor_breaks = NULL, limits=c(-20,30)) +
  #  scale_x_continuous(minor_breaks = NULL, limits=c(-30,30)) +
  #  ggtitle("RNAseq gene expression") +
  coord_fixed()


  # chemosensory PCs ====
    # only using genes annotated as chemosensory genes
legs <- read.table("2. genecount PCA/PCA_legs.txt")
  legsGr <- legs[which(grepl("Gr",rownames(legs),ignore.case = TRUE)),]
  legsIR <- legs[which(grepl("IR",rownames(legs),ignore.case = TRUE)),]
  legsOBP <- legs[which(grepl("OBP",rownames(legs),ignore.case = TRUE)),]
  legsOR <- legs[which(grepl("OR",rownames(legs),ignore.case = TRUE)),]
  legsCSP <- legs[which(grepl("CSP",rownames(legs),ignore.case = TRUE)),]
  
  legschem <- rbind(legsGr,legsIR)
  legschem <- rbind(legschem,legsOBP)
  legschem <- rbind(legschem,legsOR)
  legschem <- rbind(legschem,legsCSP)
  
coldata <- read.table("2. genecount PCA/coldata_legs.txt")
  coldata$sex <- as.factor(coldata$sex)
  coldata$group<- as.factor(coldata$group)
  coldata$sample <- as.factor(coldata$sample)

# helpful guide: RNA-seq workflow: gene-level exploratory analysis and differential expression
# https://www.bioconductor.org/packages/devel/workflows/vignettes/rnaseqGene/inst/doc/rnaseqGene.html#exploratory-analysis-and-visualization

dds <- DESeqDataSetFromMatrix(countData = legschem, colData = coldata, design = ~ sex + group)

# transform data using vsd (use rlog for smaller datasets i.e. < 30)
vst<-varianceStabilizingTransformation(dds,blind = FALSE)

# create pca data
pcaData <- plotPCA(vst, intgroup = c("group", "sex"),returnData=TRUE)

# calculate the variance        
percentVar <- round(100 * attr(pcaData, "percentVar"))

color = c("#237a01","#8ce34d","#4a0475","#ab79c9")
grouptext <- c("mature & brood 07","mature & brood 11","young & brood 07", "young & brood 11")

# ggplot
ggplot(pcaData, aes(x = PC1, y = PC2,fill=group.1,shape=sex),color="black") +
  geom_point(size = 3,alpha=0.6) +
  scale_shape_manual(values=c(24,21),name="Sex") + 
  scale_fill_manual(values=color,name = "Age & Brood", labels = grouptext) +
  xlab(paste0("PC1: ", percentVar[1], "% variance")) +
  ylab(paste0("PC2: ", percentVar[2], "% variance")) +
  guides(fill = guide_legend(override.aes = list(color=color,shape=22,size=4))) +
  #  scale_y_continuous(minor_breaks = NULL, limits=c(-20,30)) +
  #  scale_x_continuous(minor_breaks = NULL, limits=c(-30,30)) +
  #  ggtitle("RNAseq gene expression") +
  coord_fixed()

    # PC 3 and 4 ====

# change 'plotPCA' function from the DESeq2 package
pcaData <- PCA_more(vst, intgroup = c("group", "sex"),returnData=TRUE)

# calculate the variance        
percentVar <- round(100 * attr(pcaData, "percentVar"))

color = c("#237a01","#8ce34d","#4a0475","#ab79c9")
grouptext <- c("mature & brood 07","mature & brood 11","young & brood 07", "young & brood 11")
# ggplot

ggplot(pcaData, aes(x = PC3, y = PC4,fill=group.1,shape=sex),color="black") +
  geom_point(size = 3,alpha=0.6) +
  scale_shape_manual(values=c(24,21),name="Sex") + 
  scale_fill_manual(values=color,name = "Age & Brood", labels = grouptext) +
  xlab(paste0("PC3: ", percentVar[1], "% variance")) +
  ylab(paste0("PC4: ", percentVar[2], "% variance")) +
  guides(fill = guide_legend(override.aes = list(color=color,shape=22,size=4))) +
  #  scale_y_continuous(minor_breaks = NULL, limits=c(-20,30)) +
  #  scale_x_continuous(minor_breaks = NULL, limits=c(-30,30)) +
  #  ggtitle("RNAseq gene expression") +
  coord_fixed()

# ANTENNAE ====
  # data prep      
    # combine data tables ====
fem_old <- read.table("HTSeq output/merged_htseq_ant_F1_MPCP_female_old_genecounts.txt",row.names = 1)
colnames(fem_old) <- c("141","214","326")

male_old <-read.table("HTSeq output/merged_htseq_ant_F1_MPCP_genecounts.txt",row.names = 1)
colnames(male_old) <- c("274","333","352","356","365")

fem_young <-read.table("HTSeq output/merged_htseq_ant_F1_MPCP_female_young_genecounts.txt",row.names = 1)
colnames(fem_young) <- c("137","261","313","59P")

male_young <-read.table("HTSeq output/merged_htseq_ant_F1_MPCP_young_genecounts.txt",row.names = 1)
colnames(male_young) <- c("172","284","316")

antennae <- cbind(fem_old,male_old,fem_young,male_young)

write.table (antennae,"genecount PCA/PCA_antennae.txt")
    # create coldata table ====
sample <- factor(colnames(antennae))
sex <- factor(c(rep("F",3),rep("M",5),rep("F",4),rep("M",3)))
age <- factor(c(rep("old",8),rep("young",7)))
brood <- factor(c(rep("007",2),rep("011",6),rep("007",1),rep("011",2),rep("007",2),rep("011",2)))

coldata <- data.frame(row.names = colnames(antennae),sample,sex,age,brood)

coldata$group <- as.factor(paste(coldata$age,coldata$brood,sep = "."))

write.table(coldata, "genecount PCA/coldata_antennae.txt")

  # PCA ====

antennae <- read.table("2. genecount PCA/PCA_antennae.txt")
coldata <- read.table("2. genecount PCA/coldata_antennae.txt")
  coldata$sex <- as.factor(coldata$sex)
  coldata$group<- as.factor(coldata$group)
  coldata$sample <- as.factor(coldata$sample)

# helpful guide: RNA-seq workflow: gene-level exploratory analysis and differential expression
# https://www.bioconductor.org/packages/devel/workflows/vignettes/rnaseqGene/inst/doc/rnaseqGene.html#exploratory-analysis-and-visualization

dds <- DESeqDataSetFromMatrix(countData = antennae, colData = coldata, design = ~sex + group)

# transform data using vsd (use rlog for smaller datasets i.e. < 30)
vst <- vst(dds,blind = FALSE)

# create pca data
pcaData <- plotPCA(vst, intgroup = c("group", "sex"),returnData=TRUE)

# calculate the variance        
percentVar <- round(100 * attr(pcaData, "percentVar"))

color = c("#237a01","#8ce34d","#4a0475","#ab79c9")
grouptext <- c("mature & brood 07","mature & brood 11","young & brood 07", "young & brood 11")

# ggplot
ggplot(pcaData, aes(x = PC1, y = PC2,fill=group.1,shape=sex),color="black") +
  geom_point(size = 3,alpha=0.6) +
  scale_shape_manual(values=c(24,21),name="Sex") + 
  scale_fill_manual(values=color,name = "Age & Brood", labels = grouptext) +
  xlab(paste0("PC1: ", percentVar[1], "% variance")) +
  ylab(paste0("PC2: ", percentVar[2], "% variance")) +
 # labs(fill = "Age & Brood", shape ="Sex") +
  guides(fill = guide_legend(override.aes = list(color=color,shape=22,size=4))) +
#  ggtitle("RNAseq gene expression") +
  coord_fixed()

    # PC 3 and 4 ====

# change 'plotPCA' function from the DESeq2 package
pcaData <- PCA_more(vst, intgroup = c("group", "sex"),returnData=TRUE)

# calculate the variance        
percentVar <- round(100 * attr(pcaData, "percentVar"))

color = c("#237a01","#8ce34d","#4a0475","#ab79c9")
grouptext <- c("mature & brood 07","mature & brood 11","young & brood 07", "young & brood 11")

# ggplot

ggplot(pcaData, aes(x = PC3, y = PC4,fill=group.1,shape=sex),color="black") +
  geom_point(size = 3,alpha=0.6) +
  scale_shape_manual(values=c(24,21),name="Sex") + 
  scale_fill_manual(values=color,name = "Age & Brood", labels = grouptext) +
  xlab(paste0("PC3: ", percentVar[1], "% variance")) +
  ylab(paste0("PC4: ", percentVar[2], "% variance")) +
  guides(fill = guide_legend(override.aes = list(color=color,shape=22,size=4))) +
  #  scale_y_continuous(minor_breaks = NULL, limits=c(-20,30)) +
  #  scale_x_continuous(minor_breaks = NULL, limits=c(-30,30)) +
  #  ggtitle("RNAseq gene expression") +
  coord_fixed()

  # chemo PCs ====
antennae <- read.table("2. genecount PCA/PCA_antennae.txt")

antennaeGr <- antennae[which(grepl("Gr",rownames(antennae),ignore.case = TRUE)),]
antennaeIR <- antennae[which(grepl("IR",rownames(antennae),ignore.case = TRUE)),]
antennaeOBP <- antennae[which(grepl("OBP",rownames(antennae),ignore.case = TRUE)),]
antennaeOR <- antennae[which(grepl("OR",rownames(antennae),ignore.case = TRUE)),]
antennaeCSP <- antennae[which(grepl("CSP",rownames(antennae),ignore.case = TRUE)),]

antennaechem <- rbind(antennaeGr,antennaeIR)
antennaechem <- rbind(antennaechem,antennaeOBP)
antennaechem <- rbind(antennaechem,antennaeOR)
antennaechem <- rbind(antennaechem,antennaeCSP)

coldata <- read.table("2. genecount PCA/coldata_antennae.txt")
coldata$sex <- as.factor(coldata$sex)
coldata$group<- as.factor(coldata$group)
coldata$sample <- as.factor(coldata$sample)

# helpful guide: RNA-seq workflow: gene-level exploratory analysis and differential expression
# https://www.bioconductor.org/packages/devel/workflows/vignettes/rnaseqGene/inst/doc/rnaseqGene.html#exploratory-analysis-and-visualization

dds <- DESeqDataSetFromMatrix(countData = antennaechem, colData = coldata, design = ~ sex + group)

# transform data using vsd (use rlog for smaller datasets i.e. < 30)
vst<-varianceStabilizingTransformation(dds,blind = FALSE)

# create pca data
pcaData <- plotPCA(vst, intgroup = c("group", "sex"),returnData=TRUE)

# calculate the variance        
percentVar <- round(100 * attr(pcaData, "percentVar"))

color = c("#237a01","#8ce34d","#4a0475","#ab79c9")
grouptext <- c("mature & brood 07","mature & brood 11","young & brood 07", "young & brood 11")

# ggplot
ggplot(pcaData, aes(x = PC1, y = PC2,fill=group.1,shape=sex),color="black") +
  geom_point(size = 3,alpha=0.6) +
  scale_shape_manual(values=c(24,21),name="Sex") + 
  scale_fill_manual(values=color,name = "Age & Brood", labels = grouptext) +
  xlab(paste0("PC1: ", percentVar[1], "% variance")) +
  ylab(paste0("PC2: ", percentVar[2], "% variance")) +
  guides(fill = guide_legend(override.aes = list(color=color,shape=22,size=4))) +
  #  scale_y_continuous(minor_breaks = NULL, limits=c(-20,30)) +
  #  scale_x_continuous(minor_breaks = NULL, limits=c(-30,30)) +
  #  ggtitle("RNAseq gene expression") +
  coord_fixed()

    # PC 3 and 4 ====

# change 'plotPCA' function from the DESeq2 package
pcaData <- PCA_more(vst, intgroup = c("group", "sex"),returnData=TRUE)

# calculate the variance        
percentVar <- round(100 * attr(pcaData, "percentVar"))

color = c("#237a01","#8ce34d","#4a0475","#ab79c9")
grouptext <- c("mature & brood 07","mature & brood 11","young & brood 07", "young & brood 11")
# ggplot

ggplot(pcaData, aes(x = PC3, y = PC4,fill=group.1,shape=sex),color="black") +
  geom_point(size = 3,alpha=0.6) +
  scale_shape_manual(values=c(24,21),name="Sex") + 
  scale_fill_manual(values=color,name = "Age & Brood", labels = grouptext) +
  xlab(paste0("PC3: ", percentVar[1], "% variance")) +
  ylab(paste0("PC4: ", percentVar[2], "% variance")) +
  guides(fill = guide_legend(override.aes = list(color=color,shape=22,size=4))) +
  #  scale_y_continuous(minor_breaks = NULL, limits=c(-20,30)) +
  #  scale_x_continuous(minor_breaks = NULL, limits=c(-30,30)) +
  #  ggtitle("RNAseq gene expression") +
  coord_fixed()
