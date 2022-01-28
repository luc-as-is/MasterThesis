library (tidyr)
library(stringr)
library(dplyr)
library(DESeq2)
library(apeglm)
library(VennDiagram)
library(RColorBrewer)
library(readr)

setwd("~/lukeball/HostPref_R")

master<- read.table("4. DESeq2/antennae_DESeq_F1_MPCP_female_old_sig_autosomes.txt")
row.names(master)<- master$gene_id

# LEGS
## F1_MPCP_female_old ====

# read files in
F1_MPCP_female_old <- read.table("3. ASE data prep/legs/5.legs_genes_melpcyd_F1_MPCP_female_old.txt")
coldata <- read.table("3. ASE data prep/legs/6.legs_coldata_F1_MPCP_female_old.txt")

coldata$individual <- as.factor(coldata$individual)
coldata$allele <- as.factor(coldata$allele)

#separate into broods
b07 <- F1_MPCP_female_old[,1:6]
b11 <- F1_MPCP_female_old[,7:10]

coldata07 <- coldata[1:6,]
coldata11 <- coldata[7:10,]

# create dataset for DESeq2
dds <- DESeqDataSetFromMatrix(countData = b11, colData = coldata11, design = ~individual + allele)

# do not normalize counts
# Because we have sample in the design, 
# and so we are looking at ratios within sample, we do not want 
# to use DESeq2's size factor normalization, and hence set the 
# size factors all to 1.
sizeFactors(dds) <- rep(1, length(coldata11$individual))

#run DESeq
dds <- DESeq(dds,fitType = "local")

# check and set contrasts
resultsNames(dds)
res <- results(dds,contrast = c("allele","melp","cyd"))

# plots

resLFC <- lfcShrink(dds, coef="allele_melp_vs_cyd", type="apeglm")
plotMA(resLFC, ylim=c(-3,3), main = "Differential Allele Expression - mature females")


# number of genes with padj < 0.05
sum(res$padj < 0.05, na.rm=TRUE)
# b07 = 1269 genes have ASE
# b11 = 964

F1_MPCP_female_old_b11 <- res[which(res$padj < 0.05),]

write.table(F1_MPCP_female_old_b11,"4. DESeq2/legs_DESeq_F1_MPCP_female_old_b11_sig.txt")
write.table(res,"4. DESeq2/legs_DESeq_F1_MPCP_female_old_b11_all.txt")

F1_MPCP_female_old_sig_res <- read.table("4. DESeq2/legs_DESeq_F1_MPCP_female_old_b07_sig.txt")
res <- read.table("4. DESeq2/legs_DESeq_F1_MPCP_female_old_b11_all.txt")


# separate chr 21 and autosomes
chr <- read.table("chromo_gene_info.txt",header = TRUE)
res$gene_id <- row.names(res)
res <- left_join(res,chr)

resauto <- res[-which(res$scaffold=="chr21"),]
resautosig <- resauto[which(resauto$padj<0.05),]
write.table(resautosig,"4. DESeq2/legs_DESeq_F1_MPCP_female_old_b11_sig_autosomes.txt")


# proportion of genes upregulated from each allele
resautosig <- read.table("DESeq2/legs_DESeq_F1_MPCP_female_old_sig_autosomes.txt")
ressig <- read.table("DESeq2/legs_DESeq_F1_MPCP_female_old_sig.txt")

# all chromosomes
sum(ressig$log2FoldChange > 0)/nrow(ressig)
#melp= 0.3452

#autosomes
sum(resautosig$log2FoldChange > 0)/nrow(resautosig)
#melp= 0.3340381


## F1_MPCP ====

# read files in
F1_MPCP <- read.table("3. ASE data prep/legs/5.legs_genes_melpcyd_F1_MPCP.txt")
coldata <- read.table("3. ASE data prep/legs/6.legs_coldata_F1_MPCP.txt")

coldata$individual <- as.factor(coldata$individual)
coldata$allele <- as.factor(coldata$allele)

#separate into broods
b07 <- F1_MPCP[,1:2]
b11 <- F1_MPCP[,3:12]

coldata07 <- coldata[1:2,]
coldata11 <- coldata[3:12,]

# create dataset for DESeq2
dds07 <- DESeqDataSetFromMatrix(countData = b07, colData = coldata07, design = ~allele)
dds11 <- DESeqDataSetFromMatrix(countData = b11, colData = coldata11, design = ~individual + allele)

# do not normalize counts
# Because we have sample in the design, 
# and so we are looking at ratios within sample, we do not want 
# to use DESeq2's size factor normalization, and hence set the 
# size factors all to 1.
sizeFactors(dds07) <- rep(1, length(coldata07$individual))
sizeFactors(dds11) <- rep(1, length(coldata11$individual))

#run DESeq
dds07 <- DESeq(dds07,fitType = "local")
dds11 <- DESeq(dds11,fitType = "local")

# check and set contrasts
resultsNames(dds11)
res <- results(dds11,contrast = c("allele","melp","cyd"))

restable <- data.frame(row.names = rownames(res),res$log2FoldChange,res$padj)
restable[-which(grepl("HMEL",rownames(restable))),]
# plots

resLFC <- lfcShrink(dds, coef="allele_melp_vs_cyd", type="apeglm")
plotMA(resLFC, ylim=c(-3,3), main = "Differential Allele Expression - mature males")


# number of genes with padj < 0.05
sum(res$padj < 0.05, na.rm=TRUE)
# b11 = 2680 genes have ASE

F1_MPCP_sig_res <- res[which(res$padj < 0.05),]

write.table(F1_MPCP_sig_res,"4. DESeq2/legs_DESeq_F1_MPCP_b11_sig.txt")
write.table(res,"4. DESeq2/legs_DESeq_F1_MPCP_b11_all.txt")

F1_MPCP_sig_res <- read.table("4. DESeq2/legs_DESeq_F1_MPCP_b11_sig.txt")
res <- read.table("4. DESeq2/legs_DESeq_F1_MPCP_b11_all.txt")


# separate chr 21 and autosomes
chr <- read.table("chromo_gene_info.txt",header = TRUE)
res$gene_id <- row.names(res)
res <- left_join(res,chr)

resauto <- res[-which(res$scaffold=="chr21"),]
resautosig <- resauto[which(resauto$padj<0.05),]
write.table(resautosig,"4. DESeq2/legs_DESeq_F1_MPCP_b11_sig_autosomes.txt")


# proportion of genes upregulated from each allele
resautosig <- read.table("DESeq2/legs_DESeq_F1_MPCP_sig_autosomes.txt")
ressig <- read.table("DESeq2/legs_DESeq_F1_MPCP_sig.txt")

# all chromosomes
sum(ressig$log2FoldChange > 0)/nrow(ressig)
#melp= 0.4137796

#autosomes
sum(resautosig$log2FoldChange > 0)/nrow(resautosig)
#melp= 0.4067633


## F1_MPCP_female_young ====

# read files in
F1_MPCP_female_young <- read.table("3. ASE data prep/legs/5.legs_genes_melpcyd_F1_MPCP_female_young.txt")
coldata <- read.table("3. ASE data prep/legs/6.legs_coldata_F1_MPCP_female_young.txt")

coldata$individual <- as.factor(coldata$individual)
coldata$allele <- as.factor(coldata$allele)

#separate into broods
b07 <- F1_MPCP_female_young[,1:2]
b11 <- F1_MPCP_female_young[,3:8]

coldata07 <- coldata[1:2,]
coldata11 <- coldata[3:8,]

# create dataset for DESeq2
dds07 <- DESeqDataSetFromMatrix(countData = b07, colData = coldata07, design = ~allele)
dds11 <- DESeqDataSetFromMatrix(countData = b11, colData = coldata11, design = ~individual + allele)

# do not normalize counts
# Because we have sample in the design, 
# and so we are looking at ratios within sample, we do not want 
# to use DESeq2's size factor normalization, and hence set the 
# size factors all to 1.
sizeFactors(dds07) <- rep(1, length(coldata07$individual))
sizeFactors(dds11) <- rep(1, length(coldata11$individual))

#run DESeq
dds07 <- DESeq(dds07,fitType = "local")
dds11 <- DESeq(dds11,fitType = "local")

# check and set contrasts
resultsNames(dds11)
res <- results(dds11,contrast = c("allele","melp","cyd"))

restable <- data.frame(row.names = rownames(res),res$log2FoldChange,res$padj)
restable[-which(grepl("HMEL",rownames(restable))),]
# plots

resLFC <- lfcShrink(dds, coef="allele_melp_vs_cyd", type="apeglm")
plotMA(resLFC, ylim=c(-3,3), main = "Differential Allele Expression - juvenile females")


# number of genes with padj < 0.05
sum(res$padj < 0.05, na.rm=TRUE)
# b11 = 483 genes have ASE

F1_MPCP_female_young_sig_res <- res[which(res$padj < 0.05),]

write.table(F1_MPCP_female_young_sig_res,"4. DESeq2/legs_DESeq_F1_MPCP_female_young_b11_sig.txt")
write.table(res,"4. DESeq2/legs_DESeq_F1_MPCP_female_young_b11_all.txt")

F1_MPCP_female_young_sig_res <- read.table("4. DESeq2/legs_DESeq_F1_MPCP_female_young_b11_sig.txt")
res <- read.table("4. DESeq2/legs_DESeq_F1_MPCP_female_young_b11_all.txt")


# separate chr 21 and autosomes
chr <- read.table("chromo_gene_info.txt",header = TRUE)
res$gene_id <- row.names(res)
res <- left_join(res,chr)

resauto <- res[-which(res$scaffold=="chr21"),]
resautosig <- resauto[which(resauto$padj<0.05),]
write.table(resautosig,"4. DESeq2/legs_DESeq_F1_MPCP_female_young_b11_sig_autosomes.txt")


# proportion of genes upregulated from each allele
resautosig <- read.table("DESeq2/legs_DESeq_F1_MPCP_female_young_sig_autosomes.txt")
ressig <- read.table("DESeq2/legs_DESeq_F1_MPCP_female_young_sig.txt")

# all chromosomes
sum(ressig$log2FoldChange > 0)/nrow(ressig)
#melp= 0.3722488

#autosomes
sum(resautosig$log2FoldChange > 0)/nrow(resautosig)
#melp= 0.4870105


## F1_MPCP_young ====

# read files in
F1_MPCP_young <- read.table("3. ASE data prep/legs/5.legs_genes_melpcyd_F1_MPCP_young.txt")
coldata <- read.table("3. ASE data prep/legs/6.legs_coldata_F1_MPCP_young.txt")

coldata$individual <- as.factor(coldata$individual)
coldata$allele <- as.factor(coldata$allele)

#separate into broods
b07 <- F1_MPCP_young[,1:4]
b11 <- F1_MPCP_young[,5:10]

coldata07 <- coldata[1:4,]
coldata11 <- coldata[5:10,]

# create dataset for DESeq2
dds07 <- DESeqDataSetFromMatrix(countData = b07, colData = coldata07, design = ~individual + allele)
dds11 <- DESeqDataSetFromMatrix(countData = b11, colData = coldata11, design = ~individual + allele)

# do not normalize counts
# Because we have sample in the design, 
# and so we are looking at ratios within sample, we do not want 
# to use DESeq2's size factor normalization, and hence set the 
# size factors all to 1.
sizeFactors(dds07) <- rep(1, length(coldata07$individual))
sizeFactors(dds11) <- rep(1, length(coldata11$individual))

#run DESeq
dds07 <- DESeq(dds07,fitType = "local")
dds11 <- DESeq(dds11,fitType = "local")

# check and set contrasts
res07 <- results(dds07,contrast = c("allele","melp","cyd"))
res11 <- results(dds11,contrast = c("allele","melp","cyd"))

restable <- data.frame(row.names = rownames(res),res$log2FoldChange,res$padj)
restable[-which(grepl("HMEL",rownames(restable))),]

# number of genes with padj < 0.05
sum(res07$padj < 0.05, na.rm=TRUE)
# b07 = 709 genes have ASE
sum(res11$padj < 0.05, na.rm=TRUE)
# b11 = 1811 genes have ASE


F1_MPCP_young_sig_res <- res07[which(res07$padj < 0.05),]

write.table(F1_MPCP_young_sig_res,"4. DESeq2/legs_DESeq_F1_MPCP_young_b07_sig.txt")
write.table(res07,"4. DESeq2/legs_DESeq_F1_MPCP_young_b07_all.txt")

F1_MPCP_young_sig_res <- res11[which(res11$padj < 0.05),]

write.table(F1_MPCP_young_sig_res,"4. DESeq2/legs_DESeq_F1_MPCP_young_b11_sig.txt")
write.table(res11,"4. DESeq2/legs_DESeq_F1_MPCP_young_b11_all.txt")

F1_MPCP_young_sig_res <- read.table("DESeq2/legs_DESeq_F1_MPCP_young_sig.txt")
res <- read.table("4. DESeq2/legs_DESeq_F1_MPCP_young_b11_all.txt")

# separate chr 21 and autosomes
chr <- read.table("chromo_gene_info.txt",header = TRUE)

res$gene_id <- row.names(res)
res <- left_join(res,chr)

resauto <- res[-which(res$scaffold=="chr21"),]
resautosig <- resauto[which(resauto$padj<0.05),]
write.table(resautosig,"4. DESeq2/legs_DESeq_F1_MPCP_young_b11_sig_autosomes.txt")


# venn diagram ====

fem_old<- read.table("4. DESeq2/legs_DESeq_F1_MPCP_female_old_b07_sig_autosomes.txt")
male_old<- read.table("4. DESeq2/legs_DESeq_F1_MPCP_b07_sig_autosomes.txt")
fem_young<- read.table("4. DESeq2/legs_DESeq_F1_MPCP_female_young_b07_sig_autosomes.txt")
male_young<- read.table("4. DESeq2/legs_DESeq_F1_MPCP_young_b07_sig_autosomes.txt")

fem_old <- c(fem_old$gene_id)
male_old <- c(male_old$gene_id)
fem_young <- c(fem_young$gene_id)
male_young <- c(male_young$gene_id)

color = c("#a61107","#4563a8","#a1c957","#7d599c")

venn.diagram(
  x = list(male_old, fem_old, male_young, fem_young),
  category.names = c("mature males" , "mature females" , "young males", "young females"),
  fill=color,
  cat.cex = 0.9,
  filename = 'venn_diagramm_legs_b11.png',
  imagetype = "png",
  output=TRUE
)


# lists for ASE heatmaps====
  # B11====
fem_old<- read.table("4. DESeq2/legs_DESeq_F1_MPCP_female_old_b11_sig_autosomes.txt")
row.names(fem_old)<- fem_old$gene_id
male_old<- read.table("4. DESeq2/legs_DESeq_F1_MPCP_b11_sig_autosomes.txt")
row.names(male_old)<- male_old$gene_id
fem_young<- read.table("4. DESeq2/legs_DESeq_F1_MPCP_female_young_b11_sig_autosomes.txt")
row.names(fem_young)<- fem_young$gene_id
male_young<- read.table("4. DESeq2/legs_DESeq_F1_MPCP_young_b11_sig_autosomes.txt")
row.names(male_young)<- male_young$gene_id

# merge tables so they are all in the same order
list <- full_join(master,fem_old,by="gene_id")
list <- full_join(fem_old,male_old,by="gene_id")
list <- full_join(list, fem_young,by="gene_id")
list <- full_join(list,male_young,by="gene_id")

row.names(list) <- list$gene_id

# separate back into groups
fem_old_ASE <- list[,c(1:6,8:10)]
male_old_ASE <- list[,c(11:19)]
fem_young_ASE <- list[,c(20:28)]
male_young_ASE <- list[,c(29:37)]
male_young_ASE <- list[,c(38:46)]

# make lists using genes expressed in old females
male_old_ASE <- male_old_ASE[rows,]
fem_young_ASE <- fem_young_ASE[rows,]
male_young_ASE <- male_young_ASE[rows,]
fem_old_ASE <- fem_old_ASE[rows,]

write.csv(fem_old_ASE, "4. DESeq2/legs_ASE_female_old_b11_fo_genes")
write.csv(male_old_ASE, "4. DESeq2/legs_ASE_male_old_b11_fo_genes")
write.csv(fem_young_ASE, "4. DESeq2/legs_ASE_female_young_b11_fo_genes")
write.csv(male_young_ASE, "4. DESeq2/legs_ASE_male_young_b11_fo_genes")


# make list using only annotated genes present in old females
master_ASE <- read.csv("4. DESeq2/legs_ASE_female_old_fo_genes",row.names = 1)
fem_old_ASE <- read.csv("4. DESeq2/legs_ASE_female_old_b11_fo_genes",row.names = 1)
male_old_ASE <- read.csv("4. DESeq2/legs_ASE_male_old_b11_fo_genes",row.names = 1)
fem_young_ASE <- read.csv("4. DESeq2/legs_ASE_female_young_b11_fo_genes",row.names = 1)
male_young_ASE <- read.csv("4. DESeq2/legs_ASE_male_young_b11_fo_genes",row.names = 1)

male_old_ASEsmol <- male_old_ASE[-which(grepl("HMEL",rownames(fem_old_ASE),ignore.case = FALSE)==TRUE),]
fem_young_ASEsmol <- fem_young_ASE[-which(grepl("HMEL",rownames(fem_old_ASE),ignore.case = FALSE)==TRUE),]
male_young_ASEsmol <- male_young_ASE[-which(grepl("HMEL",rownames(fem_old_ASE),ignore.case = FALSE)==TRUE),]
fem_old_ASEsmol <- fem_old_ASE[-which(grepl("HMEL",rownames(fem_old_ASE),ignore.case = FALSE)==TRUE),]

write.csv(fem_old_ASEsmol, "4. DESeq2/legs_ASE_female_old_b11_smol2.csv")
write.csv(male_old_ASEsmol, "4. DESeq2/legs_ASE_male_old_b11_smol2.csv")
write.csv(fem_young_ASEsmol, "4. DESeq2/legs_ASE_female_young_b11_smol2.csv")
write.csv(male_young_ASEsmol, "4. DESeq2/legs_ASE_male_young_b11_smol2.csv")


  # b07====
fem_old<- read.table("4. DESeq2/legs_DESeq_F1_MPCP_female_old_b07_sig_autosomes.txt")
row.names(fem_old)<- fem_old$gene_id
male_young<- read.table("4. DESeq2/legs_DESeq_F1_MPCP_young_b07_sig_autosomes.txt")
row.names(male_young)<- male_young$gene_id

# merge tables so they are all in the same order
list <- full_join(master,fem_old,by="gene_id")
list <- full_join(list,male_young,by="gene_id")

row.names(list) <- list$gene_id

# separate back into groups
master_ASE <- list[,c(1:6,8:10)]
fem_old_ASE <- list[,c(11:19)]
male_young_ASE <- list[,c(20:28)]

# make lists using genes expressed in old females
male_young_ASE <- male_young_ASE[-which(is.na(fem_old_ASE)),]
fem_old_ASE <- fem_old_ASE[-which(is.na(fem_old_ASE)),]

write.csv(fem_old_ASE, "4. DESeq2/legs_ASE_female_old_b07_fo_genes")
write.csv(male_young_ASE, "4. DESeq2/legs_ASE_male_young_b07_fo_genes")

# make list using only annotated genes present in old females
master_ASE <- read.csv("4. DESeq2/legs_ASE_female_old_fo_genes",row.names = 1)
fem_old_ASE <- read.csv("4. DESeq2/legs_ASE_female_old_b07_fo_genes",row.names = 1)
male_young_ASE <- read.csv("4. DESeq2/legs_ASE_male_young_b07_fo_genes",row.names = 1)

male_young_ASEsmol <- male_young_ASE[-which(grepl("HMEL",rownames(fem_old_ASE),ignore.case = FALSE)==TRUE),]
fem_old_ASEsmol <- fem_old_ASE[-which(grepl("HMEL",rownames(fem_old_ASE),ignore.case = FALSE)==TRUE),]

write.csv(fem_old_ASEsmol, "4. DESeq2/legs_ASE_female_old_b07_smol2.csv")
write.csv(male_young_ASEsmol, "4. DESeq2/legs_ASE_male_young_b07_smol2.csv")


# ANTENNAE ====
## F1_MPCP_female_old ====

# read files in
F1_MPCP_female_old <- read.table("3. ASE data prep/antennae/5.antennae_genes_melpcyd_F1_MPCP_female_old.txt")
coldata <- read.table("3. ASE data prep/antennae/6.antennae_coldata_F1_MPCP_female_old.txt")

coldata$individual <- as.factor(coldata$individual)
coldata$allele <- as.factor(coldata$allele)

#separate into broods
b07 <- F1_MPCP_female_old[,1:4]
b11 <- F1_MPCP_female_old[,5:6]

coldata07 <- coldata[1:4,]
coldata11 <- coldata[5:6,]

# create dataset for DESeq2
dds <- DESeqDataSetFromMatrix(countData = b07, colData = coldata07, design = ~individual + allele)

# do not normalize counts
# Because we have sample in the design, 
# and so we are looking at ratios within sample, we do not want 
# to use DESeq2's size factor normalization, and hence set the 
# size factors all to 1.
sizeFactors(dds) <- rep(1, length(coldata07$individual))

#run DESeq
dds <- DESeq(dds,fitType = "local")

# check and set contrasts
resultsNames(dds)
res <- results(dds,contrast = c("allele","melp","cyd"))

# number of genes with padj < 0.05
sum(res$padj < 0.05, na.rm=TRUE)
# b07 = 869 genes have ASE

F1_MPCP_female_old_b07 <- res[which(res$padj < 0.05),]

write.table(F1_MPCP_female_old_b07,"4. DESeq2/antennae_DESeq_F1_MPCP_female_old_b07_sig.txt")
write.table(res,"4. DESeq2/antennae_DESeq_F1_MPCP_female_old_b07_all.txt")

F1_MPCP_female_old_sig_res <- read.table("4. DESeq2/antennae_DESeq_F1_MPCP_female_old_b07_sig.txt")
res <- read.table("4. DESeq2/antennae_DESeq_F1_MPCP_female_old_b07_all.txt")


# separate chr 21 and autosomes
chr <- read.table("chromo_gene_info.txt",header = TRUE)
res$gene_id <- row.names(res)
res <- left_join(res,chr)

resauto <- res[-which(res$scaffold=="chr21"),]
resautosig <- resauto[which(resauto$padj<0.05),]
write.table(resautosig,"4. DESeq2/antennae_DESeq_F1_MPCP_female_old_b07_sig_autosomes.txt")


# proportion of genes upregulated from each allele
resautosig <- read.table("DESeq2/antennae_DESeq_F1_MPCP_female_old_sig_autosomes.txt")
ressig <- read.table("DESeq2/antennae_DESeq_F1_MPCP_female_old_sig.txt")

# all chromosomes
sum(ressig$log2FoldChange > 0)/nrow(ressig)
#melp= 0.3452

#autosomes
sum(resautosig$log2FoldChange > 0)/nrow(resautosig)
#melp= 0.3340381


## F1_MPCP ====

# read files in
F1_MPCP <- read.table("3. ASE data prep/antennae/5.antennae_genes_melpcyd_F1_MPCP.txt")
coldata <- read.table("3. ASE data prep/antennae/6.antennae_coldata_F1_MPCP.txt")

coldata$individual <- as.factor(coldata$individual)
coldata$allele <- as.factor(coldata$allele)

#separate into broods

## theyre all brood 11
## ASE is same as antennae_F1_MPCP

## F1_MPCP_female_young ====

# read files in
F1_MPCP_female_young <- read.table("3. ASE data prep/antennae/5.antennae_genes_melpcyd_F1_MPCP_female_young.txt")
coldata <- read.table("3. ASE data prep/antennae/6.antennae_coldata_F1_MPCP_female_young.txt")

coldata$individual <- as.factor(coldata$individual)
coldata$allele <- as.factor(coldata$allele)

#separate into broods
b07 <- F1_MPCP_female_young[,1:4]
b11 <- F1_MPCP_female_young[,5:8]

coldata07 <- coldata[1:4,]
coldata11 <- coldata[5:8,]

# create dataset for DESeq2
dds07 <- DESeqDataSetFromMatrix(countData = b07, colData = coldata07, design = ~individual+allele)
dds11 <- DESeqDataSetFromMatrix(countData = b11, colData = coldata11, design = ~individual + allele)

# do not normalize counts
# Because we have sample in the design, 
# and so we are looking at ratios within sample, we do not want 
# to use DESeq2's size factor normalization, and hence set the 
# size factors all to 1.
sizeFactors(dds07) <- rep(1, length(coldata07$individual))
sizeFactors(dds11) <- rep(1, length(coldata11$individual))

#run DESeq
dds07 <- DESeq(dds07,fitType = "local")
dds11 <- DESeq(dds11,fitType = "local")

# check and set contrasts
res <- results(dds07,contrast = c("allele","melp","cyd"))
res <- results(dds11,contrast = c("allele","melp","cyd"))

restable <- data.frame(row.names = rownames(res),res$log2FoldChange,res$padj)
restable[-which(grepl("HMEL",rownames(restable))),]

# number of genes with padj < 0.05
sum(res$padj < 0.05, na.rm=TRUE)
# b07 = 12 genes have ASE
# b11 = 10 genes have ASE

F1_MPCP_female_young_sig_res <- res[which(res$padj < 0.05),]

write.table(F1_MPCP_female_young_sig_res,"4. DESeq2/antennae_DESeq_F1_MPCP_female_young_b11_sig.txt")
write.table(res,"4. DESeq2/antennae_DESeq_F1_MPCP_female_young_b11_all.txt")


# separate chr 21 and autosomes
F1_MPCP_female_young_sig_res <- read.table("4. DESeq2/antennae_DESeq_F1_MPCP_female_young_b11_sig.txt")
res <- read.table("4. DESeq2/antennae_DESeq_F1_MPCP_female_young_b11_all.txt")
F1_MPCP_female_young_sig_res <- read.table("4. DESeq2/antennae_DESeq_F1_MPCP_female_young_b07_sig.txt")
res <- read.table("4. DESeq2/antennae_DESeq_F1_MPCP_female_young_b07_all.txt")

chr <- read.table("chromo_gene_info.txt",header = TRUE)
res$gene_id <- row.names(res)
res <- left_join(res,chr)

resauto <- res[-which(res$scaffold=="chr21"),]
resautosig <- resauto[which(resauto$padj<0.05),]
write.table(resautosig,"4. DESeq2/antennae_DESeq_F1_MPCP_female_young_b07_sig_autosomes.txt")


# proportion of genes upregulated from each allele
resautosig <- read.table("DESeq2/antennae_DESeq_F1_MPCP_female_young_sig_autosomes.txt")
ressig <- read.table("DESeq2/antennae_DESeq_F1_MPCP_female_young_sig.txt")

# all chromosomes
sum(ressig$log2FoldChange > 0)/nrow(ressig)
#melp= 0.3722488

#autosomes
sum(resautosig$log2FoldChange > 0)/nrow(resautosig)
#melp= 0.4870105


## F1_MPCP_young ====

# read files in
F1_MPCP_young <- read.table("3. ASE data prep/antennae/5.antennae_genes_melpcyd_F1_MPCP_young.txt")
coldata <- read.table("3. ASE data prep/antennae/6.antennae_coldata_F1_MPCP_young.txt")

coldata$individual <- as.factor(coldata$individual)
coldata$allele <- as.factor(coldata$allele)

#separate into broods
b07 <- F1_MPCP_young[,1:2]
b11 <- F1_MPCP_young[,3:6]

coldata07 <- coldata[1:2,]
coldata11 <- coldata[3:6,]

# create dataset for DESeq2
dds11 <- DESeqDataSetFromMatrix(countData = b11, colData = coldata11, design = ~individual + allele)

# do not normalize counts
# Because we have sample in the design, 
# and so we are looking at ratios within sample, we do not want 
# to use DESeq2's size factor normalization, and hence set the 
# size factors all to 1.
sizeFactors(dds11) <- rep(1, length(coldata11$individual))

#run DESeq
dds11 <- DESeq(dds11,fitType = "local")

# check and set contrasts
res <- results(dds11,contrast = c("allele","melp","cyd"))

restable <- data.frame(row.names = rownames(res),res$log2FoldChange,res$padj)
restable[-which(grepl("HMEL",rownames(restable))),]

# number of genes with padj < 0.05
sum(res$padj < 0.05, na.rm=TRUE)
# b11 = 1354 genes have ASE

F1_MPCP_young_sig_res <- res[which(res$padj < 0.05),]

write.table(F1_MPCP_young_sig_res,"4. DESeq2/antennae_DESeq_F1_MPCP_young_b11_sig.txt")
write.table(res,"4. DESeq2/antennae_DESeq_F1_MPCP_young_b11_all.txt")

F1_MPCP_young_sig_res <- read.table("4. DESeq2/antennae_DESeq_F1_MPCP_young_b11_sig.txt")
res <- read.table("4. DESeq2/antennae_DESeq_F1_MPCP_young_b11_all.txt")

# separate chr 21 and autosomes
chr <- read.table("chromo_gene_info.txt",header = TRUE)

res$gene_id <- row.names(res)
res <- left_join(res,chr)

resauto <- res[-which(res$scaffold=="chr21"),]
resautosig <- resauto[which(resauto$padj<0.05),]
write.table(resautosig,"4. DESeq2/antennae_DESeq_F1_MPCP_young_b11_sig_autosomes.txt")


# venn diagram ====

fem_old<- read.table("4. DESeq2/antennae_DESeq_F1_MPCP_female_old_b07_sig_autosomes.txt")
fem_young<- read.table("4. DESeq2/antennae_DESeq_F1_MPCP_female_young_b07_sig_autosomes.txt")

fem_old <- c(fem_old$gene_id)
fem_young <- c(fem_young$gene_id)

color = c("#a61107","#4563a8","#a1c957","#7d599c")
color07 = c("#4563a8","#7d599c")

venn.diagram(
  x = list(fem_old, fem_young),
  category.names = c("mature females" ,"young females"),
  fill=color07,
  cat.cex = 0.9,
  filename = 'venn_diagramm_antennae_b07.png',
  imagetype = "png",
  output=TRUE
)

male_old<- read.table("4. DESeq2/antennae_DESeq_F1_MPCP_sig_autosomes.txt")
fem_young<- read.table("4. DESeq2/antennae_DESeq_F1_MPCP_female_young_b11_sig_autosomes.txt")
male_young<- read.table("4. DESeq2/antennae_DESeq_F1_MPCP_female_young_b11_sig_autosomes.txt")

male_old<- c(male_old$gene_id)
male_young <- c(male_young$gene_id)
fem_young <- c(fem_young$gene_id)

color = c("#a61107","#a1c957","#7d599c")

venn.diagram(
  x = list(male_old,male_young, fem_young),
  category.names = c("mature males", "young males", "young females"),
  fill=color,
  cat.cex = 0.9,
  filename = 'venn_diagramm_antennae_b11.png',
  imagetype = "png",
  output=TRUE
)

# lists for ASE heatmaps====
  #B11====
male_old<- read.table("4. DESeq2/antennae_DESeq_F1_MPCP_sig_autosomes.txt")
row.names(male_old)<- male_old$gene_id
fem_young<- read.table("4. DESeq2/antennae_DESeq_F1_MPCP_female_young_b11_sig_autosomes.txt")
row.names(fem_young)<- fem_young$gene_id
male_young<- read.table("4. DESeq2/antennae_DESeq_F1_MPCP_young_b11_sig_autosomes.txt")
row.names(male_young)<- male_young$gene_id

# merge tables so they are all in the same order
list <- full_join(master,male_old,by="gene_id")
list <- full_join(list, fem_young,by="gene_id")
list <- full_join(list,male_young,by="gene_id")

row.names(list) <- list$gene_id

# separate back into groups
master_ASE <- list[,c(1:6,8:10)]
male_old_ASE <- list[,c(11:19)]
fem_young_ASE <- list[,c(20:28)]
male_young_ASE <- list[,c(29:37)]

# make lists using genes expressed in old females
male_old_ASE <- male_old_ASE[-which(is.na(master_ASE)),]
fem_young_ASE <- fem_young_ASE[-which(is.na(master_ASE)),]
male_young_ASE <- male_young_ASE[-which(is.na(master_ASE)),]

write.csv(male_old_ASE, "4. DESeq2/antennae_ASE_male_old_b11_fo_genes")
write.csv(fem_young_ASE, "4. DESeq2/antennae_ASE_female_young_b11_fo_genes")
write.csv(male_young_ASE, "4. DESeq2/antennae_ASE_male_young_b11_fo_genes")

# make list using only annotated genes present in old females
master_ASE <- read.table("4. DESeq2/antennae_DESeq_F1_MPCP_female_old_sig_autosomes.txt")
  row.names(master_ASE)<- master_ASE$gene_id
male_old_ASE <- read.csv("4. DESeq2/antennae_ASE_male_old_b11_fo_genes",row.names = 1)
fem_young_ASE <- read.csv("4. DESeq2/antennae_ASE_female_young_b11_fo_genes",row.names = 1)
male_young_ASE <- read.csv("4. DESeq2/antennae_ASE_male_young_b11_fo_genes",row.names = 1)

male_old_ASEsmol <- male_old_ASE[-which(grepl("HMEL",rownames(master_ASE),ignore.case = FALSE)==TRUE),]
fem_young_ASEsmol <- fem_young_ASE[-which(grepl("HMEL",rownames(master_ASE),ignore.case = FALSE)==TRUE),]
male_young_ASEsmol <- male_young_ASE[-which(grepl("HMEL",rownames(master_ASE),ignore.case = FALSE)==TRUE),]

write.csv(male_old_ASEsmol, "4. DESeq2/antennae_ASE_male_old_b11_smol.csv")
write.csv(fem_young_ASEsmol, "4. DESeq2/antennae_ASE_female_young_b11_smol.csv")
write.csv(male_young_ASEsmol, "4. DESeq2/antennae_ASE_male_young_b11_smol.csv")
b11<-read.csv("4. DESeq2/antennae_ASE_female_young_b11_fo_genes")
  #B07====

fem_old<- read.table("4. DESeq2/antennae_DESeq_F1_MPCP_female_old_b07_sig_autosomes.txt")
row.names(male_old)<- male_old$gene_id
fem_young<- read.table("4. DESeq2/antennae_DESeq_F1_MPCP_female_young_b07_sig_autosomes.txt")
row.names(fem_young)<- fem_young$gene_id

# merge tables so they are all in the same order
list <- full_join(master,fem_old,by="gene_id")
list <- full_join(list, fem_young,by="gene_id")

row.names(list) <- list$gene_id

# separate back into groups
master_ASE <- list[,c(1:6,8:10)]
fem_old_ASE <- list[,c(11:19)]
fem_young_ASE <- list[,c(20:28)]

# make lists using genes expressed in old females
fem_old_ASE <- fem_old_ASE[-which(is.na(master_ASE)),]
fem_young_ASE <- fem_young_ASE[-which(is.na(master_ASE)),]

write.csv(fem_old_ASE, "4. DESeq2/antennae_ASE_female_old_b07_fo_genes")
write.csv(fem_young_ASE, "4. DESeq2/antennae_ASE_female_young_b07_fo_genes")

# make list using only annotated genes present in old females
master_ASE <- read.table("4. DESeq2/antennae_DESeq_F1_MPCP_female_old_sig_autosomes.txt")
row.names(master_ASE)<- master_ASE$gene_id
fem_old_ASE <- read.csv("4. DESeq2/antennae_ASE_female_old_b07_fo_genes",row.names = 1)
fem_young_ASE <- read.csv("4. DESeq2/antennae_ASE_female_young_b07_fo_genes",row.names = 1)

fem_old_ASEsmol <- fem_old_ASE[-which(grepl("HMEL",rownames(fem_old_ASE),ignore.case = FALSE)==TRUE),]
fem_young_ASEsmol <- fem_young_ASE[-which(grepl("HMEL",rownames(fem_old_ASE),ignore.case = FALSE)==TRUE),]

write.csv(fem_old_ASEsmol, "4. DESeq2/antennae_ASE_female_old_b07_smol2.csv")
write.csv(fem_young_ASEsmol, "4. DESeq2/antennae_ASE_female_young_b07_smol2.csv")
