library (tidyr)
library(stringr)
library(dplyr)
library(DESeq2)
library(apeglm)
library(VennDiagram)
library(RColorBrewer)
library(readr)

setwd("~/lukeball/HostPref_R")

# LEGS
    ## F1_MPCP_female_old ====

# read files in
F1_MPCP_female_old <- read.table("3. ASE data prep/legs/5.legs_genes_melpcyd_F1_MPCP_female_old.txt")
coldata <- read.table("3. ASE data prep/legs/6.legs_coldata_F1_MPCP_female_old.txt")

coldata$individual <- as.factor(coldata$individual)
coldata$allele <- as.factor(coldata$allele)

# create dataset for DESeq2
dds <- DESeqDataSetFromMatrix(countData = F1_MPCP_female_old, colData = coldata, design = ~individual + allele)

# do not normalize counts
# Because we have sample in the design, 
# and so we are looking at ratios within sample, we do not want 
# to use DESeq2's size factor normalization, and hence set the 
# size factors all to 1.
sizeFactors(dds) <- rep(1, length(coldata$individual))

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
# = 504 genes have ASE

F1_MPCP_female_old_sig_res <- res[which(res$padj < 0.05),]

write.table(F1_MPCP_female_old_sig_res,"DESeq2/legs_DESeq_F1_MPCP_female_old_sig.txt")
write.table(res,"DESeq2/legs_DESeq_F1_MPCP_female_old_all.txt")

F1_MPCP_female_old_sig_res <- read.table("DESeq2/legs_DESeq_F1_MPCP_female_old_sig.txt")
res <- read.table("DESeq2/legs_DESeq_F1_MPCP_female_old_all.txt")


# separate chr 21 and autosomes
chr <- read.table("chromo_gene_info.txt",header = TRUE)
res$gene_id <- row.names(res)
res <- left_join(res,chr)

resauto <- res[-which(res$scaffold=="chr21"),]
resautosig <- resauto[which(resauto$padj<0.05),]
write.table(resautosig,"DESeq2/legs_DESeq_F1_MPCP_female_old_sig_autosomes.txt")


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

# create dataset for DESeq2
dds <- DESeqDataSetFromMatrix(countData = F1_MPCP, colData = coldata, design = ~individual + allele)

# do not normalize counts
# Because we have sample in the design, 
# and so we are looking at ratios within sample, we do not want 
# to use DESeq2's size factor normalization, and hence set the 
# size factors all to 1.
sizeFactors(dds) <- rep(1, length(coldata$individual))

#run DESeq
dds <- DESeq(dds,fitType = "local")

# check and set contrasts
resultsNames(dds)
res <- results(dds,contrast = c("allele","melp","cyd"))

restable <- data.frame(row.names = rownames(res),res$log2FoldChange,res$padj)
restable[-which(grepl("HMEL",rownames(restable))),]
# plots

resLFC <- lfcShrink(dds, coef="allele_melp_vs_cyd", type="apeglm")
plotMA(resLFC, ylim=c(-3,3), main = "Differential Allele Expression - mature males")


# number of genes with padj < 0.05
sum(res$padj < 0.05, na.rm=TRUE)
# = 2134 genes have ASE

F1_MPCP_sig_res <- res[which(res$padj < 0.05),]

write.table(F1_MPCP_sig_res,"DESeq2/legs_DESeq_F1_MPCP_sig.txt")
write.table(res,"DESeq2/legs_DESeq_F1_MPCP_all.txt")

F1_MPCP_sig_res <- read.table("DESeq2/legs_DESeq_F1_MPCP_sig.txt")
res <- read.table("DESeq2/legs_DESeq_F1_MPCP_all.txt")


# separate chr 21 and autosomes
chr <- read.table("chromo_gene_info.txt",header = TRUE)
res$gene_id <- row.names(res)
res <- left_join(res,chr)

resauto <- res[-which(res$scaffold=="chr21"),]
resautosig <- resauto[which(resauto$padj<0.05),]
write.table(resautosig,"DESeq2/legs_DESeq_F1_MPCP_sig_autosomes.txt")


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
F1_MPCP_female_young <- read.table("ASE data prep/legs/5.legs_genes_melpcyd_F1_MPCP_female_young.txt")
coldata <- read.table("ASE data prep/legs/6.legs_coldata_F1_MPCP_female_young.txt")

coldata$individual <- as.factor(coldata$individual)
coldata$allele <- as.factor(coldata$allele)

# create dataset for DESeq2
dds <- DESeqDataSetFromMatrix(countData = F1_MPCP_female_young, colData = coldata, design = ~individual + allele)

# do not normalize counts
# Because we have sample in the design, 
# and so we are looking at ratios within sample, we do not want 
# to use DESeq2's size factor normalization, and hence set the 
# size factors all to 1.
sizeFactors(dds) <- rep(1, length(coldata$individual))

#run DESeq
dds <- DESeq(dds,fitType = "local")

# check and set contrasts
resultsNames(dds)
res <- results(dds,contrast = c("allele","melp","cyd"))

restable <- data.frame(row.names = rownames(res),res$log2FoldChange,res$padj)
restable[-which(grepl("HMEL",rownames(restable))),]
# plots

resLFC <- lfcShrink(dds, coef="allele_melp_vs_cyd", type="apeglm")
plotMA(resLFC, ylim=c(-3,3), main = "Differential Allele Expression - juvenile females")


# number of genes with padj < 0.05
sum(res$padj < 0.05, na.rm=TRUE)
# = 1045 genes have ASE

F1_MPCP_female_young_sig_res <- res[which(res$padj < 0.05),]

write.table(F1_MPCP_female_young_sig_res,"DESeq2/legs_DESeq_F1_MPCP_female_young_sig.txt")
write.table(res,"DESeq2/legs_DESeq_F1_MPCP_female_young_all.txt")

F1_MPCP_female_young_sig_res <- read.table("DESeq2/legs_DESeq_F1_MPCP_female_young_sig.txt")
res <- read.table("DESeq2/legs_DESeq_F1_MPCP_female_young_all.txt")


# separate chr 21 and autosomes
chr <- read.table("chromo_gene_info.txt",header = TRUE)
res$gene_id <- row.names(res)
res <- left_join(res,chr)

resauto <- res[-which(res$scaffold=="chr21"),]
resautosig <- resauto[which(resauto$padj<0.05),]
write.table(resautosig,"DESeq2/legs_DESeq_F1_MPCP_female_young_sig_autosomes.txt")


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
F1_MPCP_young <- read.table("ASE data prep/legs/5.legs_genes_melpcyd_F1_MPCP_young.txt")
coldata <- read.table("ASE data prep/legs/6.legs_coldata_F1_MPCP_young.txt")

coldata$individual <- as.factor(coldata$individual)
coldata$allele <- as.factor(coldata$allele)

# create dataset for DESeq2
dds <- DESeqDataSetFromMatrix(countData = F1_MPCP_young, colData = coldata, design = ~individual + allele)

# do not normalize counts
# Because we have sample in the design, 
# and so we are looking at ratios within sample, we do not want 
# to use DESeq2's size factor normalization, and hence set the 
# size factors all to 1.
sizeFactors(dds) <- rep(1, length(coldata$individual))

#run DESeq
dds <- DESeq(dds,fitType = "local")

# check and set contrasts
resultsNames(dds)
res <- results(dds,contrast = c("allele","melp","cyd"))

restable <- data.frame(row.names = rownames(res),res$log2FoldChange,res$padj)
restable[-which(grepl("HMEL",rownames(restable))),]
# plots

resLFC <- lfcShrink(dds, coef="allele_melp_vs_cyd", type="apeglm")
plotMA(resLFC, ylim=c(-3,3), main = "Differential Allele Expression - juvenile males")


# number of genes with padj < 0.05
sum(res$padj < 0.05, na.rm=TRUE)
# = 914 genes have ASE

F1_MPCP_young_sig_res <- res[which(res$padj < 0.05),]

write.table(F1_MPCP_young_sig_res,"DESeq2/legs_DESeq_F1_MPCP_young_sig.txt")
write.table(res,"DESeq2/legs_DESeq_F1_MPCP_young_all.txt")

F1_MPCP_young_sig_res <- read.table("DESeq2/legs_DESeq_F1_MPCP_young_sig.txt")
res <- read.table("DESeq2/legs_DESeq_F1_MPCP_young_all.txt")


# separate chr 21 and autosomes
chr <- read.table("chromo_gene_info.txt",header = TRUE)
res$gene_id <- row.names(res)
res <- left_join(res,chr)

resauto <- res[-which(res$scaffold=="chr21"),]
resautosig <- resauto[which(resauto$padj<0.05),]
write.table(resautosig,"DESeq2/legs_DESeq_F1_MPCP_young_sig_autosomes.txt")


# proportion of genes upregulated from each allele
resautosig <- read.table("DESeq2/legs_DESeq_F1_MPCP_young_sig_autosomes.txt")
ressig <- read.table("DESeq2/legs_DESeq_F1_MPCP_young_sig.txt")

# all chromosomes
sum(ressig$log2FoldChange > 0)/nrow(ressig)
#melp= 0.3873085

#autosomes
sum(resautosig$log2FoldChange > 0)/nrow(resautosig)
#melp= 0.4839199



  # venn diagram ====

fem_old<- read.table("4. DESeq2/legs_DESeq_F1_MPCP_female_old_sig_autosomes.txt")
male_old<- read.table("4. DESeq2/legs_DESeq_F1_MPCP_sig_autosomes.txt")
fem_young<- read.table("4. DESeq2/legs_DESeq_F1_MPCP_female_young_sig_autosomes.txt")
male_young<- read.table("4. DESeq2/legs_DESeq_F1_MPCP_young_sig_autosomes.txt")


fem_old <- c(fem_old$gene_id)
male_old <- c(male_old$gene_id)
fem_young <- c(fem_young$gene_id)
male_young <- c(male_young$gene_id)

color = c("#a61107","#4563a8","#a1c957","#7d599c")

v<-venn.diagram(
  x = list(male_old, fem_old, male_young, fem_young),
  height = 2000, width = 2250,
  category.names = c("mature males" , "mature females" , "young males", "young females"),
  fill=color,
  cat.cex = 1,
  cat.just=list(c(0.4,-1.5) , c(0.65,-1.5) , c(0.5,-1.5) , c(0.5,-1.5)),
  filename ='venn_diagramm_legs.png',
  imagetype = "png",
  output=TRUE
)

  
# lists for ASE heatmaps====

fem_old<- read.table("4. DESeq2/legs_DESeq_F1_MPCP_female_old_sig_autosomes.txt")
  row.names(fem_old)<- fem_old$gene_id
male_old<- read.table("4. DESeq2/legs_DESeq_F1_MPCP_sig_autosomes.txt")
  row.names(male_old)<- male_old$gene_id
fem_young<- read.table("4. DESeq2/legs_DESeq_F1_MPCP_female_young_sig_autosomes.txt")
  row.names(fem_young)<- fem_young$gene_id
male_young<- read.table("4. DESeq2/legs_DESeq_F1_MPCP_young_sig_autosomes.txt")
  row.names(male_young)<- male_young$gene_id

# merge tables so they are all in the same order
list <- full_join(fem_old,male_old,by="gene_id")
list <- full_join(list, fem_young,by="gene_id")
list <- full_join(list,male_young,by="gene_id")

row.names(list) <- list$gene_id

which(grepl(paste(chemogenes$ID,collapse = "|"),list$gene_id))

# separate back into groups
fem_old_ASE <- list[,c(1:6,8:10)]
male_old_ASE <- list[,c(11:19)]
fem_young_ASE <- list[,c(20:28)]
male_young_ASE <- list[,c(29:37)]

# make lists using genes expressed in old females
male_old_ASE <- male_old_ASE[-which(is.na(fem_old_ASE)),]
fem_young_ASE <- fem_young_ASE[-which(is.na(fem_old_ASE)),]
male_young_ASE <- male_young_ASE[-which(is.na(fem_old_ASE)),]
fem_old_ASE <- fem_old_ASE[-which(is.na(fem_old_ASE)),]

write.csv(fem_old_ASE, "4. DESeq2/legs_ASE_female_old_fo_genes")
write.csv(male_old_ASE, "4. DESeq2/legs_ASE_male_old_fo_genes")
write.csv(fem_young_ASE, "4. DESeq2/legs_ASE_female_young_fo_genes")
write.csv(male_young_ASE, "4. DESeq2/legs_ASE_male_young_fo_genes")

# make lists of genes unique to old females
fem_oldgenes <- rownames(fem_old)
male_oldgenes= rownames(male_old)
fem_younggenes = rownames(fem_young)
male_younggenes = rownames(male_young)

fem_old_uniqgenes = setdiff(fem_oldgenes,male_oldgenes)
fem_old_uniqgenes = setdiff(fem_old_uniqgenes,fem_younggenes)
fem_old_uniqgenes= setdiff(fem_old_uniqgenes,male_younggenes)

uniq = which(grepl(paste(fem_old_uniqgenes,collapse = "|"),row.names(fem_old)))
other<- setdiff((1:473),uniq)

fem_old_ASE_full<- fem_old_ASE[c(uniq,other),]
male_old_ASE_full<- male_old_ASE[c(uniq,other),]
fem_young_ASE_full<- fem_young_ASE[c(uniq,other),]
male_young_ASE_full<- male_young_ASE[c(uniq,other),]

write.csv(fem_old_ASE_full,"4. DESeq2/legs_ASE_female_old_full.csv")
write.csv(male_old_ASE_full,"4. DESeq2/legs_ASE_male_old_full.csv")
write.csv(fem_young_ASE_full,"4. DESeq2/legs_ASE_female_young_full.csv")
write.csv(male_young_ASE_full,"4. DESeq2/legs_ASE_male_young_full.csv")

fem_old_ASEuniq = filter(fem_old,grepl(paste(fem_old_uniqgenes,collapse = "|"),row.names(fem_old)))
  fem_old_ASEuniq<- fem_old_ASEuniq[c(41,1:40,42:43),]

write.csv(fem_old_ASEuniq,"4. DESeq2/legs_ASE_female_old_unique.csv")
write.csv(male_old_ASEuniq,"4. DESeq2/legs_ASE_male_old_unique.csv")
write.csv(fem_young_ASEuniq,"4. DESeq2/legs_ASE_female_young_unique.csv")
write.csv(male_young_ASEuniq,"4. DESeq2/legs_ASE_male_young_unique.csv")

sum(fem_old_ASEuniq$log2FoldChange > 0)/length(fem_old_ASEuniq$log2FoldChange)
#   0.3846154

# make list using only annotated genes present in old females
male_old_ASEsmol <- male_old_ASE[-which(grepl("HMEL",rownames(fem_old_ASE),ignore.case = FALSE)==TRUE),]
fem_young_ASEsmol <- fem_young_ASE[-which(grepl("HMEL",rownames(fem_old_ASE),ignore.case = FALSE)==TRUE),]
male_young_ASEsmol <- male_young_ASE[-which(grepl("HMEL",rownames(fem_old_ASE),ignore.case = FALSE)==TRUE),]
fem_old_ASEsmol <- fem_old_ASE[-which(grepl("HMEL",rownames(fem_old_ASE),ignore.case = FALSE)==TRUE),]

#write.csv(fem_old_ASEsmol, "4. DESeq2/legs_ASE_female_old_smol.csv")
#write.csv(male_old_ASEsmol, "4. DESeq2/legs_ASE_male_old_smol.csv")
#write.csv(fem_young_ASEsmol, "4. DESeq2/legs_ASE_female_young_smol.csv")
#write.csv(male_young_ASEsmol, "4. DESeq2/legs_ASE_male_young_smol.csv")

  # geneoverlap
if (!requireNamespace("BiocManager", quietly = TRUE))
install.packages("BiocManager")
BiocManager::install("GeneOverlap")
install.packages("overlapping")
library(GeneOverlap)
go.obj <- newGeneOverlap(rownames(fem_old), rownames(fem_young), genome.size = length(rownames(list)))
testgo.obj<-testGeneOverlap(go.obj)
print(testgo.obj)


# ANTENNAE ====
    ## F1_MPCP_female_old ====

# read files in
F1_MPCP_female_old <- read.table("3. ASE data prep/antennae/5.antennae_genes_melpcyd_F1_MPCP_female_old.txt")
coldata <- read.table("3. ASE data prep/antennae/6.antennae_coldata_F1_MPCP_female_old.txt")

coldata$individual <- as.factor(coldata$individual)
coldata$allele <- as.factor(coldata$allele)

# create dataset for DESeq2
dds <- DESeqDataSetFromMatrix(countData = F1_MPCP_female_old, colData = coldata, design = ~individual + allele)

# do not normalize counts
# Because we have sample in the design, 
# and so we are looking at ratios within sample, we do not want 
# to use DESeq2's size factor normalization, and hence set the 
# size factors all to 1.
sizeFactors(dds) <- rep(1, length(coldata$individual))

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
# = 181 genes have ASE

F1_MPCP_female_old_sig_res <- res[which(res$padj < 0.05),]

write.table(F1_MPCP_female_old_sig_res,"DESeq2/antennae_DESeq_F1_MPCP_female_old_sig.txt")
write.table(res,"DESeq2/antennae_DESeq_F1_MPCP_female_old_all.txt")

F1_MPCP_female_old_sig_res <- read.table("DESeq2/antennae_DESeq_F1_MPCP_female_old_sig.txt")
res <- read.table("DESeq2/antennae_DESeq_F1_MPCP_female_old_all.txt")


# separate chr 21 and autosomes
chr <- read.table("chromo_gene_info.txt",header = TRUE)
res$gene_id <- row.names(res)
res <- left_join(res,chr)

resauto <- res[-which(res$scaffold=="chr21"),]
resautosig <- resauto[which(resauto$padj<0.05),]
write.table(resautosig,"DESeq2/antennae_DESeq_F1_MPCP_female_old_sig_autosomes.txt")


# proportion of genes upregulated from each allele
resautosig <- read.table("DESeq2/antennae_DESeq_F1_MPCP_female_old_sig_autosomes.txt")
ressig <- read.table("DESeq2/antennae_DESeq_F1_MPCP_female_old_sig.txt")

# all chromosomes
sum(ressig$log2FoldChange > 0)/nrow(ressig)
#melp= 0.4198895

#autosomes
sum(resautosig$log2FoldChange > 0)/nrow(resautosig)
#melp= 0.3892216


    ## F1_MPCP ====

# read files in
F1_MPCP <- read.table("ASE data prep/antennae/5.antennae_genes_melpcyd_F1_MPCP.txt")
coldata <- read.table("ASE data prep/antennae/6.antennae_coldata_F1_MPCP.txt")

coldata$individual <- as.factor(coldata$individual)
coldata$allele <- as.factor(coldata$allele)

# create dataset for DESeq2
dds <- DESeqDataSetFromMatrix(countData = F1_MPCP, colData = coldata, design = ~individual + allele)

# do not normalize counts
# Because we have sample in the design, 
# and so we are looking at ratios within sample, we do not want 
# to use DESeq2's size factor normalization, and hence set the 
# size factors all to 1.
sizeFactors(dds) <- rep(1, length(coldata$individual))

#run DESeq
dds <- DESeq(dds,fitType = "local")

# check and set contrasts
resultsNames(dds)
res <- results(dds,contrast = c("allele","melp","cyd"))

restable <- data.frame(row.names = rownames(res),res$log2FoldChange,res$padj)
restable[-which(grepl("HMEL",rownames(restable))),]
# plots

resLFC <- lfcShrink(dds, coef="allele_melp_vs_cyd", type="apeglm")
plotMA(resLFC, ylim=c(-3,3), main = "Differential Allele Expression - mature males")


# number of genes with padj < 0.05
sum(res$padj < 0.05, na.rm=TRUE)
# = 3240 genes have ASE

F1_MPCP_sig_res <- res[which(res$padj < 0.05),]

write.table(F1_MPCP_sig_res,"DESeq2/antennae_DESeq_F1_MPCP_sig.txt")
write.table(res,"DESeq2/antennae_DESeq_F1_MPCP_all.txt")

F1_MPCP_sig_res <- read.table("DESeq2/antennae_DESeq_F1_MPCP_sig.txt")
res <- read.table("DESeq2/antennae_DESeq_F1_MPCP_all.txt")


# separate chr 21 and autosomes
chr <- read.table("chromo_gene_info.txt",header = TRUE)
res$gene_id <- row.names(res)
res <- left_join(res,chr)

resauto <- res[-which(res$scaffold=="chr21"),]
resautosig <- resauto[which(resauto$padj<0.05),]
write.table(resautosig,"DESeq2/antennae_DESeq_F1_MPCP_sig_autosomes.txt")


# proportion of genes upregulated from each allele
resautosig <- read.table("DESeq2/antennae_DESeq_F1_MPCP_sig_autosomes.txt")
ressig <- read.table("DESeq2/antennae_DESeq_F1_MPCP_sig.txt")

# all chromosomes
sum(ressig$log2FoldChange > 0)/nrow(ressig)
#melp= 0.4552469

#autosomes
sum(resautosig$log2FoldChange > 0)/nrow(resautosig)
#melp= 0.4528662


    ## F1_MPCP_female_young ====

# read files in
F1_MPCP_female_young <- read.table("ASE data prep/antennae/5.antennae_genes_melpcyd_F1_MPCP_female_young.txt")
coldata <- read.table("ASE data prep/antennae/6.antennae_coldata_F1_MPCP_female_young.txt")

coldata$individual <- as.factor(coldata$individual)
coldata$allele <- as.factor(coldata$allele)

# create dataset for DESeq2
dds <- DESeqDataSetFromMatrix(countData = F1_MPCP_female_young, colData = coldata, design = ~individual + allele)

# do not normalize counts
# Because we have sample in the design, 
# and so we are looking at ratios within sample, we do not want 
# to use DESeq2's size factor normalization, and hence set the 
# size factors all to 1.
sizeFactors(dds) <- rep(1, length(coldata$individual))

#run DESeq
dds <- DESeq(dds,fitType = "local")

# check and set contrasts
resultsNames(dds)
res <- results(dds,contrast = c("allele","melp","cyd"))

restable <- data.frame(row.names = rownames(res),res$log2FoldChange,res$padj)
restable[-which(grepl("HMEL",rownames(restable))),]
# plots

resLFC <- lfcShrink(dds, coef="allele_melp_vs_cyd", type="apeglm")
plotMA(resLFC, ylim=c(-3,3), main = "Differential Allele Expression - juvenile females")


# number of genes with padj < 0.05
sum(res$padj < 0.05, na.rm=TRUE)
# = 688 genes have ASE

F1_MPCP_female_young_sig_res <- res[which(res$padj < 0.05),]

write.table(F1_MPCP_female_young_sig_res,"DESeq2/antennae_DESeq_F1_MPCP_female_young_sig.txt")
write.table(res,"DESeq2/antennae_DESeq_F1_MPCP_female_young_all.txt")

F1_MPCP_female_young_sig_res <- read.table("DESeq2/antennae_DESeq_F1_MPCP_female_young_sig.txt")
res <- read.table("DESeq2/antennae_DESeq_F1_MPCP_female_young_all.txt")


# separate chr 21 and autosomes
chr <- read.table("chromo_gene_info.txt",header = TRUE)
res$gene_id <- row.names(res)
res <- left_join(res,chr)

resauto <- res[-which(res$scaffold=="chr21"),]
resautosig <- resauto[which(resauto$padj<0.05),]
write.table(resautosig,"DESeq2/antennae_DESeq_F1_MPCP_female_young_sig_autosomes.txt")


# proportion of genes upregulated from each allele
resautosig <- read.table("DESeq2/antennae_DESeq_F1_MPCP_female_young_sig_autosomes.txt")
ressig <- read.table("DESeq2/antennae_DESeq_F1_MPCP_female_young_sig.txt")

# all chromosomes
sum(ressig$log2FoldChange > 0)/nrow(ressig)
#melp = 0.3851744

#autosomes
sum(resautosig$log2FoldChange > 0)/nrow(resautosig)
#melp = 0.3849693


    ## F1_MPCP_young ====

# read files in
F1_MPCP_young <- read.table("ASE data prep/antennae/5.antennae_genes_melpcyd_F1_MPCP_young.txt")
coldata <- read.table("ASE data prep/antennae/6.antennae_coldata_F1_MPCP_young.txt")

coldata$individual <- as.factor(coldata$individual)
coldata$allele <- as.factor(coldata$allele)

# create dataset for DESeq2
dds <- DESeqDataSetFromMatrix(countData = F1_MPCP_young, colData = coldata, design = ~individual + allele)

# do not normalize counts
# Because we have sample in the design, 
# and so we are looking at ratios within sample, we do not want 
# to use DESeq2's size factor normalization, and hence set the 
# size factors all to 1.
sizeFactors(dds) <- rep(1, length(coldata$individual))

#run DESeq
dds <- DESeq(dds,fitType = "local")

# check and set contrasts
resultsNames(dds)
res <- results(dds,contrast = c("allele","melp","cyd"))

restable <- data.frame(row.names = rownames(res),res$log2FoldChange,res$padj)
restable[-which(grepl("HMEL",rownames(restable))),]
# plots

resLFC <- lfcShrink(dds, coef="allele_melp_vs_cyd", type="apeglm")
plotMA(resLFC, ylim=c(-3,3), main = "Differential Allele Expression - juvenile males")


# number of genes with padj < 0.05
sum(res$padj < 0.05, na.rm=TRUE)
# = 781 genes have ASE

F1_MPCP_young_sig_res <- res[which(res$padj < 0.05),]

write.table(F1_MPCP_young_sig_res,"DESeq2/antennae_DESeq_F1_MPCP_young_sig.txt")
write.table(res,"DESeq2/antennae_DESeq_F1_MPCP_young_all.txt")

F1_MPCP_young_sig_res <- read.table("DESeq2/antennae_DESeq_F1_MPCP_young_sig.txt")
res <- read.table("DESeq2/antennae_DESeq_F1_MPCP_young_all.txt")


# separate chr 21 and autosomes
chr <- read.table("chromo_gene_info.txt",header = TRUE)
res$gene_id <- row.names(res)
res <- left_join(res,chr)

resauto <- res[-which(res$scaffold=="chr21"),]
resautosig <- resauto[which(resauto$padj<0.05),]
write.table(resautosig,"DESeq2/antennae_DESeq_F1_MPCP_young_sig_autosomes.txt")


# proportion of genes upregulated from each allele
resautosig <- read.table("DESeq2/antennae_DESeq_F1_MPCP_young_sig_autosomes.txt")
ressig <- read.table("DESeq2/antennae_DESeq_F1_MPCP_young_sig.txt")

# all chromosomes
sum(ressig$log2FoldChange > 0)/nrow(ressig)
#melp= 0.4084507

#autosomes
sum(resautosig$log2FoldChange > 0)/nrow(resautosig)
#melp= 0.4050465



  # venn diagram ====

fem_old<- read.table("4. DESeq2/antennae_DESeq_F1_MPCP_female_old_sig_autosomes.txt")
male_old<- read.table("4. DESeq2/antennae_DESeq_F1_MPCP_sig_autosomes.txt")
fem_young<- read.table("4. DESeq2/antennae_DESeq_F1_MPCP_female_young_sig_autosomes.txt")
male_young<- read.table("4. DESeq2/antennae_DESeq_F1_MPCP_young_sig_autosomes.txt")

fem_old <- c(fem_old$gene_id)
male_old <- c(male_old$gene_id)
fem_young <- c(fem_young$gene_id)
male_young <- c(male_young$gene_id)

color = c("#a61107","#4563a8","#a1c957","#7d599c")

venn.diagram(
  x = list(male_old, fem_old, male_young, fem_young),
  height = 2000, width = 2250,
  category.names = c("mature males" , "mature females" , "young males", "young females"),
  fill=color,
  cat.cex = 1,
  cat.just=list(c(0.4,-1.5) , c(0.65,-1.5) , c(0.5,-1.5) , c(0.5,-1.5)),
  filename = 'venn_diagramm_antennae.png',
  imagetype = "png",
  output=TRUE
)

  #lists for ASE heatmaps ====

fem_old<- read.table("4. DESeq2/antennae_DESeq_F1_MPCP_female_old_sig_autosomes.txt")
  row.names(fem_old)<- fem_old$gene_id
male_old<- read.table("4. DESeq2/antennae_DESeq_F1_MPCP_sig_autosomes.txt")
  row.names(male_old)<- male_old$gene_id
fem_young<- read.table("4. DESeq2/antennae_DESeq_F1_MPCP_female_young_sig_autosomes.txt")
  row.names(fem_young)<- fem_young$gene_id
male_young<- read.table("4. DESeq2/antennae_DESeq_F1_MPCP_young_sig_autosomes.txt")
  row.names(male_young)<- male_young$gene_id

# merge tables so they are all in the same order
list <- full_join(fem_old,male_old,by="gene_id")
list <- full_join(list, fem_young,by="gene_id")
list <- full_join(list,male_young,by="gene_id")

row.names(list) <- list$gene_id

# separate back into groups
fem_old_ASE <- list[,c(1:6,8:10)]
male_old_ASE <- list[,c(11:19)]
fem_young_ASE <- list[,c(20:28)]
male_young_ASE <- list[,c(29:37)]

# make lists using genes expressed in old females
male_old_ASE <- male_old_ASE[-which(is.na(fem_old_ASE)),]
fem_young_ASE <- fem_young_ASE[-which(is.na(fem_old_ASE)),]
male_young_ASE <- male_young_ASE[-which(is.na(fem_old_ASE)),]
fem_old_ASE <- fem_old_ASE[-which(is.na(fem_old_ASE)),]

write.csv(fem_old_ASE, "4. DESeq2/antennae_ASE_female_old_fo_genes.csv")
write.csv(male_old_ASE, "4. DESeq2/antennae_ASE_male_old_fo_genes.csv")
write.csv(fem_young_ASE, "4. DESeq2/antennae_ASE_female_young_fo_genes.csv")
write.csv(male_young_ASE, "4. DESeq2/antennae_ASE_male_young_fo_genes.csv")


sum(fem_old$log2FoldChange.x > 0)/length(fem_old$log2FoldChange.x)
  # 0.3892216


# make list using only annotated genes
fem_old_ASEsmol <- fem_old_ASEuniq[-which(grepl("HMEL",rownames(fem_old_ASEuniq),ignore.case = FALSE)==TRUE),]
male_old_ASEsmol <- male_old_ASEuniq[-which(grepl("HMEL",rownames(fem_old_ASEuniq),ignore.case = FALSE)==TRUE),]
fem_young_ASEsmol <- fem_young_ASEuniq[-which(grepl("HMEL",rownames(fem_old_ASEuniq),ignore.case = FALSE)==TRUE),]
male_young_ASEsmol <- male_young_ASEuniq[-which(grepl("HMEL",rownames(fem_old_ASEuniq),ignore.case = FALSE)==TRUE),]


write.csv(fem_old_ASEsmol, "4. DESeq2/antennae_ASE_female_old_smol.csv")
write.csv(male_old_ASEsmol, "4. DESeq2/antennae_ASE_male_old_smol.csv")
write.csv(fem_young_ASEsmol, "4. DESeq2/antennae_ASE_female_young_smol.csv")
write.csv(male_young_ASEsmol, "4. DESeq2/antennae_ASE_male_young_smol.csv")
