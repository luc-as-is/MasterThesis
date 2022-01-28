library (tidyr)
library(dplyr)
library(DESeq2)
library(GenomicRanges)
library(apeglm)

setwd("~/lukeball/HostPref_R")

#LEGS ####
#-- preparing data ----
  #=== 1. melp reads ====

# read in melp allele ASE tables
melp_60 <- read.table("3. ASE data prep/ASE output/legs_007/60__007_melp.csv", header = TRUE)
melp_142 <- read.table("3. ASE data prep/ASE output/legs_007/142_007_melp.csv", header = TRUE)
melp_173 <- read.table("3. ASE data prep/ASE output/legs_007/173_007_melp.csv", header = TRUE)
melp_181 <- read.table("3. ASE data prep/ASE output/legs_007/181_007_melp.csv", header = TRUE)
melp_215 <- read.table("3. ASE data prep/ASE output/legs_007/215_007_melp.csv", header = TRUE)
melp_239 <- read.table("3. ASE data prep/ASE output/legs_007/239_007_melp.csv", header = TRUE)
melp_259 <- read.table("3. ASE data prep/ASE output/legs_007/259_007_melp.csv", header = TRUE)

melp_265 <- read.table("3. ASE data prep/ASE output/legs_011/265_011_melp.csv", header = TRUE)
melp_282 <- read.table("3. ASE data prep/ASE output/legs_011/282_011_melp.csv", header = TRUE)
melp_285 <- read.table("3. ASE data prep/ASE output/legs_011/285_011_melp.csv", header = TRUE)
melp_294 <- read.table("3. ASE data prep/ASE output/legs_011/294_011_melp.csv", header = TRUE)
melp_300 <- read.table("3. ASE data prep/ASE output/legs_011/300_011_melp.csv", header = TRUE)
melp_304 <- read.table("3. ASE data prep/ASE output/legs_011/304_011_melp.csv", header = TRUE)
melp_314 <- read.table("3. ASE data prep/ASE output/legs_011/314_011_melp.csv", header = TRUE)
melp_317 <- read.table("3. ASE data prep/ASE output/legs_011/317_011_melp.csv", header = TRUE)
melp_322 <- read.table("3. ASE data prep/ASE output/legs_011/322_011_melp.csv", header = TRUE)
melp_327 <- read.table("3. ASE data prep/ASE output/legs_011/327_011_melp.csv", header = TRUE)
melp_334 <- read.table("3. ASE data prep/ASE output/legs_011/334_011_melp.csv", header = TRUE)
melp_348 <- read.table("3. ASE data prep/ASE output/legs_011/348_011_melp.csv", header = TRUE)
melp_359 <- read.table("3. ASE data prep/ASE output/legs_011/359_011_melp.csv", header = TRUE)

#combine into one table
melp<-full_join(melp_142,melp_173, by=c("position","contig"))  
melp<-full_join(melp,melp_181, by=c("position","contig"))  
melp<-full_join(melp,melp_215, by=c("position","contig"))  
melp<-full_join(melp,melp_239, by=c("position","contig"))  
melp<-full_join(melp,melp_259, by=c("position","contig"))  

melp<-full_join(melp,melp_265, by=c("position","contig"))  
melp<-full_join(melp,melp_282, by=c("position","contig"))  
melp<-full_join(melp,melp_285, by=c("position","contig"))  
melp<-full_join(melp,melp_294, by=c("position","contig"))  
melp<-full_join(melp,melp_300, by=c("position","contig"))  
melp<-full_join(melp,melp_304, by=c("position","contig"))  
melp<-full_join(melp,melp_314, by=c("position","contig"))  
melp<-full_join(melp,melp_317, by=c("position","contig"))  
melp<-full_join(melp,melp_322, by=c("position","contig"))  
melp<-full_join(melp,melp_327, by=c("position","contig"))  
melp<-full_join(melp,melp_334, by=c("position","contig"))  
melp<-full_join(melp,melp_348, by=c("position","contig"))
melp<-full_join(melp,melp_359, by=c("position","contig"))
melp<-full_join(melp,melp_60, by=c("position","contig"))

# keep only scaffold, position, and ref and alt count data
melp <-melp[,c(1,2,which(grepl("Count",colnames(melp))))]
melp <- melp[,-which(grepl("total", colnames(melp)))]

# rename columns
colnames(melp)[1] <- "scaffold"
# replace "xy." in the column name with id 
id <- c(142,173,181,215,239,259,265,282,285,294,300,304,314,317,322,327,334,348,359,60)
id <- rep(id,each=2)
colnames(melp)[3:length(colnames(melp))] <- str_remove_all( paste0(id,colnames(melp)[3:length(colnames(melp))]), "[xy.]")

#write.table(melp,"ASE data prep/1.ASE_counts_legs_melp.txt")

  #=== 1. cyd reads ====

# read in cydno ASE tables
F1_60 <- read.table("3. ASE data prep/ASE output/legs_007/60__007_cyd.csv", header = TRUE)
F1_142 <- read.table("3. ASE data prep/ASE output/legs_007/142_007_cyd.csv", header = TRUE)
F1_173 <- read.table("3. ASE data prep/ASE output/legs_007/173_007_cyd.csv", header = TRUE)
F1_181 <- read.table("3. ASE data prep/ASE output/legs_007/181_007_cyd.csv", header = TRUE)
F1_215 <- read.table("3. ASE data prep/ASE output/legs_007/215_007_cyd.csv", header = TRUE)
F1_239 <- read.table("3. ASE data prep/ASE output/legs_007/239_007_cyd.csv", header = TRUE)
F1_259 <- read.table("3. ASE data prep/ASE output/legs_007/259_007_cyd.csv", header = TRUE)

F1_265 <- read.table("3. ASE data prep/ASE output/legs_011/265_011_cyd.csv", header = TRUE)
F1_282 <- read.table("3. ASE data prep/ASE output/legs_011/282_011_cyd.csv", header = TRUE)
F1_285 <- read.table("3. ASE data prep/ASE output/legs_011/285_011_cyd.csv", header = TRUE)
F1_294 <- read.table("3. ASE data prep/ASE output/legs_011/294_011_cyd.csv", header = TRUE)
F1_300 <- read.table("3. ASE data prep/ASE output/legs_011/300_011_cyd.csv", header = TRUE)
F1_304 <- read.table("3. ASE data prep/ASE output/legs_011/304_011_cyd.csv", header = TRUE)
F1_314 <- read.table("3. ASE data prep/ASE output/legs_011/314_011_cyd.csv", header = TRUE)
F1_317 <- read.table("3. ASE data prep/ASE output/legs_011/317_011_cyd.csv", header = TRUE)
F1_322 <- read.table("3. ASE data prep/ASE output/legs_011/322_011_cyd.csv", header = TRUE)
F1_327 <- read.table("3. ASE data prep/ASE output/legs_011/327_011_cyd.csv", header = TRUE)
F1_334 <- read.table("3. ASE data prep/ASE output/legs_011/334_011_cyd.csv", header = TRUE)
F1_348 <- read.table("3. ASE data prep/ASE output/legs_011/348_011_cyd.csv", header = TRUE)
F1_359 <- read.table("3. ASE data prep/ASE output/legs_011/359_011_cyd.csv", header = TRUE)


# combine into one data table, make sure they are in the same order as the melp data
cyd<-full_join(F1_142,F1_173, by = c("position", "contig"))
cyd<-full_join(cyd,F1_181, by = c("position", "contig"))  
cyd<-full_join(cyd,F1_215, by = c("position", "contig"))  
cyd<-full_join(cyd,F1_239, by = c("position", "contig") )  
cyd<-full_join(cyd,F1_259, by = c("position", "contig") )
cyd<-full_join(cyd,F1_265, by = c("position", "contig") )
cyd<-full_join(cyd,F1_282, by = c("position", "contig") )  
cyd<-full_join(cyd,F1_285, by = c("position", "contig") )  
cyd<-full_join(cyd,F1_294, by = c("position", "contig") )  
cyd<-full_join(cyd,F1_300, by = c("position", "contig") )
cyd<-full_join(cyd,F1_304, by = c("position", "contig"))
cyd<-full_join(cyd,F1_314, by = c("position", "contig"))
cyd<-full_join(cyd,F1_317, by = c("position", "contig"))
cyd<-full_join(cyd,F1_322, by = c("position", "contig"))  
cyd<-full_join(cyd,F1_327, by = c("position", "contig"))
cyd<-full_join(cyd,F1_334, by = c("position", "contig"))
cyd<-full_join(cyd,F1_348, by = c("position", "contig"))
cyd<-full_join(cyd,F1_359, by = c("position", "contig"))
cyd<-full_join(cyd,F1_60, by = c("position", "contig"))

cyd<-cyd[,c(1,2,which(grepl("Count",colnames(cyd))))]
cyd <- cyd[,-which(grepl("total", colnames(cyd)))]
colnames(cyd)[1] <- "scaffold"

id <- c(142,173,181,215,239,259,265,282,285,294,300,304,314,317,322,327,334,348,359,60)
id <- rep(id,each=2)
colnames(cyd)[3:length(colnames(cyd))] <- str_remove_all( paste0(id,colnames(cyd)[3:length(colnames(cyd))]), "[xy.]")

write.table(cyd, "ASE data prep/1.ASE_counts_legs_cyd.txt")

  #=== 1b. avg number of ASE snps per brood ====
melpb7 <- (nrow(melp_142)+nrow(melp_173)+nrow(melp_181)+nrow(melp_215)+
  nrow(melp_239)+ nrow(melp_259)+nrow(melp_60))/7
sd(c(nrow(melp_142),nrow(melp_173),nrow(melp_181),nrow(melp_215),
      nrow(melp_239), nrow(melp_259),nrow(melp_60)))

cydb7 <- (nrow(F1_142)+nrow(F1_173)+nrow(F1_181)+nrow(F1_215)+
             nrow(F1_239)+ nrow(F1_259)+nrow(F1_60))/7
sd(c(nrow(F1_142),nrow(F1_173),nrow(F1_181),nrow(F1_215),
     nrow(F1_239), nrow(F1_259),nrow(F1_60)))

melpb11 <- (nrow(melp_265)+nrow(melp_282)+nrow(melp_285)+nrow(melp_294)+
              nrow(melp_300)+nrow(melp_304)+nrow(melp_314)+nrow(melp_317)+
              nrow(melp_322)+nrow(melp_327)+nrow(melp_334)+nrow(melp_348)+
              nrow(melp_359))/13
sd(c(nrow(melp_265),nrow(melp_282),nrow(melp_285),nrow(melp_294),
    nrow(melp_300),nrow(melp_304),nrow(melp_314),nrow(melp_317),
    nrow(melp_322),nrow(melp_327),nrow(melp_334),nrow(melp_348),
    nrow(melp_359)))
  
cydb11 <- (nrow(F1_265)+nrow(F1_282)+nrow(F1_285)+nrow(F1_294)+
              nrow(F1_300)+nrow(F1_304)+nrow(F1_314)+nrow(F1_317)+
              nrow(F1_322)+nrow(F1_327)+nrow(F1_334)+nrow(F1_348)+
              nrow(F1_359))/13
sd(c(nrow(F1_265),nrow(F1_282),nrow(F1_285),nrow(F1_294),
     nrow(F1_300),nrow(F1_304),nrow(F1_314),nrow(F1_317),
     nrow(F1_322),nrow(F1_327),nrow(F1_334),nrow(F1_348),
     nrow(F1_359)))
  #=== 2. assign genes to the snps ----

# read in gene annotation
annotation <- read.table("gene_info.txt",header = TRUE)

# create ranges of where genes are
GRgenes<-makeGRangesFromDataFrame(annotation, seqnames.field=c("scaffold"))


# read in ASE count tables 
# ***do melp and cyd reads separately*** #

#choose to read in melp or cyd reads
legs <- read.table("3. ASE data prep/legs/1.ASE_counts_legs_melp.txt")
legs <- read.table("3. ASE data prep/legs/1.ASE_counts_legs_cyd.txt")

# add "Start" and "end" position (same as start because its a snp), needed by GRanges
colnames(legs)[2] <- "start"
legs$end = legs$start
GRsnps<-makeGRangesFromDataFrame(legs, seqnames.field=c("scaffold"))

#find overlaps between snps and genes
olaps <- findOverlaps(GRgenes, GRsnps)
legs_snps<-cbind(annotation[queryHits(olaps),], legs[subjectHits(olaps),])
#rename column
colnames(legs_snps)[6] <-"position"

#save tables
#write.table(legs_snps,"ASE data prep/2.legs_snps_melp.txt")
write.table(legs_snps,"ASE data prep/2.legs_snps_cyd.txt")

  #=== 3. subset into groups ====

# *** do melp and cyd separately ***

# read in melp or cyd 
legs_snps <- read.table("3. ASE data prep/legs/2.legs_snps_melp.txt")
legs_snps<- read.table("3. ASE data prep/legs/2.legs_snps_cyd.txt")

#F1_MPCP
legs_F1_MPCP <- legs_snps[,c(1:4,6,which(grepl("181",colnames(legs_snps))),which(grepl("282",colnames(legs_snps))),
                             which(grepl("294",colnames(legs_snps))),which(grepl("300",colnames(legs_snps))),
                             which(grepl("334",colnames(legs_snps))),which(grepl("348",colnames(legs_snps))))]

#F1_MPCP_female_old     
legs_F1_MPCP_female_old <- legs_snps[,c(1:4,6,which(grepl("142",colnames(legs_snps))),which(grepl("215",colnames(legs_snps))),
                                        which(grepl("259",colnames(legs_snps))),which(grepl("327",colnames(legs_snps))),which(grepl("359",colnames(legs_snps))))]

#F1_MPCP_female_young
legs_F1_MPCP_female_young <- legs_snps[,c(1:4,6,which(grepl("265",colnames(legs_snps))),which(grepl("314",colnames(legs_snps))),
                                          which(grepl("322",colnames(legs_snps))),which(grepl("60",colnames(legs_snps))))]

#F1_MPCP_young
legs_F1_MPCP_young <- legs_snps[,c(1:4,6,which(grepl("173",colnames(legs_snps))),which(grepl("239",colnames(legs_snps))),
                                   which(grepl("285",colnames(legs_snps))),which(grepl("304",colnames(legs_snps))),
                                   which(grepl("317",colnames(legs_snps))))]


#write.table(legs_F1_MPCP, file = "ASE data prep/3.legs_snps_F1_MPCP_melp.txt")
write.table(legs_F1_MPCP, file="ASE data prep/3.legs_snps_F1_MPCP_cyd.txt")

#write.table(legs_F1_MPCP_female_old, file = "ASE data prep/3.legs_snps_F1_MPCP_female_old_melp.txt")
write.table(legs_F1_MPCP_female_old, file="ASE data prep/3.legs_snps_F1_MPCP_female_old_cyd.txt")

#write.table(legs_F1_MPCP_young, file = "ASE data prep/3.legs_snps_F1_MPCP_young_melp.txt")
write.table(legs_F1_MPCP_young, file="ASE data prep/3.legs_snps_F1_MPCP_young_cyd.txt")

#write.table(legs_F1_MPCP_female_young, file = "ASE data prep/3.legs_snps_F1_MPCP_female_young_melp.txt")
write.table(legs_F1_MPCP_female_young, file="ASE data prep/3.legs_snps_F1_MPCP_female_young_cyd.txt")


  #=== 4. sum all SNPs of the same gene together ====

# F1_MPCP

#choose melp or cyd snp counts
#legs_F1_MPCP <- read.table("ASE data prep/3.legs_snps_F1_MPCP_melp.txt")
#legs_F1_MPCP <- read.table("ASE data prep/3.legs_snps_F1_MPCP_cyd.txt")

# replace NA with 0
legs_F1_MPCP[is.na(legs_F1_MPCP)] <- 0

# sum up the rows with the same gene_id
ASE_genes_legs_F1_MPCP <- rowsum(legs_F1_MPCP[6:length(colnames(legs_F1_MPCP))], legs_F1_MPCP$gene_id)

# change column names
id <- c("F1_181_ref","F1_181_alt","F1_282_ref","F1_282_alt","F1_294_ref","F1_294_alt",
        "F1_300_ref","F1_300_alt","F1_334_ref","F1_334_alt","F1_348_ref","F1_348_alt")

colnames(ASE_genes_legs_F1_MPCP)<-id


#write.table(ASE_genes_legs_F1_MPCP,"ASE data prep/4.legs_genes_F1_MPCP_melp.txt")
#write.table(ASE_genes_legs_F1_MPCP,"ASE data prep/4.legs_genes_F1_MPCP_cyd.txt")


# F1_MPCP_female_old

#read in melp/cyd snp counts
#legs_F1_MPCP_female_old <- read.table("ASE data prep/3.legs_snps_F1_MPCP_female_old_melp.txt")
legs_F1_MPCP_female_old <- read.table("ASE data prep/3.legs_snps_F1_MPCP_female_old_cyd.txt")

# replace NA with 0
legs_F1_MPCP_female_old[is.na(legs_F1_MPCP_female_old)] <- 0

# sum up the rows with the same gene_id
ASE_genes_legs_F1_MPCP_female_old <- rowsum(legs_F1_MPCP_female_old[6:length(colnames(legs_F1_MPCP_female_old))], legs_F1_MPCP_female_old$gene_id)

# change column names
id <- c("F1_142_ref","F1_142_alt","F1_215_ref","F1_215_alt","F1_259_ref","F1_259_alt",
        "F1_327_ref","F1_327_alt","F1_359_ref","F1_359_alt")

colnames(ASE_genes_legs_F1_MPCP_female_old)<-id


#write.table(ASE_genes_legs_F1_MPCP_female_old,"ASE data prep/4.legs_genes_F1_MPCP_female_old_melp.txt")
#write.table(ASE_genes_legs_F1_MPCP_female_old,"ASE data prep/4.legs_genes_F1_MPCP_female_old_cyd.txt")


# F1_MPCP_female_young

#choose melp or cyd snps
#legs_F1_MPCP_female_young <- read.table("ASE data prep/3.legs_snps_F1_MPCP_female_young_melp.txt")
legs_F1_MPCP_female_young <- read.table("ASE data prep/3.legs_snps_F1_MPCP_female_young_cyd.txt")


legs_F1_MPCP_female_young[is.na(legs_F1_MPCP_female_young)] <- 0

ASE_genes_legs_F1_MPCP_female_young <- rowsum(legs_F1_MPCP_female_young[6:length(colnames(legs_F1_MPCP_female_young))], legs_F1_MPCP_female_young$gene_id)

id <- c("F1_265_ref","F1_265_alt","F1_314_ref","F1_314_alt",
        "F1_322_ref","F1_322_alt","F1_60_ref","F1_60_alt")

colnames(ASE_genes_legs_F1_MPCP_female_young)<-id

#write.table(ASE_genes_legs_F1_MPCP_female_young, "ASE data prep/4.legs_genes_F1_MPCP_female_young_melp.txt")
write.table(ASE_genes_legs_F1_MPCP_female_young, "ASE data prep/4.legs_genes_F1_MPCP_female_young_cyd.txt")

# F1_MPCP_young

legs_F1_MPCP_young <- read.table("ASE data prep/3.legs_snps_F1_MPCP_young_melp.txt")
legs_F1_MPCP_young <- read.table("ASE data prep/3.legs_snps_F1_MPCP_young_cyd.txt")

legs_F1_MPCP_young[is.na(legs_F1_MPCP_young)] <- 0

ASE_genes_legs_F1_MPCP_young <- rowsum(legs_F1_MPCP_young[6:length(colnames(legs_F1_MPCP_young))], legs_F1_MPCP_young$gene_id)

id <- c("F1_173_ref","F1_173_alt","F1_239_ref","F1_239_alt","F1_285_ref",
        "F1_285_alt","F1_304_ref","F1_304_alt","F1_317_ref","F1_317_alt")
colnames(ASE_genes_legs_F1_MPCP_young)<-id

#write.table(ASE_genes_legs_F1_MPCP_young, "ASE data prep/4.legs_genes_F1_MPCP_young_melp.txt")
write.table(ASE_genes_legs_F1_MPCP_young, "ASE data prep/4.legs_genes_F1_MPCP_young_cyd.txt")


  #=== 5/6. combine cyd and melp tables and creat coldata table for DESeq ====
    # F1_MPCP_female_old ----

# read in files
F1_MPCP_female_old_melp <- read.table("3. ASE data prep/legs/4.legs_genes_F1_MPCP_female_old_melp.txt")
F1_MPCP_female_old_cyd <- read.table("3. ASE data prep/legs/4.legs_genes_F1_MPCP_female_old_cyd.txt")

# assign melp/cyd to ref/alt alleles
colnames(F1_MPCP_female_old_melp) <- gsub("ref", "cyd", colnames(F1_MPCP_female_old_melp))
colnames(F1_MPCP_female_old_melp) <- gsub("alt", "melp", colnames(F1_MPCP_female_old_melp))

colnames(F1_MPCP_female_old_cyd) <- gsub("ref", "melp", colnames(F1_MPCP_female_old_cyd))
colnames(F1_MPCP_female_old_cyd) <- gsub("alt", "cyd", colnames(F1_MPCP_female_old_cyd))

#create column with gene id
F1_MPCP_female_old_melp$gene <- row.names(F1_MPCP_female_old_melp)
F1_MPCP_female_old_cyd$gene <- row.names(F1_MPCP_female_old_cyd)

# combine cyd and melp snps into one table
# x = melp, y = cyd 
F1_MPCP_female_old <- full_join(F1_MPCP_female_old_melp,F1_MPCP_female_old_cyd,by="gene")
# make the gene names the row names, remove "gene" column
row.names(F1_MPCP_female_old) <- F1_MPCP_female_old$gene
F1_MPCP_female_old <- F1_MPCP_female_old[,-which(grepl("gene",colnames(F1_MPCP_female_old)))]
# sum melp/cyd counts for each individual (sum columns with the same name)
colnames(F1_MPCP_female_old)<- gsub(".x$","",colnames(F1_MPCP_female_old))
colnames(F1_MPCP_female_old)<- gsub(".y$","",colnames(F1_MPCP_female_old))
F1_MPCP_female_old[is.na(F1_MPCP_female_old)] <- 0
F1_MPCP_female_old <- t(rowsum(t(F1_MPCP_female_old), group = colnames(F1_MPCP_female_old)))

write.table(F1_MPCP_female_old,"ASE data prep/5.legs_genes_melpcyd_F1_MPCP_female_old.txt")

#create ColData table
#assign sample id to each column
individual <- factor(substr(colnames(F1_MPCP_female_old),4,6))
#assign melp or cyd alleles to each column
allele <- factor(substr(colnames(F1_MPCP_female_old),8,11))

coldata <- data.frame(row.names=colnames(F1_MPCP_female_old),individual,allele)

brood <- c(rep("007",6),rep("011",4))

coldata$brood <- factor(brood)

write.table(coldata,"ASE data prep/6.coldata_F1_MPCP_female_old.txt")

    # male old ====
# read in files
F1_MPCP_melp <- read.table("ASE data prep/4.legs_genes_F1_MPCP_melp.txt")
F1_MPCP_cyd <- read.table("ASE data prep/4.legs_genes_F1_MPCP_cyd.txt")

# assign melp/cyd to ref/alt alleles
colnames(F1_MPCP_melp) <- gsub("ref", "cyd", colnames(F1_MPCP_melp))
colnames(F1_MPCP_melp) <- gsub("alt", "melp", colnames(F1_MPCP_melp))

colnames(F1_MPCP_cyd) <- gsub("ref", "melp", colnames(F1_MPCP_cyd))
colnames(F1_MPCP_cyd) <- gsub("alt", "cyd", colnames(F1_MPCP_cyd))

#create column with gene id
F1_MPCP_melp$gene <- row.names(F1_MPCP_melp)
F1_MPCP_cyd$gene <- row.names(F1_MPCP_cyd)

# combine cyd and melp snps into one table
# x = melp, y = cyd 
F1_MPCP <- full_join(F1_MPCP_melp,F1_MPCP_cyd,by="gene")
# make the gene names the row names, remove "gene" column
row.names(F1_MPCP) <- F1_MPCP$gene
F1_MPCP <- F1_MPCP[,-which(grepl("gene",colnames(F1_MPCP)))]
# sum melp/cyd counts for each individual (sum columns with the same name)
colnames(F1_MPCP)<- gsub(".x$","",colnames(F1_MPCP))
colnames(F1_MPCP)<- gsub(".y$","",colnames(F1_MPCP))
F1_MPCP[is.na(F1_MPCP)] <- 0
F1_MPCP <- t(rowsum(t(F1_MPCP), group = colnames(F1_MPCP)))

write.table(F1_MPCP,"ASE data prep/5.legs_genes_melpcyd_F1_MPCP.txt")

#create ColData table
#assign sample id to each column
individual <- factor(substr(colnames(F1_MPCP),4,6))
#assign melp or cyd alleles to each column
allele <- factor(substr(colnames(F1_MPCP),8,11))

brood <- factor(c(rep("007",2),rep("011",10)))

coldata <- data.frame(row.names=colnames(F1_MPCP),individual,allele,brood)

write.table(coldata,"ASE data prep/6.legs_coldata_F1_MPCP.txt")


    # female young ====
# read in files
F1_MPCP_female_young_melp <- read.table("ASE data prep/4.legs_genes_F1_MPCP_female_young_melp.txt")
F1_MPCP_female_young_cyd <- read.table("ASE data prep/4.legs_genes_F1_MPCP_female_young_cyd.txt")

# assign melp/cyd to ref/alt alleles
colnames(F1_MPCP_female_young_melp) <- gsub("ref", "cyd", colnames(F1_MPCP_female_young_melp))
colnames(F1_MPCP_female_young_melp) <- gsub("alt", "melp", colnames(F1_MPCP_female_young_melp))

colnames(F1_MPCP_female_young_cyd) <- gsub("ref", "melp", colnames(F1_MPCP_female_young_cyd))
colnames(F1_MPCP_female_young_cyd) <- gsub("alt", "cyd", colnames(F1_MPCP_female_young_cyd))

#create column with gene id
F1_MPCP_female_young_melp$gene <- row.names(F1_MPCP_female_young_melp)
F1_MPCP_female_young_cyd$gene <- row.names(F1_MPCP_female_young_cyd)

# combine cyd and melp snps into one table
# x = melp, y = cyd 
F1_MPCP_female_young <- full_join(F1_MPCP_female_young_melp,F1_MPCP_female_young_cyd,by="gene")
# make the gene names the row names, remove "gene" column
row.names(F1_MPCP_female_young) <- F1_MPCP_female_young$gene
F1_MPCP_female_young <- F1_MPCP_female_young[,-which(grepl("gene",colnames(F1_MPCP_female_young)))]
# sum melp/cyd counts for each individual (sum columns with the same name)
colnames(F1_MPCP_female_young)<- gsub(".x$","",colnames(F1_MPCP_female_young))
colnames(F1_MPCP_female_young)<- gsub(".y$","",colnames(F1_MPCP_female_young))
F1_MPCP_female_young[is.na(F1_MPCP_female_young)] <- 0
F1_MPCP_female_young <- t(rowsum(t(F1_MPCP_female_young), group = colnames(F1_MPCP_female_young)))

write.table(F1_MPCP_female_young,"ASE data prep/5.legs_genes_melpcyd_F1_MPCP_female_young.txt")

#create ColData table
#assign sample id to each column
individual <- factor(substr(colnames(F1_MPCP_female_young),4,6))
#assign melp or cyd alleles to each column
allele <- factor(substr(colnames(F1_MPCP_female_young),8,11))

brood <- factor(c(rep("011",6),rep("007",2)))

coldata <- data.frame(row.names=colnames(F1_MPCP_female_young),individual,allele,brood)

coldata$allele[7] <- "cyd"
coldata$allele[8] <- "melp"

write.table(coldata,"ASE data prep/6.legs_coldata_F1_MPCP_female_young.txt")
    # male young ====
# read in files
F1_MPCP_young_melp <- read.table("ASE data prep/4.legs_genes_F1_MPCP_young_melp.txt")
F1_MPCP_young_cyd <- read.table("ASE data prep/4.legs_genes_F1_MPCP_young_cyd.txt")

# assign melp/cyd to ref/alt alleles
colnames(F1_MPCP_young_melp) <- gsub("ref", "cyd", colnames(F1_MPCP_young_melp))
colnames(F1_MPCP_young_melp) <- gsub("alt", "melp", colnames(F1_MPCP_young_melp))

colnames(F1_MPCP_young_cyd) <- gsub("ref", "melp", colnames(F1_MPCP_young_cyd))
colnames(F1_MPCP_young_cyd) <- gsub("alt", "cyd", colnames(F1_MPCP_young_cyd))

#create column with gene id
F1_MPCP_young_melp$gene <- row.names(F1_MPCP_young_melp)
F1_MPCP_young_cyd$gene <- row.names(F1_MPCP_young_cyd)

# combine cyd and melp snps into one table
# x = melp, y = cyd 
F1_MPCP_young <- full_join(F1_MPCP_young_melp,F1_MPCP_young_cyd,by="gene")
# make the gene names the row names, remove "gene" column
row.names(F1_MPCP_young) <- F1_MPCP_young$gene
F1_MPCP_young <- F1_MPCP_young[,-which(grepl("gene",colnames(F1_MPCP_young)))]
# sum melp/cyd counts for each individual (sum columns with the same name)
colnames(F1_MPCP_young)<- gsub(".x$","",colnames(F1_MPCP_young))
colnames(F1_MPCP_young)<- gsub(".y$","",colnames(F1_MPCP_young))
F1_MPCP_young[is.na(F1_MPCP_young)] <- 0
F1_MPCP_young <- t(rowsum(t(F1_MPCP_young), group = colnames(F1_MPCP_young)))

write.table(F1_MPCP_young,"ASE data prep/5.legs_genes_melpcyd_F1_MPCP_young.txt")

#create ColData table
#assign sample id to each column
individual <- factor(substr(colnames(F1_MPCP_young),4,6))
#assign melp or cyd alleles to each column
allele <- factor(substr(colnames(F1_MPCP_young),8,11))

brood <- factor(brood <- c(rep("007",4),rep("011",6)))

coldata <- data.frame(row.names=colnames(F1_MPCP_young),individual,allele)

write.table(coldata,"ASE data prep/6.legs_coldata_F1_MPCP_young.txt")


#ANTENNAE ####
#-- preparing data ----
  #=== 1. melp reads ====

# read in melp allele ASE tables
melp_59P <- read.table("3. ASE data prep/ASE output//antennae_007/59P_007_melp.csv", header = TRUE)
melp_137 <- read.table("3. ASE data prep/ASE output//antennae_007/137_007_melp.csv", header = TRUE)
melp_141 <- read.table("3. ASE data prep/ASE output//antennae_007/141_007_melp.csv", header = TRUE)
melp_172 <- read.table("3. ASE data prep/ASE output//antennae_007/172_007_melp.csv", header = TRUE)
#melp_210 <- read.table("3. ASE data prep/ASE output//antennae_007/210_007_melp.csv", header = TRUE)
  # 210 is corrupted
melp_214 <- read.table("3. ASE data prep/ASE output//antennae_007/214_007_melp.csv", header = TRUE)

melp_261 <- read.table("3. ASE data prep/ASE output//antennae_011/261_011_melp.csv", header = TRUE)
melp_274 <- read.table("3. ASE data prep/ASE output//antennae_011/274_011_melp.csv", header = TRUE)
melp_284 <- read.table("3. ASE data prep/ASE output//antennae_011/284_011_melp.csv", header = TRUE)
melp_313 <- read.table("3. ASE data prep/ASE output//antennae_011/313_011_melp.csv", header = TRUE)
melp_316 <- read.table("3. ASE data prep/ASE output//antennae_011/316_011_melp.csv", header = TRUE)
melp_326 <- read.table("3. ASE data prep/ASE output//antennae_011/326_011_melp.csv", header = TRUE)
melp_333 <- read.table("3. ASE data prep/ASE output//antennae_011/333_011_melp.csv", header = TRUE)
melp_352 <- read.table("3. ASE data prep/ASE output//antennae_011/352_011_melp.csv", header = TRUE)
melp_356 <- read.table("3. ASE data prep/ASE output//antennae_011/356_011_melp.csv", header = TRUE)
melp_365 <- read.table("3. ASE data prep/ASE output//antennae_011/365_011_melp.csv", header = TRUE)


#combine into one table
melp<-full_join(melp_137,melp_141, by=c("position","contig"))  
melp<-full_join(melp,melp_172, by=c("position","contig"))  
melp<-full_join(melp,melp_214, by=c("position","contig"))  
melp<-full_join(melp,melp_261, by=c("position","contig"))  
melp<-full_join(melp,melp_274, by=c("position","contig"))  
melp<-full_join(melp,melp_284, by=c("position","contig"))  
melp<-full_join(melp,melp_313, by=c("position","contig"))  
melp<-full_join(melp,melp_316, by=c("position","contig"))  
melp<-full_join(melp,melp_326, by=c("position","contig"))  
melp<-full_join(melp,melp_333, by=c("position","contig"))  
melp<-full_join(melp,melp_352, by=c("position","contig"))  
melp<-full_join(melp,melp_356, by=c("position","contig"))  
melp<-full_join(melp,melp_365, by=c("position","contig"))  
melp<-full_join(melp,melp_59P, by=c("position","contig"))  

# keep only scaffold, position, and ref and alt count data
melp <-melp[,c(1,2,which(grepl("Count",colnames(melp))))]
melp <- melp[,-which(grepl("total", colnames(melp)))]

# rename columns
colnames(melp)[1] <- "scaffold"
# replace "xy." in the column name with id 
id <- c(137,141,172,214,261,274,284,313,316,326,333,352,356,365,59)
id <- rep(id,each=2)
colnames(melp)[3:length(colnames(melp))] <- str_remove_all( paste0(id,colnames(melp)[3:length(colnames(melp))]), "[xy.]")

write.table(melp,"ASE data prep/1.ASE_counts_antennae_melp.txt")

  #=== 1. cyd reads ====

# read in cydno ASE tables
cyd_59P <- read.table("3. ASE data prep/ASE output/antennae_007/59P_007_cyd.csv", header = TRUE)
cyd_137 <- read.table("3. ASE data prep/ASE output/antennae_007/137_007_cyd.csv", header = TRUE)
cyd_141 <- read.table("3. ASE data prep/ASE output/antennae_007/141_007_cyd.csv", header = TRUE)
cyd_172 <- read.table("3. ASE data prep/ASE output/antennae_007/172_007_cyd.csv", header = TRUE)
#cyd_210 <- read.table("3. ASE data prep/ASE output/antennae_007/210_007_cyd.csv", header = TRUE)
# 210 is corrupted
cyd_214 <- read.table("3. ASE data prep/ASE output/antennae_007/214_007_cyd.csv", header = TRUE)

cyd_261 <- read.table("3. ASE data prep/ASE output/antennae_011/261_011_cyd.csv", header = TRUE)
cyd_274 <- read.table("3. ASE data prep/ASE output/antennae_011/274_011_cyd.csv", header = TRUE)
cyd_284 <- read.table("3. ASE data prep/ASE output/antennae_011/284_011_cyd.csv", header = TRUE)
cyd_313 <- read.table("3. ASE data prep/ASE output/antennae_011/313_011_cyd.csv", header = TRUE)
cyd_316 <- read.table("3. ASE data prep/ASE output/antennae_011/316_011_cyd.csv", header = TRUE)
cyd_326 <- read.table("3. ASE data prep/ASE output/antennae_011/326_011_cyd.csv", header = TRUE)
cyd_333 <- read.table("3. ASE data prep/ASE output/antennae_011/333_011_cyd.csv", header = TRUE)
cyd_352 <- read.table("3. ASE data prep/ASE output/antennae_011/352_011_cyd.csv", header = TRUE)
cyd_356 <- read.table("3. ASE data prep/ASE output/antennae_011/356_011_cyd.csv", header = TRUE)
cyd_365 <- read.table("3. ASE data prep/ASE output/antennae_011/365_011_cyd.csv", header = TRUE)


#combine into one table
cyd<-full_join(cyd_137,cyd_141, by=c("position","contig"))  
cyd<-full_join(cyd,cyd_172, by=c("position","contig"))  
cyd<-full_join(cyd,cyd_214, by=c("position","contig"))  
cyd<-full_join(cyd,cyd_261, by=c("position","contig"))  
cyd<-full_join(cyd,cyd_274, by=c("position","contig"))  
cyd<-full_join(cyd,cyd_284, by=c("position","contig"))  
cyd<-full_join(cyd,cyd_313, by=c("position","contig"))  
cyd<-full_join(cyd,cyd_316, by=c("position","contig"))  
cyd<-full_join(cyd,cyd_326, by=c("position","contig"))  
cyd<-full_join(cyd,cyd_333, by=c("position","contig"))  
cyd<-full_join(cyd,cyd_352, by=c("position","contig"))  
cyd<-full_join(cyd,cyd_356, by=c("position","contig"))  
cyd<-full_join(cyd,cyd_365, by=c("position","contig"))  
cyd<-full_join(cyd,cyd_59P, by=c("position","contig"))  

# keep only scaffold, position, and ref and alt count data

cyd<-cyd[,c(1,2,which(grepl("Count",colnames(cyd))))]
cyd <- cyd[,-which(grepl("total", colnames(cyd)))]
colnames(cyd)[1] <- "scaffold"

id <- c(137,141,172,214,261,274,284,313,316,326,333,352,356,365,59)
id <- rep(id,each=2)
colnames(cyd)[3:length(colnames(cyd))] <- str_remove_all( paste0(id,colnames(cyd)[3:length(colnames(cyd))]), "[xy.]")

write.table(cyd, "ASE data prep/1.ASE_counts_antennae_cyd.txt")

  #=== 1b. avg number of ASE snps per brood ====
(nrow(melp_137)+nrow(melp_141)+nrow(melp_172)+nrow(melp_214)+
             nrow(melp_59P))/5
sd(c(nrow(melp_137),nrow(melp_141),nrow(melp_172),nrow(melp_214),
     nrow(melp_59P)))

(nrow(cyd_137)+nrow(cyd_141)+nrow(cyd_172)+nrow(cyd_214)+
    nrow(cyd_59P))/5
sd(c(nrow(cyd_137),nrow(cyd_141),nrow(cyd_172),nrow(cyd_214),
     nrow(cyd_59P)))

(nrow(melp_261)+nrow(melp_274)+nrow(melp_284)+nrow(melp_313)+
              nrow(melp_316)+nrow(melp_326)+nrow(melp_333)+nrow(melp_352)+
              nrow(melp_356)+nrow(melp_365))/10
sd(c(nrow(melp_261),nrow(melp_274),nrow(melp_284),nrow(melp_313),
       nrow(melp_316),nrow(melp_326),nrow(melp_333),nrow(melp_352),
       nrow(melp_356),nrow(melp_365)))

(nrow(cyd_261)+nrow(cyd_274)+nrow(cyd_284)+nrow(cyd_313)+
  nrow(cyd_316)+nrow(cyd_326)+nrow(cyd_333)+nrow(cyd_352)+
  nrow(cyd_356)+nrow(cyd_365))/10
sd(c(nrow(cyd_261),nrow(cyd_274),nrow(cyd_284),nrow(cyd_313),
     nrow(cyd_316),nrow(cyd_326),nrow(cyd_333),nrow(cyd_352),
     nrow(cyd_356),nrow(cyd_365)))
  #=== 2. assign genes to the snps ----

# read in gene annotation
annotation <- read.table("gene_info.txt",header = TRUE)

# create ranges of where genes are
GRgenes<-makeGRangesFromDataFrame(annotation, seqnames.field=c("scaffold"))


# read in ASE count tables 
# ***do melp and cyd reads separately*** #

#choose to read in melp or cyd reads
  #antennae <- read.table("ASE data prep/1.ASE_counts_antennae_melp.txt")
  #antennae <- read.table("ASE data prep/1.ASE_counts_antennae_cyd.txt")

# add "Start" and "end" position (same as start because its a snp), needed by GRanges
colnames(antennae)[2] <- "start"
antennae$end = antennae$start
GRsnps<-makeGRangesFromDataFrame(antennae, seqnames.field=c("scaffold"))

#find overlaps between snps and genes
olaps <- findOverlaps(GRgenes, GRsnps)
antennae_snps<-cbind(annotation[queryHits(olaps),], antennae[subjectHits(olaps),])
#rename column
colnames(antennae_snps)[6] <-"position"

#save tables
#write.table(antennae_snps,"ASE data prep/2.antennae_snps_melp.txt")
#write.table(antennae_snps,"ASE data prep/2.antennae_snps_cyd.txt")

  #=== 3. subset into groups ====

# *** do melp and cyd separately ***

# read in melp or cyd 
antennae_snps <- read.table("3. ASE data prep/antennae/2.antennae_snps_melp.txt")
antennae_snps<- read.table("3. ASE data prep/antennae/2.antennae_snps_cyd.txt")

#F1_MPCP
antennae_F1_MPCP <- antennae_snps[,c(1:4,6,which(grepl("274",colnames(antennae_snps))),which(grepl("333",colnames(antennae_snps))),
                                     which(grepl("352",colnames(antennae_snps))),which(grepl("356",colnames(antennae_snps))),
                                     which(grepl("365",colnames(antennae_snps))))]

#F1_MPCP_female_old     
antennae_F1_MPCP_female_old <- antennae_snps[,c(1:4,6,which(grepl("141",colnames(antennae_snps))),which(grepl("214",colnames(antennae_snps)))
                                                ,which(grepl("326",colnames(antennae_snps))))]

#F1_MPCP_female_young
antennae_F1_MPCP_female_young <- antennae_snps[,c(1:4,6,which(grepl("137",colnames(antennae_snps))),which(grepl("261",colnames(antennae_snps)))
                                                  ,which(grepl("313",colnames(antennae_snps))),which(grepl("59",colnames(antennae_snps))))]

#F1_MPCP_young
antennae_F1_MPCP_young <- antennae_snps[,c(1:4,6,which(grepl("172",colnames(antennae_snps))),which(grepl("284",colnames(antennae_snps))),
                                           which(grepl("316",colnames(antennae_snps))))]


#write.table(antennae_F1_MPCP, file = "ASE data prep/3.antennae_snps_F1_MPCP_melp.txt")
#write.table(antennae_F1_MPCP, file="ASE data prep/3.antennae_snps_F1_MPCP_cyd.txt")

#write.table(antennae_F1_MPCP_female_old, file = "ASE data prep/3.antennae_snps_F1_MPCP_female_old_melp.txt")
#write.table(antennae_F1_MPCP_female_old, file="ASE data prep/3.antennae_snps_F1_MPCP_female_old_cyd.txt")

#write.table(antennae_F1_MPCP_young, file = "ASE data prep/3.antennae_snps_F1_MPCP_young_melp.txt")
#write.table(antennae_F1_MPCP_young, file="ASE data prep/3.antennae_snps_F1_MPCP_young_cyd.txt")

#write.table(antennae_F1_MPCP_female_young, file = "ASE data prep/3.antennae_snps_F1_MPCP_female_young_melp.txt")
#write.table(antennae_F1_MPCP_female_young, file="ASE data prep/3.antennae_snps_F1_MPCP_female_young_cyd.txt")


  #=== 4. sum all SNPs of the same gene together ====

# F1_MPCP

#choose melp or cyd snp counts
  #antennae_F1_MPCP <- read.table("ASE data prep/3.antennae_snps_F1_MPCP_melp.txt")
  antennae_F1_MPCP <- read.table("ASE data prep/3.antennae_snps_F1_MPCP_cyd.txt")

# replace NA with 0
antennae_F1_MPCP[is.na(antennae_F1_MPCP)] <- 0

# sum up the rows with the same gene_id
ASE_genes_antennae_F1_MPCP <- rowsum(antennae_F1_MPCP[6:length(colnames(antennae_F1_MPCP))], antennae_F1_MPCP$gene_id)

# change column names
id <- c("F1_274_ref","F1_274_alt","F1_333_ref","F1_333_alt","F1_352_ref","F1_352_alt","F1_356_ref","F1_356_alt",
        "F1_365_ref","F1_365_alt")

colnames(ASE_genes_antennae_F1_MPCP)<-id


#write.table(ASE_genes_antennae_F1_MPCP,"ASE data prep/4.antennae_genes_F1_MPCP_melp.txt")
#write.table(ASE_genes_antennae_F1_MPCP,"ASE data prep/4.antennae_genes_F1_MPCP_cyd.txt")


# F1_MPCP_female_old

#read in melp/cyd snp counts
#antennae_F1_MPCP_female_old <- read.table("ASE data prep/3.antennae_snps_F1_MPCP_female_old_melp.txt")
#antennae_F1_MPCP_female_old <- read.table("ASE data prep/3.antennae_snps_F1_MPCP_female_old_cyd.txt")

# replace NA with 0
antennae_F1_MPCP_female_old[is.na(antennae_F1_MPCP_female_old)] <- 0

# sum up the rows with the same gene_id
ASE_genes_antennae_F1_MPCP_female_old <- rowsum(antennae_F1_MPCP_female_old[6:length(colnames(antennae_F1_MPCP_female_old))], antennae_F1_MPCP_female_old$gene_id)

# change column names
id <- c("F1_141_ref","F1_141_alt","F1_214_ref","F1_214_alt","F1_326_ref","F1_326_alt")

colnames(ASE_genes_antennae_F1_MPCP_female_old)<-id


#write.table(ASE_genes_antennae_F1_MPCP_female_old,"ASE data prep/4.antennae_genes_F1_MPCP_female_old_melp.txt")
#write.table(ASE_genes_antennae_F1_MPCP_female_old,"ASE data prep/4.antennae_genes_F1_MPCP_female_old_cyd.txt")


# F1_MPCP_female_young

#choose melp or cyd snps
antennae_F1_MPCP_female_young <- read.table("ASE data prep/3.antennae_snps_F1_MPCP_female_young_melp.txt")
antennae_F1_MPCP_female_young <- read.table("ASE data prep/3.antennae_snps_F1_MPCP_female_young_cyd.txt")


antennae_F1_MPCP_female_young[is.na(antennae_F1_MPCP_female_young)] <- 0

ASE_genes_antennae_F1_MPCP_female_young <- rowsum(antennae_F1_MPCP_female_young[6:length(colnames(antennae_F1_MPCP_female_young))], antennae_F1_MPCP_female_young$gene_id)

id <- c("F1_137_ref","F1_137_alt","F1_261_ref","F1_261_alt","F1_313_ref","F1_313_alt","F1_59_ref","F1_59_alt")

colnames(ASE_genes_antennae_F1_MPCP_female_young)<-id

#write.table(ASE_genes_antennae_F1_MPCP_female_young, "ASE data prep/4.antennae_genes_F1_MPCP_female_young_melp.txt")
#write.table(ASE_genes_antennae_F1_MPCP_female_young, "ASE data prep/4.antennae_genes_F1_MPCP_female_young_cyd.txt")

# F1_MPCP_young

antennae_F1_MPCP_young <- read.table("ASE data prep/3.antennae_snps_F1_MPCP_young_melp.txt")
antennae_F1_MPCP_young <- read.table("ASE data prep/3.antennae_snps_F1_MPCP_young_cyd.txt")

antennae_F1_MPCP_young[is.na(antennae_F1_MPCP_young)] <- 0

ASE_genes_antennae_F1_MPCP_young <- rowsum(antennae_F1_MPCP_young[6:length(colnames(antennae_F1_MPCP_young))], antennae_F1_MPCP_young$gene_id)

id <- c("F1_172_ref","F1_172_alt","F1_284_ref","F1_284_alt","F1_316_ref","F1_316_alt")
colnames(ASE_genes_antennae_F1_MPCP_young)<-id

#write.table(ASE_genes_antennae_F1_MPCP_young, "ASE data prep/4.antennae_genes_F1_MPCP_young_melp.txt")
#write.table(ASE_genes_antennae_F1_MPCP_young, "ASE data prep/4.antennae_genes_F1_MPCP_young_cyd.txt")


  #=== 5/6. combine cyd and melp tables and creat coldata table for DESeq ====
    # F1_MPCP_female_old ----

# read in files
F1_MPCP_female_old_melp <- read.table("3. ASE data prep/antennae/4.antennae_genes_F1_MPCP_female_old_melp.txt")
F1_MPCP_female_old_cyd <- read.table("3. ASE data prep/antennae/4.antennae_genes_F1_MPCP_female_old_cyd.txt")

# assign melp/cyd to ref/alt alleles
colnames(F1_MPCP_female_old_melp) <- gsub("ref", "cyd", colnames(F1_MPCP_female_old_melp))
colnames(F1_MPCP_female_old_melp) <- gsub("alt", "melp", colnames(F1_MPCP_female_old_melp))

colnames(F1_MPCP_female_old_cyd) <- gsub("ref", "melp", colnames(F1_MPCP_female_old_cyd))
colnames(F1_MPCP_female_old_cyd) <- gsub("alt", "cyd", colnames(F1_MPCP_female_old_cyd))

#create column with gene id
F1_MPCP_female_old_melp$gene <- row.names(F1_MPCP_female_old_melp)
F1_MPCP_female_old_cyd$gene <- row.names(F1_MPCP_female_old_cyd)

# combine cyd and melp snps into one table
# x = melp, y = cyd 
F1_MPCP_female_old <- full_join(F1_MPCP_female_old_melp,F1_MPCP_female_old_cyd,by="gene")
# make the gene names the row names, remove "gene" column
row.names(F1_MPCP_female_old) <- F1_MPCP_female_old$gene
F1_MPCP_female_old <- F1_MPCP_female_old[,-which(grepl("gene",colnames(F1_MPCP_female_old)))]
# sum melp/cyd counts for each individual (sum columns with the same name)
colnames(F1_MPCP_female_old)<- gsub(".x$","",colnames(F1_MPCP_female_old))
colnames(F1_MPCP_female_old)<- gsub(".y$","",colnames(F1_MPCP_female_old))
F1_MPCP_female_old[is.na(F1_MPCP_female_old)] <- 0
F1_MPCP_female_old <- t(rowsum(t(F1_MPCP_female_old), group = colnames(F1_MPCP_female_old)))

write.table(F1_MPCP_female_old,"ASE data prep/5.antennae_genes_melpcyd_F1_MPCP_female_old.txt")

#create ColData table
#assign sample id to each column
individual <- factor(substr(colnames(F1_MPCP_female_old),4,6))
#assign melp or cyd alleles to each column
allele <- factor(substr(colnames(F1_MPCP_female_old),8,11))

brood <- factor(c(rep("007",4),rep("011",2)))

coldata <- data.frame(row.names=colnames(F1_MPCP_female_old),individual,allele,brood)

write.table(coldata,"ASE data prep/6.antennae_coldata_F1_MPCP_female_old.txt")

    # male old ====
# read in files
F1_MPCP_melp <- read.table("ASE data prep/4.antennae_genes_F1_MPCP_melp.txt")
F1_MPCP_cyd <- read.table("ASE data prep/4.antennae_genes_F1_MPCP_cyd.txt")

# assign melp/cyd to ref/alt alleles
colnames(F1_MPCP_melp) <- gsub("ref", "cyd", colnames(F1_MPCP_melp))
colnames(F1_MPCP_melp) <- gsub("alt", "melp", colnames(F1_MPCP_melp))

colnames(F1_MPCP_cyd) <- gsub("ref", "melp", colnames(F1_MPCP_cyd))
colnames(F1_MPCP_cyd) <- gsub("alt", "cyd", colnames(F1_MPCP_cyd))

#create column with gene id
F1_MPCP_melp$gene <- row.names(F1_MPCP_melp)
F1_MPCP_cyd$gene <- row.names(F1_MPCP_cyd)

# combine cyd and melp snps into one table
# x = melp, y = cyd 
F1_MPCP <- full_join(F1_MPCP_melp,F1_MPCP_cyd,by="gene")
# make the gene names the row names, remove "gene" column
row.names(F1_MPCP) <- F1_MPCP$gene
F1_MPCP <- F1_MPCP[,-which(grepl("gene",colnames(F1_MPCP)))]
# sum melp/cyd counts for each individual (sum columns with the same name)
colnames(F1_MPCP)<- gsub(".x$","",colnames(F1_MPCP))
colnames(F1_MPCP)<- gsub(".y$","",colnames(F1_MPCP))
F1_MPCP[is.na(F1_MPCP)] <- 0
F1_MPCP <- t(rowsum(t(F1_MPCP), group = colnames(F1_MPCP)))

write.table(F1_MPCP,"ASE data prep/5.antennae_genes_melpcyd_F1_MPCP.txt")

#create ColData table
#assign sample id to each column
individual <- factor(substr(colnames(F1_MPCP),4,6))
#assign melp or cyd alleles to each column
allele <- factor(substr(colnames(F1_MPCP),8,11))

brood <- factor(c(rep("011",10)))

coldata <- data.frame(row.names=colnames(F1_MPCP),individual,allele,brood)

write.table(coldata,"ASE data prep/6.antennae_coldata_F1_MPCP.txt")


    # female young ====
# read in files
F1_MPCP_female_young_melp <- read.table("ASE data prep/4.antennae_genes_F1_MPCP_female_young_melp.txt")
F1_MPCP_female_young_cyd <- read.table("ASE data prep/4.antennae_genes_F1_MPCP_female_young_cyd.txt")

# assign melp/cyd to ref/alt alleles
colnames(F1_MPCP_female_young_melp) <- gsub("ref", "cyd", colnames(F1_MPCP_female_young_melp))
colnames(F1_MPCP_female_young_melp) <- gsub("alt", "melp", colnames(F1_MPCP_female_young_melp))

colnames(F1_MPCP_female_young_cyd) <- gsub("ref", "melp", colnames(F1_MPCP_female_young_cyd))
colnames(F1_MPCP_female_young_cyd) <- gsub("alt", "cyd", colnames(F1_MPCP_female_young_cyd))

#create column with gene id
F1_MPCP_female_young_melp$gene <- row.names(F1_MPCP_female_young_melp)
F1_MPCP_female_young_cyd$gene <- row.names(F1_MPCP_female_young_cyd)

# combine cyd and melp snps into one table
# x = melp, y = cyd 
F1_MPCP_female_young <- full_join(F1_MPCP_female_young_melp,F1_MPCP_female_young_cyd,by="gene")
# make the gene names the row names, remove "gene" column
row.names(F1_MPCP_female_young) <- F1_MPCP_female_young$gene
F1_MPCP_female_young <- F1_MPCP_female_young[,-which(grepl("gene",colnames(F1_MPCP_female_young)))]
# sum melp/cyd counts for each individual (sum columns with the same name)
colnames(F1_MPCP_female_young)<- gsub(".x$","",colnames(F1_MPCP_female_young))
colnames(F1_MPCP_female_young)<- gsub(".y$","",colnames(F1_MPCP_female_young))
F1_MPCP_female_young[is.na(F1_MPCP_female_young)] <- 0
F1_MPCP_female_young <- t(rowsum(t(F1_MPCP_female_young), group = colnames(F1_MPCP_female_young)))

write.table(F1_MPCP_female_young,"ASE data prep/5.antennae_genes_melpcyd_F1_MPCP_female_young.txt")

#create ColData table
#assign sample id to each column
individual <- factor(substr(colnames(F1_MPCP_female_young),4,6))
#assign melp or cyd alleles to each column
allele <- factor(substr(colnames(F1_MPCP_female_young),8,11))

brood <- factor(c(rep("007",2),rep("011",4),rep("007",2)))

coldata <- data.frame(row.names=colnames(F1_MPCP_female_young),individual,allele)

coldata$allele[7] <- "cyd"
coldata$allele[8] <- "melp"

write.table(coldata,"ASE data prep/6.antennae_coldata_F1_MPCP_female_young.txt")
    # male young ====
# read in files
F1_MPCP_young_melp <- read.table("ASE data prep/4.antennae_genes_F1_MPCP_young_melp.txt")
F1_MPCP_young_cyd <- read.table("ASE data prep/4.antennae_genes_F1_MPCP_young_cyd.txt")

# assign melp/cyd to ref/alt alleles
colnames(F1_MPCP_young_melp) <- gsub("ref", "cyd", colnames(F1_MPCP_young_melp))
colnames(F1_MPCP_young_melp) <- gsub("alt", "melp", colnames(F1_MPCP_young_melp))

colnames(F1_MPCP_young_cyd) <- gsub("ref", "melp", colnames(F1_MPCP_young_cyd))
colnames(F1_MPCP_young_cyd) <- gsub("alt", "cyd", colnames(F1_MPCP_young_cyd))

#create column with gene id
F1_MPCP_young_melp$gene <- row.names(F1_MPCP_young_melp)
F1_MPCP_young_cyd$gene <- row.names(F1_MPCP_young_cyd)

# combine cyd and melp snps into one table
# x = melp, y = cyd 
F1_MPCP_young <- full_join(F1_MPCP_young_melp,F1_MPCP_young_cyd,by="gene")
# make the gene names the row names, remove "gene" column
row.names(F1_MPCP_young) <- F1_MPCP_young$gene
F1_MPCP_young <- F1_MPCP_young[,-which(grepl("gene",colnames(F1_MPCP_young)))]
# sum melp/cyd counts for each individual (sum columns with the same name)
colnames(F1_MPCP_young)<- gsub(".x$","",colnames(F1_MPCP_young))
colnames(F1_MPCP_young)<- gsub(".y$","",colnames(F1_MPCP_young))
F1_MPCP_young[is.na(F1_MPCP_young)] <- 0
F1_MPCP_young <- t(rowsum(t(F1_MPCP_young), group = colnames(F1_MPCP_young)))

write.table(F1_MPCP_young,"ASE data prep/5.antennae_genes_melpcyd_F1_MPCP_young.txt")

#create ColData table
#assign sample id to each column
individual <- factor(substr(colnames(F1_MPCP_young),4,6))
#assign melp or cyd alleles to each column
allele <- factor(substr(colnames(F1_MPCP_young),8,11))

brood <- factor(c(rep("007",2),rep("011",4)))

coldata <- data.frame(row.names=colnames(F1_MPCP_young),individual,allele,brood)

write.table(coldata,"ASE data prep/6.antennae_coldata_F1_MPCP_young.txt")
