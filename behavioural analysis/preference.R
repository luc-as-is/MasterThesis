library(lme4)
library(emmeans)
library(ggplot2)
library(AER)

setwd("~/lukeball/HostPref_R/1. behaviour")

data <- read.csv("hp_data_clean.csv", header = T, stringsAsFactors = F)

data$newID<-1:nrow(data)
data$hp_acceptance <- data$meni/data$total_eggs
data$type[!data$type%in%c("BC","BM","CP","MP")]<-"F1"

#remove sterile individuals
data$sterile <- as.numeric(is.na(data$hp_acceptance))
plotdata <- data[data$sterile==0,]

#### by female type ====
#reorder x axis
plotdata$type <- factor(plotdata$type, levels = c ("CP", "BC", "F1", "BM", "MP"))
#ggplot
ggplot(plotdata, aes(type, hp_acceptance, fill=type)) + 
  geom_jitter(size = 0.1 * plotdata$total_eggs) 
cbind(plotdata$meni,plotdata$viti)

###
#stats

library(MASS)

glm <- glm.nb(meni/viti ~ type, data = plotdata)


#glmer and anova
mod1 <- glmer(cbind(meni,viti) ~ type + (1|newID), family = "binomial", data = plotdata)
mod2 <- glmer(cbind(meni,viti) ~ (1|newID), family = "binomial", data = plotdata)
mod3 <- glmer(cbind(meni,viti) ~ type + (1|data_set) + (1|newID), family = "binomial", data = plotdata)
mod4 <- glmer(cbind(meni,viti) ~ type + (1|brood), family = "binomial", data = plotdata)
mod5 <- glmer(cbind(meni,viti) ~ type + (1|data_set) + (1|newID) +(1|brood), family = "binomial", data = plotdata)

anova(mod1, mod2, mod3,mod4,mod5)


## mod 3 is best
emmeans(mod3, pairwise ~type)

#check for overdispersion
summary(mod3)

sqrt(sum(c(resid(mod3),mod3@u)^2)/length(resid(mod3)))
# = ~0.881 

#means and confidence intervals
stats <- summary(emmeans(mod3, specs = "type"))
stats$emmean <- plogis(stats$emmean)
stats$asymp.LCL <- plogis(stats$asymp.LCL)
stats$asymp.UCL <- plogis(stats$asymp.UCL)
  # fix MP cofidence intervals
stats[5,5] <- 1.0

###

# plot ====
## by female type ====
par(mar=c(3.2,3.5,0.5,4))

plot(1,1, type = "n", xaxt="n",yaxt="n",xlab="",ylab="",
     xlim=c(0.5,5.5),xaxs="i",ylim=0:1)
abline(v=c(1.5,2.5,3.5,4.5))

points(jitter(rep(1,sum(plotdata$type=="CP")),7),
       plotdata$hp_acceptance[plotdata$type=="CP"],
       pch=21,bg=adjustcolor("gray",0.5),
       cex=0.2*sqrt(plotdata$total_eggs[plotdata$type=="CP"]))

points(jitter(rep(2,sum(plotdata$type=="BC")),7),
       plotdata$hp_acceptance[plotdata$type=="BC"],
       pch=21,bg=adjustcolor("gray",0.5),
       cex=0.2*sqrt(plotdata$total_eggs[plotdata$type=="BC"]))

points(jitter(rep(3,sum(plotdata$type=="F1")),6),
       plotdata$hp_acceptance[plotdata$type=="F1"],
       pch=21,bg=adjustcolor("gray",0.5),
       cex=0.2*sqrt(plotdata$total_eggs[plotdata$type=="F1"]))

points(jitter(rep(4,sum(plotdata$type=="BM")),5),
       plotdata$hp_acceptance[plotdata$type=="BM"],
       pch=21,bg=adjustcolor("gray",0.5),
       cex=0.2*sqrt(plotdata$total_eggs[plotdata$type=="BM"]))

points(jitter(rep(5,sum(plotdata$type=="MP")),4),
       plotdata$hp_acceptance[plotdata$type=="MP"],
       pch=21,bg=adjustcolor("gray",0.5),
       cex=0.2*sqrt(plotdata$total_eggs[plotdata$type=="MP"]))


axis(1,at=1:5,c(expression(italic("H. cydno")),"BCC","F1","BCM",expression(italic("H. melpomene"))),
     lwd=0,line=-0.5, cex.axis=0.8)

axis(2,las=2)

mtext("Female Type",1,line=2)
mtext(substitute(paste("Acceptance for ",italic(' P. menispermifolia'))),2,line=2.5)

segments((1:5)-0.2, stats$asymp.LCL, (1:5)-0.2, stats$asymp.UCL, col = "#f52c2c",lwd = 2)
segments((1:5)-0.3, stats$emmean, (1:5)-0.1, stats$emmean, col = "#f52c2c",lwd = 3)
         
#check range of data points
range(plotdata$total_eggs[plotdata$type%in%c("CP", "BC", "F1", "BM", "MP")])

#Create legend: Take some values that seem more or less adequate
dots_for<-c(10,30,50,70)

par(xpd=T)

points(rep(5.85, 4),
       c(0.4,0.51,0.62,0.73),
       cex=0.2*sqrt(dots_for),pch=21, 
       bg=adjustcolor("gray",0.5))
text(rep(6.25, 4),
     c(0.4,0.51,0.62,0.73),
     as.character(dots_for))
par(xpd=F)

dev.copy(png, "behaviour.png",res=300,width=1800, height=1100)
dev.off()

## by brood====
# plot each female type separately, 
# grouping individuals by brood to look for brood effects

# group by female type
F1 <- plotdata[which(plotdata$type=="F1"),]
  F1 <- F1[-which(is.na(F1$brood)),]
BCC <- plotdata[which(plotdata$type=="BC"),]
BCM <- plotdata[which(plotdata$type=="BM"),]
  BCM <- BCM[-which(is.na(BCM$brood)),]

  ggplot(BCM, aes(brood, hp_acceptance)) + 
    geom_point() 
  cbind(plotdata$meni,plotdata$viti)
  
## stats


# GLMM
  #F1 ====
lm1 <- glmer(cbind(meni,viti) ~ (1|brood), family = "binomial", data = F1)
lm2 <- glmer(cbind(meni,viti) ~ (1|newID), family = "binomial", data = F1)
lm3 <- glmer(cbind(meni,viti) ~ brood + (1|newID), family = "binomial", data = F1)
lm4 <- glmer(cbind(meni,viti) ~ (1|brood) + (1|newID), family = "binomial", data = F1)

anova(lm1,lm2,lm3,lm4)


emmeans(lm3, pairwise ~brood)

#check for overdispersion again
summary(lm3)

sqrt(sum(c(resid(lm3),lm3@u)^2)/length(resid(lm3)))

# = 0.7988674

#means and confidence intervals
stats.F1 <- summary(emmeans(lm3, specs = "brood"))
stats.F1$emmean <- plogis(stats.F1$emmean)
stats.F1$asymp.LCL <- plogis(stats.F1$asymp.LCL)
stats.F1$asymp.UCL <- plogis(stats.F1$asymp.UCL)
stats.F1

  stats.F1[2,5]<- 1
    ## plot ====
par(mar=c(3.2,3.5,0.5,4.5))

plot(1,1, type = "n", xaxt="n",yaxt="n",xlab="",ylab="",
     xlim=c(0.5,3.5),xaxs="i",ylim=0:1)
abline(v=c(1.5,2.5))

points(rep(1,sum(F1$brood=="c1")),
       F1$hp_acceptance[F1$brood=="c1"],
       pch=21,bg=adjustcolor("gray",0.5),
       cex=0.5*sqrt(F1$total_eggs[F1$brood=="c1"]))

points(jitter(rep(2,sum(F1$brood=="c6")),6),
       F1$hp_acceptance[F1$brood=="c6"],
       pch=21,bg=adjustcolor("gray",0.5),
       cex=0.5*sqrt(F1$total_eggs[F1$brood=="c6"]))

points(rep(3,sum(F1$brood=="N32")),
       F1$hp_acceptance[F1$brood=="N32"],
       pch=21,bg=adjustcolor("gray",0.5),
       cex=0.5*sqrt(F1$total_eggs[F1$brood=="N32"]))


axis(1,at=1:3,c("C1","C6","N32"),
     lwd=0,line=-0.5)

axis(2,las=2)

mtext("F1 Brood",1,line=2)
mtext(substitute(paste("Acceptance for ",italic(' P. menispermifolia'))),2,line=2.5)

segments((1:3)-0.2, stats.F1$asymp.LCL, (1:3)-0.2, stats.F1$asymp.UCL, col = "#f52c2c",lwd = 2)
segments((1:3)-0.3, stats.F1$emmean, (1:3)-0.1, stats.F1$emmean, col = "#f52c2c",lwd=3)
#check range of data points
range(F1$total_eggs[F1$brood%in%c("c1", "c6", "N32")])

#Take some values that seem more or less adequate
dots_for<-c(10,25,40,55)

par(xpd=T)

points(rep(3.75, 4),
       c(0.4,0.51,0.62,0.73),
       cex=0.5*sqrt(dots_for),pch=21, 
       bg=adjustcolor("gray",0.5))
text(rep(4, 4),
     c(0.4,0.51,0.62,0.73),
     as.character(dots_for))
par(xpd=F)

dev.copy(png, "behaviour broods - F1.png",res=300,width=1800, height=1100)
dev.off()
  # BCC ====
lm1 <- glmer(cbind(meni,viti) ~ (1|brood), family = "binomial", data = BCC)
lm2 <- glmer(cbind(meni,viti) ~ (1|newID), family = "binomial", data = BCC)
lm3 <- glmer(cbind(meni,viti) ~ brood + (1|newID), family = "binomial", data = BCC)
lm4 <- glmer(cbind(meni,viti) ~ (1|brood) + (1|newID), family = "binomial", data = BCC)

anova(lm1,lm2,lm3)

emmeans(lm3, pairwise ~brood)

#check for overdispersion again
summary(lm3)

sqrt(sum(c(resid(lm3),lm3@u)^2)/length(resid(lm3)))

# = 1.056174

#means and confidence intervals
stats.BCC <- summary(emmeans(lm3, specs = "brood"))
stats.BCC$emmean <- plogis(stats.BCC$emmean)
stats.BCC$asymp.LCL <- plogis(stats.BCC$asymp.LCL)
stats.BCC$asymp.UCL <- plogis(stats.BCC$asymp.UCL)
stats.BCC

stats.BCC[which(stats.BCC$asymp.LCL==0),5]<- 1
    ## plot ====
par(mar=c(3.2,3.5,0.5,4))

plot(1,1, type = "n", xaxt="n",yaxt="n",xlab="",ylab="",
     xlim=c(0.5,19.5),xaxs="i",ylim=0:1)
abline(v=c((1:19)+0.5))

points(rep(1,sum(BCC$brood=="C10")),
       BCC$hp_acceptance[BCC$brood=="C10"],
       pch=21,bg=adjustcolor("gray",0.5),
       cex=0.3*sqrt(BCC$total_eggs[BCC$brood=="C10"]))
points(rep(2,sum(BCC$brood=="C11")),
       BCC$hp_acceptance[BCC$brood=="C11"],
       pch=21,bg=adjustcolor("gray",0.5),
       cex=0.3*sqrt(BCC$total_eggs[BCC$brood=="C11"]))
points(rep(3,sum(BCC$brood=="C14")),
       BCC$hp_acceptance[BCC$brood=="C14"],
       pch=21,bg=adjustcolor("gray",0.5),
       cex=0.3*sqrt(BCC$total_eggs[BCC$brood=="C14"]))
points(rep(4,sum(BCC$brood=="C18")),
       BCC$hp_acceptance[BCC$brood=="C18"],
       pch=21,bg=adjustcolor("gray",0.5),
       cex=0.3*sqrt(BCC$total_eggs[BCC$brood=="C18"]))
points(rep(5,sum(BCC$brood=="C21")),
       BCC$hp_acceptance[BCC$brood=="C21"],
       pch=21,bg=adjustcolor("gray",0.5),
       cex=0.3*sqrt(BCC$total_eggs[BCC$brood=="C21"]))
points(rep(6,sum(BCC$brood=="C24")),
       BCC$hp_acceptance[BCC$brood=="C24"],
       pch=21,bg=adjustcolor("gray",0.5),
       cex=0.3*sqrt(BCC$total_eggs[BCC$brood=="C24"]))
points(rep(7,sum(BCC$brood=="C8")),
       BCC$hp_acceptance[BCC$brood=="C8"],
       pch=21,bg=adjustcolor("gray",0.5),
       cex=0.3*sqrt(BCC$total_eggs[BCC$brood=="C8"]))
points(rep(8,sum(BCC$brood=="H5")),
       BCC$hp_acceptance[BCC$brood=="H5"],
       pch=21,bg=adjustcolor("gray",0.5),
       cex=0.3*sqrt(BCC$total_eggs[BCC$brood=="H5"]))
points(rep(9,sum(BCC$brood=="N10")),
       BCC$hp_acceptance[BCC$brood=="N10"],
       pch=21,bg=adjustcolor("gray",0.5),
       cex=0.3*sqrt(BCC$total_eggs[BCC$brood=="N10"]))
points(rep(10,sum(BCC$brood=="N12")),
       BCC$hp_acceptance[BCC$brood=="N12"],
       pch=21,bg=adjustcolor("gray",0.5),
       cex=0.3*sqrt(BCC$total_eggs[BCC$brood=="N12"]))
points(rep(11,sum(BCC$brood=="N14")),
       BCC$hp_acceptance[BCC$brood=="N14"],
       pch=21,bg=adjustcolor("gray",0.5),
       cex=0.3*sqrt(BCC$total_eggs[BCC$brood=="N14"]))
points(rep(12,sum(BCC$brood=="N15")),
       BCC$hp_acceptance[BCC$brood=="N15"],
       pch=21,bg=adjustcolor("gray",0.5),
       cex=0.3*sqrt(BCC$total_eggs[BCC$brood=="N15"]))
points(rep(13,sum(BCC$brood=="N17")),
       BCC$hp_acceptance[BCC$brood=="N17"],
       pch=21,bg=adjustcolor("gray",0.5),
       cex=0.3*sqrt(BCC$total_eggs[BCC$brood=="N17"]))
points(rep(14,sum(BCC$brood=="N2")),
       BCC$hp_acceptance[BCC$brood=="N2"],
       pch=21,bg=adjustcolor("gray",0.5),
       cex=0.3*sqrt(BCC$total_eggs[BCC$brood=="N2"]))
points(rep(15,sum(BCC$brood=="N35")),
       BCC$hp_acceptance[BCC$brood=="N35"],
       pch=21,bg=adjustcolor("gray",0.5),
       cex=0.3*sqrt(BCC$total_eggs[BCC$brood=="N35"]))
points(rep(16,sum(BCC$brood=="N49")),
       BCC$hp_acceptance[BCC$brood=="N49"],
       pch=21,bg=adjustcolor("gray",0.5),
       cex=0.3*sqrt(BCC$total_eggs[BCC$brood=="N49"]))
points(rep(17,sum(BCC$brood=="N62")),
       BCC$hp_acceptance[BCC$brood=="N62"],
       pch=21,bg=adjustcolor("gray",0.5),
       cex=0.3*sqrt(BCC$total_eggs[BCC$brood=="N62"]))
points(rep(18,sum(BCC$brood=="N64")),
       BCC$hp_acceptance[BCC$brood=="N64"],
       pch=21,bg=adjustcolor("gray",0.5),
       cex=0.3*sqrt(BCC$total_eggs[BCC$brood=="N64"]))
points(rep(19,sum(BCC$brood=="N66")),
       BCC$hp_acceptance[BCC$brood=="N66"],
       pch=21,bg=adjustcolor("gray",0.5),
       cex=0.3*sqrt(BCC$total_eggs[BCC$brood=="N66"]))


axis(1,at=1:19,unique(BCC$brood),
     lwd=0,line=-0.5,las=2)

axis(2,las=2)

mtext("BCC Brood",1,line=2)
mtext(substitute(paste("Acceptance for ",italic(' P. menispermifolia'))),2,line=2.5)

remove <- c(2,5,6,16,17,19)
keep <- c(1:19)
keep<-keep[-remove]

segments(keep-0.2, stats.BCC$asymp.LCL[keep], keep-0.2, stats.BCC$asymp.UCL[keep], col = "#f52c2c",lwd = 2)
segments((1:19)-0.3, stats.BCC$emmean, (1:19)-0.1, stats.BCC$emmean, col = "#f52c2c",lwd=3)

#check range of data points
range(BCC$total_eggs[BCC$brood%in%unique(BCC$brood)])

#Take some values that seem more or less adequate
dots_for<-c(10,30,50,70)

par(xpd=T)

points(rep(20.75, 4),
       c(0.4,0.51,0.62,0.73),
       cex=0.3*sqrt(dots_for),pch=21, 
       bg=adjustcolor("gray",0.5))
text(rep(22, 4),
     c(0.4,0.51,0.62,0.73),
     as.character(dots_for))
par(xpd=F)
dev.copy(png, "behaviour broods - BCC.png",res=300,width=1800, height=1100)
dev.off()
  
#BCM ====
lm1 <- glmer(cbind(meni,viti) ~ (1|brood), family = "binomial", data = BCM)
lm2 <- glmer(cbind(meni,viti) ~ (1|newID), family = "binomial", data = BCM)
lm3 <- glmer(cbind(meni,viti) ~ (1|data_set), family = "binomial", data = BCM)

lm3 <- glmer(cbind(meni,viti) ~ brood + (1|newID), family = "binomial", data = BCM)
lm4 <- glmer(cbind(meni,viti) ~ brood + (1|newID) + (1|data_set), family = "binomial", data = BCM)

anova(lm1,lm2,lm3,lm4)


emmeans(lm3, pairwise ~brood)

#check for overdispersion again
summary(lm3)

sqrt(sum(c(resid(lm3),lm3@u)^2)/length(resid(lm3)))

# = 0.7988674

#means and confidence intervals
stats.BCM <- summary(emmeans(lm3, specs = "brood"))
stats.BCM$emmean <- plogis(stats.BCM$emmean)
stats.BCM$asymp.LCL <- plogis(stats.BCM$asymp.LCL)
stats.BCM$asymp.UCL <- plogis(stats.BCM$asymp.UCL)
stats.BCM

N41 = c(1,1,1,1,0.7837838)
# calculate mean and CIs for brood N41
stats.BCM [3,2] <- mean(N41)

error <- qnorm(0.975)*sd(N41)/sqrt(5)
stats.BCM [3,5] <- mean(N41)-error
stats.BCM [3,6] <- mean(N41)+error


stats.BCM[(1:2),5]<- 1
    ## plot ====
par(mar=c(3.2,3.5,0.5,4.5))

plot(1,1, type = "n", xaxt="n",yaxt="n",xlab="",ylab="",
     xlim=c(0.5,3.5),xaxs="i",ylim=0:1)
abline(v=c(1.5,2.5))

points(jitter(rep(1,sum(BCM$brood=="c4")),10),
       BCM$hp_acceptance[BCM$brood=="c4"],
       pch=21,bg=adjustcolor("gray",0.5),
       cex=0.7*sqrt(BCM$total_eggs[BCM$brood=="c4"]))

points(jitter(rep(2,sum(BCM$brood=="c5")),10),
       BCM$hp_acceptance[BCM$brood=="c5"],
       pch=21,bg=adjustcolor("gray",0.5),
       cex=0.7*sqrt(BCM$total_eggs[BCM$brood=="c5"]))

points(jitter(rep(3,sum(BCM$brood=="N41")),7),
       BCM$hp_acceptance[BCM$brood=="N41"],
       pch=21,bg=adjustcolor("gray",0.5),
       cex=0.7*sqrt(BCM$total_eggs[BCM$brood=="N41"]))


axis(1,at=1:3,c("C4","C5","N41"),
     lwd=0,line=-0.5)

axis(2,las=2)

mtext("BCM Brood",1,line=2)
mtext(substitute(paste("Acceptance for ",italic(' P. menispermifolia'))),2,line=2.5)

segments((1:3)-0.2, stats.BCM$asymp.LCL, (1:3)-0.2, stats.BCM$asymp.UCL, col = "#f52c2c",lwd = 2)
segments((1:3)-0.3, stats.BCM$emmean, (1:3)-0.1, stats.BCM$emmean, col = "#f52c2c",lwd=3)

#check range of data points
range(BCM$total_eggs[BCM$brood%in%unique(BCM$brood)])

#Take some values that seem more or less adequate
dots_for<-c(10,20,30,40)

par(xpd=T)

points(rep(3.75, 4),
       c(0.38,0.5,0.62,0.74),
       cex=0.7*sqrt(dots_for),pch=21, 
       bg=adjustcolor("gray",0.5))
text(rep(4, 4),
     c(0.38,0.5,0.62,0.74),
     as.character(dots_for))
par(xpd=F)
dev.copy(png, "behaviour broods - BCM.png",res=300,width=1800, height=1100)
dev.off()
