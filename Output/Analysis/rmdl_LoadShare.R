rm(list=ls())
library(readr)
library(ggplot2)
library(gridExtra)
library(plyr)


########################################################################################
### Read data from remodelling of both bone and interface
########################################################################################

Rmdl_11 <- read_csv("OutputLoadCentrality_RMDL11.txt")

Rmdl_11$subjectCode <- as.factor(Rmdl_11$subjectCode)
Rmdl_11$LegSide <- as.factor(Rmdl_11$LegSide)
Rmdl_11$ImplantSize <- as.factor(Rmdl_11$ImplantSize)
Rmdl_11$alpha <- as.factor(Rmdl_11$alpha)
Rmdl_11$Lshare100 <- 100*Rmdl_11$Lshare
Rmdl_11$Rmdl.TBCMT <- rep("Yes",length(Rmdl_11$Lshare))
Rmdl_11$Rmdl.BONE <- rep("Yes",length(Rmdl_11$Lshare))
Rmdl_11$Rmdl <- rep("Both",length(Rmdl_11$Lshare))

########################################################################################
### Read data from remodelling of both bone only
########################################################################################

Rmdl_10 <- read_csv("OutputLoadCentrality_RMDL10.txt")

Rmdl_10$subjectCode <- as.factor(Rmdl_10$subjectCode)
Rmdl_10$LegSide <- as.factor(Rmdl_10$LegSide)
Rmdl_10$ImplantSize <- as.factor(Rmdl_10$ImplantSize)
Rmdl_10$alpha <- as.factor(Rmdl_10$alpha)
Rmdl_10$Lshare100 <- 100*Rmdl_10$Lshare
Rmdl_10$Rmdl.TBCMT <- rep("No",length(Rmdl_10$Lshare))
Rmdl_10$Rmdl.BONE <- rep("Yes",length(Rmdl_10$Lshare))
Rmdl_10$Rmdl <- rep("Both",length(Rmdl_10$Lshare))


########################################################################################
### Read data from remodelling of both bone only
########################################################################################
# To be added ?

########################################################################################
### Merge DataSet
########################################################################################

DF = rbind(Rmdl_11,Rmdl_10)
DF$Rmdl.TBCMT <- as.factor(DF$Rmdl.TBCMT)
DF$Rmdl.BONE <- as.factor(DF$Rmdl.BONE)
DF$LoadStep <- as.factor(DF$LoadStep)
DF$Mois <- DF$Epoch
DF$Alignement <- as.factor(DF$Alignement)
revalue(DF$Alignement, c("Kine"="Cinématique", "Mech"="Mécanique")) -> DF$Alignement
DF$Epoch <- as.numeric(DF$Epoch)
DF$Align.Rmdl.TBCMT <- paste(DF$Alignement, "-",  DF$Rmdl.TBCMT)
DF$Align.Rmdl.Bone <- paste(DF$Alignement, "-",  DF$Rmdl.BONE)
DF <- subset(DF, Epoch<72 )


DF_mean <- aggregate(Lshare100 ~ Mois + Alignement + Align.Rmdl.TBCMT, DF, FUN = 'mean')
########################################################################################
### Bar PLOTS  , color = Step , color = SubjectCode
########################################################################################
bar1 <- ggplot(DF, aes(x = Mois, y = Lshare100 , fill = Alignement)) + 
  geom_bar(stat="identity", position="dodge")+
  theme_classic() + xlab('Mois') + ylab('Load Bypass (%)') + 
  theme(axis.text=element_text(size=14), axis.title = element_text(size=20))

bar1 + facet_grid(Rmdl.TBMCT ~ LoadStep)

bar_all <- ggplot(DF, aes(x = Mois, y = Lshare100 , fill = Align.Rmdl.TBCMT)) + 
  geom_bar(stat="identity", position="dodge")+
  theme_classic() + xlab('Mois') + ylab('Load Bypass (%)') + 
  theme(axis.text=element_text(size=14), axis.title = element_text(size=20))
bar_all + facet_grid(. ~ LoadStep)

bar_all <- ggplot(DF_mean, aes(x = Mois, y = Lshare100 , fill = Align.Rmdl.TBCMT)) + 
  geom_bar(stat="identity", position="dodge")+
  theme_classic() + xlab('Mois') + ylab('Load Bypass (%)') + 
  theme(axis.text=element_text(size=14), axis.title = element_text(size=20))
bar_all


########################################################################################
### Point PLOTS  , color = Step , color = SubjectCode
########################################################################################
DF <- subset(DF, Epoch<72 & Epoch>-0.5   )
DF$Mois <- as.numeric(DF$Mois) 
DF_mean <- aggregate(Lshare100 ~ Mois + Alignement + Align.Rmdl.TBCMT + Rmdl.TBCMT, DF, FUN = 'mean')
#Varus
pt <- ggplot(DF_mean, aes(Mois, Lshare100 , color = Alignement, shape = Rmdl.TBCMT)) + 
  geom_point(size=2) + theme_classic() +
  geom_smooth(method = 'loess', se=FALSE) +
  xlab('Mois') +
  ylab('Load Bypass (%)') + 
  theme(axis.text=element_text(size=14), axis.title = element_text(size=20))
pt

point_all <- ggplot(DF_mean, aes(x = Mois, y = Lshare100 , color = Align.Rmdl.TBCMT)) + 
  geom_point(size=2) + theme_classic() +
  geom_smooth(method = 'loess', se=FALSE) + 
  theme_classic() + xlab('Mois') + ylab('Load Bypass (%)') + 
  theme(axis.text=element_text(size=14), axis.title = element_text(size=20))
point_all
