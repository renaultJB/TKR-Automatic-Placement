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
Rmdl_11$alpha.F <- as.factor(Rmdl_11$alpha)
Rmdl_11$Epoch.F <- as.factor(Rmdl_11$Epoch)
Rmdl_11$Lshare100 <- 100*Rmdl_11$Lshare
Rmdl_11$Lshare.Top <- 100*Rmdl_11$Lsup/Rmdl_11$Ltot
Rmdl_11$Min.Mean <- 2*Rmdl_11$Min.Mean
Rmdl_11$Rmdl.TBCMT <- rep("Yes",length(Rmdl_11$Lshare))
Rmdl_11$Rmdl.BONE <- rep("Yes",length(Rmdl_11$Lshare))
Rmdl_11$Rmdl <- rep("Both",length(Rmdl_11$Lshare))
revalue(Rmdl_11$Alignement, c("Kine"="Cinématique", "Mech"="Mécanique")) -> Rmdl_11$Alignement

Rmdl_11_init <- Rmdl_11





# Rmdl_11 <- subset(Rmdl_11, subjectCode == 'GUI_R' )

########################################################################################
### Plot Lshare vs centrality
########################################################################################
# & LoadStep == "GC_50"
Rmdl_11 <- subset(Rmdl_11, Epoch == 0 | Epoch == 12 | Epoch == 36 | Epoch == 72 & LoadStep != "CU"  )
Rmdl_11 <- aggregate(Lshare100 ~ Epoch.F + alpha + Min.Mean + Alignement + CV , Rmdl_11, FUN = 'mean')
ggplot(Rmdl_11, aes(Min.Mean, Lshare100 , color = Epoch.F, fill = Epoch.F )) + 
  geom_point(size=3) + theme_classic() +
  stat_smooth(method = "lm", formula = y ~ x + I(x^2), size = 1) +
  xlab('Centrality.MM') + ylab('Fprop (%)') +
  theme(axis.text=element_text(size=16), axis.title = element_text(size=16)) +
  theme(legend.text=element_text(size=16)) +
  labs(fill = "Mois") + labs(color = "Mois")

ggplot(Rmdl_11, aes(Min.Mean, Lshare100 , color = Epoch.F, fill = Epoch.F )) + 
  geom_point(size=3) + theme_classic() +
  stat_smooth(method = "lm") +
  xlab('Centrality.MM') + ylab('Fprop (%)') +
  theme(axis.text=element_text(size=16), axis.title = element_text(size=16)) +
  theme(legend.text=element_text(size=16)) +
  labs(fill = "Mois") + labs(color = "Mois")

ggplot(Rmdl_11, aes(alpha, Lshare100 , color = Epoch.F, fill = Epoch.F )) + 
  geom_point(size=3) + theme_classic() +
  stat_smooth(method = "lm", formula = y ~ x + I(x^2), size = 1) +
  xlab('HKA (°)') + ylab('Fprop (%)') +
  theme(axis.text=element_text(size=16), axis.title = element_text(size=16)) +
  theme(legend.text=element_text(size=16)) +
  labs(fill = "Mois") + labs(color = "Mois")

Rmdl_11 <- Rmdl_11_init
Rmdl_11 <- subset(Rmdl_11, Epoch == 0 | Epoch == 12 | Epoch == 36 | Epoch == 72 & LoadStep != "CU"  )
Rmdl_11 <- aggregate(Lshare.Top ~ Epoch.F + Min.Mean + Alignement + CV , Rmdl_11, FUN = 'mean')
ggplot(Rmdl_11, aes(Min.Mean, Lshare.Top , color = Epoch.F, fill = Epoch.F )) + 
  geom_point(size=3) + theme_classic() +
  stat_smooth(method = "lm", formula = y ~ x + I(x^2), size = 1) +
  xlab('Centrality.MM') + ylab('Fprop Top (%)') +
  theme(axis.text=element_text(size=16), axis.title = element_text(size=16)) +
  theme(legend.text=element_text(size=16)) +
  labs(fill = "Mois") + labs(color = "Mois")


########################################################################################
### Plot Lshare vs centrality
########################################################################################
Rmdl_11 <- Rmdl_11_init
Rmdl_11 <- subset(Rmdl_11, (alpha = alpha != 3.5= 0 | alpha == -3.48)  & subjectCode == 'TRO_L'  )


bar_all <- ggplot(Rmdl_11, aes(x = Epoch.F, y = Lshare100 , fill = Alignement)) + 
  geom_bar(stat="identity", position="dodge")+
  theme_classic() + xlab('Mois') + ylab('Fprop (%)') + 
  theme(axis.text=element_text(size=20), axis.title = element_text(size=20))
bar_all + facet_grid(. ~ LoadStep)

Rmdl_11 <- Rmdl_11_init
Rmdl_11 <- subset(Rmdl_11, (alpha != -1.75 | alpha != 1.75 | alpha != 3.5))
bar_all <- ggplot(Rmdl_11, aes(x = Epoch.F, y = Lshare100 , fill = Alignement)) + 
  geom_bar(stat="identity", position="dodge")+
  theme_classic() + xlab('Mois') + ylab('Fprop (%)') + 
  theme(axis.text=element_text(size=12), axis.title = element_text(size=16))
bar_all + facet_grid(subjectCode ~ LoadStep) +
  theme(strip.text.x = element_text(size = 16)) +
  theme(strip.text.y = element_text(size = 16))

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
revalue(Rmdl_10$Alignement, c("Kine"="Cinématique", "Mech"="Mécanique")) -> Rmdl_10$Alignement



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

########################################################################################
### Bar PLOTS  , color = Step , color = SubjectCode
########################################################################################
bar1 <- ggplot(DF, aes(x = Mois, y = Lshare100 , fill = Alignement)) + 
  geom_bar(stat="identity", position="dodge")+
  theme_classic() + xlab('Mois') + ylab('Load Bypass (%)') + 
  theme(axis.text=element_text(size=14), axis.title = element_text(size=20))

bar1 + facet_grid(Rmdl.TBCMT ~ LoadStep)


DF2 <- subset(DF, Rmdl.TBCMT = "Yes" )
DF2$Epoch <- as.factor(DF2$Epoch)
bar_all <- ggplot(DF2, aes(x = Epoch, y = Lshare100 , fill = Alignement)) + 
  geom_bar(stat="identity", position="dodge")+
  theme_classic() + xlab('Mois') + ylab('Fprop (%)') + 
  theme(axis.text=element_text(size=20), axis.title = element_text(size=20))
bar_all + facet_grid(. ~ LoadStep)

DF_mean <- aggregate(Lshare100 ~ Mois + Alignement + Align.Rmdl.TBCMT + Rmdl.TBCMT, DF, FUN = 'mean')
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
pt <- ggplot(DF_mean, aes(Mois, Lshare100 , color = Rmdl.TBCMT, shape = Alignement)) + 
  geom_point(size=4) + theme_classic() +
  geom_smooth(method = 'loess', se=FALSE) +
  xlab('Mois') +
  ylab('Fprop (%)') + 
  theme(axis.text=element_text(size=16), axis.title = element_text(size=16)) +
  theme(legend.text=element_text(size=16), legend.position= c(0.8, 0.2))
pt

point_all <- ggplot(DF_mean, aes(x = Mois, y = Lshare100 , color = Align.Rmdl.TBCMT)) + 
  geom_point(size=2) + theme_classic() +
  geom_smooth(method = 'loess', se=FALSE) + 
  theme_classic() + xlab('Mois') + ylab('Load Bypass (%)') + 
  theme(axis.text=element_text(size=14), axis.title = element_text(size=20))
point_all
