rm(list=ls())
library(readr)
library(ggplot2)
library(gridExtra)

DF <- read_csv("OutputLoadCentrality.txt")

DF$subjectCode <- as.factor(DF$subjectCode)
DF$LegSide <- as.factor(DF$LegSide)
DF$ImplantSize <- as.factor(DF$ImplantSize)
DF$alpha <- as.factor(DF$alpha)
DF$Lshare100 <- 100*DF$Lshare


########################################################################################
### GENERAL PLOTS  , color = Step , color = SubjectCode
########################################################################################

#Varus
ggplot(DF, aes(Epoch, Lshare100 , color = alpha)) + 
  geom_point(size=2) + theme_classic() +
  stat_smooth(method = "lm", formula = y ~ x + I(x^2), size = 1, se=FALSE) +
  xlab('Epoch') +
  ylab('Load Bypass (%)') + 
  theme(axis.text=element_text(size=14), axis.title = element_text(size=20))

ggplot(DF, aes(x = Epoch, y = Lshare100 , fill = alpha)) + 
  geom_bar(stat="identity", position=position_dodge())+
  theme_classic() + xlab('Epoch') + ylab('Load Bypass (%)') + 
  theme(axis.text=element_text(size=14), axis.title = element_text(size=20))
