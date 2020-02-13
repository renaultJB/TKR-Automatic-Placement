rm(list=ls())

library(readr)
library(relaimpo)
library(ggplot2)
library(gridExtra)


DF <- read_csv("~/These/EuroMech/analyse/AnalyseCentrality/OutputLoadCentrality.txt")
DF <- read_csv("~/These/EuroMech/analyse/AnalyseCentrality/OutputLoadCentralityWCCM.txt")

DF <- DF[ DF$Min.Mean > 0.05, ]

DF$subjectCode <- as.factor(DF$subjectCode)
DF$LegSide <- as.factor(DF$LegSide)
DF$ImplantSize <- as.factor(DF$ImplantSize)
DF$PCAGaitN <- as.factor(DF$PCAGaitN)

DF$CV <- 100*DF$CV
DF$Lshare <- 100*DF$Lshare
DF$Ratio.MM <- 100*DF$Min.Mean


########################################################################################
# Fit 1 
lmRMMCV <- lm(Lshare ~ Ratio.MM + CV, data = DF )
interceptFit <- summary(lmRMMCV)$coefficients[1]
RMM_C <- summary(lmRMMCV)$coefficients[2]
CV_C <- summary(lmRMMCV)$coefficients[3]
summary(lmRMMCV)
summary(lmRMMCV)$r.squared

# New Fit variable : 
DF$Fit.RMM.CV <- 100-6*(interceptFit + RMM_C*DF$Ratio.MM + CV_C*DF$CV - 5.5)
lmFitRMMCV <- lm(Lshare ~ Fit.RMM.CV, data = DF )
summary(lmFitRMMCV)$r.squared
########################################################################################

########################################################################################
# Fit quadratic alpha
lmAA2 <- lm(Lshare ~ alpha + I(alpha^2), data = DF )
summary(lmAA2)$r.squared

# Fit simple MM 
lmMM <- lm(Lshare ~ Ratio.MM, data = DF )
summary(lmMM)
summary(lmMM)$r.squared

# Fit simple CV 
lmCV <- lm(Lshare ~ CV, data = DF )
summary(lmCV)
summary(lmCV)$r.squared

########################################################################################

########################################################################################
### GENERAL PLOTS  , color = Step , color = SubjectCode
########################################################################################

#Varus
ggplot(DF, aes(alpha, Lshare , color = Step, shape = PCAGaitN)) + 
  geom_point(size=2) + theme_classic() +
  stat_smooth(method = "lm", formula = y ~ x + I(x^2), size = 1, se=FALSE) +
  xlab('Valgus angle (°)') +
  ylab('Load Bypass (%)') + 
  theme(axis.text=element_text(size=14), axis.title = element_text(size=20))

#Ratio.MM
ggplot(DF, aes(Ratio.MM, Lshare , color = Step, shape = PCAGaitN)) + 
  geom_point(size=2) + theme_classic() +
  geom_smooth(method='lm',se=FALSE) +
  xlab('Centrality (M/M) (%)') +
  ylab('Load Bypass (%)') + 
  theme(axis.text=element_text(size=14), axis.title = element_text(size=20))

#CV
ggplot(DF, aes(CV, Lshare , color = Step, shape = PCAGaitN)) + 
  geom_point(size=2) + theme_classic() +
  geom_smooth(method='lm',se=FALSE) +
  xlab('Centrality (CV)(%)') +
  ylab('Load Bypass (%)') + 
  theme(axis.text=element_text(size=14), axis.title = element_text(size=20))

#Fit
ggplot(DF, aes(Fit.RMM.CV, Lshare , color = Step, shape = PCAGaitN)) + 
  geom_point(size=2) + theme_classic() +
  geom_smooth(method='lm',se=FALSE) +
  xlab('Centrality (Fitted from M/M and CV) (%)') +
  ylab('Load Bypass (%)') + 
  theme(axis.text=element_text(size=14), axis.title = element_text(size=20))



########################################################################################
### GaitCycle Plots
########################################################################################

#Varus
ggplot(DF, aes(alpha, Lshare, color = PCAGaitN)) + 
  geom_point(size=2) + theme_classic() +
  stat_smooth(method = "lm", formula = y ~ x + I(x^2), size = 1, se=FALSE) +
  xlab('Valgus angle (°)') +
  ylab('Load Bypass (%)') + 
  theme(axis.text=element_text(size=14), axis.title = element_text(size=20))

#Ratio.MM
ggplot(DF, aes(Ratio.MM, Lshare, color = PCAGaitN)) + 
  geom_point(size=2) + theme_classic() +
  geom_smooth(method='lm',se=FALSE) +
  xlab('Centrality (M/M) (%)') +
  ylab('Load Bypass (%)') + 
  theme(axis.text=element_text(size=14), axis.title = element_text(size=20))

#CV
ggplot(DF, aes(CV, Lshare, color = PCAGaitN)) + 
  geom_point(size=2) + theme_classic() +
  geom_smooth(method='lm',se=FALSE) +
  xlab('Centrality (CV)(%)') +
  ylab('Load Bypass (%)') + 
  theme(axis.text=element_text(size=14), axis.title = element_text(size=20))

#Fit
ggplot(DF, aes(Fit.RMM.CV, Lshare, color = PCAGaitN)) + 
  geom_point(size=2) + theme_classic() +
  geom_smooth(method='lm',se=FALSE) +
  xlab('Centrality (Fitted from M/M and CV) (%)') +
  ylab('Load Bypass (%)') + 
  theme(axis.text=element_text(size=14), axis.title = element_text(size=20))

########################################################################################
### Subject Plots
########################################################################################

#Varus
ggplot(DF, aes(alpha, Lshare, color = subjectCode)) + 
  geom_point(size=2) + theme_classic() +
  stat_smooth(method = "lm", formula = y ~ x + I(x^2), size = 1, se=FALSE) +
  xlab('Valgus angle (°)') +
  ylab('Load Bypass (%)') + 
  theme(axis.text=element_text(size=14), axis.title = element_text(size=20))

#Ratio.MM
ggplot(DF, aes(Ratio.MM, Lshare, color = subjectCode)) + 
  geom_point(size=2) + theme_classic() +
  geom_smooth(method='lm',se=FALSE) +
  xlab('Centrality (M/M) (%)') +
  ylab('Load Bypass (%)') + 
  theme(axis.text=element_text(size=14), axis.title = element_text(size=20))

#CV
ggplot(DF, aes(CV, Lshare, color = subjectCode)) + 
  geom_point(size=2) + theme_classic() +
  geom_smooth(method='lm',se=FALSE) +
  xlab('Centrality (CV)(%)') +
  ylab('Load Bypass (%)') + 
  theme(axis.text=element_text(size=14), axis.title = element_text(size=20))

#Fit
ggplot(DF, aes(Fit.RMM.CV, Lshare, color = subjectCode)) + 
  geom_point(size=2) + theme_classic() +
  geom_smooth(method='lm',se=FALSE) +
  xlab('Centrality (Fitted from M/M and CV) (%)') +
  ylab('Load Bypass (%)') + 
  theme(axis.text=element_text(size=14), axis.title = element_text(size=20))


########################################################################################
##### CORRELATION BETWEEN ALPHA AND CENTRALATIES
########################################################################################

#Ratio.MM vs alpha
ggplot(DF, aes(alpha, Ratio.MM, color = subjectCode)) + 
  geom_point(size=2) + theme_classic() +
  stat_smooth(method = "lm", formula = y ~ x + I(x^2), size = 1, se=FALSE) +
  ylab('Centrality (M/M) (%)') +
  xlab('Valgus angle (°)') + 
  theme(axis.text=element_text(size=14), axis.title = element_text(size=20))

#CV vs alpha
ggplot(DF, aes(alpha, CV, color = subjectCode)) + 
  geom_point(size=2) + theme_classic() +
  stat_smooth(method = "lm", formula = y ~ x + I(x^2), size = 1, se=FALSE) +
  ylab('Centrality (CV)(%)') +
  xlab('Valgus angle (°)') + 
  theme(axis.text=element_text(size=14), axis.title = element_text(size=20))

#Fit vs alpha
ggplot(DF, aes(alpha, Ratio.MM, color = subjectCode)) + 
  geom_point(size=2) + theme_classic() +
  stat_smooth(method = "lm", formula = y ~ x + I(x^2), size = 1, se=FALSE) +
  ylab('Centrality (Fitted from M/M and CV) (%)') +
  xlab('Valgus angle (°)') + 
  theme(axis.text=element_text(size=14), axis.title = element_text(size=20))