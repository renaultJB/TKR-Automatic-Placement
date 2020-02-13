rm(list=ls())
library(readr)
library(ggplot2)
library(gridExtra)

DF <- read_csv("OutputLoadCentrality.txt")
# DF <- DF[ DF$Min.Mean > 0.05, ]


DF$subjectCode <- as.factor(DF$subjectCode)
DF$LegSide <- as.factor(DF$LegSide)
DF$ImplantSize <- as.factor(DF$ImplantSize)
DF$
DF$Corctalpha <-  -DF$alpha -DF$Varus
DF$CV <- 100*DF$CV
DF$Ctr <- DF$CV + 1/DF$Mean.Distance
DF$Min.Min <- DF$Min.Mean*DF$Mean.Distance

DF$Lshare <- 100*DF$Lshare

DF$Ratio.Poutre <- DF$Min.Min/(2*DF$Mean.Distance)

# LIM_L <- subset(x = DF, subjectCode == 'DES_R')

DF$Th.Poutre <- 12*(DF$Min.Min/(2*DF$Mean.Distance))^2 - 12*(DF$Min.Min/(2*DF$Mean.Distance)) + 4

DF$Th.Poutre2 <- 1/16*(1-3*DF$Ratio.Poutre+3*DF$Ratio.Poutre^2)/(DF$Ratio.Poutre^3*(1-DF$Ratio.Poutre)^3)

DF$Th.Poutre3 <- 3*1000/DF$Mean.Distance^3*DF$Th.Poutre2

DF$alpha.squared <- DF$alpha^2





DF$Ratio.MM <- DF$Min.Min/(2*DF$Mean.Distance)

DF$Ratio.MM <- 100*DF$Min.Mean

DF$Ratio.MM.squared <- DF$Min.Min/(2*DF$Mean.Distance)^2

DF$Ratio.MM.inv <- 1/(DF$Min.Min/(2*DF$Mean.Distance))



#DF$Th.Plaque <- (DF$Min.Min/(DF$Mean.Distance*2))^-0.41

DF$Fit.PoutreCV <- 0.2142 + 0.2450*DF$Th.Poutre + 0.2679*DF$CV 

DF$Fit.Ratio.poly2 <- 0.2944 + 4.4275*DF$Ratio.MM.squared - 0.7074*DF$Ratio.MM

DF$Fit.Ratio.Power <- 0.5 -0.5013*DF$Ratio.MM^0.2251

DF$Fit.Poutre.inv <- 0.2423073 + 9.5718569*DF$Ratio.MM.squared - 0.6794561*DF$Ratio.MM + 0.0006922*DF$Ratio.MM.inv

# DF$Fit.RatioMM.CV <- 0.1730 + -0.4101*DF$Ratio.MM + 0.3015*DF$CV


DF$Fit.RatioMM.CV <- 0.1534 + -0.2794*DF$Ratio.MM + 0.2542*DF$CV

DF$Fit.alpha.squared <- 0.1667430 + 0.0004004*DF$alpha.squared + 0.0034765*DF$alpha

DF$Fit.MinMax.CV <- 0.1556 + 0.2087*DF$CV -0.1847*DF$Min.Max


DF$Global<-1/(1/10000+1/(DF$Th.Poutre3)+1/20000)


lmAlphaS <- lm(Lshare ~ subjectCode + alpha + I(alpha^2), data = DF )
summary(lmAlphaS)
summary(lmAlphaS)$r.squared

lme.AlphaS <- lme(Lshare ~ alpha+ I(alpha^2), data = DF, random = ~1|subjectCode)
summary(lme.AlphaS)

lmRMM <- lm(Lshare ~ Ratio.MM , data = DF )
summary(lmRMM)
summary(lmRMM)$r.squared

lme.RMM <- lme(Lshare ~ Ratio.MM , data = DF , random = ~1|subjectCode)
summary(lme.RMM)

lmRMMS <- lm(Lshare ~ Ratio.MM + I(Ratio.MM^2), data = DF )
summary(lmRMMS)
summary(lmRMMS)$r.squared


######################################
###             Raw Data           ###
######################################

# plot the result for the raw data

ggplot(DF, aes(CV, Lshare)) + 
  geom_point() + theme_classic() +
  geom_smooth(method='lm')

ggplot(DF, aes(Th.Poutre, Lshare)) + 
  geom_point() + theme_classic() +
  geom_smooth(method='lm')

ggplot(DF, aes(alpha, Lshare,)) + 
  geom_point(size=2) + theme_classic() +
  stat_smooth(method = "lm", formula = y ~ x + I(x^2), size = 1) +
  xlab('Valgus angle (°)') +
  ylab('Load Bypass (%') + 
  theme(axis.text=element_text(size=14), axis.title = element_text(size=20))

ggplot(DF, aes(CV, Lshare, color = Step)) + 
  geom_point(size=2) + theme_classic() +
  geom_smooth(method='lm') +
  xlab('Centrality (%)') +
  ylab('Load Bypass (%)') + 
  theme(axis.text=element_text(size=14), axis.title = element_text(size=20))

fit <- lm(Lshare ~ DF$alpha + I(DF$alpha^2), data = DF)
anova(fit)

fit <- lm(Lshare ~ DF$Ratio.MM + DF$Ratio.MM.squared , data = DF)
anova(fit)

fit <- lm(Lshare ~ DF$Ratio.MM + DF$Ratio.MM.squared + DF$Ratio.MM.inv , data = DF)
anova(fit)


fit <- lm(Lshare ~ CV + Mean.Distance, data = DF)
anova(fit)

fit <- lm(Lshare ~ Ratio.MM + alpha + CV, data = DF)
cor(y=DF$Lshare, x=DF$Th.Poutre)
anova(fit)

relImportance <- calc.relimp(fit, type = "lmg", rela = TRUE)
sort(relImportance$lmg, decreasing=TRUE)

relImportance <- calc.relimp(fit, type = "pmvd", rela = TRUE)
sort(relImportance$pmvd, decreasing=TRUE)

# 

cor(y=DF$Lshare, x=DF$Ratio.MM)
cor(y=DF$Lshare, x=DF$CV)

fit <- lm(Lshare ~ CV + Ratio.MM + Min.Min + Mean.Distance + Min.Max + Min.Mean + Ratio.AP + Ratio.ML, data = DF)
anova(fit)

relImportance <- calc.relimp(fit, type = "lmg", rela = TRUE)
sort(relImportance$lmg, decreasing=TRUE)

relImportance <- calc.relimp(fit, type = "pmvd", rela = TRUE)
sort(relImportance$pmvd, decreasing=TRUE)




######################################
###          Scaled Data           ###
######################################

# Scale data

scaled_DF<-data.frame(scale(DF[, 4:length(DF)]))
scaled_DF$subjectCode<-DF$subjectCode
scaled_DF$LegSide <- DF$LegSide
scaled_DF$ImplantSize <- DF$ImplantSize


# plot the result for the raw data




######################################
## Scaled Data only for S5 implant  ##
######################################

# Subste of the data data
scaled_DF_S5 <- subset(scaled_DF, ImplantSize=='S5' )

# Plots
ggplot(scaled_DF_S5, aes(CV, Lshare)) + 
  geom_point() + theme_classic() +
  geom_smooth(method='lm')

ggplot(scaled_DF_S5, aes(Th.Poutre, Lshare)) + 
  geom_point() + theme_classic() +
  geom_smooth(method='lm')

ggplot(scaled_DF_S5, aes(alpha, Lshare)) + 
  geom_point() + theme_classic() +
  geom_smooth(method='lm')




##########################################
## Scaled Data only for alpha -5 to +5  ##
##########################################

# Subste of the data data for Low centered alpha
DF_LA <- subset(DF, alpha > -6 )
scaled_DF_LA<-data.frame(scale(DF_LA[, 4:length(DF_LA)]))
scaled_DF_LA$subjectCode<-DF_LA$subjectCode
scaled_DF_LA$LegSide <- DF_LA$LegSide
scaled_DF_LA$ImplantSize <- DF_LA$ImplantSize




# Plots
ggplot(scaled_DF_LA, aes(CV, Lshare)) + 
  geom_point() + theme_classic() +
  geom_smooth(method='lm')

ggplot(scaled_DF_LA, aes(Th.Poutre, Lshare)) + 
  geom_point() + theme_classic() +
  geom_smooth(method='lm')

ggplot(scaled_DF_LA, aes(alpha, Lshare)) + 
  geom_point() + theme_classic() +
  geom_smooth(method='lm')


##########################################
## Scaled Data only for varus only  ##
##########################################

# Subste of the data data for Low centered alpha
DF_Ovarus <- subset(DF, alpha < 1.0 ) #& alpha > -6.0
scaled_DF_Ovarus<-data.frame(scale(DF_Ovarus[, 4:length(DF_Ovarus)]))
scaled_DF_Ovarus$subjectCode<-DF_LA$subjectCode
scaled_DF_Ovarus$LegSide <- DF_LA$LegSide
scaled_DF_Ovarus$ImplantSize <- DF_LA$ImplantSize




# Plots
ggplot(scaled_DF_Ovarus, aes(CV, Lshare)) + 
  geom_point() + theme_classic() +
  geom_smooth(method='lm')

ggplot(scaled_DF_Ovarus, aes(Th.Poutre, Lshare)) + 
  geom_point() + theme_classic() +
  geom_smooth(method='lm')

ggplot(scaled_DF_Ovarus, aes(alpha, Lshare)) + 
  geom_point() + theme_classic() +
  geom_smooth(method='lm')


##########################################
##        Well aligned Data only        ##
##########################################

# Subste of the data data for Low centered alpha
DF_WA <- subset(DF, alpha < 3.0 & alpha > -3.0 ) #& alpha > -6.0

# Pearsons correlation coefficient
cor(y=scaled_DF_WA$Lshare, x=scaled_DF_WA$Th.Poutre)
cor(y=scaled_DF_WA$Lshare, x=scaled_DF_WA$CV)

# Plots
ggplot(scaled_DF_WA, aes(alpha, Lshare)) + 
  geom_point() + theme_classic() +
  geom_smooth(method='lm')

ggplot(scaled_DF_WA, aes(CV, Lshare)) + 
  geom_point() + theme_classic() +
  geom_smooth(method='lm')

ggplot(scaled_DF_WA, aes(Th.Poutre, Lshare)) + 
  geom_point() + theme_classic() +
  geom_smooth(method='lm')

ggplot(scaled_DF_WA, aes(Th.Poutre, Lshare)) + 
  geom_point() + theme_classic() +
  geom_smooth(method='lm')


##########################################
##              Ok Data only            ##
##########################################

# Subste of the data data for Low centered alpha
DF_Ok <- subset(DF, Ratio.MM > 0.075 ) #& alpha > -6.0

# 
# DF_Ok$Fit.RatioMM.CV <- 0.1534 + -0.2794*DF_Ok$Ratio.MM + 0.2542*DF_Ok$CV
# 
# fit <- lm(Lshare ~ Ratio.MM + CV , data = DF_Ok)
# 

# Pearson correlation coefficient
cor(y=DF_Ok$Lshare, x=DF_Ok$Ratio.MM)
cor(y=DF_Ok$Lshare, x=DF_Ok$CV)
cor(y=DF_Ok$Lshare, x=DF_Ok$Fit.RatioMM.CV)

AAA <- count(DF_Ok, vars = subjectCode)

# Plots
ggplot(DF_Ok, aes(CV, Lshare)) + 
  geom_point() + theme_classic() +
  geom_smooth(method='lm')

ggplot(DF_Ok, aes(Ratio.MM, Lshare)) + 
  geom_point() + theme_classic() +
  geom_smooth(method='lm')

ggplot(DF_Ok, aes(alpha, Lshare)) + 
  geom_point() + theme_classic() +
  geom_smooth(method='lm')

ggplot(DF_Ok, aes(Fit.RatioMM.CV, Lshare)) + 
  geom_point() + theme_classic() +
  geom_smooth(method='lm')

ggplot(DF_Ok, aes(alpha, Lshare,colour=subjectCode)) + 
  geom_point() + theme_classic() +
  geom_smooth(method='lm',se=FALSE)

ggplot(DF_Ok, aes(CV, Lshare,colour=subjectCode)) + 
  geom_point() + theme_classic() +
  geom_smooth(method='lm',se=FALSE)

ggplot(DF_Ok, aes(Ratio.MM, Lshare,colour=subjectCode)) + 
  geom_point() + theme_classic() +
  geom_smooth(method='lm',se=FALSE)

ggplot(DF_Ok, aes(Fit.RatioMM.CV, Lshare,colour=subjectCode)) + 
  geom_point() + theme_classic() +
  geom_smooth(method='lm',se=FALSE)

##########################################
##              Ok Data only  Scaled          ##
##########################################

# Subste of the data data for Low centered alpha
DF_Ok <- subset(DF, Min.Mean > 0.2 ) #& alpha > -6.0
scaled_DF_Ok<-data.frame(scale(DF_Ok[, 4:length(DF_Ok)]))
scaled_DF_Ok$subjectCode<-DF_Ok$subjectCode
scaled_DF_Ok$LegSide <- DF_Ok$LegSide
scaled_DF_Ok$ImplantSize <- DF_Ok$ImplantSize


scaled_DF_Ok$Fit.RatioMM.CV <- -4.300e-01*scaled_DF_Ok$Ratio.MM + 4.482e-01*scaled_DF_Ok$CV

fit <- lm(Lshare ~ Ratio.MM + CV , data = scaled_DF_Ok)


# Pearson correlation coefficient
cor(y=scaled_DF_Ok$Lshare, x=scaled_DF_Ok$Th.Poutre)
cor(y=scaled_DF_Ok$Lshare, x=scaled_DF_Ok$CV)

AAA <- count(scaled_DF_Ok, vars = subjectCode)

# Plots
ggplot(scaled_DF_Ok, aes(CV, Lshare)) + 
  geom_point() + theme_classic() +
  geom_smooth(method='lm')

ggplot(scaled_DF_Ok, aes(Th.Poutre, Lshare)) + 
  geom_point() + theme_classic() +
  geom_smooth(method='lm')

ggplot(scaled_DF_Ok, aes(alpha, Lshare)) + 
  geom_point() + theme_classic() +
  geom_smooth(method='lm')

ggplot(scaled_DF_Ok, aes(Th.Poutre, Lshare,colour=subjectCode)) + 
  geom_point() + theme_classic() +
  geom_smooth(method='lm',se=FALSE)

ggplot(scaled_DF_Ok, aes(alpha, Lshare,colour=subjectCode)) + 
  geom_point() + theme_classic() +
  geom_smooth(method='lm',se=FALSE)

ggplot(scaled_DF_Ok, aes(CV, Lshare,colour=subjectCode)) + 
  geom_point() + theme_classic() +
  geom_smooth(method='lm',se=FALSE)

ggplot(scaled_DF_Ok, aes(Ratio.MM, Lshare,colour=subjectCode)) + 
  geom_point() + theme_classic() +
  geom_smooth(method='lm',se=FALSE)

ggplot(scaled_DF_Ok, aes(Th.Poutre, Lshare,colour=subjectCode)) + 
  geom_point() + theme_classic() +
  geom_smooth(method='lm',se=FALSE)

ggplot(scaled_DF_Ok, aes(Fit.RatioMM.CV, Lshare,colour=subjectCode)) + 
  geom_point() + theme_classic() +
  geom_smooth(method='lm',se=FALSE)


##########################################
##         T.Test and boxplots          ##
##########################################


DF_Ok$Outlier.alpha <- c('Yes', 'No')[(DF_Ok$alpha>-3 & DF_Ok$alpha<3)+1L]
DF_Ok$Outlier.alpha <- as.factor(DF_Ok$Outlier.alpha)

DF_Ok$Outlier.CV <- c('Yes', 'No')[(DF_Ok$CV<0.25)+1L]
DF_Ok$Outlier.CV <- as.factor(DF_Ok$Outlier.CV)

DF_Ok$Outlier.Ratio.MM <- c('Yes', 'No')[(DF_Ok$Ratio.MM>0.20)+1L]
DF_Ok$Outlier.Ratio.MM <- as.factor(DF_Ok$Outlier.Ratio.MM)

DF_Ok$Outlier.Fit.RatioMM.CV <- c('Yes', 'No')[(DF_Ok$Fit.RatioMM.CV>0.16)+1L]
DF_Ok$Outlier.Fit.RatioMM.CV <- as.factor(DF_Ok$Outlier.Fit.RatioMM.CV)


# t test
t.test(Lshare~Outlier.alpha, data = DF_Ok)
t.test(Lshare~Outlier.CV, data = DF_Ok)
t.test(Lshare~Outlier.Ratio.MM, data = DF_Ok)
t.test(Lshare~Outlier.Fit.RatioMM.CV, data = DF_Ok)

boxplot(Lshare~Outlier.alpha, data = DF_Ok)
boxplot(Lshare~Outlier.CV, data = DF_Ok)
t.test(Lshare~Outlier.Ratio.MM, data = DF_Ok)
boxplot(Lshare~Outlier.Fit.RatioMM.CV, data = DF_Ok)




##############################
##        Functions         ##
##############################


get_legend<-function(myggplot){
  tmp <- ggplot_gtable(ggplot_build(myggplot))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]
  return(legend)
}

write_Tiff_SubjectCode <- function(DF,name){
  p1 <- ggplot(DF, aes(alpha, Lshare,colour=subjectCode)) + 
    geom_point() + theme_classic() +
    geom_smooth(method='lm',se=FALSE)
  
  p2 <- ggplot(DF, aes(CV, Lshare,colour=subjectCode)) + 
    geom_point() + theme_classic() +
    geom_smooth(method='lm',se=FALSE) + theme(legend.position="none")
  
  p3 <- ggplot(DF, aes(Ratio.MM, Lshare,colour=subjectCode)) + 
    geom_point() + theme_classic() +
    geom_smooth(method='lm',se=FALSE) + theme(legend.position="none")
  
  p4 <- ggplot(DF, aes(Fit.RatioMM.CV, Lshare,colour=subjectCode)) + 
    geom_point() + theme_classic() + 
    geom_smooth(method='lm',se=FALSE) + theme(legend.position="none")
  
  legend <- get_legend(p1)
  
  p1 <- p1 + theme(legend.position="none")
  
  tiff(paste(name, ".tiff", sep=""), width = 20, height = 18, units = 'cm', res = 300)
  grid.arrange(p1, p2, legend, p3, p4 , ncol=3, nrow=2, widths=c(2, 2, 0.5))
  dev.off()
  
}

write_Tiff_Pooled <- function(DF,name){
  p1 <- ggplot(DF, aes(alpha, Lshare)) + 
    geom_point() + theme_classic() +
    geom_smooth(method='lm')
  
  p2 <- ggplot(DF, aes(CV, Lshare)) + 
    geom_point() + theme_classic() +
    geom_smooth(method='lm') 
  
  p3 <- ggplot(DF, aes(Ratio.MM, Lshare)) + 
    geom_point() + theme_classic() +
    geom_smooth(method='lm')
  
  p4 <- ggplot(DF, aes(Fit.RatioMM.CV, Lshare)) + 
    geom_point() + theme_classic() + 
    geom_smooth(method='lm')

  
  tiff(paste(name, ".tiff", sep=""), width = 20, height = 20, units = 'cm', res = 300)
  grid.arrange(p1, p2, p3, p4 , ncol=2, nrow=2, widths=c(2, 2))
  dev.off()
  
}

