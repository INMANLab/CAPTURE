########### Testing RWN Models #################
########################### Initialization ##########################
rm(list=ls(all=TRUE))

x = readline()
D:\CAPTURE Project\GUI\Data
setwd(x)
getwd()

library(pacman)
p_load(data.table,
       reshape2,
       ez,
       lme4,
       lmerTest,
       ggplot2,
       grid,
       tidyr,
       plyr,
       dplyr,
       effects,
       gridExtra,
       DescTools,
       Cairo, #alternate image writing package with superior performance.
       corrplot,
       knitr,
       PerformanceAnalytics,
       afex,
       ggpubr,
       readxl,
       officer,
       psych,
       rstatix,
       emmeans,
       ggformula,
       export,
       scales,
       correlation,
       sjPlot,
       interactions,
       jtools,
       MuMIn)


RD = "D:\\CAPTURE Project\\GUI\\Data\\"
WD = "D:\\CAPTURE Project\\GUI\\Data\\"
###################################################### Load data -----
datSeg = read.csv(paste(RD,"MultSegData.csv",sep = ""))
datAH = datSeg[datSeg$RegionLabel == "AntHipp",]
datAH = datAH %>% group_by(Pt,Chan,Walk) %>%
  mutate(Epoch = seq(1,n(),1)) %>%
  as.data.frame()

A = datAH[,c("Pt", "Chan", "Walk","Epoch")]
names(datAH)
predictorVars = c("Vel", "Fix", "Amb", "KDE", "EDA")
idCols = c("Pt", "Chan", "Walk","Epoch")
outcomeVars = c("wvTheta", "wvAlpha", "wvBeta", "wvGamma", "wvHG", 
                "mtTheta", "mtAlpha", "mtBeta", "mtGamma", "mtHG")


######################## Data Uniqueness ----------
datModel = datAH
datModel$DV =  datModel$wvTheta
summaryDat = datModel %>% group_by(Pt,Chan,Walk) %>%
  summarise(N_DV = n(), NUnique_DV = length(unique(DV)),
            N_Vel = sum(!is.na(Vel)), NUnique_Vel = length(unique(Vel)),
            N_Fix = sum(!is.na(Fix)), NUnique_Fix = length(unique(Fix)),
            N_Amb = sum(!is.na(Amb)), NUnique_Amb = length(unique(Amb)),
            N_KDE = sum(!is.na(KDE)), NUnique_KDE = length(unique(KDE)),
            N_EDA = sum(!is.na(EDA)), NUnique_EDA = length(unique(EDA))) %>%
  as.data.frame()

  #   group_by(SID,sessID,testName,timeStamp) %>%
  #   summarize(CorrL = mean(CorrL,na.rm = T),
  #             diff = mean(diff,na.rm = T),
  #             IncorrL = mean(IncorrL,na.rm = T))

######################## Mixed Level Model fit ----------
datModel = datAH
datModel$DV =  datModel$wvTheta
# datModel$DV =  datModel$mtTheta

mNULL1 = lmer(DV ~ 1+(1|Pt), data=datModel, na.action=na.exclude)
mNULL2 = lmer(DV ~ 1+(1|Pt/Chan), data=datModel, na.action=na.exclude)
mNULL3 = lmer(DV ~ 1+(1|Pt/Chan/Walk), data=datModel, na.action=na.exclude)
# mNULL4 = lmer(DV ~ 1+(1|Pt/Walk)+(1|Chan), data=datModel, na.action=na.exclude)
# mNULL5 = lmer(DV ~ 1+(1|Pt/Walk)+(1|Pt/Chan), data=datModel, na.action=na.exclude)
# anova(mNULL1,mNULL2,mNULL3,mNULL4,mNULL5)
anova(mNULL1,mNULL2,mNULL3)

M = mNULL3
summary(M)
tab_model(M,show.se = T,show.stat = T,show.est = T)
performance::performance(M)
hist(resid(M))
VIF(M)


m0 = lmer(DV ~ 1+(1|Pt/Chan/Walk), data=datModel, na.action=na.exclude)
m1 = lmer(DV ~ 
            Vel+ 
            (1|Pt/Chan/Walk), data=datModel, na.action=na.exclude)
m2 = lmer(DV ~ 
            Vel+ 
            Fix+ 
            (1|Pt/Chan/Walk), data=datModel, na.action=na.exclude)
m3 = lmer(DV ~ 
            Vel+ 
            Fix+ 
            Amb+ 
            (1|Pt/Chan/Walk), data=datModel, na.action=na.exclude)
m4 = lmer(DV ~ 
            Vel+ 
            Fix+ 
            Amb+ 
            KDE+ 
            (1|Pt/Chan/Walk), data=datModel, na.action=na.exclude)
m5 = lmer(DV ~ 
            Vel+ 
            Fix+ 
            Amb+ 
            KDE+ 
            EDA+
            (1|Pt/Chan/Walk), data=datModel, na.action=na.exclude)
m6 = lmer(DV ~ 
            Vel+ 
            Fix+ 
            Amb+ 
            KDE+ 
            EDA+
            (1|Pt/Chan/Walk)+(Vel+Fix+Amb+KDE+EDA|Walk), data=datModel, na.action=na.exclude)
m7 = lmer(DV ~ 
            Vel+ 
            Fix+ 
            Amb+ 
            KDE+ 
            EDA+
            (1|Pt:Chan)+(1|Pt:Chan:Walk)+(Vel+Fix+Amb+KDE+EDA|Pt), data=datModel, na.action=na.exclude)
anova(m0,m1,m2,m3,m4,m5)

M = m7
summary(M)
tab_model(M,show.se = T,show.stat = T,show.est = T)
performance::performance(M)
hist(resid(M))
VIF(M)
######################## Mixed Level Model fit Separate Walks iterative ----------
m = NULL
resTab = NULL
for (walkID in seq(1,8)){
  print(walkID)
  datModel = datAH[datAH$Walk == walkID,]
  datModel$DV =  datModel$wvTheta
  
  mTemp = lmer(DV ~
                 Vel+ 
                 Fix+ 
                 Amb+ 
                 KDE+ 
                 EDA+
                 (1|Pt/Chan), data=datModel, na.action=na.exclude)
  m[[walkID]] = mTemp
  s = summary(m[[walkID]])
  s = s$coefficients
  s = data.frame(s)
  s$R2m = r.squaredGLMM(m[[walkID]])[1]
  s$R2c = r.squaredGLMM(m[[walkID]])[2]
  s$walkID = walkID 
  s$Sig = " "
  names(s) = c("Estimate", "StdErr", "df", "t","p", "R2m","R2c", "walkID", "Sig")
  s$Sig[s$p<0.05] = "*"
  resTab = rbind(resTab,s)
}
# summaryDat = datModel %>% group_by(Pt,Chan,Walk) %>%
#   summarise(N_DV = n(), NUnique_DV = length(unique(DV)),
#             N_Vel = sum(!is.na(Vel)), NUnique_Vel = length(unique(Vel)),
#             N_Fix = sum(!is.na(Fix)), NUnique_Fix = length(unique(Fix)),
#             N_Amb = sum(!is.na(Amb)), NUnique_Amb = length(unique(Amb)),
#             N_KDE = sum(!is.na(KDE)), NUnique_KDE = length(unique(KDE)),
#             N_EDA = sum(!is.na(EDA)), NUnique_EDA = length(unique(EDA))) %>%
#   as.data.frame()
M = m[[1]]
A = tab_model(M,show.se = T,show.stat = T,show.est = T)
summary(M)
performance::performance(M)
hist(resid(M))
VIF(M)

######################## Mixed Level Model fit Separate Walks Averaged Over Channels iterative ----------
m = NULL
resTab = NULL
for (walkID in seq(1,8)){
  print(walkID)
  datModel = datAH[datAH$Walk == walkID,]
  datModel$DV =  datModel$wvTheta
  datModel = datModel %>% group_by(Pt,Epoch) %>%
    summarise(DV = mean(DV,na.rm=T),
              Vel = mean(DV,na.rm=T),
              Fix = mean(Fix,na.rm=T),
              Amb = mean(Amb,na.rm=T),
              KDE = mean(KDE,na.rm=T),
              EDA = mean(EDA,na.rm=T)) %>%
    as.data.frame()
  
  mTemp = lmer(DV ~
                 Vel+ 
                 Fix+ 
                 Amb+ 
                 KDE+ 
                 EDA+
                 (1|Pt), data=datModel, na.action=na.exclude)
  m[[walkID]] = mTemp
  s = summary(m[[walkID]])
  s = s$coefficients
  s = data.frame(s)
  s$R2m = r.squaredGLMM(m[[walkID]])[1]
  s$R2c = r.squaredGLMM(m[[walkID]])[2]
  s$walkID = walkID 
  s$Sig = " "
  names(s) = c("Estimate", "StdErr", "df", "t","p", "R2m","R2c", "walkID", "Sig")
  s$Sig[s$p<0.05] = "*"
  resTab = rbind(resTab,s)
}
# summaryDat = datModel %>% group_by(Pt,Chan,Walk) %>%
#   summarise(N_DV = n(), NUnique_DV = length(unique(DV)),
#             N_Vel = sum(!is.na(Vel)), NUnique_Vel = length(unique(Vel)),
#             N_Fix = sum(!is.na(Fix)), NUnique_Fix = length(unique(Fix)),
#             N_Amb = sum(!is.na(Amb)), NUnique_Amb = length(unique(Amb)),
#             N_KDE = sum(!is.na(KDE)), NUnique_KDE = length(unique(KDE)),
#             N_EDA = sum(!is.na(EDA)), NUnique_EDA = length(unique(EDA))) %>%
#   as.data.frame()
M = m[[3]]
tab_model(M,show.se = T,show.stat = T,show.est = T)
summary(M)
performance::performance(M)
hist(resid(M))
VIF(M)
######################## Mixed Level Model fit Separate Walks test ----------
datModel = datAH[datAH$Walk == 1,]
datModel$DV =  datModel$wvTheta
# datModel$DV =  datModel$mtTheta

mNULL1 = lmer(DV ~ 1+(1|Pt), data=datModel, na.action=na.exclude)
mNULL2 = lmer(DV ~ 1+(1|Pt/Chan), data=datModel, na.action=na.exclude)
anova(mNULL1,mNULL2,mNULL3)

M = mNULL2
summary(M)
tab_model(M,show.se = T,show.stat = T,show.est = T)
performance::performance(M)
hist(resid(M))
VIF(M)


m0 = lmer(DV ~ 1+(1|Pt/Chan), data=datModel, na.action=na.exclude)
m1 = lmer(DV ~ 
            Vel+ 
            (1|Pt/Chan), data=datModel, na.action=na.exclude)
m2 = lmer(DV ~ 
            Vel+ 
            Fix+ 
            (1|Pt/Chan), data=datModel, na.action=na.exclude)
m3 = lmer(DV ~ 
            Vel+ 
            Fix+ 
            Amb+ 
            (1|Pt/Chan), data=datModel, na.action=na.exclude)
m4 = lmer(DV ~ 
            Vel+ 
            Fix+ 
            Amb+ 
            KDE+ 
            (1|Pt/Chan), data=datModel, na.action=na.exclude)
m5 = lmer(DV ~ 
            Vel+ 
            Fix+ 
            Amb+ 
            KDE+ 
            EDA+
            (1|Pt/Chan), data=datModel, na.action=na.exclude)
anova(m0,m1,m2,m3,m4,m5)

M = m1
tab_model(M,show.se = T,show.stat = T,show.est = T)
summary(M)
performance::performance(M)
hist(resid(M))
VIF(M)

######################## Mixed Level Model fit Separate Walks test Averages over Channels----------
datModel = datAH[datAH$Walk == 1,]
datModel$DV =  datModel$wvTheta
datModel = datModel %>% group_by(Pt,Epoch,Walk) %>%
  summarise(DV = mean(DV,na.rm=T),
            Vel = mean(DV,na.rm=T),
            Fix = mean(Fix,na.rm=T),
            Amb = mean(Amb,na.rm=T),
            KDE = mean(KDE,na.rm=T),
            EDA = mean(EDA,na.rm=T)) %>%
  as.data.frame()



m0 = lmer(DV ~ 1+(1|Pt), data=datModel, na.action=na.exclude)
m1 = lmer(DV ~ 
            Vel+ 
            (1|Pt), data=datModel, na.action=na.exclude)
m2 = lmer(DV ~ 
            Vel+ 
            Fix+ 
            (1|Pt), data=datModel, na.action=na.exclude)
m3 = lmer(DV ~ 
            Vel+ 
            Fix+ 
            Amb+ 
            (1|Pt), data=datModel, na.action=na.exclude)
m4 = lmer(DV ~ 
            Vel+ 
            Fix+ 
            Amb+ 
            KDE+ 
            (1|Pt), data=datModel, na.action=na.exclude)
m5 = lmer(DV ~ 
            Vel+ 
            Fix+ 
            Amb+ 
            KDE+ 
            EDA+
            (1|Pt), data=datModel, na.action=na.exclude)
anova(m0,m1,m2,m3,m4,m5)

M = m5
tab_model(M,show.se = T,show.stat = T,show.est = T)
summary(M)
performance::performance(M)
hist(resid(M))
VIF(M)
