########### Testing RWN Models #################
########################### Initialization ##########################
rm(list=ls(all=TRUE))

x = readline()
D:\CAPTURE Project\GUI\Data
setwd(x)
getwd()

library(pacman)
p_load(ez,
       lme4,
       ggplot2,
       tidyr,
       plyr,
       dplyr,
       effects,
       rstatix,
       sjPlot,
       flexplot,
       MuMIn)
# devtools::install_github("dustinfife/flexplot", ref="development")

RD = "D:\\CAPTURE Project\\GUI\\Data\\"
WD = "D:\\CAPTURE Project\\GUI\\Data\\"
###################################################### Load data -----
datSeg = read.csv(paste(RD,"MultSegData.csv",sep = ""))
datAH = datSeg[datSeg$RegionLabel == "AntHipp",]
datAH = datAH %>% group_by(Pt,Chan,Walk) %>%
  mutate(Epoch = seq(1,n(),1)) %>%
  as.data.frame()

# A = datAH[,c("Pt", "Chan", "Walk","Epoch")]
names(datAH)
predictorVars = c("Vel", "Fix", "Amb", "KDE", "EDA")
idCols = c("Pt", "Chan", "Walk","Epoch")
outcomeVars = c("wvTheta", "wvAlpha", "wvBeta", "wvGamma", "wvHG", 
                "mtTheta", "mtAlpha", "mtBeta", "mtGamma", "mtHG")

datAH = as.data.frame(lapply(datAH, FUN = function(x) {x[is.nan(x)]=NA;return(x)}))
######################## Data Uniqueness ----------
datModel = datAH
# datModel = as.data.frame(lapply(datModel, FUN = function(x) {x[is.nan(x)]=NA;return(x)}))
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

######################## Plot tests for Random Effects ----------
datModel = datAH
datModel$DV = datModel$wvTheta
datModel$DV =  as.numeric(scale(datModel$DV,center = T))
datModel$Vel =  as.numeric(scale(datModel$Vel,center = T))
datModel$Fix =  as.numeric(scale(datModel$Fix,center = T))
datModel$EDA =  as.numeric(scale(datModel$EDA,center = T))
datModel$KDE =  as.numeric(scale(datModel$KDE,center = T))
datModel$Amb =  as.numeric(scale(datModel$Amb,center = T))

datFit = datModel[complete.cases(datModel[,c("Pt","DV","Vel", "Fix", "Amb", "KDE", "EDA")]),]


datModel$IV = datModel$KDE
names(datModel)


m1 = lme4::lmer(DV ~ KDE+ 
                  (KDE|uPtChan), data=datFit, na.action=na.exclude)


fixed_effects <- fixef(m1)
random_effects <- ranef(m1)$uPtChan

# Create a data frame for regression lines
line_data <- datFit %>%
  dplyr::select(uPtChan, Walk, KDE) %>%
  unique() %>%
  dplyr::mutate(
    Intercept = fixed_effects["(Intercept)"] + random_effects[as.character(uPtChan), "(Intercept)"],
    Slope = fixed_effects["KDE"] + random_effects[as.character(uPtChan), "KDE"],
    Fitted = Intercept + Slope * Walk
  )

datFit$uPtChan = as.factor(datFit$uPtChan)
datFit$Walk = as.factor(datFit$KDE)
line_data$uPtChan = as.factor(line_data$uPtChan)
line_data$Walk = as.factor(line_data$KDE)
# Plot actual data with regression lines for each group
ggplot(datFit, aes(x = KDE, y = DV,colour = uPtChan)) +
  geom_point(alpha = 0.6) +  # Scatterplot of data
  geom_abline(data = line_data,
              aes(intercept = Intercept, slope = Slope, colour = uPtChan),
              linewidth = 1, show.legend = FALSE) +
  geom_abline(data = data.frame(intercept = fixed_effects["(Intercept)"],KDE = fixed_effects["KDE"]),
              aes(intercept = intercept, slope = KDE),colour = "black",
              linewidth = 2)

######################## Progress Report ----------
datModel = datAH
datModel$DV = datModel$wvTheta
datModel$DV =  as.numeric(scale(datModel$DV,center = T))
datModel$Vel =  as.numeric(scale(datModel$Vel,center = T))
datModel$Fix =  as.numeric(scale(datModel$Fix,center = T))
datModel$EDA =  as.numeric(scale(datModel$EDA,center = T))
datModel$KDE =  as.numeric(scale(datModel$KDE,center = T))
datModel$Amb =  as.numeric(scale(datModel$Amb,center = T))

datFit = datModel[complete.cases(datModel[,c("Pt","DV","Vel", "Fix", "Amb", "KDE", "EDA")]),]
predictorVars = c("Vel", "Fix", "Amb", "KDE", "EDA")
outcomeVars = c("wvTheta", "wvAlpha", "wvBeta", "wvGamma", "wvHG")
Models = NULL




Models[["NULL"]] = m0
statDat = NULL
for (dvName in outcomeVars){
  datFit$DV = datFit[[dvName]]
  m0 = lme4::lmer(DV ~1+
                    (1|Pt/Chan)+(1|Walk), data=datFit, na.action=na.exclude)
  for (ivName in predictorVars) {
    print(paste(dvName,ivName,sep = "_"))
    datFit$IV = datFit[[ivName]]
    m = lme4::lmer(DV ~ IV+
                     (1|Pt/Chan)+(IV|Walk), data=datFit, na.action=na.exclude)
    Models[[paste(dvName,ivName,sep = "_")]] = m
    A = model.comparison(Models[["NULL"]], Models[[ivName]])
    temp = data.frame(DV = dvName,IV = ivName)
    temp$R_Squared_Change = A$r_squared_change[4]
    A = anova(Models[["NULL"]], Models[[ivName]])
    temp$Chisq = A$Chisq[2]
    temp$DF = A$Df[2]
    temp$p = A$`Pr(>Chisq)`[2]
    statDat = rbind(statDat,temp)
  }
}
statDat$Sig = ifelse(statDat$p<0.05,"*","")

tab_model(m0,show.se = T,show.stat = T,show.est = T, df.method = "satterthwaite", show.df = T)


m7 = lmer(DV ~ 
            Vel+ 
            Fix+ 
            Amb+ 
            KDE+ 
            EDA+
            (1|Pt/Chan)+(1|Walk), data=datModel, na.action=na.exclude)
summary(m7)
fixef(m7)
ranef(m7)
# This is for Consistency across walks
m8 = lmer(DV ~ 
            Vel+ 
            Fix+ 
            Amb+ 
            KDE+ 
            EDA+
            (1|Pt:Chan:Walk)+(Vel+Fix+Amb+KDE+EDA|Walk), data=datModel, na.action=na.exclude)



anova(m0,m1,m2,m3,m4,m5)


######################## Mixed Level Model fit ----------
datModel = datAH
datModel$DV =  as.numeric(scale(datModel$VT,center = T))
# datModel$DV =  datModel$mtTheta


# lme(ACC ~ ageBase+ageChange+Condition+Vicinity+
#       ageBase:Vicinity+
#       ageChange:Vicinity+Condition:Vicinity,
#     random = list(~1+ageChange|SID,
#                   ~1|Order), data=datModel)


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
m8 = lmer(DV ~ 
            Vel+ 
            Fix+ 
            Amb+ 
            KDE+ 
            EDA+
            (1|Pt:Chan:Walk)+(Vel+Fix+Amb+KDE+EDA|Pt), data=datModel, na.action=na.exclude)

m9 = lme(DV ~ 
           Vel+ 
           Fix+ 
           Amb+ 
           KDE+ 
           EDA+,
          random = list(~1|Pt:Chan:Walk,~Vel+Fix+Amb+KDE+EDA|Pt), data=datModel)

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
