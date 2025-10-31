########################### Initialization ##########################
rm(list=ls(all=TRUE))

x = readline()
D:\CAPTURE Project\LenskyData\First Fixation Data
setwd(x)
getwd()

# install.packages('pacman')
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
       jtools)


RD = r"(D:\CAPTURE Project\LenskyData\First Fixation Data)"
WD = r"(D:\CAPTURE Project\LenskyData\First Fixation Data)"
################################# 
fileName = "MergedData_NoPartialSaccades.csv"
varNames <- c("pID","wID","AOI","Response","AOIviewed","TimeToFirstFixation","RecTask")      
dat = read.csv(file.path(RD,fileName))
names(dat)
dat = dat[,varNames]
dat$TimeToFirstFixation = dat$TimeToFirstFixation/1000 # convert to seconds
dat = dat[!is.na(dat$TimeToFirstFixation),]
dat$TimetoFirstFrame = dat$TimeToFirstFixation
dat$Event = "LMFirstFixation"
dat$Memory = ifelse(dat$Response>=4,"Remembered","Forgotten")
dat$Confidence = paste(dat$Response,sep = "")
dat$Threshold = case_when(dat$Response>=5~"more than 5",
                           dat$Response<=2~"Less than 2")
dat$RecTask = paste(dat$RecTask," Rec Task",sep = "")

dat$Description = paste(dat$Memory,dat$Confidence,dat$Threshold,dat$RecTask,dat$AOI,sep = "];[")
dat$Description = paste("[",dat$Description,"]",sep = "")
varNames <- c("pID","wID","TimetoFirstFrame","Event","Description","Memory","Confidence","Threshold","RecTask","AOI") 
dat = dat[,varNames] 
write.csv(dat,"TimetoFirstFixation.csv",row.names = F)
############################### test ----
fRate = 29.9664
pID = 1
wID = 6
AOI = "LM001"
dP = dat[dat$pID == pID & dat$wID == wID & dat$AOI == AOI,]
dP$TimeToFirstFixation
dP$TimeToFirstFixation/fRate