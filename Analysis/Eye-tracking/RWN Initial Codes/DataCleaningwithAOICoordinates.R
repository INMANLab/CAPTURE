########################### Initialization ##########################
rm(list=ls(all=TRUE))

x = readline()
Z:\Data\RealWorldNavigationCory\Blickshift Analytics Dynamic AOIs\Test_Data1
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

################################# 
ETfileList = list.files(pattern = "R.*_W.*csv")

  

fileNum = 1
dat = NULL
for (pID in 1:5){
  for (wID in 1:8){
    ETfileList = list.files(pattern = paste("R",pID,".W",wID,".csv",sep = ""))
    if(length(ETfileList)==1){
      print(paste("loading file:",fileNum," Patient:",pID," Walk:",wID,sep = ""))
      fileNum = fileNum+1
      datTemp = read.csv(ETfileList,sep = ";")
      if (ncol(datTemp)==1){
        datTemp = read.csv(ETfileList,sep = ",")
      }
      datTemp$pID = pID
      datTemp$wID = wID
      dat = rbind(dat,datTemp)
    }
  }
}

varnames <- c("pID","wID","AOI","TotalGazeDuration","NormalizedGazeDuration","AverageGazeDuration",
              "MaximumGazeDuration","MinimumGazeDuration","GazeCount","TimeToFirstFixation",
              "GazedAtBy","AOITransitionRate","FixationCount","SaccadeCount",
              "GazePointValidity","AverageFixationDuration","AverageSaccadeDuration","AverageSaccadeLength",
              "FixationRate","FixationSaccadeTimeRatio","ScanPathLength","ScanPathDuration",
              "ScanPathArea","MeanPupilDilation")    
unique(dat$AOI)

ConcurrentAOIs = dat[grepl(".*#.*",dat$AOI),c("pID","wID","AOI")]
write.csv(ConcurrentAOIs,"ConcurrentAOIs.csv",row.names = F)


dat$AOI = gsub(" ","",dat$AOI)
unique(dat$TotalGazeDuration)


######################################
pDat = dat[dat$pID==4 & dat$wID==1,]

unique(pDat$AOI)
