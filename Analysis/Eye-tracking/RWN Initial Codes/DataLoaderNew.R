########################### Initialization ##########################
rm(list=ls(all=TRUE))

x = readline()
Z:\Data\RealWorldNavigationCory\Blickshift Analytics Dynamic AOIs\LoadandMergeETandRecogTaskData
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

varnames <- c("AOI","TotalGazeDuration","NormalizedGazeDuration","AverageGazeDuration",
"MaximumGazeDuration","MinimumGazeDuration","GazeCount","TimeToFirstFixation",
"GazedAtBy","AOITransitionRate","FixationCount","SaccadeCount",
"GazePointValidity","AverageFixationDuration","AverageSaccadeDuration","AverageSaccadeLength",
"FixationRate","FixationSaccadeTimeRatio","ScanPathLength","ScanPathDuration",
"ScanPathArea","MeanPupilDilation")      

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
      datTemp = datTemp[,varnames]
      datTemp$pID = pID
      datTemp$wID = wID
      dat = rbind(dat,datTemp)
    }
  }
}
dat$AOI = gsub(" ","",dat$AOI)
dat = dat[!(dat$AOI %in% c("LM098,LM099","LM025,LM023","LM082,LM081",
                           "LM071,LM072","LM031,LM030","LM013n,LM014","LM095,LM098",
                           "LM009,LM008","LM145,LM135","LM140a,LM135","LM137,LM136",
                           "LM109,LM110")),]

unique(dat$AOI)

respDat = NULL
ETfileList = list.files(path = paste(getwd(),"/ResponseData",sep = ""),pattern = ".*xlsx")
for (pID in 1:5){
  datTemp = read_excel(path = paste(getwd(), "/ResponseData/","RW",pID,".xlsx",sep = ""))
  datTemp = datTemp[,c("Response(1-6)","ImageFile")]
  datTemp$ImageFile = gsub(".JPEG","",datTemp$ImageFile)
  datTemp = datTemp[grepl("LM",datTemp$ImageFile),]
  datTemp$pID = pID
  names(datTemp) = c("Response","RecTaskLab","pID")
  respDat = rbind(respDat,datTemp)
}
unique(respDat$RecTaskLab)


labelKey = read_excel("~/Test_Data1/RecognitionTask_LandmarkLabel_KEY.xlsx")
names(labelKey) = c("a","AOI","RecTaskLab")
labelKey$AOI = gsub(".JPEG","",labelKey$AOI)
labelKey$AOI = gsub(".jpg","",labelKey$AOI)
labelKey = labelKey[,c("AOI","RecTaskLab")]

dat = merge(dat,labelKey,all.x = T, by="AOI")
A = dat[is.na(dat$RecTaskLab),]

dat = merge(dat,respDat,all.x = T, by=c("pID","RecTaskLab"))
dat$binaryResp = ifelse(dat$Response>3,"Remembered","NotRemembered")
dat$Hit = ifelse(dat$Response>3,1,0)
dat = dat[,c("pID","wID" , "RecTaskLab", "AOI","Response", "binaryResp", "Hit", names(dat)[ !(names(dat) %in% c("pID","wID" , "RecTaskLab", "AOI","Response", "binaryResp"))])]
write.csv(dat,file="MergedData.csv",row.names = F)

###############################

dat = read.csv(file="MergedData.csv")

############################### Graphs -----
names(dat)
graphDat = dat[!is.na(dat$Response),]
nrow(dat) - nrow(graphDat)

ggplot(graphDat,aes(x=binaryResp,y=FixationCount,fill=binaryResp)) + 
  geom_bar(stat="summary",fun="mean",position="dodge")+
  geom_point(position = position_jitterdodge(jitter.width = .4,
                                             jitter.height = 0,
                                             dodge.width = .9),shape = 21)+
  stat_summary(fun.data = "mean_se", geom="errorbar",position="dodge",color = "Black")+
  theme_bw(base_family = "serif")

ggplot(graphDat,aes(x=binaryResp,y=FixationCount,fill=binaryResp)) + 
  geom_bar(stat="summary",fun="mean",position="dodge")+
  stat_summary(fun.data = "mean_se", geom="errorbar",position="dodge",color = "Black")

graphDat$wID = as.factor(graphDat$wID)
ggplot(graphDat,aes(x=binaryResp,y=FixationCount,fill=wID)) + 
  geom_bar(stat="summary",fun="mean",position="dodge")+
  geom_point(position = position_jitterdodge(jitter.width = .4,
                                           jitter.height = 0,
                                           dodge.width = .9),shape = 21)+
  stat_summary(fun.data = "mean_se", geom="errorbar",position="dodge",color = "Black")
  
graphDat$wID = as.factor(graphDat$wID)
ggplot(graphDat,aes(x=binaryResp,y=FixationCount,fill=wID)) + 
  geom_bar(stat="summary",fun="mean",position="dodge")+
  geom_point(position = position_jitterdodge(jitter.width = .4,
                                             jitter.height = 0,
                                             dodge.width = .9),shape = 21)+
  stat_summary(fun.data = "mean_se", geom="errorbar",position="dodge",color = "Black") +
  facet_wrap(~pID)

ggplot(graphDat,aes(x=binaryResp,y=FixationCount,fill=wID)) + 
  geom_boxplot()+
  facet_wrap(~pID)

datOutlier = dat[dat$FixationCount > 10, ]
View(datOutlier)

docDat = unique(dat[,c("pID","wID" , "RecTaskLab", "AOI","Response", "binaryResp", "Hit","FixationCount")])
View(docDat)

docDat = unique(dat[,c("pID","wID", "AOI","FixationCount")])

docDat = docDat[is.numeric(docDat$FixationCount), ]
docDat = docDat[(docDat$FixationCount == 0) | (docDat$FixationCount >= 1), ]
docDat$FixationCount = ifelse(docDat$FixationCount == 0, 0, 1)
docDat$FixationCount[docDat$pID == 5 & docDat$wID == 8 & docDat$AOI == "LM069"] = NA
docDat = unique(docDat[,c("pID","wID", "AOI","FixationCount")])
docDat = pivot_wider(docDat,id_cols = c("pID","wID"), names_from = "AOI", values_from = "FixationCount")

docDat %>%
  dplyr::summarise(n = dplyr::n(), .by = c(pID, wID, AOI)) %>%
  dplyr::filter(n > 1) %>%
  as.data.frame()
write.csv(docDat,file="DocumentationFile.csv",row.names = F)
