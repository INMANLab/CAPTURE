########################### Initialization ##########################
rm(list=ls(all=TRUE))

x = readline()
Z:\Data\RealWorldNavigationCory\Blickshift Analytics Dynamic AOIs\Test_Data1
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
       jtools)

#################################
# ETfileList = list.files(pattern = "R.*_W.*csv")

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
dat$AOI = gsub(" ","",dat$AOI)
dat = dat[!(dat$AOI %in% c("N/A","0","LM098,LM099","LM025,LM023","LM082,LM081",
                         "LM071,LM072","LM031,LM030","LM013n,LM014","LM095,LM098",
                         "LM009,LM008","LM145,LM135","LM140a,LM135","LM137,LM136",
                         "LM109,LM110")),]

unique(dat$AOI)

respDat = NULL
# ETfileList = list.files(path = paste(getwd(),"/ResponseData",sep = ""),pattern = ".*xlsx")
for (pID in 1:5){
  datTemp = read_excel(path = paste(getwd(),"/ResponseData/","RW",pID,".xlsx",sep = ""))
  datTemp = datTemp[,c("Response(1-6)","ImageFile")]
  datTemp$ImageFile = gsub(".JPEG","",datTemp$ImageFile)
  datTemp = datTemp[grepl("LM",datTemp$ImageFile),]
  datTemp$pID = pID
  names(datTemp) = c("Response","RecTaskLab","pID")
  respDat = rbind(respDat,datTemp)
}
unique(respDat$RecTaskLab)


labelKey = read_excel("RecognitionTask_LandmarkLabel_KEY.xlsx")
names(labelKey) = c("a","AOI","RecTaskLab")
labelKey$AOI = gsub(".JPEG","",labelKey$AOI)
labelKey$AOI = gsub(".jpg","",labelKey$AOI)
labelKey = labelKey[,c("AOI","RecTaskLab")]

dat = merge(dat,labelKey,all.x = T, by="AOI")
A = dat[is.na(dat$RecTaskLab),]

dat = merge(dat,respDat,all.x = T, by=c("pID","RecTaskLab"))
dat$binaryResp = ifelse(dat$Response>3,"Remembered","NotRemembered")

dat = dat[,c("pID","wID" , "RecTaskLab", "AOI","Response", "binaryResp", names(dat)[ !(names(dat) %in% c("pID","wID" , "RecTaskLab", "AOI","Response", "binaryResp"))])]
write.csv(dat,file="MergedData.csv",row.names = F)

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
  stat_summary(fun.data = "mean_se", geom="errorbar",position="dodge",color = "Black")

ggplot(graphDat[graphDat$wID!=3,],aes(x=binaryResp,y=FixationCount,fill=wID)) + 
  geom_bar(stat="summary",fun="mean",position="dodge")+
  stat_summary(fun.data = "mean_se", geom="errorbar",position="dodge",color = "Black")


ggplot(graphDat[graphDat$wID!=3,],aes(x=binaryResp,y=FixationCount,fill=binaryResp)) + 
  geom_bar(stat="summary",fun="mean",position="dodge")+
  stat_summary(fun.data = "mean_se", geom="errorbar",position="dodge",color = "Black")+
  facet_wrap(~wID)
  
