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

#--------> TODO Compute the number of times each AOI were observed across walks within each participants
#--------> TODO Compute the proportional size of the AOI relative to screen

dat = read.csv(file="MergedData_NoPartialSaccades.csv")

############################### Graphs -----
names(dat)
#Check to see if we have any NA or missing
graphDat = dat[!is.na(dat$Response),]
nrow(dat) - nrow(graphDat)

##############################  FixationCount --------
graphDat = dat
graphDat$ACC = graphDat$Acc1_3vs4_6

# TotalGazeDuration
# SaccadeCount
# FixationCount
# A = graphDat[graphDat$FixationCount>50,]

statData = dat
statData = statData %>% group_by(pID, AOI,Acc1_3vs4_6) %>%
                summarise(FixationCount = mean(FixationCount,na.rm=T)) %>% as.data.frame()

statData = statData %>% group_by(pID,Acc1_3vs4_6) %>%
  summarise(FixationCount = mean(FixationCount,na.rm=T)) %>% as.data.frame()

t.test(FixationCount~Acc1_3vs4_6, data = statData)

ggplot(statData,aes(x=Acc1_3vs4_6,y=FixationCount,fill=Acc1_3vs4_6)) + 
  geom_bar(stat="summary",fun="mean",position="dodge")+
  geom_point(position = position_jitterdodge(jitter.width = .4,
                                             jitter.height = 0,
                                             dodge.width = .9),shape = 21)+
  stat_summary(fun.data = "mean_se", geom="errorbar",position="dodge",color = "Black")+
  theme_bw(base_family = "serif")+
  stat_compare_means(method = "t.test",label = "p.format", paired = F)



statData = dat
ggplot(statData,aes(x=Acc1_3vs4_6, y=FixationCount, fill=Acc1_3vs4_6)) + 
  geom_bar(stat="summary",fun="mean",position="dodge")+
  # geom_point(position = position_jitterdodge(jitter.width = .4,
  #                                            jitter.height = 0,
  #                                            dodge.width = .9),shape = 21)+
  stat_summary(fun.data = "mean_se", geom="errorbar",position="dodge",color = "Black")+
  theme_bw(base_family = "serif")+
  stat_compare_means(method = "t.test",label = "p.format", paired = F,label.y = 0)


statData = dat
statData = statData[statData$AOIviewed==1,]
statData = statData %>% group_by(pID, AOI,Acc1_3vs4_6) %>%
  summarise(FixationCount = mean(FixationCount,na.rm=T)) %>% as.data.frame()

ggplot(statData,aes(x=Acc1_3vs4_6,y=FixationCount,fill=Acc1_3vs4_6)) + 
  geom_bar(stat="summary",fun="mean",position="dodge")+
  # geom_point(position = position_jitterdodge(jitter.width = .4,
  #                                            jitter.height = 0,
  #                                            dodge.width = .9),shape = 21)+
  stat_summary(fun.data = "mean_se", geom="errorbar",position="dodge",color = "Black")+
  theme_bw(base_family = "serif")+
  facet_wrap(~pID)+
  stat_compare_means(method = "t.test",label = "p.format", paired = F, label.y = 5)


##############################  TotalGazeDuration --------
statData = dat
statData = statData %>% group_by(pID, AOI,Acc1_3vs4_6) %>%
  summarise(TotalGazeDuration = mean(TotalGazeDuration,na.rm=T)) %>% as.data.frame()

statData = statData %>% group_by(pID,Acc1_3vs4_6) %>%
  summarise(TotalGazeDuration = mean(TotalGazeDuration,na.rm=T)) %>% as.data.frame()

t.test(TotalGazeDuration~Acc1_3vs4_6, data = statData)

ggplot(statData,aes(x=Acc1_3vs4_6,y=TotalGazeDuration,fill=Acc1_3vs4_6)) + 
  geom_bar(stat="summary",fun="mean",position="dodge")+
  geom_point(position = position_jitterdodge(jitter.width = .4,
                                             jitter.height = 0,
                                             dodge.width = .9),shape = 21)+
  stat_summary(fun.data = "mean_se", geom="errorbar",position="dodge",color = "Black")+
  theme_bw(base_family = "serif")+
  stat_compare_means(method = "t.test",label = "p.format", paired = F)


statData = dat
ggplot(statData,aes(x=Acc1_3vs4_6,y=TotalGazeDuration,fill=Acc1_3vs4_6)) + 
  geom_bar(stat="summary",fun="mean",position="dodge")+
  # geom_point(position = position_jitterdodge(jitter.width = .4,
  #                                            jitter.height = 0,
  #                                            dodge.width = .9),shape = 21)+
  stat_summary(fun.data = "mean_se", geom="errorbar",position="dodge",color = "Black")+
  theme_bw(base_family = "serif")+
  stat_compare_means(method = "t.test",label = "p.format", paired = F,label.y = 0)



statData = dat
statData = statData[statData$AOIviewed==1,]
statData = statData %>% group_by(pID, AOI,Acc1_3vs4_6) %>%
  summarise(TotalGazeDuration = mean(TotalGazeDuration,na.rm=T)) %>% as.data.frame()

ggplot(statData,aes(x=Acc1_3vs4_6,y=TotalGazeDuration,fill=Acc1_3vs4_6)) + 
  geom_bar(stat="summary",fun="mean",position="dodge")+
  # geom_point(position = position_jitterdodge(jitter.width = .4,
  #                                            jitter.height = 0,
  #                                            dodge.width = .9),shape = 21)+
  stat_summary(fun.data = "mean_se", geom="errorbar",position="dodge",color = "Black")+
  theme_bw(base_family = "serif")+
  facet_wrap(~pID)+
  stat_compare_means(method = "t.test",label = "p.format", paired = F, label.y = 5)

##############################  Saccade Count --------
statData = dat
statData = statData %>% group_by(pID, AOI,Acc1_3vs4_6) %>%
  summarise(SaccadeCount = mean(SaccadeCount,na.rm=T)) %>% as.data.frame()

statData = statData %>% group_by(pID,Acc1_3vs4_6) %>%
  summarise(SaccadeCount = mean(SaccadeCount,na.rm=T)) %>% as.data.frame()

t.test(SaccadeCount~Acc1_3vs4_6, data = statData)

ggplot(statData,aes(x=Acc1_3vs4_6,y=SaccadeCount,fill=Acc1_3vs4_6)) + 
  geom_bar(stat="summary",fun="mean",position="dodge")+
  geom_point(position = position_jitterdodge(jitter.width = .4,
                                             jitter.height = 0,
                                             dodge.width = .9),shape = 21)+
  stat_summary(fun.data = "mean_se", geom="errorbar",position="dodge",color = "Black")+
  theme_bw(base_family = "serif")+
  stat_compare_means(method = "t.test",label = "p.format", paired = F)


statData = dat
ggplot(statData,aes(x=Acc1_3vs4_6,y=SaccadeCount,fill=Acc1_3vs4_6)) + 
  geom_bar(stat="summary",fun="mean",position="dodge")+
  # geom_point(position = position_jitterdodge(jitter.width = .4,
  #                                            jitter.height = 0,
  #                                            dodge.width = .9),shape = 21)+
  stat_summary(fun.data = "mean_se", geom="errorbar",position="dodge",color = "Black")+
  theme_bw(base_family = "serif")+
  stat_compare_means(method = "t.test",label = "p.format", paired = F,label.y = 0)


statData = dat
statData = statData[statData$AOIviewed==1,]
statData = statData %>% group_by(pID, AOI,Acc1_3vs4_6) %>%
  summarise(SaccadeCount = mean(SaccadeCount,na.rm=T)) %>% as.data.frame()

ggplot(statData,aes(x=Acc1_3vs4_6,y=SaccadeCount,fill=Acc1_3vs4_6)) + 
  geom_bar(stat="summary",fun="mean",position="dodge")+
  # geom_point(position = position_jitterdodge(jitter.width = .4,
  #                                            jitter.height = 0,
  #                                            dodge.width = .9),shape = 21)+
  stat_summary(fun.data = "mean_se", geom="errorbar",position="dodge",color = "Black")+
  theme_bw(base_family = "serif")+
  facet_wrap(~pID)+
  stat_compare_means(method = "t.test",label = "p.format", paired = F, label.y = 5)


##############################  Average Gaze Duration --------
statData = dat
statData = statData %>% group_by(pID, AOI,Acc1_3vs4_6) %>%
  summarise(AverageFixationDuration = mean(AverageFixationDuration,na.rm=T)) %>% as.data.frame()

statData = statData %>% group_by(pID,Acc1_3vs4_6) %>%
  summarise(AverageFixationDuration = mean(AverageFixationDuration,na.rm=T)) %>% as.data.frame()

t.test(AverageFixationDuration~Acc1_3vs4_6, data = statData)
ggplot(statData,aes(x=Acc1_3vs4_6,y=AverageFixationDuration,fill=Acc1_3vs4_6)) + 
  geom_bar(stat="summary",fun="mean",position="dodge")+
  geom_point(position = position_jitterdodge(jitter.width = .4,
                                             jitter.height = 0,
                                             dodge.width = .9),shape = 21)+
  stat_summary(fun.data = "mean_se", geom="errorbar",position="dodge",color = "Black")+
  theme_bw(base_family = "serif")+
  stat_compare_means(method = "t.test",label = "p.format", paired = F)


statData = dat
ggplot(statData,aes(x=Acc1_3vs4_6,y=AverageFixationDuration,fill=Acc1_3vs4_6)) + 
  geom_bar(stat="summary",fun="mean",position="dodge")+
  # geom_point(position = position_jitterdodge(jitter.width = .4,
  #                                            jitter.height = 0,
  #                                            dodge.width = .9),shape = 21)+
  stat_summary(fun.data = "mean_se", geom="errorbar",position="dodge",color = "Black")+
  theme_bw(base_family = "serif")+
  stat_compare_means(method = "t.test",label = "p.format", paired = F,label.y = 0)


statData = dat
statData = statData[statData$AOIviewed==1,]
statData = statData %>% group_by(pID, AOI,Acc1_3vs4_6) %>%
  summarise(AverageFixationDuration = mean(AverageFixationDuration,na.rm=T)) %>% as.data.frame()

ggplot(statData,aes(x=Acc1_3vs4_6,y=AverageFixationDuration,fill=Acc1_3vs4_6)) + 
  geom_bar(stat="summary",fun="mean",position="dodge")+
  # geom_point(position = position_jitterdodge(jitter.width = .4,
  #                                            jitter.height = 0,
  #                                            dodge.width = .9),shape = 21)+
  stat_summary(fun.data = "mean_se", geom="errorbar",position="dodge",color = "Black")+
  theme_bw(base_family = "serif")+
  facet_wrap(~pID)+
  stat_compare_means(method = "t.test",label = "p.format", paired = F, label.y = 5)


##############################  FixationCount --------
statData = dat
statData = statData %>% group_by(pID, AOI, RecTask, Acc1_3vs4_6) %>%
  summarise(FixationCount = mean(FixationCount,na.rm=T)) %>% as.data.frame()

statData = statData %>% group_by(pID, RecTask, Acc1_3vs4_6) %>%
  summarise(FixationCount = mean(FixationCount,na.rm=T)) %>% as.data.frame()

# t.test(FixationCount~Acc1_3vs4_6, data = statData)
statData$RecTask = factor(statData$RecTask, levels = c("Before","After"))

ggplot(statData,aes(x=RecTask,y=FixationCount,fill=Acc1_3vs4_6)) + 
  geom_bar(stat="summary",fun="mean",position="dodge")+
  geom_point(position = position_jitterdodge(jitter.width = .4,
                                             jitter.height = 0,
                                             dodge.width = .9),shape = 21)+
  stat_summary(fun.data = "mean_se", geom="errorbar",position="dodge",color = "Black")+
  theme_bw(base_family = "serif")+
  stat_compare_means(method = "t.test",label = "p.format", paired = F, label.y = 0)


statData = statData %>% group_by(pID, RecTask) %>%
  summarise(FixationCount = mean(FixationCount,na.rm=T)) %>% as.data.frame()

statData$RecTask = factor(statData$RecTask, levels = c("Before","After"))

ggplot(statData,aes(x=RecTask,y=FixationCount,fill=RecTask)) + 
  geom_bar(stat="summary",fun="mean",position="dodge")+
  geom_point(position = position_jitterdodge(jitter.width = .4,
                                             jitter.height = 0,
                                             dodge.width = .9),shape = 21)+
  stat_summary(fun.data = "mean_se", geom="errorbar",position="dodge",color = "Black")+
  theme_bw(base_family = "serif")+
  facet_wrap(~RecTask)+
  stat_compare_means(method = "t.test",label = "p.format", paired = F, label.y = 0)

##############################  FixationCount --------
# TotalGazeDuration
# SaccadeCount
# FixationCount
# A = graphDat[graphDat$FixationCount>50,]

statData = dat
statData$ACC = statData$Acc1_3vs4_6
statData$IV = statData$TotalGazeDuration
statData = statData %>% group_by(pID, AOI,ACC) %>%
  summarise(IV = mean(IV,na.rm=T)) %>% as.data.frame()

statData = statData %>% group_by(pID,ACC) %>%
  summarise(IV = mean(IV,na.rm=T)) %>% as.data.frame()

t.test(IV~ACC, data = statData)

ggplot(statData,aes(x=ACC,y=IV,fill=ACC)) + 
  geom_bar(stat="summary",fun="mean",position="dodge")+
  geom_line(aes(group = pID))+
  geom_point()+
  # geom_point(position = position_jitterdodge(jitter.width = .4,
  #                                            jitter.height = 0,
  #                                            dodge.width = .9),shape = 21)+
  stat_summary(fun.data = "mean_se", geom="errorbar",position="dodge",color = "Black")+
  theme_bw(base_family = "serif")+
  stat_compare_means(method = "t.test",label = "p.format", paired = F)



statData = dat
statData$ACC = statData$Acc1_3vs4_6
statData$IV = statData$FixationCount
ggplot(statData,aes(x=ACC, y=IV, fill=ACC)) + 
  geom_bar(stat="summary",fun="mean",position="dodge")+
  geom_point(position = position_jitterdodge(jitter.width = .4,
                                             jitter.height = 0,
                                             dodge.width = .9),shape = 21)+
  stat_summary(fun.data = "mean_se", geom="errorbar",position="dodge",color = "Black")+
  theme_bw(base_family = "serif")+
  stat_compare_means(method = "t.test",label = "p.format", paired = F,label.y = 0)


statData = dat
statData$ACC = statData$Acc1_3vs4_6
statData$IV = statData$FixationCount
statData = statData[statData$AOIviewed==1,]
statData = statData %>% group_by(pID, AOI,Acc1_3vs4_6) %>%
  summarise(FixationCount = mean(FixationCount,na.rm=T)) %>% as.data.frame()

ggplot(statData,aes(x=Acc1_3vs4_6,y=FixationCount,fill=Acc1_3vs4_6)) + 
  geom_bar(stat="summary",fun="mean",position="dodge")+
  # geom_point(position = position_jitterdodge(jitter.width = .4,
  #                                            jitter.height = 0,
  #                                            dodge.width = .9),shape = 21)+
  stat_summary(fun.data = "mean_se", geom="errorbar",position="dodge",color = "Black")+
  theme_bw(base_family = "serif")+
  facet_wrap(~pID)+
  stat_compare_means(method = "t.test",label = "p.format", paired = F, label.y = 5)


###################### Old backup-------------
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
