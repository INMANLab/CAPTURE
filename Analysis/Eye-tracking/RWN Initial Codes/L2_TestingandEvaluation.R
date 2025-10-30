########################### Initialization ##########################
rm(list=ls(all=TRUE))

x = readline()
C:/Users/Lensky Augustin/Downloads/LoadandMergeETandRecogTaskData/LoadandMergeETandRecogTaskData
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
       jtools,
       pROC)

#--------> TODO Compute the number of times each AOI were observed across walks within each participants
#--------> TODO Compute the proportional size of the AOI relative to screen
############################### Load Data -------
dat = read.csv(file="MergedData_NoPartialSaccades.csv")
behavioralDat = read.csv("BehavioralPerformance.csv")
trialLevelPerformanceDat = read.csv("TrialLevelBehaviorData.csv")
############################### Variables -------
gDefaultPars = theme(
  axis.title.x = element_text(size = 15, face = "bold"), 
  axis.title.y = element_text(size = 15, face = "bold"),
  panel.border = element_blank(), 
  panel.grid = element_blank(),
  axis.line = element_line(), 
  axis.text = element_text(colour = "black", size = rel(1.3), face = "bold"),
  legend.title = element_text(colour = "black", size = rel(1), face = "bold"),
  legend.text = element_text(colour = "black", size = rel(1), face = "bold"), 
  legend.position = "top")

############################### Functions ----
BarGraphRemvsForgotten <- function(graphDat, colorCodes = c("seagreen2", "purple3")){ 
  p = ggplot(graphDat,aes(x=X,y=Y,fill=Fill)) +
    geom_bar(stat="summary",fun="mean",position="dodge") +
    geom_point(position = position_jitterdodge(jitter.width = .4,
                                               jitter.height = 0,
                                               dodge.width = .9),
               shape = 21, size = 3, alpha = .8, colour = "black")+
    stat_summary(fun.data = "mean_se", geom="errorbar", position = position_dodge(width = 0.9),
                 width = 0.2, linewidth = 1, color = "Black")+
    theme_bw(base_family = "serif")+
    #stat_compare_means(method = "t.test",label = "p.format", paired = F) +
    scale_fill_manual(values = colorCodes) +
    labs(fill = fillLab) +
    xlab(xLab) +
    ylab(yLab)
  
  p = p + scale_y_continuous(expand = expansion(mult = c(0,.05))) + gDefaultPars
  plot(p)
  return(p)
}
############################### Graphs -----
names(dat)
#Check to see if we have any NA or missing
graphDat = dat[!is.na(dat$Response),]
nrow(dat) - nrow(graphDat)


####################### Multilevel Model Sample Code ############################

datModel = dat
datModel$ACC = ifelse(datModel$Acc1_3vs4_6=="Scored<=3",0,1)
datModel$RecTask = factor(datModel$RecTask,levels=c("Before","After"))


#datModel = RemoveOutliers(datModel,factorNames = c("sessID","testName"),varNames = "value",Criteria=3)
#datModel = datModel[!is.na(datModel$value),]

m = lmer(SaccadeCount ~ RecTask+ACC+
           (1|pID/wID), data=datModel, na.action=na.exclude)
# m = glmer(ACC ~ RecTask+FixationCount+
#                  (1|pID/wID), data=datModel, family = binomial,na.action=na.exclude)
tab_model(m,show.se = T,show.stat = T,show.est = T)
summary(m)
plot(allEffects(m, default.levels=4))

# mLogGP = lme(fixed = ACC ~ testName+ageChange+value+ageBase,
#              random = 1+testName|SID,
#              data=datModel)

m = glmer(ACC ~ RecTask+FixationCount+
            (1|pID:wID), data=datModel, family = binomial,na.action=na.exclude)
tab_model(m,show.se = T,show.stat = T,show.est = T)
summary(m)
plot(allEffects(m, default.levels=4))

####################### Multilevel Model Sample Code (Data before recognition task) ############################

datModel = dat
datModel$ACC = ifelse(datModel$Acc1_3vs4_6=="Scored<=3",0,1)
datModel = datModel[datModel$RecTask=="Before",]

#datModel = RemoveOutliers(datModel,factorNames = c("sessID","testName"),varNames = "value",Criteria=3)
#datModel = datModel[!is.na(datModel$value),]

m = lmer(SaccadeCount ~ ACC+
           (1|pID/wID), data=datModel, na.action=na.exclude)
# m = glmer(ACC ~ RecTask+FixationCount+
#                  (1|pID/wID), data=datModel, family = binomial,na.action=na.exclude)
tab_model(m,show.se = T,show.stat = T,show.est = T)
summary(m)
plot(allEffects(m, default.levels=4))

m = lmer(FixationCount ~ ACC+
           (1|pID/wID), data=datModel, na.action=na.exclude)
# m = glmer(ACC ~ RecTask+FixationCount+
#                  (1|pID/wID), data=datModel, family = binomial,na.action=na.exclude)
tab_model(m,show.se = T,show.stat = T,show.est = T)
summary(m)
plot(allEffects(m, default.levels=4))

m = lmer(TotalGazeDuration ~ ACC+
           (1|pID/wID), data=datModel, na.action=na.exclude)
# m = glmer(ACC ~ RecTask+FixationCount+
#                  (1|pID/wID), data=datModel, family = binomial,na.action=na.exclude)
tab_model(m,show.se = T,show.stat = T,show.est = T)
summary(m)
plot(allEffects(m, default.levels=4))

m = lmer(ACC ~ (1|pID/wID), data=datModel, na.action=na.exclude)
# m = glmer(ACC ~ RecTask+FixationCount+
#                  (1|pID/wID), data=datModel, family = binomial,na.action=na.exclude)
tab_model(m,show.se = T,show.stat = T,show.est = T)
summary(m)
plot(allEffects(m, default.levels=4))
########################### Recognition Task D Prime Score ##################################
statData = dat
statData = statData %>% group_by(pID, AOI,Acc1_3vs4_6,RecTask) %>%
  summarise(ACC = mean(ACC,na.rm=T)) %>% as.data.frame()

statData = statData %>% group_by(pID,Acc1_3vs4_6,RecTask) %>%
  summarise(ACC = mean(ACC,na.rm=T)) %>% as.data.frame()


t.test(TotalGazeDuration~Acc1_3vs4_6, data = statData)
fSize = 16
gDefault =   theme(
  text = element_text(family = "serif"),
  strip.text = element_text(size=fSize-2, face="bold"),
  plot.title = element_text(size = fSize),
  axis.title = element_text(size = fSize),
  legend.title = element_text(size=fSize),
  legend.text=element_text(size=fSize-4),
  axis.text = element_text(size = fSize-4),
  panel.background = element_rect(fill = "transparent", color = "black", 
                                  linewidth = 0.7, linetype =  "solid"),
  plot.background = element_rect(fill = "transparent"),
  axis.line = element_line(linewidth = .7, linetype =  "solid"),
  strip.background = element_rect(fill="transparent"),
  panel.grid.major.x = element_blank(),
  panel.grid.minor = element_blank()
)
ggplot(datModel,aes(x=Acc1_3vs4_6,y=ACC, fill='aes')) + 
  geom_bar(stat="summary",fun="mean",position="dodge")+
  geom_point(position = position_jitterdodge(jitter.width = .4,
                                             jitter.height = 0,
                                             dodge.width = .9),shape = 21)+
  gDefault +
  stat_summary(fun.data = "mean_se", geom="errorbar",position="dodge",color = "Black")+
  theme_bw(base_family = "serif")+
  stat_compare_means(method = "t.test",label = "p.format", paired = F)
graph2ppt(file= "TotalGazeDurationGraphAfter.pptx", width = 8, height = 5)

############################### Updated Graphs (By walk)
graphDat = dat
graphDat$ACC = graphDat$Acc1_3vs4_6

statData = dat
statData = statData %>% group_by(pID, wID, AOI,Acc1_3vs4_6,RecTask) %>%
  summarise(FixationCount = mean(FixationCount,na.rm=T)) %>% as.data.frame()
statData = statData %>% group_by(pID, wID, AOI,Acc1_3vs4_6,RecTask) %>%
  summarise(FixationCount = mean(FixationCount,na.rm=T)) %>% as.data.frame()

statData = statData %>% group_by(pID,wID,Acc1_3vs4_6,RecTask) %>%
  summarise(FixationCount = mean(FixationCount,na.rm=T)) %>% as.data.frame()
statData$RecTask = factor(statData$RecTask, levels = c("Before", "After"))
statData$Acc1_3vs4_6 = factor(statData$Acc1_3vs4_6, levels = c("Scored>=4", "Scored<=3"))

t.test(FixationCount~Acc1_3vs4_6, data = statData)

p = ggplot(statData,aes(x=wID,y=FixationCount,fill=Acc1_3vs4_6)) +
  geom_bar(stat="summary",fun="mean",position="dodge") +
  geom_point(position = position_jitterdodge(jitter.width = .4,
                                             jitter.height = 0,
                                             dodge.width = .9),shape = 21, size = 3)+
  stat_summary(fun.data = "mean_se", geom="errorbar", position = position_dodge(width = 0.9),
               width = 0.4, linewidth = 0.7, color = "Black")+
  theme_bw(base_family = "serif")+
  #stat_compare_means(method = "t.test",label = "p.format", paired = F) +
  scale_fill_manual(values = c("Scored>=4" = "seagreen2", "Scored<=3" = "purple3"), 
                    labels = c("Remembered", "Forgotten")) +
  labs(fill = "Task Results") +
  xlab("Walk") +
  ylab("Fixation Count")

p + scale_y_continuous(expand = c(0,0)) + theme(
  axis.title.x = element_text(size = 15, face = "bold"), axis.title.y = element_text(size = 15, face = "bold"),
  panel.border = element_blank(), panel.grid = element_blank(),
  axis.line = element_line(), axis.text = element_text(colour = "black", size = rel(1.3), face = "bold"),
  legend.title = element_text(colour = "black", size = rel(1), face = "bold"),
  legend.text = element_text(colour = "black", size = rel(1), face = "bold"), legend.position = "top"
)
graph2ppt(file= "FixCountByWalk.pptx", width = 8, height = 5)




statData = dat
statData = statData %>% group_by(pID, wID, AOI,Acc1_3vs4_6,RecTask) %>%
  summarise(SaccadeCount = mean(SaccadeCount,na.rm=T)) %>% as.data.frame()
statData = statData %>% group_by(pID, wID, AOI,Acc1_3vs4_6,RecTask) %>%
  summarise(SaccadeCount = mean(SaccadeCount,na.rm=T)) %>% as.data.frame()

statData = statData %>% group_by(pID,wID,Acc1_3vs4_6,RecTask) %>%
  summarise(SaccadeCount = mean(SaccadeCount,na.rm=T)) %>% as.data.frame()
statData$RecTask = factor(statData$RecTask, levels = c("Before", "After"))
statData$Acc1_3vs4_6 = factor(statData$Acc1_3vs4_6, levels = c("Scored>=4", "Scored<=3"))

t.test(SaccadeCount~Acc1_3vs4_6, data = statData)

p = ggplot(statData,aes(x=wID,y=SaccadeCount,fill=Acc1_3vs4_6)) +
  geom_bar(stat="summary",fun="mean",position="dodge") +
  geom_point(position = position_jitterdodge(jitter.width = .4,
                                             jitter.height = 0,
                                             dodge.width = .9),shape = 21, size = 3)+
  stat_summary(fun.data = "mean_se", geom="errorbar", position = position_dodge(width = 0.9),
               width = 0.4, linewidth = 0.7, color = "Black")+
  theme_bw(base_family = "serif")+
  #stat_compare_means(method = "t.test",label = "p.format", paired = F) +
  scale_fill_manual(values = c("Scored>=4" = "seagreen2", "Scored<=3" = "purple3"), 
                    labels = c("Remembered", "Forgotten")) +
  labs(fill = "Task Results") +
  xlab("Walk") +
  ylab("Saccade Count")

p + scale_y_continuous(expand = c(0,0)) + theme(
  axis.title.x = element_text(size = 15, face = "bold"), axis.title.y = element_text(size = 15, face = "bold"),
  panel.border = element_blank(), panel.grid = element_blank(),
  axis.line = element_line(), axis.text = element_text(colour = "black", size = rel(1.3), face = "bold"),
  legend.title = element_text(colour = "black", size = rel(1), face = "bold"),
  legend.text = element_text(colour = "black", size = rel(1), face = "bold"), legend.position = "top"
)
graph2ppt(file= "SacCountByWalk.pptx", width = 8, height = 5)



statData = dat
statData = statData %>% group_by(pID, wID, AOI,Acc1_3vs4_6,RecTask) %>%
  summarise(TotalGazeDuration = mean(TotalGazeDuration,na.rm=T)) %>% as.data.frame()
statData = statData %>% group_by(pID, wID, AOI,Acc1_3vs4_6,RecTask) %>%
  summarise(TotalGazeDuration = mean(TotalGazeDuration,na.rm=T)) %>% as.data.frame()

statData = statData %>% group_by(pID,wID,Acc1_3vs4_6,RecTask) %>%
  summarise(TotalGazeDuration = mean(TotalGazeDuration,na.rm=T)) %>% as.data.frame()
statData$RecTask = factor(statData$RecTask, levels = c("Before", "After"))
statData$Acc1_3vs4_6 = factor(statData$Acc1_3vs4_6, levels = c("Scored>=4", "Scored<=3"))

t.test(TotalGazeDuration~Acc1_3vs4_6, data = statData)

p = ggplot(statData,aes(x=wID,y=TotalGazeDuration,fill=Acc1_3vs4_6)) +
  geom_bar(stat="summary",fun="mean",position="dodge") +
  geom_point(position = position_jitterdodge(jitter.width = .4,
                                             jitter.height = 0,
                                             dodge.width = .9),shape = 21, size = 3)+
  stat_summary(fun.data = "mean_se", geom="errorbar", position = position_dodge(width = 0.9),
               width = 0.4, linewidth = 0.7, color = "Black")+
  theme_bw(base_family = "serif")+
  #stat_compare_means(method = "t.test",label = "p.format", paired = F) +
  scale_fill_manual(values = c("Scored>=4" = "seagreen2", "Scored<=3" = "purple3"), 
                    labels = c("Remembered", "Forgotten")) +
  labs(fill = "Task Results") +
  xlab("Walk") +
  ylab("Total Fixation Duration (ms)")

p + scale_y_continuous(expand = c(0,0)) + theme(
  axis.title.x = element_text(size = 15, face = "bold"), axis.title.y = element_text(size = 15, face = "bold"),
  panel.border = element_blank(), panel.grid = element_blank(),
  axis.line = element_line(), axis.text = element_text(colour = "black", size = rel(1.3), face = "bold"),
  legend.title = element_text(colour = "black", size = rel(1), face = "bold"),
  legend.text = element_text(colour = "black", size = rel(1), face = "bold"), legend.position = "top"
)
graph2ppt(file= "TotalGazeDurByWalk.pptx", width = 8, height = 5)
############################################# Updated Graphs (By B&A)
statData = dat
statData = statData %>% group_by(pID, AOI,Acc1_3vs4_6,RecTask) %>%
  summarise(FixationCount = mean(FixationCount,na.rm=T)) %>% as.data.frame()
statData = statData %>% group_by(pID, AOI,Acc1_3vs4_6,RecTask) %>%
  summarise(FixationCount = mean(FixationCount,na.rm=T)) %>% as.data.frame()

statData = statData %>% group_by(pID,Acc1_3vs4_6,RecTask) %>%
  summarise(FixationCount = mean(FixationCount,na.rm=T)) %>% as.data.frame()
statData$RecTask = factor(statData$RecTask, levels = c("Before", "After"))
statData$Acc1_3vs4_6 = factor(statData$Acc1_3vs4_6, levels = c("Scored>=4", "Scored<=3"))

p = ggplot(statData,aes(x=RecTask,y=FixationCount,fill=Acc1_3vs4_6)) +
  geom_bar(stat="summary",fun="mean",position="dodge") +
  geom_point(position = position_jitterdodge(jitter.width = .4,
                                             jitter.height = 0,
                                             dodge.width = .9),shape = 21, size = 3)+
  stat_summary(fun.data = "mean_se", geom="errorbar", position = position_dodge(width = 0.9),
               width = 0.4, linewidth = 0.7, color = "Black")+
  theme_bw(base_family = "serif")+
  #stat_compare_means(method = "t.test",label = "p.format", paired = F) +
  scale_fill_manual(values = c("Scored>=4" = "seagreen2", "Scored<=3" = "purple3"), 
                    labels = c("Remembered", "Forgotten")) +
  labs(fill = "Task Results") +
  xlab("Recognition Task") +
  ylab("Fixation Count")

p + scale_y_continuous(expand = c(0,0)) + theme(
  axis.title.x = element_text(size = 15, face = "bold"), axis.title.y = element_text(size = 15, face = "bold"),
  panel.border = element_blank(), panel.grid = element_blank(),
  axis.line = element_line(), axis.text = element_text(colour = "black", size = rel(1.3), face = "bold"),
  legend.title = element_text(colour = "black", size = rel(1), face = "bold"),
  legend.text = element_text(colour = "black", size = rel(1), face = "bold"), legend.position = "top"
)
graph2ppt(file= "FixCountByB&AResults.pptx", width = 8, height = 5)



statData = dat
statData = statData %>% group_by(pID, AOI,Acc1_3vs4_6,RecTask) %>%
  summarise(SaccadeCount = mean(SaccadeCount,na.rm=T)) %>% as.data.frame()
statData = statData %>% group_by(pID, AOI,Acc1_3vs4_6,RecTask) %>%
  summarise(SaccadeCount = mean(SaccadeCount,na.rm=T)) %>% as.data.frame()

statData = statData %>% group_by(pID,Acc1_3vs4_6,RecTask) %>%
  summarise(SaccadeCount = mean(SaccadeCount,na.rm=T)) %>% as.data.frame()
statData$RecTask = factor(statData$RecTask, levels = c("Before", "After"))
statData$Acc1_3vs4_6 = factor(statData$Acc1_3vs4_6, levels = c("Scored>=4", "Scored<=3"))

p = ggplot(statData,aes(x=RecTask,y=SaccadeCount,fill=Acc1_3vs4_6)) +
  geom_bar(stat="summary",fun="mean",position="dodge") +
  geom_point(position = position_jitterdodge(jitter.width = .4,
                                             jitter.height = 0,
                                             dodge.width = .9),shape = 21, size = 3)+
  stat_summary(fun.data = "mean_se", geom="errorbar", position = position_dodge(width = 0.9),
               width = 0.4, linewidth = 0.7, color = "Black")+
  theme_bw(base_family = "serif")+
  #stat_compare_means(method = "t.test",label = "p.format", paired = F) +
  scale_fill_manual(values = c("Scored>=4" = "seagreen2", "Scored<=3" = "purple3"), 
                    labels = c("Remembered", "Forgotten")) +
  labs(fill = "Task Results") +
  xlab("Recognition Task") +
  ylab("Saccade Count")

p + scale_y_continuous(expand = c(0,0)) + theme(
  axis.title.x = element_text(size = 15, face = "bold"), axis.title.y = element_text(size = 15, face = "bold"),
  panel.border = element_blank(), panel.grid = element_blank(),
  axis.line = element_line(), axis.text = element_text(colour = "black", size = rel(1.3), face = "bold"),
  legend.title = element_text(colour = "black", size = rel(1), face = "bold"),
  legend.text = element_text(colour = "black", size = rel(1), face = "bold"), legend.position = "top"
)
graph2ppt(file= "SacCountByB&AResults.pptx", width = 8, height = 5)



statData = dat
statData = statData %>% group_by(pID, AOI,Acc1_3vs4_6,RecTask) %>%
  summarise(TotalGazeDuration = mean(TotalGazeDuration,na.rm=T)) %>% as.data.frame()
statData = statData %>% group_by(pID, AOI,Acc1_3vs4_6,RecTask) %>%
  summarise(TotalGazeDuration = mean(TotalGazeDuration,na.rm=T)) %>% as.data.frame()

statData = statData %>% group_by(pID,Acc1_3vs4_6,RecTask) %>%
  summarise(TotalGazeDuration = mean(TotalGazeDuration,na.rm=T)) %>% as.data.frame()
statData$RecTask = factor(statData$RecTask, levels = c("Before", "After"))
statData$Acc1_3vs4_6 = factor(statData$Acc1_3vs4_6, levels = c("Scored>=4", "Scored<=3"))

p = ggplot(statData,aes(x=RecTask,y=TotalGazeDuration,fill=Acc1_3vs4_6)) +
  geom_bar(stat="summary",fun="mean",position="dodge") +
  geom_point(position = position_jitterdodge(jitter.width = .4,
                                             jitter.height = 0,
                                             dodge.width = .9),shape = 21, size = 3)+
  stat_summary(fun.data = "mean_se", geom="errorbar", position = position_dodge(width = 0.9),
               width = 0.4, linewidth = 0.7, color = "Black")+
  theme_bw(base_family = "serif")+
  #stat_compare_means(method = "t.test",label = "p.format", paired = F) +
  scale_fill_manual(values = c("Scored>=4" = "seagreen2", "Scored<=3" = "purple3"), 
                    labels = c("Remembered", "Forgotten")) +
  labs(fill = "Task Results") +
  xlab("Recognition Task") +
  ylab("Total Fixation Duration (ms)")

p + scale_y_continuous(expand = c(0,0)) + theme(
  axis.title.x = element_text(size = 15, face = "bold"), axis.title.y = element_text(size = 15, face = "bold"),
  panel.border = element_blank(), panel.grid = element_blank(),
  axis.line = element_line(), axis.text = element_text(colour = "black", size = rel(1.3), face = "bold"),
  legend.title = element_text(colour = "black", size = rel(1), face = "bold"),
  legend.text = element_text(colour = "black", size = rel(1), face = "bold"), legend.position = "top"
)
graph2ppt(file= "TotalGazeDurByB&AResults.pptx", width = 8, height = 5)
###################################### Updated Graphs --------
########################################### TotalGazeDuration (Results before recognition task) ####
statData = dat
statData$DV = statData$TotalGazeDuration
statData$Accuracy = statData$Acc1_3vs4_6

statData = statData[statData$RecTask=="Before",]
statData = statData %>% group_by(pID, AOI,Accuracy,RecTask) %>%
  summarise(DV = mean(DV,na.rm=T),N=n()) %>% as.data.frame()
statData = statData %>% group_by(pID,Accuracy,RecTask) %>%
  summarise(DV = mean(DV,na.rm=T),N=n()) %>% as.data.frame()
statData$Accuracy = factor(statData$Accuracy, levels = c("Scored>=4", "Scored<=3"),
                           labels = c("Remembered", "Forgotten"))

A = statData %>% pivot_wider(id_cols = c(pID),names_from = Accuracy,
                             values_from = DV) %>% as.data.frame()
t.test(A$Remembered,
       A$Forgotten, paired = T)  # Fix p-Values for all other tests to be paired.

graphDat = statData
graphDat$X = graphDat$Accuracy
graphDat$Y = graphDat$DV
graphDat$Fill = graphDat$Accuracy
xLab = "Recognition Task Results"
yLab = "Total Fixation Duration (ms)"
fillLab = "Task Results"
BarGraphRemvsForgotten(graphDat)
graph2ppt(file= "TotalGazeDurBefore.pptx", width = 8, height = 5)


########################################### FixationCount (Results before recognition task) ####
statData = dat
statData$DV = statData$FixationCount
statData$Accuracy = statData$Acc1_3vs4_6

statData = statData[statData$RecTask=="Before",]
statData = statData %>% group_by(pID, AOI,Accuracy,RecTask) %>%
  summarise(DV = mean(DV,na.rm=T),N=n()) %>% as.data.frame()
statData = statData %>% group_by(pID,Accuracy,RecTask) %>%
  summarise(DV = mean(DV,na.rm=T),N=n()) %>% as.data.frame()
statData$Accuracy = factor(statData$Accuracy, levels = c("Scored>=4", "Scored<=3"),
                           labels = c("Remembered", "Forgotten"))

t.test(DV~Accuracy, data = statData)

graphDat = statData
graphDat$X = graphDat$Accuracy
graphDat$Y = graphDat$DV
graphDat$Fill = graphDat$Accuracy
xLab = "Recognition Task Results"
yLab = "Fixation Count"
fillLab = "Task Results"
BarGraphRemvsForgotten(graphDat)
graph2ppt(file= "FixCountBefore.pptx", width = 8, height = 5)

########################################### SaccadeCount (Results before recognition task) ####
statData = dat
statData$DV = statData$SaccadeCount
statData$Accuracy = statData$Acc1_3vs4_6

statData = statData[statData$RecTask=="Before",]
statData = statData %>% group_by(pID, AOI,Accuracy,RecTask) %>%
  summarise(DV = mean(DV,na.rm=T),N=n()) %>% as.data.frame()
statData = statData %>% group_by(pID,Accuracy,RecTask) %>%
  summarise(DV = mean(DV,na.rm=T),N=n()) %>% as.data.frame()
statData$Accuracy = factor(statData$Accuracy, levels = c("Scored>=4", "Scored<=3"),
                           labels = c("Remembered", "Forgotten"))

t.test(DV~Accuracy, data = statData)

graphDat = statData
graphDat$X = graphDat$Accuracy
graphDat$Y = graphDat$DV
graphDat$Fill = graphDat$Accuracy
xLab = "Recognition Task Results"
yLab = "Saccade Count"
fillLab = "Task Results"
BarGraphRemvsForgotten(graphDat)
graph2ppt(file= "SacCountBefore.pptx", width = 8, height = 5)

########################################### TotalGazeDuration (before vs. after recognition task) ####
statData = dat
statData$DV = statData$TotalGazeDuration
statData$Accuracy = statData$Acc1_3vs4_6

statData = statData %>% group_by(pID, AOI,Accuracy,RecTask) %>%
  summarise(DV = mean(DV,na.rm=T),N=n()) %>% as.data.frame()
statData = statData %>% group_by(pID,Accuracy,RecTask) %>%
  summarise(DV = mean(DV,na.rm=T),N=n()) %>% as.data.frame()
statData$Accuracy = factor(statData$Accuracy, levels = c("Scored>=4", "Scored<=3"),
                           labels = c("Remembered", "Forgotten"))
statData$RecTask = factor(statData$RecTask, levels = c("Before", "After"))

t.test(DV~Accuracy, data = statData)

graphDat = statData
graphDat$X = graphDat$Accuracy
graphDat$Y = graphDat$DV
graphDat$Fill = graphDat$RecTask
xLab = "Recognition Task Results"
yLab = "Total Gaze Duration (ms)"
fillLab = "Walk Time"
BarGraphRemvsForgotten(graphDat, colorCodes = c("#DB532C","#555D79"))
graph2ppt(file= "TotalGazeDurByB&AResults.pptx", width = 8, height = 5)


########################################### FixationCount (before vs. after recognition task) ####
statData = dat
statData$DV = statData$FixationCount
statData$Accuracy = statData$Acc1_3vs4_6

statData = statData %>% group_by(pID, AOI,Accuracy,RecTask) %>%
  summarise(DV = mean(DV,na.rm=T),N=n()) %>% as.data.frame()
statData = statData %>% group_by(pID,Accuracy,RecTask) %>%
  summarise(DV = mean(DV,na.rm=T),N=n()) %>% as.data.frame()
statData$Accuracy = factor(statData$Accuracy, levels = c("Scored>=4", "Scored<=3"),
                           labels = c("Remembered", "Forgotten"))
statData$RecTask = factor(statData$RecTask, levels = c("Before", "After"))

t.test(DV~Accuracy, data = statData)

graphDat = statData
graphDat$X = graphDat$Accuracy
graphDat$Y = graphDat$DV
graphDat$Fill = graphDat$RecTask
xLab = "Recognition Task Results"
yLab = "Fixation Count"
fillLab = "Walk Time"
BarGraphRemvsForgotten(graphDat, colorCodes = c("#DB532C","#555D79"))
graph2ppt(file= "FixCountByB&AResults.pptx", width = 8, height = 5)

########################################### SaccadeCount (before vs. after recognition task) ####
statData = dat
statData$DV = statData$SaccadeCount
statData$Accuracy = statData$Acc1_3vs4_6

statData = statData %>% group_by(pID, AOI,Accuracy,RecTask) %>%
  summarise(DV = mean(DV,na.rm=T),N=n()) %>% as.data.frame()
statData = statData %>% group_by(pID,Accuracy,RecTask) %>%
  summarise(DV = mean(DV,na.rm=T),N=n()) %>% as.data.frame()
statData$Accuracy = factor(statData$Accuracy, levels = c("Scored>=4", "Scored<=3"),
                           labels = c("Remembered", "Forgotten"))
statData$RecTask = factor(statData$RecTask, levels = c("Before", "After"))

t.test(DV~Accuracy, data = statData)

graphDat = statData
graphDat$X = graphDat$Accuracy
graphDat$Y = graphDat$DV
graphDat$Fill = graphDat$RecTask
xLab = "Recognition Task Results"
yLab = "Saccade Count"
fillLab = "Walk Time"
BarGraphRemvsForgotten(graphDat, colorCodes = c("#DB532C","#555D79"))
graph2ppt(file= "SacCountByB&AResults.pptx", width = 8, height = 5)

########################################### TotalGazeDuration (across Walks) ####
statData = dat
statData$DV = statData$TotalGazeDuration
statData$Accuracy = statData$Acc1_3vs4_6
statData = statData %>% group_by(pID, wID,Accuracy) %>%
  summarise(DV = mean(DV,na.rm=T),N=n()) %>% as.data.frame()
statData$Accuracy = factor(statData$Accuracy, levels = c("Scored>=4", "Scored<=3"),
                           labels = c("Remembered", "Forgotten"))

A = statData %>% pivot_wider(id_cols = c(pID,wID),names_from = Accuracy,
                             values_from = DV) %>% as.data.frame()

t.test(A$Remembered[statData$wID==1],
       A$Forgotten[statData$wID==1], paired = T)
t.test(A$Remembered[statData$wID==2],
       A$Forgotten[statData$wID==2], paired = T)
t.test(A$Remembered[statData$wID==3],
       A$Forgotten[statData$wID==3], paired = T)
t.test(A$Remembered[statData$wID==4],
       A$Forgotten[statData$wID==4], paired = T)
t.test(A$Remembered[statData$wID==5],
       A$Forgotten[statData$wID==5], paired = T)
t.test(A$Remembered[statData$wID==6],
       A$Forgotten[statData$wID==6], paired = T)
t.test(A$Remembered[statData$wID==7],
       A$Forgotten[statData$wID==7], paired = T)
t.test(A$Remembered[statData$wID==8],
       A$Forgotten[statData$wID==8], paired = T)

graphDat = statData
graphDat$X = graphDat$Accuracy
graphDat$Y = graphDat$DV
graphDat$Fill = as.factor(graphDat$wID)
xLab = "Accuracy"
yLab = "Total Gaze Duration (ms)"
fillLab = "Walk Index"
BarGraphRemvsForgotten(graphDat,colorCodes=palette.colors(8,palette = "Set 3"))

graphDat = statData
graphDat$X = factor(graphDat$wID)
graphDat$Y = graphDat$DV
graphDat$Fill = graphDat$Accuracy
xLab = "Walk Index"
yLab = "Total Gaze Duration (ms)"
fillLab = "Accuracy"
BarGraphRemvsForgotten(graphDat)

graph2ppt(file= "TotalGazeDurByWalk.pptx", width = 8, height = 5)


########################################### FixationCount (across Walks) ####
statData = dat
statData$DV = statData$FixationCount
statData$Accuracy = statData$Acc1_3vs4_6
statData = statData %>% group_by(pID, wID,Accuracy) %>%
  summarise(DV = mean(DV,na.rm=T),N=n()) %>% as.data.frame()
statData$Accuracy = factor(statData$Accuracy, levels = c("Scored>=4", "Scored<=3"),
                           labels = c("Remembered", "Forgotten"))

A = statData %>% pivot_wider(id_cols = c(pID,wID),names_from = Accuracy,
                             values_from = DV) %>% as.data.frame()

t.test(A$Remembered[statData$wID==1],
       A$Forgotten[statData$wID==1], paired = T)
t.test(A$Remembered[statData$wID==2],
       A$Forgotten[statData$wID==2], paired = T)
t.test(A$Remembered[statData$wID==3],
       A$Forgotten[statData$wID==3], paired = T)
t.test(A$Remembered[statData$wID==4],
       A$Forgotten[statData$wID==4], paired = T)
t.test(A$Remembered[statData$wID==5],
       A$Forgotten[statData$wID==5], paired = T)
t.test(A$Remembered[statData$wID==6],
       A$Forgotten[statData$wID==6], paired = T)
t.test(A$Remembered[statData$wID==7],
       A$Forgotten[statData$wID==7], paired = T)
t.test(A$Remembered[statData$wID==8],
       A$Forgotten[statData$wID==8], paired = T)

graphDat = statData
graphDat$X = graphDat$Accuracy
graphDat$Y = graphDat$DV
graphDat$Fill = as.factor(graphDat$wID)
xLab = "Accuracy"
yLab = "Fixation Count"
fillLab = "Walk Index"
BarGraphRemvsForgotten(graphDat,colorCodes=palette.colors(8,palette = "Set 3"))

graphDat = statData
graphDat$X = factor(graphDat$wID)
graphDat$Y = graphDat$DV
graphDat$Fill = graphDat$Accuracy
xLab = "Walk Index"
yLab = "Fixation Count"
fillLab = "Accuracy"
BarGraphRemvsForgotten(graphDat)

graph2ppt(file= "FixCountByWalk.pptx", width = 8, height = 5)

########################################### SaccadeCount (across Walks) ####
statData = dat
statData$DV = statData$SaccadeCount
statData$Accuracy = statData$Acc1_3vs4_6
statData = statData %>% group_by(pID, wID,Accuracy) %>%
  summarise(DV = mean(DV,na.rm=T),N=n()) %>% as.data.frame()
statData$Accuracy = factor(statData$Accuracy, levels = c("Scored>=4", "Scored<=3"),
                           labels = c("Remembered", "Forgotten"))

A = statData %>% pivot_wider(id_cols = c(pID,wID),names_from = Accuracy,
                             values_from = DV) %>% as.data.frame()

t.test(A$Remembered[statData$wID==1],
       A$Forgotten[statData$wID==1], paired = T)
t.test(A$Remembered[statData$wID==2],
       A$Forgotten[statData$wID==2], paired = T)
t.test(A$Remembered[statData$wID==3],
       A$Forgotten[statData$wID==3], paired = T)
t.test(A$Remembered[statData$wID==4],
       A$Forgotten[statData$wID==4], paired = T)
t.test(A$Remembered[statData$wID==5],
       A$Forgotten[statData$wID==5], paired = T)
t.test(A$Remembered[statData$wID==6],
       A$Forgotten[statData$wID==6], paired = T)
t.test(A$Remembered[statData$wID==7],
       A$Forgotten[statData$wID==7], paired = T)
t.test(A$Remembered[statData$wID==8],
       A$Forgotten[statData$wID==8], paired = T)

graphDat = statData
graphDat$X = graphDat$Accuracy
graphDat$Y = graphDat$DV
graphDat$Fill = as.factor(graphDat$wID)
xLab = "Accuracy"
yLab = "Saccade Count"
fillLab = "Walk Index"
BarGraphRemvsForgotten(graphDat,colorCodes=palette.colors(8,palette = "Set 3"))

graphDat = statData
graphDat$X = factor(graphDat$wID)
graphDat$Y = graphDat$DV
graphDat$Fill = graphDat$Accuracy
xLab = "Walk Index"
yLab = "Saccade Count"
fillLab = "Accuracy"
BarGraphRemvsForgotten(graphDat)

graph2ppt(file= "SacCountByWalk.pptx", width = 8, height = 5)

#################################### D Prime graph ----
dpri = drime 
dpri = dpri[dpri$pID != "Average",]
dpri$pID = as.factor(dpri$pID)

ggplot(dpri,aes(x=" ",y=dPrime,fill = 1)) +
  geom_bar(stat="summary",fun="mean") +
  geom_point(position = position_jitterdodge(jitter.width = .4,
                                             jitter.height = 0,
                                             dodge.width = .9),shape = 21, size = 3)+
  stat_summary(fun.data = "mean_se", geom="errorbar",  width = 0.2, linewidth = 1, color = "black")+
  xlab("d Prime")+
  scale_y_continuous(expand = expansion(mult = c(0,.05)))+
  gDefault

ggplot(dpri,aes(x=" ",y=dPrime,fill = "1")) +
  geom_bar(stat="summary",fun="mean") +
  geom_point(position = position_jitterdodge(jitter.width = .4,
                                             jitter.height = 0,
                                             dodge.width = .9),shape = 21, size = 3, aes(fill = pID))+
  stat_summary(fun.data = "mean_se", geom="errorbar",  width = 0.2, linewidth = 1, color = "black")+
  xlab("d Prime")+
  scale_y_continuous(expand = expansion(mult = c(0,.05)))+
  gDefault


#################################### D Prime graph ----
dpri = drime 
dpri = dpri[dpri$pID != "Average",]
dpri$pID = as.factor(dpri$pID)
dpri = dpri[order(dpri$HR),]

ggplot(dpri,aes(x=FAR,y=HR)) +
  geom_point(shape = 21, size = 3)+
  # geom_line()+
  stat_smooth(method = "lm", formula = y~splines::bs(x,1),col = "red")+
  xlab("False Positive Rate")+
  ylab("Hit Rate")+
  expand_limits(x=c(0,1))+
  xlim(c(0,1))+
  ylim(c(0,1))
scale_y_continuous(expand = expansion(mult = c(0,.05)))+
  gDefault





##############################  FixationCount --------
graphDat = dat
graphDat$ACC = graphDat$Acc1_3vs4_6

# TotalGazeDuration
# SaccadeCount
# FixationCount
# A = graphDat[graphDat$FixationCount>50,]

statData = dat
statData = statData %>% group_by(pID, AOI,Acc1_3vs4_6,RecTask) %>%
  summarise(FixationCount = mean(FixationCount,na.rm=T)) %>% as.data.frame()
statData = statData %>% group_by(pID, AOI,Acc1_3vs4_6,RecTask) %>%
  summarise(FixationCount = mean(FixationCount,na.rm=T)) %>% as.data.frame()

statData = statData %>% group_by(pID,Acc1_3vs4_6,RecTask) %>%
  summarise(FixationCount = mean(FixationCount,na.rm=T)) %>% as.data.frame()

t.test(FixationCount~Acc1_3vs4_6, data = statData)

statData$RecTask = factor(statData$RecTask, levels = c("Before", "After"))
statData$Acc1_3vs4_6 = factor(statData$Acc1_3vs4_6, levels = c("Scored>=4", "Scored<=3"))
ggplot(statData,aes(x=RecTask,y=FixationCount,fill=Acc1_3vs4_6)) + 
  geom_bar(stat="summary",fun="mean",position="dodge")+
  geom_point(position = position_jitterdodge(jitter.width = .4,
                                             jitter.height = 0,
                                             dodge.width = .9),shape = 21)+
  stat_summary(fun.data = "mean_se", geom="errorbar",position="dodge",color = "Black")+
  theme_bw(base_family = "serif")+
  stat_compare_means(method = "t.test",label = "p.format", paired = F)
graph2ppt(file= "FixationCountAfterO.pptx", width = 8, height = 5)


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
statData = statData %>% group_by(pID, AOI,Acc1_3vs4_6,RecTask) %>%
  summarise(TotalGazeDuration = mean(TotalGazeDuration,na.rm=T)) %>% as.data.frame()

statData = statData %>% group_by(pID,Acc1_3vs4_6,RecTask) %>%
  summarise(TotalGazeDuration = mean(TotalGazeDuration,na.rm=T)) %>% as.data.frame()

t.test(TotalGazeDuration~Acc1_3vs4_6, data = statData)

fSize = 16
gDefault =   theme(
  text = element_text(family = "serif"),
  strip.text = element_text(size=fSize-2, face="bold"),
  plot.title = element_text(size = fSize),
  axis.title = element_text(size = fSize),
  legend.title = element_text(size=fSize),
  legend.text=element_text(size=fSize-4),
  axis.text = element_text(size = fSize-4),
  panel.background = element_rect(fill = "transparent", color = "black", 
                                  linewidth = 0.7, linetype =  "solid"),
  plot.background = element_rect(fill = "transparent"),
  axis.line = element_line(linewidth = .7, linetype =  "solid"),
  strip.background = element_rect(fill="transparent"),
  panel.grid.major.x = element_blank(),
  panel.grid.minor = element_blank()
)
ggplot(statData,aes(x=RecTask,y=TotalGazeDuration,fill=Acc1_3vs4_6)) + 
  geom_bar(stat="summary",fun="mean",position="dodge")+
  geom_point(position = position_jitterdodge(jitter.width = .4,
                                             jitter.height = 0,
                                             dodge.width = .9),shape = 21)+
  gDefault +
  stat_summary(fun.data = "mean_se", geom="errorbar",position="dodge",color = "Black")+
  theme_bw(base_family = "serif")+
  stat_compare_means(method = "t.test",label = "p.format", paired = F)
graph2ppt(file= "TotalGazeDurationGraphAfter.pptx", width = 8, height = 5)

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
statData = statData %>% group_by(pID, AOI,Acc1_3vs4_6,RecTask) %>%
  summarise(SaccadeCount = mean(SaccadeCount,na.rm=T)) %>% as.data.frame()

statData = statData %>% group_by(pID,Acc1_3vs4_6,RecTask) %>%
  summarise(SaccadeCount = mean(SaccadeCount,na.rm=T)) %>% as.data.frame()


t.test(SaccadeCount~Acc1_3vs4_6, data = statData)

fSize = 16
gDefault =   theme(
  text = element_text(family = "serif"),
  strip.text = element_text(size=fSize-2, face="bold"),
  plot.title = element_text(size = fSize),
  axis.title = element_text(size = fSize),
  legend.title = element_text(size=fSize),
  legend.text=element_text(size=fSize-4),
  axis.text = element_text(size = fSize-4),
  panel.background = element_rect(fill = "transparent", color = "black", 
                                  linewidth = 0.7, linetype =  "solid"),
  plot.background = element_rect(fill = "transparent"),
  axis.line = element_line(linewidth = .7, linetype =  "solid"),
  strip.background = element_rect(fill="transparent"),
  panel.grid.major.x = element_blank(),
  panel.grid.minor = element_blank()
)

ggplot(statData,aes(x=RecTask,y=SaccadeCount,fill=Acc1_3vs4_6)) + 
  geom_bar(stat="summary",fun="mean",position="dodge")+
  geom_point(position = position_jitterdodge(jitter.width = .4,
                                             jitter.height = 0,
                                             dodge.width = .9),shape = 21)+
  stat_summary(fun.data = "mean_se", geom="errorbar",position="dodge",color = "Black")+
  theme_bw(base_family = "serif")+
  stat_compare_means(method = "t.test",label = "p.format", paired = F)
graph2ppt(file= "SaccadeCountAfter.pptx", width = 8, height = 5)


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

########################### Initialization ##########################
rm(list=ls(all=TRUE))

x = readline()
C:/Users/Lensky Augustin/Downloads/LoadandMergeETandRecogTaskData/LoadandMergeETandRecogTaskData
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
statData = statData %>% group_by(pID, AOI,Acc1_3vs4_6,RecTask) %>%
                summarise(FixationCount = mean(FixationCount,na.rm=T)) %>% as.data.frame()

statData = statData %>% group_by(pID,Acc1_3vs4_6,RecTask) %>%
  summarise(FixationCount = mean(FixationCount,na.rm=T)) %>% as.data.frame()

t.test(FixationCount~Acc1_3vs4_6, data = statData)

ggplot(statData,aes(x=RecTask,y=FixationCount,fill=Acc1_3vs4_6)) + 
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
statData = statData %>% group_by(pID, AOI,Acc1_3vs4_6,RecTask) %>%
  summarise(TotalGazeDuration = mean(TotalGazeDuration,na.rm=T)) %>% as.data.frame()

statData = statData %>% group_by(pID,Acc1_3vs4_6,RecTask) %>%
  summarise(TotalGazeDuration = mean(TotalGazeDuration,na.rm=T)) %>% as.data.frame()

t.test(TotalGazeDuration~Acc1_3vs4_6, data = statData)

ggplot(statData,aes(x=RecTask,y=TotalGazeDuration,fill=Acc1_3vs4_6)) + 
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
statData = statData %>% group_by(pID, AOI,Acc1_3vs4_6,RecTask) %>%
  summarise(SaccadeCount = mean(SaccadeCount,na.rm=T)) %>% as.data.frame()

statData = statData %>% group_by(pID,Acc1_3vs4_6,RecTask) %>%
  summarise(SaccadeCount = mean(SaccadeCount,na.rm=T)) %>% as.data.frame()

t.test(SaccadeCount~Acc1_3vs4_6, data = statData)

ggplot(statData,aes(x=RecTask,y=SaccadeCount,fill=Acc1_3vs4_6)) + 
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

####################### Multilevel Model Sample Code ############################

datModel = dat
datModel$ACC = ifelse(datModel$Acc1_3vs4_6=="Scored<=3",0,1)
datModel$RecTask = factor(datModel$RecTask,levels=c("Before","After"))

#datModel = RemoveOutliers(datModel,factorNames = c("sessID","testName"),varNames = "value",Criteria=3)
#datModel = datModel[!is.na(datModel$value),]

m = lmer(SaccadeCount ~ RecTask+ACC+
            (1|pID/wID), data=datModel, na.action=na.exclude)
# m = glmer(ACC ~ RecTask+FixationCount+
#                  (1|pID/wID), data=datModel, family = binomial,na.action=na.exclude)
tab_model(m,show.se = T,show.stat = T,show.est = T)
summary(m)
plot(allEffects(m, default.levels=4))

# mLogGP = lme(fixed = ACC ~ testName+ageChange+value+ageBase,
#              random = 1+testName|SID,
#              data=datModel)

m = glmer(ACC ~ RecTask+FixationCount+
            (1|pID:wID), data=datModel, family = binomial,na.action=na.exclude)
tab_model(m,show.se = T,show.stat = T,show.est = T)
summary(m)
plot(allEffects(m, default.levels=4))