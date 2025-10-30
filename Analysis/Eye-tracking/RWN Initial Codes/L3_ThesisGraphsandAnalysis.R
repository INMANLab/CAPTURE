########################### Initialization ##########################
rm(list=ls(all=TRUE))

x = readline()
D:\CAPTURE Project\LenskyData\LoadandMergeETandRecogTaskData
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
       pROC,
       ROCR,
       RColorBrewer,
       viridis)

#--------> TODO Compute the number of times each AOI were observed across walks within each participants
#--------> TODO Compute the proportional size of the AOI relative to screen
#--------> Figure out the problem with Response level 3 in the MLM
#--------> Fix the Scaling of the logistic model
############################### Load Data -------
dat = read.csv(file="MergedData_NoPartialSaccades.csv")
behavioralDat = read.csv("BehavioralPerformance.csv")
trialLevelPerformanceDat = read.csv("TrialLevelBehaviorData.csv")
############################### Variables -------
gDefaultPars = theme(
  text = element_text(family = "sans"),
  axis.title.x = element_text(size = 15, face = "bold"), 
  axis.title.y = element_text(size = 15, face = "bold"),
  panel.border = element_blank(), 
  panel.grid = element_blank(),
  axis.line = element_line(),
  panel.background = element_blank(),
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
###################################### MultiLevel Model --------
#------------------------ Multilevel Model with Gaze duration ----
datModel = dat
datModel$ACC = ifelse(datModel$Acc1_3vs4_6=="Scored<=3",0,1)
datModel$RecTask = factor(datModel$RecTask,levels=c("Before","After"))
datModel$Response = as.numeric(datModel$Response)
datModel$DV = datModel$TotalGazeDuration
DVString = "Gaze Duration"

m = lmer(DV ~ RecTask+ACC+
           (1|pID/wID), data=datModel, na.action=na.exclude)

m = lmer(DV ~ RecTask+Response+
           (1|pID/wID), data=datModel, na.action=na.exclude)

tab_model(m,show.se = T,show.stat = T,show.est = T)
summary(m)
plot(allEffects(m, default.levels=4))

allEffects(m)

plot(effect('RecTask',m),
     rug=F,multiline=T,grid=F, 
     ci.style = "bars",
     main="",xlab = "Task Time",ylab = DVString,
     xaxs="i",yaxs = "i")

plot(effect('Response',m, xlevels = list(c(1,2,3,4,5,6)) ),
     rug=F,multiline=T,grid=F, 
     ci.style = "bands",
     main="",xlab = "Response",ylab = DVString,
     xlim =c(1,6) ,xaxs="i",yaxs = "i")

#------------------------ Multilevel Model with Gaze duration ----
datModel = dat
datModel$ACC = ifelse(datModel$Acc1_3vs4_6=="Scored<=3",0,1)
datModel$RecTask = factor(datModel$RecTask,levels=c("Before","After"))
datModel$Response = as.numeric(datModel$Response)
datModel$DV = datModel$TotalGazeDuration
DVString = "Gaze Duration"

m = lmer(DV ~ RecTask+ACC+
           (1|pID/wID), data=datModel, na.action=na.exclude)

m = lmer(DV ~ RecTask+Response+
           (1|pID/wID), data=datModel, na.action=na.exclude)

tab_model(m,show.se = T,show.stat = T,show.est = T)
summary(m)
plot(allEffects(m, default.levels=4))

allEffects(m)

plot(effect('RecTask',m),
     rug=F,multiline=T,grid=F, 
     ci.style = "bars",
     main="",xlab = "Task Time",ylab = DVString,
     xaxs="i",yaxs = "i")

plot(effect('Response',m, xlevels = list(c(1,2,3,4,5,6)) ),
     rug=F,multiline=T,grid=F, 
     ci.style = "bands",
     main="",xlab = "Response",ylab = DVString,
     xlim =c(1,6) ,xaxs="i",yaxs = "i")

#------------------------ MLM within Each walk ------
datModel = dat
datModel$ACC = ifelse(datModel$Acc1_3vs4_6=="Scored<=3",0,1)
datModel$RecTask = factor(datModel$RecTask,levels=c("Before","After"))
datModel$Response = as.numeric(datModel$Response)
datModel$DV = datModel$TotalGazeDuration
DVString = "Gaze Duration"
datModel$DV = scale(datModel$DV,center = F, scale = T)
datModel$wID = as.factor(datModel$wID)
datModel$pID = as.factor(datModel$pID)


m = glmer(ACC ~ DV:wID+
            (1|pID), data=datModel, family = binomial,na.action=na.exclude)
tab_model(m,show.se = T,show.stat = T,show.est = T)
summary(m)
plot(allEffects(m, default.levels=4))

allEffects(m)

plot(effect('DV*wID',mod = m, x.var = "DV"),
     rug=T,multiline=T,grid=F, 
     ci.style = "bands",
     main="",xlab = DVString,ylab ="Accuracy" ,
     xaxs="i",yaxs = "i",xlim = c(0,4))


p =  interact_plot(m, pred = DV, modx = wID, 
                   # modx.values = c(8,10,12),
                   plot.points = F, point.alpha = 1,point.size = .75,point.shape = F,
                   line.thickness = 1.2,vary.lty=T, interval = F,int.type = "confidence",
                   linearity.check =F, rug =F,partial.residuals = F,centered = "all",
                   # y.label = "Item-Location Memory Accuracy",
                   x.label = DVString,
                   colors =palette.colors(8,palette = "Set 3") 
)
p = p+gDefaultPars
plot(p)


p =  interact_plot(m, pred = DV, modx = wID, 
                   # modx.values = c(8,10,12),
                   plot.points = F, point.alpha = 1,point.size = .75,point.shape = F,
                   line.thickness = 1.2,vary.lty=T, interval = T,int.type = "confidence",
                   linearity.check =F, rug =F,partial.residuals = F,centered = "all",
                   # y.label = "Item-Location Memory Accuracy",
                   x.label = DVString,
                   colors =palette.colors(8,palette = "Set 3") 
                   )
p = p+facet_wrap(~wID)+gDefaultPars
plot(p)

#------------------------ Multilevel Model with Saccade Counts ----
datModel = dat
datModel$ACC = ifelse(datModel$Acc1_3vs4_6=="Scored<=3",0,1)
datModel$RecTask = factor(datModel$RecTask,levels=c("Before","After"))
datModel$Response = as.numeric(datModel$Response)
datModel$DV = datModel$SaccadeCount
DVString = "Saccade Counts"

m = lmer(DV ~ RecTask+ACC+
            (1|pID/wID), data=datModel, na.action=na.exclude)

m = lmer(DV ~ RecTask+Response+
           (1|pID/wID), data=datModel, na.action=na.exclude)

tab_model(m,show.se = T,show.stat = T,show.est = T)
summary(m)
plot(allEffects(m, default.levels=4))

allEffects(m)

# A = datModel$SaccadeCount[datModel$Response==4]

# hist(datModel$SaccadeCount[datModel$Response==3 & datModel$SaccadeCount<10],10)

plot(effect('RecTask',m),
     rug=F,multiline=T,grid=F, 
     ci.style = "bars",
     main="",xlab = "Task Time",ylab = DVString,
     xaxs="i",yaxs = "i")

plot(effect('Response',m, xlevels = list(c(1,2,3,4,5,6)) ),
     rug=F,multiline=T,grid=F, 
     ci.style = "bands",
     main="",xlab = "Response",ylab = DVString,
     xlim =c(1,6) ,xaxs="i",yaxs = "i")





#------------------------ Multilevel Model Sample Code (Data before recognition task) ############################

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
###################################### Bar Graphs Forgotten vs.Remembered --------
#------------------------ TotalGazeDuration (Results before recognition task) ----
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


#------------------------ FixationCount (Results before recognition task) ####
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

#------------------------ SaccadeCount (Results before recognition task) ####
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

#------------------------ TotalGazeDuration (before vs. after recognition task) ####
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


#------------------------ FixationCount (before vs. after recognition task) ####
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

#------------------------ SaccadeCount (before vs. after recognition task) ####
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

#------------------------ TotalGazeDuration (across Walks) ####
statData = dat
statData$DV = statData$TotalGazeDuration
statData$Accuracy = statData$Acc1_3vs4_6
statData = statData %>% group_by(pID, wID,Accuracy) %>%
  summarise(DV = mean(DV,na.rm=T),N=n()) %>% as.data.frame()
statData$Accuracy = factor(statData$Accuracy, levels = c("Scored>=4", "Scored<=3"),
                           labels = c("Remembered", "Forgotten"))

A = statData %>% pivot_wider(id_cols = c(pID,wID),names_from = Accuracy,
                             values_from = DV) %>% as.data.frame()

t.test(A$Remembered[A$wID==1],
       A$Forgotten[A$wID==1], paired = T)
t.test(A$Remembered[A$wID==2],
       A$Forgotten[A$wID==2], paired = T)
t.test(A$Remembered[A$wID==3],
       A$Forgotten[A$wID==3], paired = T)
t.test(A$Remembered[A$wID==4],
       A$Forgotten[A$wID==4], paired = T)
t.test(A$Remembered[A$wID==5],
       A$Forgotten[A$wID==5], paired = T)
t.test(A$Remembered[A$wID==6],
       A$Forgotten[A$wID==6], paired = T)
t.test(A$Remembered[A$wID==7],
       A$Forgotten[A$wID==7], paired = T)
t.test(A$Remembered[A$wID==8],
       A$Forgotten[A$wID==8], paired = T)

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



# ttests  and total Gaze Duration vs fixation duration
#------------------------ FixationCount (across Walks) ####
statData = dat
statData$DV = statData$FixationCount
statData$Accuracy = statData$Acc1_3vs4_6
statData = statData %>% group_by(pID, wID,Accuracy) %>%
  summarise(DV = mean(DV,na.rm=T),N=n()) %>% as.data.frame()
statData$Accuracy = factor(statData$Accuracy, levels = c("Scored>=4", "Scored<=3"),
                           labels = c("Remembered", "Forgotten"))

A = statData %>% pivot_wider(id_cols = c(pID,wID),names_from = Accuracy,
                             values_from = DV) %>% as.data.frame()

t.test(A$Remembered[A$wID==1],
       A$Forgotten[A$wID==1], paired = T)
t.test(A$Remembered[A$wID==2],
       A$Forgotten[A$wID==2], paired = T)
t.test(A$Remembered[A$wID==3],
       A$Forgotten[A$wID==3], paired = T)
t.test(A$Remembered[A$wID==4],
       A$Forgotten[A$wID==4], paired = T)
t.test(A$Remembered[A$wID==5],
       A$Forgotten[A$wID==5], paired = T)
t.test(A$Remembered[A$wID==6],
       A$Forgotten[A$wID==6], paired = T)
t.test(A$Remembered[A$wID==7],
       A$Forgotten[A$wID==7], paired = T)
t.test(A$Remembered[A$wID==8],
       A$Forgotten[A$wID==8], paired = T)

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

#------------------------ SaccadeCount (across Walks) ####
statData = dat
statData$DV = statData$SaccadeCount
statData$Accuracy = statData$Acc1_3vs4_6
statData = statData %>% group_by(pID, wID,Accuracy) %>%
  summarise(DV = mean(DV,na.rm=T),N=n()) %>% as.data.frame()
statData$Accuracy = factor(statData$Accuracy, levels = c("Scored>=4", "Scored<=3"),
                           labels = c("Remembered", "Forgotten"))

A = statData %>% pivot_wider(id_cols = c(pID,wID),names_from = Accuracy,
                             values_from = DV) %>% as.data.frame()

t.test(A$Remembered[A$wID==1],
       A$Forgotten[A$wID==1], paired = T)
t.test(A$Remembered[A$wID==2],
       A$Forgotten[A$wID==2], paired = T)
t.test(A$Remembered[A$wID==3],
       A$Forgotten[A$wID==3], paired = T)
t.test(A$Remembered[A$wID==4],
       A$Forgotten[A$wID==4], paired = T)
t.test(A$Remembered[A$wID==5],
       A$Forgotten[A$wID==5], paired = T)
t.test(A$Remembered[A$wID==6],
       A$Forgotten[A$wID==6], paired = T)
t.test(A$Remembered[A$wID==7],
       A$Forgotten[A$wID==7], paired = T)
t.test(A$Remembered[A$wID==8],
       A$Forgotten[A$wID==8], paired = T)

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


###################################### Behavioral Data Analysis --------
#------------------------ D Prime graph ----
dpri = behavioralDat 
dpri$pID = as.factor(dpri$pID)

ggplot(dpri,aes(x=" ",y=dPrime,fill = 1)) +
  geom_bar(stat="summary",fun="mean") +
  geom_point(position = position_jitterdodge(jitter.width = .4,
                                             jitter.height = 0,
                                             dodge.width = .9),shape = 21, size = 3)+
  stat_summary(fun.data = "mean_se", geom="errorbar",  width = 0.2, linewidth = 1, color = "black")+
  xlab("d Prime")+
  scale_y_continuous(expand = expansion(mult = c(0,.05)))+
  gDefaultPars

ggplot(dpri,aes(x=" ",y=dPrime,fill = "1")) +
  geom_bar(stat="summary",fun="mean") +
  geom_point(position = position_jitterdodge(jitter.width = .4,
                                             jitter.height = 0,
                                             dodge.width = .9),shape = 21, size = 3, aes(fill = pID))+
  stat_summary(fun.data = "mean_se", geom="errorbar",  width = 0.2, linewidth = 1, color = "black")+
  xlab("d Prime")+
  scale_y_continuous(expand = expansion(mult = c(0,.05)))+
  gDefaultPars


#------------------------ ROC graph ----
bhDat = trialLevelPerformanceDat[,c("pID", "AOI","Condition", "Response", "ResponseBinary", "Performance")]
A = bhDat %>% group_by(pID,Condition,Response) %>%
  summarize(N=n()) %>% as.data.frame()

bhDat$trialLabels = as.numeric(ifelse(bhDat$Condition=="Old", 1, 0))
bhDat$trialID = factor(bhDat$AOI,levels = unique(bhDat$AOI),
                       labels = 1:length(unique(bhDat$AOI))) %>% as.numeric()
bhDat$Response = (bhDat$Response-min(bhDat$Response))/(max(bhDat$Response)-min(bhDat$Response)) 

# patientColors = palette.colors(n=7,palette = "Okabe-Ito")
# patientColors = patientColors[c(2:4,6,7)]
patientColors = viridis(5, alpha = 1, begin = 0, end = .9, direction = -1, option = "D")
for (pIdx in unique(bhDat$pID)){
  roc_obj <- roc(response = bhDat$trialLabels[bhDat$pID == pIdx], 
                 predictor = bhDat$Response[bhDat$pID == pIdx], direction = "<",
                 smooth=TRUE, smooth.n=10, percent = T)
  
  if(pIdx == 1){
    aucInfo = NULL
    A = plot(roc_obj, legacy.axes = TRUE, xlab = "False Alarm Rate (%)", ylab = "Hit Rate (%)", col = patientColors[pIdx],
         asp =NULL,lwd=3,
         xaxs="i",yaxs = "i",
         )
  }
  else{
    A = plot(roc_obj, legacy.axes = TRUE, xlab = "False Alarm Rate (%)", ylab = "Hit Rate (%)", col = patientColors[pIdx], add = TRUE,
         asp =NULL,lwd=3,
         xaxs="i",yaxs = "i",
         )
  }
  aucInfo = c(aucInfo,A$auc)
}
legend("bottomright", legend=paste("Patient ",unique(bhDat$pID),", AUC = ",round(aucInfo),sep = ""),
       col=patientColors, lwd=3)


patientColors = viridis(5, alpha = 1, begin = 0, end = .9, direction = -1, option = "D")
for (pIdx in unique(bhDat$pID)){
  roc_obj <- roc(response = bhDat$trialLabels[bhDat$pID == pIdx], 
                 predictor = bhDat$Response[bhDat$pID == pIdx], direction = "<",
                 smooth=F, percent = T)
  
  if(pIdx == 1){
    aucInfo = NULL
    A = plot(roc_obj, legacy.axes = TRUE, xlab = "False Alarm Rate (%)", ylab = "Hit Rate (%)", col = patientColors[pIdx],
             asp =NULL,lwd=3,
             xaxs="i",yaxs = "i",
    )
  }
  else{
    A = plot(roc_obj, legacy.axes = TRUE, xlab = "False Alarm Rate (%)", ylab = "Hit Rate (%)", col = patientColors[pIdx], add = TRUE,
             asp =NULL,lwd=3,
             xaxs="i",yaxs = "i",
    )
  }
  aucInfo = c(aucInfo,A$auc)
}
legend("bottomright", legend=paste("Patient ",unique(bhDat$pID),", AUC = ",round(aucInfo),sep = ""),
       col=patientColors, lwd=3)


# Sanity check with another package

patientColors = palette.colors(n=6,palette = "Okabe-Ito")
patientColors = patientColors[2:6]
for (pIdx in unique(bhDat$pID)){
  pred <- prediction(predictions = bhDat$Response[bhDat$pID == pIdx], 
                     labels = bhDat$trialLabels[bhDat$pID == pIdx])
  perf <- performance(pred,"tpr","fpr")
  if(pIdx == 1){
    plot(perf,colorize=FALSE,add = FALSE,col = patientColors[pIdx],lwd=2,xaxs="i",yaxs = "i")
  }
  else{
    plot(perf,colorize=FALSE,add = TRUE,col = patientColors[pIdx],lwd=2,xaxs="i",yaxs = "i")
  }
}
legend("bottomright", legend=paste("P_",unique(bhDat$pID),sep = ""),
       col=patientColors, lwd=2)









