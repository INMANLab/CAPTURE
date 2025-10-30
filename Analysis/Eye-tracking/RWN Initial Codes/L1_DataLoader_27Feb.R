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

################################# Load ET Data ----
#------> Load and merge ET csv files    
fileNum = 1
dat = NULL
for (pID in 1:5){
  for (wID in 1:8){
    ETfileList = list.files(pattern = paste("R",pID,".W",wID,".csv",sep = ""))
    if(length(ETfileList) > 0){
      print(paste("loading file:",fileNum," Patient:",pID," Walk:",wID,sep = ""))
      fileNum = fileNum+1
      datTemp = read.csv(ETfileList,sep = ",")
      if (ncol(datTemp)==1){
        datTemp = read.csv(ETfileList,sep = ";")
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
dat = dat[,varnames]
#------> remove extra space in AOI names
dat$AOI = gsub(" ","",dat$AOI)
#------> remove N/A and None AOIs as they are extra rows reported by blikshift for non-AOI gaze stats
dat = dat[grepl("LM.*",dat$AOI),]
#---------------------> Checks on the number of Participants, Walks, and AOIs
unique(dat[,c("pID","wID","AOI")]) %>% group_by(pID,wID) %>%
  summarise(N=n()) %>% as.data.frame()
summary(dat[,c("pID","wID")])

#-------------------> Temporarily for this analysis we removed the concurrent AOIs
unique(dat$AOI[grepl(".*#.*",dat$AOI)])
dat %>% group_by(pID) %>%
  summarise(N = n(),Concurrent = sum(grepl(".*#.*",AOI)), Percent = sum(grepl(".*#.*",AOI))/n()) %>% as.data.frame()
conCurrentAOIs = dat[grepl(".*#.*",dat$AOI),]
dat = dat[!(dat$AOI %in% unique(dat$AOI[grepl(".*#.*",dat$AOI)])),]
unique(dat$AOI)
################################# Load Recognition Task data ----
respDat = NULL
#ETfileList = list.files(path = paste(getwd(),"/ResponseData",sep = ""),pattern = ".*xlsx")
for (pID in 1:5){
  datTemp = read_excel(path = paste(getwd(), "/ResponseData/","RW",pID,".xlsx",sep = ""))
  datTemp = datTemp[,c("ReactionTime","Condition","Response(1-6)","ImageFile")]
  datTemp$pID = pID
  names(datTemp) = c("RT","Condition","Response","RecTaskLab","pID")
  respDat = rbind(respDat,datTemp)
}
unique(respDat$RecTaskLab)
#-----> Remove extra rows 
respDat = respDat[grepl(".*JPEG",respDat$RecTaskLab),]
#-----> Select only old items
behDat = respDat
respDat = respDat[grepl("LM",respDat$RecTaskLab),]
unique(respDat$RecTaskLab)
#------> Remove JPEG extension
respDat$RecTaskLab = gsub(".JPEG","",respDat$RecTaskLab)
behDat$RecTaskLab = gsub(".JPEG","",behDat$RecTaskLab)

#################################  Load and Check the Key File ----
labelKey = as.data.frame(read_excel("RecognitionTask_LandmarkLabel_KEY.xlsx"))
names(labelKey) = c("Order","AOI","RecTaskLab","R4Changed","R5Changed")
labelKey[grepl(".*a.JPEG",labelKey$AOI),]
labelKey[grepl(".*jpg",labelKey$AOI),]

labelKey$AOI = gsub(".JPEG","",labelKey$AOI)
labelKey$AOI = gsub(".jpg","",labelKey$AOI)
labelKey = labelKey[,c("AOI","RecTaskLab","R4Changed", "R5Changed")]

#################################  Checking the labels ET with Rec Data ----
#-----------------> Check the consistency of AOIs in the labelKey and Rec task
RecAOIs = unique(respDat$RecTaskLab)
KeyRecAOIs = unique(labelKey$RecTaskLab)
setequal(RecAOIs, KeyRecAOIs) # The answer has to be TRUE

#-----------------> Check the consistency of AOIs in the labelKey/Rec and ETData

#-----> Number of unique AOI within Participants
dat %>% group_by(pID,wID) %>%
  summarise(AOI_Num = length(unique(AOI))) %>% as.data.frame()


ETAOIs = unique(dat$AOI)
RecAOIs = unique(respDat$RecTaskLab)
KeyRecAOIs = unique(labelKey$RecTaskLab)
KeyETAOIs = unique(labelKey$AOI)

differentLabs = setdiff(ETAOIs,KeyETAOIs)
#-----> These has to fixed in the blikshift 
# Are they fixed?
FixedLabels = c("LM1127nb",
          "LM1128",
          "LM07n",
          "LM013N",
          "LM090a")
A = dat[dat$AOI %in% FixedLabels,]

#################################  Merging ET with AOI Label Key File to fix the AOI column values ----
#-------> Apply Key file to remap AOI names for: P2:W1-5  P5:W1-4 P5:W6-8
dat$PWID = paste("PW",dat$pID,dat$wID,sep = "_")
missLabeled = c("PW_2_1","PW_2_2","PW_2_3","PW_2_4","PW_2_5",
                "PW_5_1","PW_5_2","PW_5_3","PW_5_4","PW_5_6",
                "PW_5_7","PW_5_8")

dat = merge(dat,labelKey[,c("AOI","RecTaskLab")],all.x = T, by="AOI")
dat0 = dat

dat = dat0
A = dat[is.na(dat$RecTaskLab) & dat$PWID %in% missLabeled,c("PWID","AOI","RecTaskLab")]
dat$AOI[dat$PWID %in% missLabeled] = dat$RecTaskLab[dat$PWID %in% missLabeled]
# dat$AOI[dat$AOI %in% FixedLabels] = dat$RecTaskLab[dat$AOI %in% FixedLabels]

# Checking to see if all discrepencies are resolved
ETAOIs = unique(dat$AOI)
RecAOIs = unique(respDat$RecTaskLab)
KeyRecAOIs = unique(labelKey$RecTaskLab)

setequal(RecAOIs, KeyRecAOIs) # The answer has to be TRUE
setequal(ETAOIs, KeyRecAOIs) 
setdiff(KeyRecAOIs,ETAOIs)
# "LM123" was never in the FOV of participants

################################## Remove Changed AOIs for R4 R5 ----
patient4Changed = labelKey$RecTaskLab[labelKey$R4Changed == "Changed"]
patient5Changed = labelKey$RecTaskLab[labelKey$R5Changed == "Changed"]
A = dat[((dat$pID == 4) & (dat$AOI %in% patient4Changed)) |
          ((dat$pID == 5) & (dat$AOI %in% patient5Changed)),c("pID","wID","AOI")]
dat = dat[!(((dat$pID == 4) & (dat$AOI %in% patient4Changed)) |
            ((dat$pID == 5) & (dat$AOI %in% patient5Changed))),]
# write.csv(A, "Patient4and5Changeditems.csv", row.names = F)



################################## Merge with Response Data from Rec Task ----
dat$RecTaskLab = dat$AOI
dat = merge(dat,respDat,all.x = T, by=c("pID","RecTaskLab"))
#--------> Check NAs
A = dat[is.na(dat$Response),] #Expected to have no NA

################################## Compute behavioral performance ----
dat$Acc1_3vs4_6 = case_when(dat$Response<=3 ~ "Scored<=3",
                            dat$Response>=4 ~ "Scored>=4",
                            TRUE~NA)
dat$Acc1_3vs6 = case_when(dat$Response<=3 ~ "Scored<=3",
                          dat$Response==6 ~ "Scored=6",
                          TRUE~NA)
dat$Acc1_2vs5_6 = case_when(dat$Response<=2 ~ "Scored<=2",
                          dat$Response>=5 ~ "Scored>=5",
                          TRUE~NA)
dat$Acc1_4vs6 = case_when(dat$Response<=4 ~ "Scored<=4",
                          dat$Response==6 ~ "Scored=6",
                            TRUE~NA)

dat$RecTask = ifelse(dat$wID>5,"After","Before")
dat$RecTask[dat$pID == 4 & dat$wID == 6] = unique("Before")

dat$AOIviewed = ifelse(dat$TotalGazeDuration>0,1,0)

dat = dat[,c("pID","wID" , "AOI","Response", "Acc1_3vs4_6","Acc1_3vs6","Acc1_2vs5_6","Acc1_4vs6", "RecTask","AOIviewed", names(dat)[ !(names(dat) %in% c("pID","wID" , "AOI","Response", "Acc1_3vs4_6","Acc1_3vs6","Acc1_2vs5_6","Acc1_4vs6", "RecTask","AOIviewed"))])]
write.csv(dat,file="MergedData_NoPartialSaccades.csv",row.names = F)


################################## Create Behavioral Data ----
dat$wID_TaskStat = paste("W",dat$wID,"_",dat$RecTask,sep = "")
datWide = dat[, c("pID","AOI","Response","wID_TaskStat","TotalGazeDuration")]
datWide$Condition = unique("Old")
##----> Turn long dat to wide: NA values are those that had appeared in the videos and 0 values are those AOIs that appeared on the walk video but had no gaze
datWide = pivot_wider(data = datWide,
                      names_from ="wID_TaskStat", values_from = "TotalGazeDuration") 

names(behDat) = c("RT", "OldCondition", "Response", "AOI", "pID")
behDat = merge(behDat,datWide, by =c("pID","AOI","Response"), all.x = T )
behDat$Condition[is.na(behDat$Condition)] = "New"


behDat$ResponseBinary = ifelse(behDat$Response>=4,"Old","New")
behDat$Performance = case_when(behDat$Condition == "Old" & behDat$ResponseBinary == "Old" ~ "Hit",
                               behDat$Condition == "Old" & behDat$ResponseBinary == "New" ~ "Miss",
                               behDat$Condition == "New" & behDat$ResponseBinary == "Old" ~ "FA",
                               behDat$Condition == "New" & behDat$ResponseBinary == "New" ~ "CR")

behDat %>% group_by(pID,Condition) %>%
  summarise(N=n()) %>% as.data.frame()

# pID Condition   N
# 1       New    155
# 1       Old    145
# 2       New    157
# 2       Old    143
# 3       New    158
# 3       Old    142
# 4       New    166
# 4       Old    134
# 5       New    179
# 5       Old    121

behDat %>% group_by(pID) %>%
  summarise(N=n()) %>% as.data.frame()


behDatTrialBase = behDat
#---------> Save Trial Level behavioral Data
write.csv(behDatTrialBase,"TrialLevelBehaviorData.csv",row.names = F)

############## Compute dPrime ---------
behDat = behDatTrialBase
behDat = behDat %>% group_by(pID,Performance) %>%
          summarise(N=n()) %>% as.data.frame()
behDat = pivot_wider(data = behDat,
                      names_from ="Performance", values_from = "N")

#Correct Location Rate
behDat$HR = behDat$Hit/(behDat$Hit+behDat$Miss)
behDat$HR[behDat$HR==0] = .5/(behDat$Hit[behDat$HR==0]+behDat$Miss[behDat$HR==0])
behDat$HR[behDat$HR==1] = ((behDat$Hit[behDat$HR==1]+behDat$Miss[behDat$HR==1])-0.5)/(behDat$Hit[behDat$HR==1]+behDat$Miss[behDat$HR==1])

behDat$FAR = behDat$FA/(behDat$FA+behDat$CR)
behDat$FAR[behDat$FAR==0] = .5/(behDat$FA[behDat$FAR==0]+behDat$CR[behDat$FAR==0])
behDat$FAR[behDat$FAR==1] = ((behDat$FA[behDat$FAR==1]+behDat$CR[behDat$FAR==1])-0.5)/(behDat$FA[behDat$FAR==1]+behDat$CR[behDat$FAR==1])

behDat$dPrime = qnorm(behDat$HR)-qnorm(behDat$FAR)
write.csv(behDat,"BehavioralPerformance.csv",row.names = F)

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

behDat$pID = as.factor(behDat$pID)
ggplot(behDat,aes(x=pID,y=dPrime,fill = pID)) + 
  geom_bar(stat = "identity")+
  gDefault

graph2ppt(file= "dPrimeGraph.pptx", width = 8, height = 5)
