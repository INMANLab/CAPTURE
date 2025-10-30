########################### Initialization ##########################
rm(list=ls(all=TRUE))

x = readline()
Z:\Data\RealWorldNavigationCory\Blickshift Analytics Dynamic AOIs\AOI_Coordinates
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




datTemp = read.csv("RW4_4.csv",sep = ",")
# names(datTemp)[grepl("LM110.*",names(datTemp))]

dat = datTemp[,c("gaze.x..px.","gaze.y..px.","fixation.x..px..fixation.","fixation.y..px..fixation.",
                 "LM109.X", "LM109.Y", "LM109.Width", "LM109.Height",
                 "LM110.X", "LM110.Y", "LM110.Width", "LM110.Height")]
dat = dat[dat$LM109.X>0 | dat$LM110.X>0,  ]
write.csv(dat,"ExampleRW4_4.csv",row.names = F)
