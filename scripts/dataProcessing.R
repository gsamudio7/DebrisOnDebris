
library(data.table)
library(dplyr)
library(plotly)
library(RColorBrewer)

# Read in trimmed data
events <- fread("data/events.csv")

# Summarize conjunction events of concern for each FragNum
concernSummary <- rbindlist(list(
  
  # PcBest
  events[!is.na(PcBest),
    .(bool=time2TCA==min(time2TCA) & time2TCA < 1,
      fragNum="PcBest",
      Pc_min=PcBest),
    by=eventNumber][bool==TRUE,!"bool"],
  
  # PcFrag10
  events[!is.na(PcFrag10),
   .(bool=time2TCA==min(time2TCA) & time2TCA < 1,
     fragNum="PcFrag10",
     Pc_min=PcFrag10),
   by=eventNumber][bool==TRUE,!"bool"],
  
  # PcFrag100
  events[!is.na(PcFrag100),
   .(bool=time2TCA==min(time2TCA) & time2TCA < 1,
     fragNum="PcFrag100",
     Pc_min=PcFrag100),
   by=eventNumber][bool==TRUE,!"bool"],
  
  # PcFrag1000
  events[!is.na(PcFrag1000),
   .(bool=time2TCA==min(time2TCA) & time2TCA < 1,
     fragNum="PcFrag1000",
     Pc_min=PcFrag1000),
   by=eventNumber][bool==TRUE,!"bool"],
  
  # PcFrag10000
  events[!is.na(PcFrag10000),
   .(bool=time2TCA==min(time2TCA) & time2TCA < 1,
     fragNum="PcFrag10000",
     Pc_min=PcFrag10000),
   by=eventNumber][bool==TRUE,!"bool"]
))

concernSummary[,"fragNum" := as.factor(fragNum)]
concernSummary[,"Pc_min" := as.double(Pc_min)]

# How many are 0?
zeroCount <- merge(
  concernSummary[Pc_min==0,.(Zero_count=.N),by=fragNum],
  concernSummary[,.(Total=.N),by=fragNum],
  by="fragNum")
zeroCount$Proportion <- round(zeroCount$Zero_count/zeroCount$Total,3)
x <- zeroCount$fragNum
zeroCount[,"Fragment Size" := case_when(
  x == "PcBest" ~ ">= 1",
  x == "PcFrag10" ~ ">= 10",
  x == "PcFrag100" ~ ">= 100",
  x == "PcFrag1000" ~ ">= 1000",
  x == "PcFrag10000" ~ ">= 10000"
)]
zeroCount[,"Prob=0 Events" := Zero_count]
zeroCount <- zeroCount[,c("Fragment Size","Prob=0 Events","Proportion")]
save(zeroCount,file="products/zeroCount.RData")

# How many concern events at varying frag numbers values?
concernCount <- merge(
  concernSummary[Pc_min >= 1e-5,.(`Concern Events`=.N),by=fragNum],
  concernSummary[,.(Total=.N),by=fragNum]
)
concernCount$Proportion <- formatC(concernCount$`Concern Events`/concernCount$Total,
                                   format="e",digits=2)
x <- concernCount$fragNum
concernCount[,"Fragment Size" := case_when(
  x == "PcBest" ~ ">= 1",
  x == "PcFrag10" ~ ">= 10",
  x == "PcFrag100" ~ ">= 100",
  x == "PcFrag1000" ~ ">= 1000",
  x == "PcFrag10000" ~ ">= 10000"
)]

concernCount <- concernCount[,c("Fragment Size","Concern Events","Proportion")] 
save(concernCount,file="products/concernCount.RData")



# Explore distribution of Collision probabilities, for different fragment numbers
x <- concernSummary$fragNum
concernSummary[,"Fragment Size" := case_when(
  x == "PcBest" ~ ">= 1",
  x == "PcFrag10" ~ ">= 10",
  x == "PcFrag100" ~ ">= 100",
  x == "PcFrag1000" ~ ">= 1000",
  x == "PcFrag10000" ~ ">= 10000"
)]

Pc_at_1_day <- 
  concernSummary %>%
  plot_ly(
    type="histogram",
    x=~log(Pc_min),
    color=~`Fragment Size`,
    nbinsx=35
  ) %>%
  layout(xaxis = list(title="<b>Collision Probability</b>",
                      tickvals = seq(-720,-20,100),
                      ticktext = seq(-720,-20,100) %>% exp() %>% formatC(format="e",digits=2),
                      tickfont = list(size = 10),
                      tickangle = 90,
                      gridcolor="#333333"),
         yaxis = list(title="<b>Frequency</b>",
                      gridcolor="#333333"),
         legend = list(title = list(text = "<b>Fragment Size</b>")),
         plot_bgcolor  = "#444444",
         paper_bgcolor = "#444444",
         font = list(color = '#FFFFFF'),
         shapes = list(
           list(type= "line",
                line = list(color='#FFFFFF'),
                x0 = quantile(log(concernSummary$Pc_min),.95,na.rm=TRUE),
                x1 = quantile(log(concernSummary$Pc_min),.95,na.rm=TRUE),
                y0 = 0, y1 = 2695),
           
           list(type ="line",
                line = list(color='#FFFFFF'),
                x0 = log(1e-5), x1 = log(1e-5),
                y0 = 0, y1 = 2150),
           list(type ="line",
                line = list(color='#FFFFFF'),
                x0 = log(1e-100), x1 = log(1e-100),
                y0 = 0, y1 = 1900),
           list(type ="line",
                line = list(color='#FFFFFF'),
                x0 = log(1e-200), x1 = log(1e-200),
                y0 = 0, y1 = 1650)),
           
         annotations = list(
           list(text=paste("<b>.95 Quartile:</b><br>",
                           quantile(concernSummary$Pc_min,.95,na.rm=TRUE) %>%
                             formatC(format="e",digits=2)),
                x = log(1e-45),
                y = 2600, showarrow=FALSE),
           list(text = paste("<b>1e-5</b>"),
                x = log(1e-5), y = 2250, showarrow=FALSE),
           list(text = paste("<b>1e-100</b>"),
                x = log(1e-100), y = 2000, showarrow=FALSE),
           list(text = paste("<b>1e-200</b>"),
                x = log(1e-200), y = 1750, showarrow=FALSE))
  )
  

save(Pc_at_1_day,file="products/Pc_at_1_day.RData")

# How does the Pc change as we approach TCA
events[,"TCA_Bin" := as.factor(ceiling(time2TCA))]

# Include the Pc values at given days to TCA for each fragNum
# Summarize conjunction events of concern for each FragNum
concernSummary_at_TOI <- rbindlist(list(
  
  # PcBest
  events[!is.na(PcBest),
     .(bool=time2TCA==min(time2TCA),
       fragNum="PcBest",
       Pc_at_TOI=PcBest),
     by=c("TCA_Bin","eventNumber")][bool==TRUE,!"bool"],
  
  # PcFrag10
  events[!is.na(PcFrag10),
     .(bool=time2TCA==min(time2TCA),
       fragNum="PcFrag10",
       Pc_at_TOI=PcFrag10),
     by=c("TCA_Bin","eventNumber")][bool==TRUE,!"bool"],
  
  # PcFrag100
  events[!is.na(PcFrag100),
     .(bool=time2TCA==min(time2TCA),
       fragNum="PcFrag100",
       Pc_at_TOI=PcFrag100),
     by=c("TCA_Bin","eventNumber")][bool==TRUE,!"bool"],
  
  # PcFrag1000
  events[!is.na(PcFrag1000),
     .(bool=time2TCA==min(time2TCA),
       fragNum="PcFrag1000",
       Pc_at_TOI=PcFrag1000),
     by=c("TCA_Bin","eventNumber")][bool==TRUE,!"bool"],
  
  # PcFrag10000
  events[!is.na(PcFrag10000),
     .(bool=time2TCA==min(time2TCA),
       fragNum="PcFrag10000",
       Pc_at_TOI=PcFrag10000),
     by=c("TCA_Bin","eventNumber")][bool==TRUE,!"bool"]
))

concernSummary_at_TOI[,"fragNum" := as.factor(fragNum)]
concernSummary_at_TOI[,"Pc_at_TOI" := as.double(Pc_at_TOI)]

# Remove all data with Pr Collision <= 1e-10
concernSummary[Pc_min <= 1e-10,.N] # There's 237543
concernSummary <- concernSummary[Pc_min > 1e-10]
concernSummary_at_TOI <- concernSummary_at_TOI[Pc_at_TOI > 1e-10]

# Save and push to Git
save(concernSummary,concernSummary_at_TOI,file="data/concernData.RData")

# Only focus on PcBest, PbFrag10, and PcFrag100
concernSummary_at_TOI <- concernSummary_at_TOI[fragNum %in% c("PcBest","PcFrag10","PcFrag100")]

# Save and push to Git
save(concernSummary,concernSummary_at_TOI,file="data/concernData.RData")


