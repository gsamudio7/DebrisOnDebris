
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
zeroCount[,"Proportion" := Zero_count/Total]
zeroCount <- zeroCount[,!"Total"]
zeroCount

# How many concern events at varying frag numbers values?
concernCount <- merge(
  concernSummary[Pc_min >= 1e-5,.(Concern_Event_Count=.N),by=fragNum],
  concernSummary[,.(Total=.N),by=fragNum]
)
concernCount[,"Proportion" := Concern_Event_Count/Total]
concernCount # Frag 10000 has no concern events


# Explore distribution of Collision probabilities, for different fragment numbers
Pc_at_1_day <- concernSummary %>%
  plot_ly(
    type="histogram",
    x=~log(Pc_min),
    color=~fragNum,
    nbinsx=35
  ) %>%
  layout(title = '<b>Collision Probability Distribution\nwhen TCA < 1 day</b>',
         xaxis = list(title="<b>Collision Probability</b>",
                      tickvals = seq(-720,-20,100),
                      ticktext = seq(-720,-20,100) %>% exp() %>% formatC(format="e",digits=2),
                      tickfont = list(size = 10),
                      tickangle = 90,
                      gridcolor="#333333"),
         yaxis = list(title="<b>Frequency</b>",
                      gridcolor="#333333"),
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
     .(bool=time2TCA==max(time2TCA),
       fragNum="PcBest",
       Pc_at_TOI=PcBest),
     by=c("TCA_Bin","eventNumber")][bool==TRUE,!"bool"],
  
  # PcFrag10
  events[!is.na(PcFrag10),
     .(bool=time2TCA==max(time2TCA),
       fragNum="PcFrag10",
       Pc_at_TOI=PcFrag10),
     by=c("TCA_Bin","eventNumber")][bool==TRUE,!"bool"],
  
  # PcFrag100
  events[!is.na(PcFrag100),
     .(bool=time2TCA==max(time2TCA),
       fragNum="PcFrag100",
       Pc_at_TOI=PcFrag100),
     by=c("TCA_Bin","eventNumber")][bool==TRUE,!"bool"],
  
  # PcFrag1000
  events[!is.na(PcFrag1000),
     .(bool=time2TCA==max(time2TCA),
       fragNum="PcFrag1000",
       Pc_at_TOI=PcFrag1000),
     by=c("TCA_Bin","eventNumber")][bool==TRUE,!"bool"],
  
  # PcFrag10000
  events[!is.na(PcFrag10000),
     .(bool=time2TCA==max(time2TCA),
       fragNum="PcFrag10000",
       Pc_at_TOI=PcFrag10000),
     by=c("TCA_Bin","eventNumber")][bool==TRUE,!"bool"]
))

concernSummary_at_TOI[,"fragNum" := as.factor(fragNum)]
concernSummary_at_TOI[,"Pc_at_TOI" := as.double(Pc_at_TOI)]

# Save and push to Git
save(concernSummary,concernSummary_at_TOI,file="data/concernData.RData")

# Start here ####

library(data.table)
library(dplyr)
library(plotly)
library(RColorBrewer)

load("data/concernData.RData")

# How many zeros?
zeroCount <- merge(
  concernSummary_at_TOI[Pc_at_TOI==0,.(Zero_count=.N),by=c("TCA_Bin","fragNum")],
  concernSummary_at_TOI[,.(Total=.N),by=c("TCA_Bin","fragNum")],
  by=c("TCA_Bin","fragNum")
)
zeroCount[,"Proportion" := Zero_count/Total]
zeroCount <- zeroCount[,!"Total"]
zeroCount # Make a pretty plot

# How does quantity change over time ####
concernSummary_at_TOI[fragNum=="PcBest" & Pc_at_TOI >= 1e-5,.N,by=TCA_Bin] %>% 
  plot_ly(
    type="bar",
    x=~TCA_Bin,
    y=~N,
    color=I("royalblue4")
  )




# How does quantity of NEW concern events change over time
concernEvents <- concernSummary_at_TOI[fragNum=="PcBest" & Pc_at_TOI >= 1e-5]

newConcernEvents <- list("10"=concernEvents[TCA_Bin==10,unique(eventNumber)])
for (i in (9:1)) {
  oldEvents <- concernEvents[TCA_Bin==i + 1,unique(eventNumber)] 
  newEvents <- concernEvents[TCA_Bin==i,unique(eventNumber)] 
  newConcernEvents[[as.character(i)]] = setdiff(newEvents,oldEvents)
}
newConcernSummary <- merge(
  data.table("TCA_Bin"=as.factor(10:1),
             "newCount"=unlist(lapply(newConcernEvents,length))),
  concernEvents[,.(concernCount=.N),by=TCA_Bin],
  by="TCA_Bin"
) 
  
newConcernSummary_plot <- newConcernSummary %>% 
  plot_ly(
  type="bar",
  x=~TCA_Bin,
  y=~concernCount,
  color=I("#222222"),
  name="Total"
) %>%
  
  add_trace(
    y=~newCount,
    name="New",
    color=I("royalblue3")
  ) %>%
  
  layout(
    xaxis=list(title="<b>Days to TCA Bin</b>",
               gridcolor="#333333"),
    yaxis=list(title="<b>Count</b>",
               gridcolor="#333333"),
    plot_bgcolor  = "#444444",
    paper_bgcolor = "#444444",
    font = list(color = '#FFFFFF')
  )



save(newConcernSummary_plot,file="products/newConcernSummary_plot.RData")


# Count the number of new events plus events that are carried through
'%ni%' <- Negate('%in%')
events_at_5 <- concernSummary_at_TOI[TCA_Bin==5,unique(eventNumber)]
length(events_at_5)
events_at_5_and_1 <- concernSummary_at_TOI[eventNumber %in% events_at_5 & TCA_Bin==1,unique(eventNumber)]
length(events_at_5_and_1)
events_at_1_and_not_at_5 <- concernSummary_at_TOI[eventNumber %ni% events_at_5 & TCA_Bin==1,
                                                  unique(eventNumber)]
length(events_at_1_and_not_at_5)

# Verify
length(events_at_5_and_1) + length(events_at_1_and_not_at_5) == concernSummary_at_TOI[TCA_Bin==1,uniqueN(eventNumber)]

# What is the distribution of Pc at 5 for events that make it to 1
concernSummary_at_TOI[fragNum=="PcBest" & 
                      eventNumber %in% events_at_5_and_1 & 
                      TCA_Bin == 5 |
                        
                      fragNum=="PcBest" & 
                      eventNumber %in% events_at_5_and_1 & 
                      TCA_Bin == 1] %>%
  plot_ly(
    type="box",
    x=~TCA_Bin,
    y=~log(Pc_at_TOI)
  ) 

# What is the distribution of the delta?
concernSummary_Pc_at_1_from_5 <- concernSummary_at_TOI[
  fragNum=="PcBest" & 
  eventNumber %in% events_at_5_and_1 & 
  TCA_Bin == 1 |
  
  fragNum=="PcBest" & 
  eventNumber %in% events_at_5_and_1 & 
  TCA_Bin == 5] %>% setorder(eventNumber,-TCA_Bin)

concernSummary_diff_at_1_from_5 <- concernSummary_Pc_at_1_from_5[
  ,.(Delta=diff(Pc_at_TOI)),by=eventNumber] 

# How many 0
concernSummary_diff_at_1_from_5[Delta==0,.N]/concernSummary_diff_at_1_from_5[,.N] # 0.64

# How many increased
concernSummary_diff_at_1_from_5[Delta > 0,.N]/concernSummary_diff_at_1_from_5[,.N] # 0.085

# How many decreased
concernSummary_diff_at_1_from_5[Delta < 0,.N]/concernSummary_diff_at_1_from_5[,.N] 

# Box Plot
concernSummary_Pc_at_1_from_5 %>%
  plot_ly(
    type="box",
    y=~log(Pc_at_TOI),
    x=~as.integer(TCA_Bin) - 1
  ) %>%
  layout(
    yaxis = list(title="<b>Probability of Collision</b>",
                 tickvals = seq(-700,0,100),
                 ticktext = seq(-700,0,100) %>% exp() %>% formatC(format="e",digits=2),
                 tickfont = list(size = 10)),
    xaxis = list(title="<b>Days to TCA</b>")
  )

concernSummary_Pc_at_1_from_5 %>%
  plot_ly(
    type="scatter",
    mode="lines",
    y=~log(Pc_at_TOI),
    x=~as.integer(TCA_Bin) - 1
  ) 

# Conclusion
## No obvious advantage in having a higher Pc_warn threshold at 5 days than the Pc_concern rate

# Question remains, what Pc_warn value at 5 days is the ideal to use to 
## minimize False negatives and minimize False Positives




# Plot
Pc_by_TCA_plot <- concernSummary_at_TOI[fragNum=="PcBest"] %>%
  plot_ly(
    type="box",
    y=~log(Pc_at_TOI),
    x=~TCA_Bin,
    hoverinfo="none"
  ) %>%
  layout(
    yaxis = list(title="<b>Probability of Collision</b>",
                 tickvals = seq(-700,0,100),
                 ticktext = seq(-700,0,100) %>% exp() %>% formatC(format="e",digits=2),
                 tickfont = list(size = 10),
                 gridcolor="#333333"),
    xaxis = list(title="<b>Days to TCA</b>",
                 gridcolor="#333333"),
    plot_bgcolor  = "#444444",
    paper_bgcolor = "#444444",
    font = list(color = '#FFFFFF')
  )
  
save(Pc_by_TCA_plot,file="products/Pc_by_TCA_plot.RData")
  


# Count total positives (events of concern) that we have visibility on at 5 days out
Pc_concern <- 1e-5
(P_at_5 <- concernSummary_at_TOI[TCA_Bin==5 & PcBest_at_TOI >= Pc_concern,.N]) # 21
(N_at_5 <- concernSummary_at_TOI[TCA_Bin==5 & PcBest_at_TOI < Pc_concern,.N]) # 29323

(P <- concernSummary[fragNum=="PcBest" & Pc_min >= Pc_concern,.N]) # 58
(N <- concernSummary[fragNum=="PcBest" & Pc_min < Pc_concern,.N]) # 48193

# Produce confusion matrix
debrisConfusion <- function(concernData_at_TOI,concernData,
                            Pc_warn,
                            days_to_TCA,
                            Pc_concern=1e-5,
                            frag_number_Pc="PcBest",
                            verbose=FALSE,
                            rate=TRUE) {
  
  df <- merge(
    concernData_at_TOI[TCA_Bin==days_to_TCA & fragNum==frag_number_Pc,c("eventNumber","Pc_at_TOI")],
    concernData[fragNum==frag_number_Pc,c("eventNumber","Pc_min")],
    by="eventNumber",
    all=FALSE
  )
  
  P <- concernData[Pc_min >= Pc_concern,.N] 
  N <- concernData[Pc_min < Pc_concern,.N]
  
  if (rate==TRUE) {
    #TP <- df[Pc_at_TOI >= Pc_warn & Pc_min >= Pc_concern,.N]
    FP <- df[Pc_at_TOI >= Pc_warn & Pc_min < Pc_concern,.N]
    #TN <- df[Pc_at_TOI < Pc_warn & Pc_min < Pc_concern,.N]
    FN <- df[Pc_at_TOI < Pc_warn & Pc_min >= Pc_concern,.N]
    warningCount <- df[Pc_at_TOI >= Pc_warn,.N]
      
    result_data_table <- data.table("Warn_Threshold"=Pc_warn %>% formatC(format="e",digits=2),
                                    "False Negative"=FN/P,
                                    "False Positive"=FP/N,
                                    "Missed"=FN,
                                    "FalseAlarms"=FP,
                                    "Warning_Count"=warningCount,
                                    "Concern_Events"=P,
                                    "Events"=nrow(df)
                                    
    )
    
  } else {
    TP <- df[Pc_at_TOI >= Pc_warn & Pc_min >= Pc_concern,.N]
    FP <- df[Pc_at_TOI >= Pc_warn & Pc_min < Pc_concern,.N]
    TN <- df[Pc_at_TOI < Pc_warn & Pc_min < Pc_concern,.N]
    FN <- df[Pc_at_TOI < Pc_warn & Pc_min >= Pc_concern,.N]
    
    result_data_table <- data.table("Warn_Threshold"=c(Pc_warn,Pc_warn,Pc_warn,Pc_warn),
                                    "Count"=c(FNR,FPR,TPR,TNR),
                                    "Type"=c("False Negative",
                                             "False Positive",
                                             "True Positive",
                                             "True Negative")
    )
  }
  if (verbose==TRUE) {
    cat(sprintf(
      "Pc_warn: %g\nPc_concern: %g\n\nTPR: %g | FPR: %g \nTNR: %g | FNR: %g",
      Pc_warn,Pc_concern,TPR,FPR,TNR,FNR)
    )
  }
  return(result_data_table)
}

# Get data for a range of values
concernRates <- lapply(
  seq(from=1.94e-22,
      to=1e-5,
      by=1e-8),
  debrisConfusion,
  concernData_at_TOI=concernSummary_at_TOI,
  concernData=concernSummary,
  days_to_TCA=5) %>% rbindlist() 

# Plot
concernRates_plot <- 

concernRates %>%
  plot_ly(type="scatter",
          mode="lines",
          x=~Warn_Threshold,
          y=~Rate,
          text=~paste0("<b>Warning Threshold: </b>",Warn_Threshold,"<br>",
                       "<b>Missed Concern Events: </b>",Missed,"<br>",
                       "<b>False Alarms: </b>",FalseAlarms,"<br>",
                       "<b>Warning Count: </b>",Warning_Count,"<br>",
                       "<b>Concern Events: </b>",Concern_Events,"<br>",
                       "<b>Events: </b>",Events),
          hoverinfo="text",
          color=~`Rate Type`) %>%
  
  layout(xaxis = list(title="<b>Warn Threshold</b>",
                      tickvals = seq(0,10e-6,2e-6),
                      ticktext = seq(0,10e-6,2e-6) %>% formatC(format="e",digits=2),
                      tickfont = list(size = 13),
                      tickangle=45,
                      gridcolor="#333333"),
         yaxis = list(title="<b>Confusion Rate</b>",
                      gridcolor = "#333333"),
         plot_bgcolor  = "#444444",
         paper_bgcolor = "#444444",
         font = list(color = '#FFFFFF')
  )





save(concernRates_plot,file="products/concernRates_plot.RData")





# FP against FN plot
concernRates_plot <- 
concernRates %>% 
  plot_ly(type="scatter",
          mode="lines+markers",
          x=~log(`False Positive`),
          y=~`False Negative`,
          text=~paste0("<b>Warning Threshold: </b>",Warn_Threshold,"<br>",
                       "<b>Missed Concern Events: </b>",Missed,"<br>",
                       "<b>False Alarms: </b>",FalseAlarms,"<br>",
                       "<b>Warning Count: </b>",Warning_Count,"<br>",
                       "<b>Concern Events: </b>",Concern_Events,"<br>",
                       "<b>Events: </b>",Events),
          hoverinfo="text",
          marker=list(size=3.5),
          line=list(width=1)
    ) %>%
  
  layout(
    xaxis = list(title="<b>False Alarm Rate (FP)</b>",
                 tickvals = seq(-10,-4,1),
                 ticktext = seq(-10,-4,1) %>% exp() %>% formatC(format="e",digits=2),
                 gridcolor="#333333"
            ),
    yaxis = list(title="<b>Missed Concern Event Rate (FN)</b>",
                 gridcolor="#333333",
                 zeroline = TRUE,
                 zerolinecolor="#333333"),
    plot_bgcolor  = "#444444",
    paper_bgcolor = "#444444",
    font = list(color = '#FFFFFF')
  )


# Function Demo ####

# Test 
rm(list=ls())

# Read in required data 
load("data/concernData.RData")

# Read in debrisCounfuction function
source("scripts/debrisConfusion.R")

# Execute function for the following inputs
results <- debrisConfusion(concernData_at_TOI=concernSummary_at_TOI,
                           concernData=concernSummary,
                           Pc_warn=2.06e-9,
                           days_to_TCA=5,
                           Pc_concern=1e-5,
                           frag_number_Pc="PcBest",
                           verbose=FALSE,
                           rate=FALSE)
