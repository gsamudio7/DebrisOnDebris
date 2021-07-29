
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
                      tickangle = 90),
         yaxis = list(title="<b>Frequency</b>"),
         hovermode = "x unified",
         
         shapes = list(
           list(type ="line",
                line = list(color="black"),
                x0 = log(1e-5), x1 = log(1e-5),
                y0 = 0, y1 = 2400),
           list(type ="line",
                line = list(color="black"),
                x0 = log(1e-100), x1 = log(1e-100),
                y0 = 0, y1 = 2150),
           list(type ="line",
                line = list(color="black"),
                x0 = log(1e-200), x1 = log(1e-200),
                y0 = 0, y1 = 1900)),
           
         annotations = list(
           list(text = paste("<b>1e-5</b>"),
                x = log(1e-5), y = 2500, showarrow=FALSE),
           list(text = paste("<b>1e-100</b>"),
                x = log(1e-100), y = 2250, showarrow=FALSE),
           list(text = paste("<b>1e-200</b>"),
                x = log(1e-200), y = 2000, showarrow=FALSE))
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






# How many zeros?
zeroCount <- merge(
  concernSummary_at_TOI[Pc_at_TOI==0,.(Zero_count=.N),by=c("TCA_Bin","fragNum")],
  concernSummary_at_TOI[,.(Total=.N),by=c("TCA_Bin","fragNum")],
  by=c("TCA_Bin","fragNum")
)
zeroCount[,"Proportion" := Zero_count/Total]
zeroCount <- zeroCount[,!"Total"]
zeroCount # Make a pretty plot


# Plot
Pc_by_TCA_plot <- concernSummary_at_TOI %>%
  plot_ly(
    type="box",
    x=~TCA_Bin,
    y=~log(PcBest_at_TOI)
  ) %>%
  layout(
    yaxis = list(title="<b>Probability of Collision</b>",
                 tickvals = seq(-700,0,100),
                 ticktext = seq(-700,0,100) %>% exp() %>% formatC(format="e",digits=2),
                 tickfont = list(size = 10))
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
    by="eventNumber"
  )
  
  P <- concernData[fragNum==frag_number_Pc & Pc_min >= Pc_concern,.N] 
  N <- concernData[fragNum==frag_number_Pc & Pc_min < Pc_concern,.N]
  
  if (rate==TRUE) {
    TPR <- df[Pc_at_TOI >= Pc_warn & Pc_min >= Pc_concern,.N]/P
    FPR <- df[Pc_at_TOI >= Pc_warn & Pc_min < Pc_concern,.N]/N
    TNR <- df[Pc_at_TOI < Pc_warn & Pc_min < Pc_concern,.N]/N
    FNR <- df[Pc_at_TOI < Pc_warn & Pc_min >= Pc_concern,.N]/P
    
    result_data_table <- data.table("Warn_Threshold"=c(Pc_warn,Pc_warn,Pc_warn,Pc_warn),
                                    "Rate"=c(FNR,FPR,TPR,TNR),
                                    "Rate Type"=c("False Negative Rate",
                                                  "False Positive Rate",
                                                  "True Positive Rate",
                                                  "True Negative Rate")
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
  seq(from=1e-9,
      to=1e-4,
      by=1e-6),
  debrisConfusion,
  concernData_at_TOI=concernSummary_at_TOI,
  concernData=concernSummary,
  days_to_TCA=5) %>% rbindlist() %>%

# Plot
  plot_ly(type="scatter",
          mode="lines",
          x=~log(Warn_Threshold),
          y=~Rate,
          color=~`Rate Type`) %>%
  layout(title = '<b>5 days to TCA and 1e-5 Concern Probability\nConcern Rates</b>',
         xaxis = list(title="<b>Warn Threshold</b>",
                      tickvals = c(-20,-17.5,-15,-12.5,-10),
                      ticktext = c(-20,-17.5,-15,-12.5,-10) %>% exp() %>% formatC(format="e",digits=2),
                      tickfont = list(size = 10),
                      tickangle=45),
         yaxis = list(title="<b>Rate</b>")
  )


save(concernRates,file="products/concernRates.RData")


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
