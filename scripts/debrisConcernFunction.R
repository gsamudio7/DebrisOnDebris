
library(data.table)
library(dplyr)
library(plotly)
library(RColorBrewer)

# Read in trimmed data
events <- fread("data/events.csv")

# Get the events of concern with Pc_concern = 1e-5
Pc_concernVector <- c(1e-3,1e-4,1e-5)

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

# Include the PcBest values at given days to TCA
concernSummary_at_TOI <- events[!is.na(PcBest),
  .(PcBest_at_TOI=PcBest,
    bool=time2TCA==max(time2TCA)),
  by=c("TCA_Bin","eventNumber")][bool==TRUE,!"bool"] %>%

  merge.data.table(concernSummary[fragNum=="PcBest"],by="eventNumber")

# How many zeros?
zeroCount <- merge(
  concernSummary_at_TOI[PcBest_at_TOI==0,.(Zero_count=.N),by=TCA_Bin],
  concernSummary_at_TOI[,.(Total=.N),by=TCA_Bin],
  by="TCA_Bin")
zeroCount[,"Proportion" := Zero_count/Total]
zeroCount <- zeroCount[,!"Total"]
zeroCount

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
debrisConfusion <- function(Pc_warn,days_to_TCA,Pc_concern=1e-5,verbose=FALSE) {
  
  df <- concernSummary_at_TOI[TCA_Bin==days_to_TCA]
  
  TPR <- df[PcBest_at_TOI >= Pc_warn & Pc_min >= Pc_concern,.N]/P
  FPR <- df[PcBest_at_TOI >= Pc_warn & Pc_min < Pc_concern,.N]/N
  TNR <- df[PcBest_at_TOI < Pc_warn & Pc_min < Pc_concern,.N]/N
  FNR <- df[PcBest_at_TOI < Pc_warn & Pc_min >= Pc_concern,.N]/P
  if (verbose==TRUE) {
    cat(sprintf(
      "Pc_warn: %g\nPc_concern: %g\n\nTPR: %g | FPR: %g \nTNR: %g | FNR: %g",
      Pc_warn,Pc_concern,TPR,FPR,TNR,FNR)
    )
  }
  return(data.table("Warn_Threshold"=c(Pc_warn,Pc_warn,Pc_warn,Pc_warn),
                    "Rate"=c(FNR,FPR,TPR,TNR),
                    "Rate Type"=c("False Negative Rate",
                                  "False Positive Rate",
                                  "True Positive Rate",
                                  "True Negative Rate"))
  )
}

# Get data for a range of values
concernRates <- lapply(
  seq(from=1e-9,
      to=1e-4,
      by=1e-6),
  debrisConfusion,
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







# Define function
debrisConcern <- function(Pc_concern_vector=c(1e-3,1e-4,1e-5),
                          days_to_TCA=5,
                          threshold_Pc=1e-7,
                          number_of_fragments=1,
                          verbose=TRUE) {
  
  # Function that takes as input:
  ## Pc threshold of concern
  ## days to TCA
  ## Number of fragments
  ## Threshold Prob of Collision 
  
  # and gives output:
  ## number of events warned that resulted in a collision
  ## number of events warned that did not result in a collision
  ## number of events not warned that resulted in a collision
  ## number of events not warned that did not result in a collision
  
  
  # Look up table for Pc types
  fragList <- list("1"="PcBest",
                   "10"="PcFrag10",
                   "100"="PcFrag100",
                   "1000"="PcFrag1000",
                   "10000"="PcFrag10000")
  
  # Initiate results list
  resultsList <- list()
  
  # Iterate through Pc_concern values
  for (Pc in Pc_concern) {
  
    # Make boolean column of collision or not
    # First get the eventNumbers that result in a collision:
    concernEvents <- 
      events[min(time2TCA) < 1 & 
      events[[fragList[[number_of_fragments %>% as.character()]]]] <= Pc,
      eventNumber]
    
    nonConcernEvents <- 
      events[min(time2TCA) < 1 & 
      events[[fragList[[number_of_fragments %>% as.character()]]]] > Pc,
      eventNumber]
    
    # Print update
    if (verbose==TRUE) {
      cat(sprintf("%g events of concern\n",length(concernEvents)))}
    
    # Create temporary boolean based on events of concern, so we know for all recorded conjunctions, 
    # for all days2TCA when that event results in a collision
    events[,"Concern" := as.factor(eventNumber %in% concernEvents)]
    
    # Get the event numbers that have at least the input days to TCA
    criticalEvents <- events[time2TCA > days_to_TCA,
                             unique(eventNumber)]
    
    # Screen data to only these event numbers
    criticalData <- events[eventNumber %in% criticalEvents]
    
    # Get the data for events with a Pc above the given threshold on the given days to TCA
    warned <- 
      # On the given days to TCA
      criticalData[time2TCA > days_to_TCA & time2TCA <= round(days_to_TCA) + 1 &
                     
                     # Get the conjunctions above the threshold Pc
                     criticalData[[fragList[[number_of_fragments %>% as.character()]]]] > threshold_Pc]
    
    # Here get the opposite
    notWarned <- 
      criticalData[time2TCA > days_to_TCA & time2TCA <= round(days_to_TCA) + 1 &
                     criticalData[[fragList[[number_of_fragments %>% as.character()]]]] <= threshold_Pc]
    
    # Print results
    if (verbose==TRUE) {
      cat(sprintf(
        "\n\n\nConfusion Matrix:\ntrue positive rate: %g  false positive rate: %g\nfalse negative rate: %g  true negative rate: %g\n",
        warned[Concern==TRUE,.N]/totalEventCount,warned[Concern==FALSE,.N]/totalEventCount,
        notWarned[Concern==TRUE,.N]/totalEventCount,notWarned[Concern==FALSE,.N]/totalEventCount))
    }
    
    resultsList[[as.character(Pc)]] <- list("Pc_concern"=Pc,
                                            "threshold_Pc"=threshold_Pc,
                                            "number_of_fragments"=number_of_fragments,
                                            "truePositiveRate"=warned[Concern==TRUE,.N]/totalEventCount,
                                            "falsePositiveRate"=warned[Concern==FALSE,.N]/totalEventCount,
                                            "falseNegativeRate"=notWarned[Concern==TRUE,.N]/totalEventCount,
                                            "trueNegativeRate"=notWarned[Concern==FALSE,.N]/totalEventCount)
  }
    
  return(resultsList)
}

results <- debrisConcern(threshold_Pc=1e-9)

# Visualize false negative rate changing as you vary the warning threshold
warningThresholds <- c(1e-5,1e-10,1e-15)
falseNegData <- data.table(
  "PcConcern"=c(),
  "WarningTHreshold"=c(),
  "FragmentNum"=c(),
  "falseNegRate"=c()
)

for (warning in warningThresholds) {
  
  warningResult <- debrisConcern(threshold_Pc=warning)
  falseNegData
}

