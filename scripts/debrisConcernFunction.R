
library(data.table)
library(dplyr)
library(plotly)
library(RColorBrewer)

# Read in trimmed data
events <- fread("data/events.csv")

# Get the events of concern with Pc_concern = 1e-3
Pc_concern <- 1e-5

concernSummary <- events[!is.na(PcBest),
  .(bool=time2TCA==min(time2TCA) & time2TCA < 1,
    PcBest=PcBest),
  by=eventNumber][bool==TRUE,!"bool"] 

# Include the PcBest values at 5 days to TCA
concernSummary_at_5 <- events[!is.na(PcBest) &
  eventNumber %in% concernSummary$eventNumber & 
  time2TCA <= 5 & time2TCA > 4,.(PcBest_at_5=PcBest,
                                 bool=time2TCA==max(time2TCA)),by=eventNumber][bool==TRUE,!"bool"] %>%

  merge.data.table(concernSummary,by="eventNumber")


# Count total positives (events of concern) that we have visibility on at 5 days out
P_at_5 <- concernSummary_at_5[PcBest >= Pc_concern,.N] # 22
N_at_5 <- concernSummary_at_5[PcBest < Pc_concern,.N] # 29322

P <- concernSummary[PcBest >= Pc_concern,.N] # 58
N <- concernSummary[PcBest < Pc_concern,.N] # 48193

# Count false negatives (PcBest_at_5 < Pc_warn and PcBest >= Pc_concern) if we use naive approach (risk averse)
Pc_warn <- concernSummary_at_5[PcBest >= Pc_concern,min(PcBest_at_5)]

# Produce confusion matrix
debrisConfusion <- function(Pc_warn,Pc_concern=1e-5,verbose=FALSE) {
  TPR <- concernSummary_at_5[PcBest >= Pc_concern & PcBest_at_5 >= Pc_warn,.N]/P
  FPR <- concernSummary_at_5[PcBest < Pc_concern & PcBest_at_5 >= Pc_warn,.N]/N
  TNR <- concernSummary_at_5[PcBest < Pc_concern & PcBest_at_5 < Pc_warn,.N]/N
  FNR <- concernSummary_at_5[PcBest >= Pc_concern & PcBest_at_5 < Pc_warn,.N]/P
  if (verbose==TRUE) {
    cat(sprintf(
      "Pc_warn: %g\nPc_concern: %g\n\nTPR: %g | FPR: %g \nTNR: %g | FNR: %g",
      Pc_warn,Pc_concern,TPR,FPR,TNR,FNR)
    )
  }
  return(data.table("Warn_Threshold"=c(Pc_warn,Pc_warn),
                    "Rate"=c(FNR,FPR),
                    "Rate Type"=c("False Negative Rate",
                                  "False Positive Rate"))
  )
}

# Use min (risk averse)
debrisConfusion(Pc_warn=Pc_warn)

# Use max (risk tolerant)
debrisConfusion(Pc_warn=concernSummary_at_5[PcBest >= Pc_concern,max(PcBest_at_5)])

# Get data for a range of values
rates <- lapply(
  seq(from=1e-9,
      to=1e-4,
      by=1e-7),
  debrisConfusion) %>% rbindlist()


# Plot
rates %>%
  plot_ly(
    type="scatter",
    mode="lines",
    x=~log(Warn_Threshold),
    y=~Rate,
    color=~`Rate Type`
  )



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

