
library(data.table)
library(dplyr)
library(plotly)
library(RColorBrewer)

# Read in trimmed data
events <- fread("data/events.csv")

# Define function
debrisConcern <- function(Pc_concern=c(1e-3,1e-4,1e-5),
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
      
      events[time2TCA < 1 & 
               
      events[[fragList[[number_of_fragments %>% as.character()]]]] > Pc,
             
      unique(eventNumber)]
    
    nonConcernEvents <- 
      events[time2TCA < 1 & 
               
      events[[fragList[[number_of_fragments %>% as.character()]]]] <= Pc,
             
      unique(eventNumber)]
      
    totalEventCount <- length(concernEvents) + length(nonConcernEvents)
    
    # Verify:
    events[time2TCA < 1 & is.na(events[[fragList[[number_of_fragments %>% as.character()]]]]),
           uniqueN(eventNumber)] + totalEventCount == events[time2TCA < 1,uniqueN(eventNumber)]
    
    
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

