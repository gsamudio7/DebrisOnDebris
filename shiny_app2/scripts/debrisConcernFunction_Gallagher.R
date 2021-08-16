# Define function
library(tidyverse)
library(data.table)
# Read in trimmed data


debrisConcern <- function(Pc_concern=c(1e-3,1e-4,1e-5),
                          days_to_TCA=5,
                          threshold_Pc=1e-7,
                          number_of_fragments=1,
                          verbose=TRUE,
                          events = df) {
  
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
  
  # Initiate dataframe
  df <- data.frame(
    name = c("Cat_A", "Cat_B", "Cat_C"),
    true_neg = c(0,0,0),
    false_neg = c(0,0,0),
    false_pos = c(0,0,0),
    true_pos = c(0,0,0))
  
  # Initiate counter
  index = 0
  
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
    # if (verbose==TRUE) {
    #   cat(sprintf("%g events of concern\n",length(concernEvents)))}
    
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
    
    # # Print results
    # if (verbose==TRUE) {
    #   cat(sprintf(
    #     "\n\n\nConfusion Matrix:\ntrue positive rate: %g  false positive rate: %g\nfalse negative rate: %g  true negative rate: %g\n",
    #     warned[Concern==TRUE,.N]/totalEventCount,warned[Concern==FALSE,.N]/totalEventCount,
    #     notWarned[Concern==TRUE,.N]/totalEventCount,notWarned[Concern==FALSE,.N]/totalEventCount))
    # }
    
   
    
    
    
    
    # resultsList[[as.character(Pc)]] <- list("Pc_concern"=Pc,
    #                                         "threshold_Pc"=threshold_Pc,
    #                                         "number_of_fragments"=number_of_fragments,
    #                                         "truePositiveRate"=warned[Concern==TRUE,.N]/totalEventCount,
    #                                         "falsePositiveRate"=warned[Concern==FALSE,.N]/totalEventCount,
    #                                         "falseNegativeRate"=notWarned[Concern==TRUE,.N]/totalEventCount,
    #                                         "trueNegativeRate"=notWarned[Concern==FALSE,.N]/totalEventCount)
    index <- index + 1 
    df[index,2] <- notWarned[Concern==FALSE,.N] # true neg
    df[index,3] <- notWarned[Concern==TRUE,.N] # false neg
    df[index,4] <- warned[Concern==FALSE,.N] # false pos
    df[index,5] <- warned[Concern==TRUE,.N] #true pos
    
    cat(paste0(df[index,1], " Concern Category", "\n"))
    cat(paste0(df[index,2], " True Neg", "\n"))
    cat(paste0(df[index,3], " false Neg", "\n"))
    cat(paste0(df[index,4], " false pos", "\n"))
    cat(paste0(df[index,5], " True pos", "\n"))
    cat("\n")
  }
  
  return(df)
}

#debrisConcern()
