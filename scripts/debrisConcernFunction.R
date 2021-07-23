# Read in trimmed data
events <- fread("data/events.csv")

# Define function
debrisConcern <- function(Pc_concern=1e-5,
                          days_to_TCA=5,
                          threshold_Pc=1e-7,
                          number_of_fragments=1) {
  # Function that takes as input:
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
  
  # Make boolean column of collision or not
  # First get the eventNumbers that result in a collision:
  concernEvents <- 
    
    events[time2TCA < 1 & 
             
             events[[fragList[[number_of_fragments %>% as.character()]]]] > Pc_concern,
           
           unique(eventNumber)]
  
  # Print update
  cat(sprintf("Found %g events of concern\n",length(concernEvents)))
  
  # Create temporary boolean based on events of concern, so we know for all recorded conjunctions, 
  # for all days2TCA when that event results in a collision
  events[,"Concern" := as.factor(eventNumber %in% concernEvents)]
  
  # Get the event numbers that have at least the input days to TCA
  criticalEvents <- events[time2TCA > days_to_TCA,unique(eventNumber)]
  
  # Screen data to only these event numbers
  criticalData <- events[eventNumber %in% criticalEvents]
  
  # Print update
  cat(sprintf("Found %g events with at least %g days to TCA\n",
              length(criticalEvents), 
              days_to_TCA))
  
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
  cat(sprintf(
    "\n\n\nConfusion Matrix:\ntrue positive: %g  false positive: %g\nfalse negative: %g  true negative: %g\n",
    warned[Concern==TRUE,.N],warned[Concern==FALSE,.N],
    notWarned[Concern==TRUE,.N],notWarned[Concern==FALSE,.N]))
  
  
  return(list("truePositive"=warned[Concern==TRUE,.N],
              "falsePositive"=warned[Concern==FALSE,.N],
              "falseNegative"=notWarned[Concern==TRUE,.N],
              "trueNegative"=notWarned[Concern==FALSE,.N]))
}