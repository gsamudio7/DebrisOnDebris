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
                                    "Count"=c(FN,FP,TP,TN),
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



# Function Demo ####
# Working directory must be DebrisOnDebris project

# # Test 
# rm(list=ls())
# 
# # Read in required data 
# load("data/concernData.RData")
# 
# # Read in debrisCounfuction function
# source("scripts/debrisConfusion.R")
# 
# # Execute function for the following inputs
# results <- debrisConfusion(concernData_at_TOI=concernSummary_at_TOI,
#                            concernData=concernSummary,
#                            Pc_warn=2.06e-9,
#                            days_to_TCA=5,
#                            Pc_concern=1e-5,
#                            frag_number_Pc="PcBest",
#                            verbose=FALSE,
#                            rate=FALSE)