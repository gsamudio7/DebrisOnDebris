debrisConfusion <- function(df,
                            Pc_warn,
                            concernProbData=concernProbs,
                            fragVector=debrisData[,unique(fragLabel)],
                            verbose=FALSE
) {
  
  # Initiate result list
  result <- list()
  
  # Loop over frag numbers
  for (frag in fragVector) {
    
    # Get the Pc concern value for the fragment size
    Pc_concern <- concernProbs[`Fragment Size`==frag,as.numeric(`.75 Super Quantile`)]
    
    print("hi charles ")
    print(class(Pc_concern))
   
    
    # Count missed events (FN) and false alarms (FP) and warnings
    # Accumulate in a data.table
    result[[paste0(frag)]] <- data.table(
      "WarnThreshold"=Pc_warn %>% formatC(format="e",digits=2),
      "Missed"=round(df[fragLabel==frag & Pc_at_TOI < Pc_warn & Pc_min >= Pc_concern,.N]*scale),
      "FalseAlarms"=round(df[fragLabel==frag & Pc_at_TOI >= Pc_warn & Pc_min < Pc_concern,.N]*scale),
      "WarningCount"=round(df[fragLabel==frag & Pc_at_TOI >= Pc_warn,.N]*scale),
      "fragLabel"=frag,
      "TCA_Bin"=df[,unique(TCA_Bin)]
    )
    
    print("did it work?????????????????????????????????")
  }
  return(rbindlist(result))
}


# the old function

# debrisConfusion <- function(df,
#                             Pc_warn,
#                             Pc_concern=1e-5,
#                             fragVector=debrisData[,unique(fragLabel)],
#                             verbose=FALSE
# ) {
#   print("hello charles")
#   print(fragVector)
#   # Initiate result list
#   result <- list()
#   
#   # Loop over frag numbers
#   for (frag in fragVector) {
#     
#     # Count missed events (FN) and false alarms (FP) and warnings
#     # Accumulate in a data.table
#     result[[paste0(frag)]] <- data.table(
#       "WarnThreshold"=Pc_warn %>% formatC(format="e",digits=2),
#       "Missed"=round(df[fragLabel==frag & Pc_at_TOI < Pc_warn & Pc_min >= Pc_concern,.N]*scale),
#       "FalseAlarms"=round(df[fragLabel==frag & Pc_at_TOI >= Pc_warn & Pc_min < Pc_concern,.N]*scale),
#       "WarningCount"=round(df[fragLabel==frag & Pc_at_TOI >= Pc_warn,.N]*scale),
#       "fragLabel"=frag,
#       "TCA_Bin"=df[,unique(TCA_Bin)]
#     )
#   }
#   return(rbindlist(result))
# }
