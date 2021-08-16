debrisConfusion <- function(concernData_at_TOI,concernData,
                            Pc_warn,
                            days_to_TCA,
                            Pc_concern=1e-5,
                            fragVector=c("PcBest","PcFrag10","PcFrag100"),
                            verbose=FALSE,
                            run_monte_carlo=FALSE
) {
  
  # Initiate result list
  result <- list()
  
  # Prepare data
  df <- merge(
    concernData_at_TOI[TCA_Bin==days_to_TCA,c("eventNumber","Pc_at_TOI","fragNum")],
    concernData[,c("eventNumber","Pc_min","fragNum")],
    by=c("eventNumber","fragNum")
  )
  
  if (run_monte_carlo==TRUE) {
    df <- df[sample(.N, 180,replace=TRUE)]}
  
  # Loop over frag numbers
  for (frag in fragVector) {
    
    # Count missed events (FN) and false alarms (FP) and warnings
    # Accumulate in a data.table
    result[[paste0(frag)]] <- data.table(
      "WarnThreshold"=Pc_warn %>% formatC(format="e",digits=2),
      "Missed"=df[fragNum==frag & Pc_at_TOI < Pc_warn & Pc_min >= Pc_concern,.N],
      "FalseAlarms"=df[fragNum==frag & Pc_at_TOI >= Pc_warn & Pc_min < Pc_concern,.N],
      "WarningCount"=df[fragNum==frag & Pc_at_TOI >= Pc_warn,.N],
      "FragmentSize"=frag
    )
  }
  return(rbindlist(result))
}