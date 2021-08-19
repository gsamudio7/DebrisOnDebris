evalPerformance <- function(
  numSamples=1000,
  optimalThresholdList=list(
    concernCount_5$optThresholds,
    concernCount_4$optThresholds,
    concernCount_3$optThresholds
  )) {
  
  # Aggregate point estimates over 3,4,5
  cat("\nPreparing data\n")
  agg <- rbindlist(optimalThresholdList)
  aggData <- agg[,.(aggMissed=sum(Missed),
                    aggFalseAlarms=sum(FalseAlarms),
                    aggWarnings=sum(WarningCount)),by=fragLabel]
  
  # Run monte carlo using the optimal thresholds ####
  cat("Sampling\n")
  aggSim <- data.table()
  for (i in 1:numSamples) {
    
    # Sample and compute false alarm and missed with optimal warning thresholds 
    aggMonte <- rbindlist(list(
      lapply(1:5,function(j) {
        debrisConfusion(
          df=debrisData[TCA_Bin == 5][sample(.N,nrow(debrisData[TCA_Bin == 5]),replace=TRUE)],
          Pc_warn=concernCount_5$optThresholds$WarnThreshold[j] %>% as.numeric(),
          fragVector=concernCount_5$optThresholds$fragLabel[j])
      }) %>% rbindlist(),
      
      lapply(1:5,function(j) {
        debrisConfusion(
          df=debrisData[TCA_Bin == 4][sample(.N,nrow(debrisData[TCA_Bin == 4]),replace=TRUE)],
          Pc_warn=concernCount_4$optThresholds$WarnThreshold[j] %>% as.numeric(),
          fragVector=concernCount_4$optThresholds$fragLabel[j])
      }) %>% rbindlist(),
      
      lapply(1:5,function(j) {
        debrisConfusion(
          df=debrisData[TCA_Bin == 3][sample(.N,nrow(debrisData[TCA_Bin == 3]),replace=TRUE)],
          Pc_warn=concernCount_3$optThresholds$WarnThreshold[j] %>% as.numeric(),
          fragVector=concernCount_3$optThresholds$fragLabel[j])
      }) %>% rbindlist()
    ))
    
    # Organize results
    aggSim <- rbindlist(list(aggSim,
                             aggMonte[,.(aggMissed=sum(Missed),
                                         aggFalseAlarms=sum(FalseAlarms),
                                         aggWarnings=sum(WarningCount)),by=fragLabel]))
    
  }
  
  # Calculate 95% CI
  cat("Computing 95% CI\n")
  aggSim_and_Data <- merge(
    aggSim[,.(missedHigh=quantile(aggMissed,.975),
              missedLow=quantile(aggMissed,.025),
              falseAlarmHigh=quantile(aggFalseAlarms,.975),
              falseAlarmLow=quantile(aggFalseAlarms,.025),
              warningsHigh=quantile(aggWarnings,.975),
              warningsLow=quantile(aggWarnings,.025)),by=fragLabel],
    aggData,
    "fragLabel")
  
  # Organize results
  finalRec <- merge(aggSim_and_Data,agg,by="fragLabel")[
    ,.(`Fragment Size`=fragLabel,
       `Days to TCA`=TCA_Bin,
       `Warn Threshold`=WarnThreshold),by=fragLabel]
  
  finalRec <- finalRec[,!"fragLabel"]
  
  finalPerformance <- aggSim_and_Data[,.(
    `Fragment Size`=fragLabel,
    `Missed`=paste0(aggMissed," [",missedLow," - ",missedHigh,"]"),
    `False Alarms`=paste0(aggFalseAlarms," [",falseAlarmLow," - ",falseAlarmHigh,"]"),
    `Warnings`=paste0(aggWarnings," [",warningsLow," - ",warningsHigh,"]")
  )]
  
  # Plot 
  cat("Plotting\n")
  final <- 
    aggSim_and_Data %>%
    plot_ly(
      type="scatter",
      mode="markers",
      x=~aggFalseAlarms,
      y=~aggMissed,
      marker=list(size=15,
                  line=list(color="#FFFFFF",
                            width=3.25)),
      color=~fragLabel,
      colors=c("blue","dodgerblue","gray"),
      marker=list(size=12),
      hoverinfo="text",
      text=~paste("<b>Total Missed:</b> ",aggMissed,"<br>",
                  "<b>95% CI:</b> [",round(missedLow)," - ",round(missedHigh),"]<br><br>", 
                  "<b>Total False Alarms:</b> ",aggFalseAlarms,"<br>",
                  "<b>95% CI:</b> [",round(falseAlarmLow)," - ",round(falseAlarmHigh),"]<br><br>", 
                  "<b>Total Warnings:</b>",aggWarnings,"<br>",
                  "<b>95% CI:</b> [",round(warningsLow)," - ",round(warningsHigh),"]")
    ) %>%
    
    layout(
      shapes=list(
        list(
          type="rect",
          fillcolor = "blue", line=list(color="#FFFFFF"),opacity = 0.3,
          x0 = aggSim_and_Data[fragLabel==">= 1",falseAlarmLow], 
          x1 = aggSim_and_Data[fragLabel==">= 1",falseAlarmHigh], xref = "x",
          y0 = aggSim_and_Data[fragLabel==">= 1",missedLow], 
          y1 = aggSim_and_Data[fragLabel==">= 1",missedHigh], yref = "y"),
        list(
          type="rect",
          fillcolor = "dodgerblue",line=list(color="#FFFFFF"),opacity = 0.3,
          x0 = aggSim_and_Data[fragLabel==">= 10",falseAlarmLow], 
          x1 = aggSim_and_Data[fragLabel==">= 10",falseAlarmHigh], xref = "x",
          y0 = aggSim_and_Data[fragLabel==">= 10",missedLow], 
          y1 = aggSim_and_Data[fragLabel==">= 10",missedHigh], yref = "y"),
        list(
          type="rect",
          fillcolor = "gray", line=list(color="#FFFFFF"),opacity = 0.3,
          x0 = aggSim_and_Data[fragLabel==">= 100",falseAlarmLow], 
          x1 = aggSim_and_Data[fragLabel==">= 100",falseAlarmHigh], xref = "x",
          y0 = aggSim_and_Data[fragLabel==">= 100",missedLow], 
          y1 = aggSim_and_Data[fragLabel==">= 100",missedHigh], yref = "y")
      ),
      legend = list(title = list(text="<b>Fragment Size</b>")),
      xaxis = list(title="<b>False Alarms (FP)</b>",
                   gridcolor="#333333"),
      yaxis = list(title="<b>Missed Concern Events (FN)</b>",
                   gridcolor="#333333",
                   zeroline = TRUE,
                   zerolinecolor="#333333"),
      plot_bgcolor  = "#444444",
      paper_bgcolor = "#444444",
      font = list(color = '#FFFFFF')
    ) %>% plotly_build()
  
  return(list("plot"=final,
              "finalRec"=finalRec,
              "finalPerformance"=finalPerformance))
}