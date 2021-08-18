trade_off_plot <- function(tcaBins=c(5),
                           test_thresholds=c(seq(from=1e-8,to=1e-5,by=1e-8),1e-5),
                           runMonte=FALSE,
                           numSamples=NULL,
                           missTolerance=list(
                             ">= 1"=15,
                             ">= 10"=10,
                             ">= 100"=5,
                             ">= 1000"=1,
                             ">= 10000"=0)) {
  
  # Compute point estimate
  cat("\nComputing point estimates")
  concernCount <- purrr::map(
    test_thresholds,
    debrisConfusion,
    df=debrisData[TCA_Bin %in% tcaBins]) %>% rbindlist() 
  
  # Compute optimal thresholds 
  cat("\nComputing optimal thresholds")
  optimal_thresholds <- lapply(
    concernCount[,unique(fragLabel)],
    function(x) {concernCount[fragLabel==x & scale*Missed <= missTolerance[[x]]][
      FalseAlarms==min(FalseAlarms)][
        WarnThreshold==min(WarnThreshold),
        c("WarnThreshold","Missed","FalseAlarms","fragLabel","WarningCount","TCA_Bin")] %>%
        unique()}
  ) %>% rbindlist()
  
  # FP against FN plot ####
  cat("\nPlotting\n")
  trade_off <- 
    concernCount %>% 
    plot_ly(type="scatter",
            mode="lines",
            x=~FalseAlarms,
            y=~Missed,
            color=~fragLabel,
            legendgroup=~fragLabel,
            text=~paste0("<b>Warning Threshold: </b>",WarnThreshold,"<br>",
                         "<b>Warning Count: </b>",round(WarningCount),"<br>",
                         "<b>Missed: </b>",round(Missed),"<br>",
                         "<b>False Alarms: </b>",round(FalseAlarms),"<br>"),
            hoverinfo="text",
            colors=c("blue4","dodgerblue","lightgray"),
            line=list(width=3.25)
    ) %>%
    
    layout(
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
    ) %>%
    
    add_trace(
      data=optimal_thresholds,
      type="scatter",
      mode="markers+lines",
      marker=list(size=15,
                  color = 'rgba(17, 157, 255,0)',
                  line=list(color="mediumspringgreen",
                            width=3.25)),
      x=~FalseAlarms,
      y=~Missed,
      color='rgba(17, 157, 255,0)',
      line=list(color='rgba(17, 157, 255,0)'),
      name=paste("<b>Miss Tolerances:</b><br>",
                 paste0(names(missTolerance)," | ",missTolerance,collapse="\n")),
      showlegend=TRUE
    ) 
  
  return(list("plot"=trade_off,
              "optThresholds"=optimal_thresholds))}




