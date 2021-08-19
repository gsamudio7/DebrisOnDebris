library(plotly)
library(tidyverse)
library(data.table)


# This code runs for ever!!
concernRates <- lapply(
  seq(from=1.94e-22,
      to=1e-5,
      by=1e-8),
  debrisConfusion,
  concernData_at_TOI=concernSummary_at_TOI,
  concernData=concernSummary,
  days_to_TCA=5) %>% rbindlist() 


concernRates_plot <- 
  concernRates %>% 
  plot_ly(type="scatter",
          mode="lines",
          x=~log(`False Positive`),
          y=~`False Negative`,
          text=~paste0("<b>Warning Threshold: </b>",Warn_Threshold,"<br>",
                       "<b>Missed Concern Events: </b>",Missed,"<br>",
                       "<b>False Alarms: </b>",FalseAlarms,"<br>",
                       "<b>Warning Count: </b>",Warning_Count,"<br>",
                       "<b>Concern Events: </b>",Concern_Events,"<br>",
                       "<b>Events: </b>",Events),
          hoverinfo="text",
          color=I("#2359c4"),
          #marker=list(size=3.5),
          line=list(width=3.5)
  ) %>%
  
  layout(
    xaxis = list(title="<b>False Alarm Rate (FP)</b>",
                 tickvals = seq(-10,-4,1),
                 ticktext = seq(-10,-4,1) %>% exp() %>% formatC(format="e",digits=2),
                 gridcolor="#333333"
    ),
    yaxis = list(title="<b>Missed Concern Event Rate (FN)</b>",
                 gridcolor="#333333",
                 zeroline = TRUE,
                 zerolinecolor="#333333"),
    plot_bgcolor  = "#444444",
    paper_bgcolor = "#444444",
    font = list(color = '#FFFFFF')
  )