
# Adpated Gabe's script to make a function for the shiny app
FN_FP_Plot <-function(pc_concern =1e-5,
                      TCA_days = 5,
                      frags=c("PcBest","PcFrag10","PcFrag100")
                      
                      ) {
  
concernCount <- purrr::map(
  seq(from=1e-7,to=1e-5,by=1e-6),
  debrisConfusion,
  concernData_at_TOI=concernSummary_at_TOI,
  concernData=concernSummary,
  Pc_concern=pc_concern,
  fragVector=frags,
  days_to_TCA=TCA_days,
  run_monte_carlo=FALSE) %>% rbindlist() 

# Make Concern Prob and Frag size factors
fragLabel <- case_when(
  concernCount[,FragmentSize] == "PcBest" ~ ">= 1",
  concernCount[,FragmentSize] == "PcFrag10" ~ ">= 10",
  concernCount[,FragmentSize] == "PcFrag100" ~ ">= 100"
)
concernCount[,"FragLabel" := fragLabel]

# FP against FN plot ####
plot <- concernCount %>% 
  plot_ly(type="scatter",
          mode="lines",
          x=~FalseAlarms,
          y=~Missed,
          color=~FragLabel,
          text=~paste0("<b>Warning Threshold: </b>",WarnThreshold,"<br>",
                       "<b>Missed Concern Events: </b>",Missed,"<br>",
                       "<b>False Alarms: </b>",FalseAlarms,"<br>",
                       "<b>Warning Count: </b>",WarningCount,"<br>"),
          hoverinfo="text",
          colors=c("yellow","orange","#990000"),
          line=list(width=3.5)
  ) %>%
  
  layout(
    xaxis = list(title="<b>False Alarms (FP)</b>",
                 gridcolor="#333333"
    ),
    yaxis = list(title="<b>Missed Concern Events (FN)</b>",
                 gridcolor="#333333",
                 zeroline = TRUE,
                 zerolinecolor="#333333"),
    plot_bgcolor  = "#444444",
    paper_bgcolor = "#444444",
    font = list(color = '#FFFFFF')
  )

return(plot)}

