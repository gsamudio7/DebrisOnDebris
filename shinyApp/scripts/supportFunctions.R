
library(data.table)
library(dplyr)
library(repmis)
library(stringr)
library(plotly)
library(RColorBrewer)
library(tidyr)


# Function to download from dropbox ####
downloadFromDropBox <- function(dropbox_csv_url) {
  cat("downloading from dropbox\n")
  str_sub(dropbox_csv_url, nchar(dropbox_csv_url)) <- "1" # Make sure the last digit is 1
  return(repmis::source_data(dropbox_csv_url,header=FALSE) %>% data.table())
}


# Function to process data downloaded from dropbox ####
processDebris <- function(raw) {
  
  # Keep relevant columns
  cat("processing data\n")
  raw <- raw[,c(1,2,3,11,14,15,16,17,45,46,58)]
  
  # Name columns
  raw <- raw %>%
    setnames(c("V1","V2","V3",
               "V11","V14","V15",
               "V16","V17","V45",
               "V46","V58"),
             c("Primary","Secondary","MissDistance",
               "PcBest","PcFrag10","PcFrag100",
               "PcFrag1000","PcFrag10000",
               "time_of_screening","TCA",
               "eventNumber"))
  
  # Generate time to TCA column and remove redundant columns
  raw[,"time2TCA" := TCA - time_of_screening]
  
  # Determine scalar to transform quantities to yearly counts
  scale <- 365/raw[,uniqueN(round(time_of_screening))]
  
  # Gather observations and summarize 
  processedData <- raw[,c("PcBest","PcFrag10","PcFrag100",
              "PcFrag1000","PcFrag10000","eventNumber","time2TCA")] %>%
    gather("fragLabel","Pc",-eventNumber,-time2TCA) %>% data.table()
  suppressWarnings(processedData[,"Pc" := as.numeric(Pc)])
  
  # Remove observations with less than 1e-10 Collision Probability
  processedData <- processedData[Pc > 1e-10]
  
  # Bin days to TCA
  processedData[,"TCA_Bin" := as.factor(ceiling(time2TCA))]
  
  # Get the latest recorded Pc at each days to TCA bin
  processedData <- processedData[
    ,.(bool=time2TCA==max(time2TCA),Pc=Pc)
    ,by=c("TCA_Bin","fragLabel","eventNumber")][bool==TRUE][,!"bool"]
  
  # Aesthetic fragment size label
  processedData[,"fragLabel" := as.factor(case_when(
    fragLabel == "PcBest" ~ ">= 1",
    fragLabel == "PcFrag10" ~ ">= 10",
    fragLabel == "PcFrag100" ~ ">= 100",
    fragLabel == "PcFrag1000" ~ ">= 1000",
    fragLabel == "PcFrag10000" ~ ">= 10000",
  ))]

  # Clean up RAM and return
  rm(raw)
  return(list("data"=processedData %>% spread(TCA_Bin, Pc),
              "concernProbs"=processedData[Pc >= quantile(Pc,.75),
                                           .(`.75 Superquantile`=mean(Pc)),
                                           by=fragLabel],
              "scale"=scale))
}

# Debris Confusion Function ####
debrisConfusion <- function(Pc_warn,
                            tca_of_interest,
                            debrisInfo) {
  
  # Filter to days to TCA of interest
  df <- cbind(
    debrisInfo$data[,c("fragLabel",
                       "eventNumber",
                       "1")],
    debrisInfo$data[[as.character(tca_of_interest)]]) 
  df[,"Pc_at_TCA" := V2]
  
  return(purrr::map_dfr(
    .x=debrisInfo$data[,unique(fragLabel)],
    .f=function(frag) {
      
      # Get the Pc concern value for the fragLabel
      Pc_concern <- debrisInfo$concernProbs[
        fragLabel==frag,`.75 Superquantile`]

      # Count missed events (FN) and false alarms (FP) and warnings
      # Accumulate in a data.table
      return(data.table(
        "WarnThreshold"=Pc_warn %>% formatC(format="e",digits=2),
        "Missed"=round(df[fragLabel==frag & Pc_at_TCA < Pc_warn & 
                                            `1` >= Pc_concern,.N]*debrisInfo$scale),
        "FalseAlarms"=round(df[fragLabel==frag & Pc_at_TCA >= Pc_warn & 
                                                 `1` < Pc_concern,.N]*debrisInfo$scale),
        "WarningCount"=round(df[fragLabel==frag & Pc_at_TCA >= Pc_warn,.N]*debrisInfo$scale),
        "fragLabel"=frag,
        "TCA_Bin"=tca_of_interest
    ))
  }))
}


# Trade off plot function
trade_off_plot <- function(tcaBin=5,
                           test_thresholds=seq(from=1e-7,to=1e-5,by=1e-7),
                           missTolerance,
                           debrisInfo) {
  
  # Compute point estimate
  cat("\nComputing point estimates")
  concernCount <- purrr::map_dfr(
    test_thresholds,
    debrisConfusion,
    tca_of_interest=tcaBin,
    debrisInfo=debrisInfo) 
  
  # Compute optimal thresholds 
  cat("\nComputing optimal thresholds")
  optimal_thresholds <- purrr::map_dfr(
    .x=concernCount[,unique(fragLabel)],
    .f=function(x) {concernCount[fragLabel==x & Missed <= missTolerance[[x]]][
      FalseAlarms==min(FalseAlarms)][
        WarnThreshold==min(WarnThreshold),
        c("WarnThreshold","Missed","FalseAlarms","fragLabel","WarningCount","TCA_Bin")]}
  ) 
  
  # FP against FN plot ####
  cat("\nPlotting\n")
  trade_off <- suppressWarnings(
    concernCount %>% 
    plot_ly(type="scatter",
            mode="lines",
            colors=colorRampPalette(c("dodgerblue","orange"))(5)) %>%
    add_trace(
            x=~FalseAlarms,
            y=~Missed,
            color=~fragLabel,
            legendgroup=~fragLabel,
            text=~paste0("<b>Warning Threshold: </b>",WarnThreshold,"<br>",
                         "<b>Warning Count: </b>",round(WarningCount),"<br>",
                         "<b>Missed: </b>",round(Missed),"<br>",
                         "<b>False Alarms: </b>",round(FalseAlarms),"<br>"),
            hoverinfo="text",
            line=list(width=3.75)
    ) %>%
    
    layout(
      legend = list(title = list(text="<b>Fragment Size</b>")),
      xaxis = list(title="<b>False Alarms (FP)</b>",
                   gridcolor="#222222"),
      yaxis = list(title="<b>Missed Concern Events (FN)</b>",
                   gridcolor="#222222",
                   zeroline = TRUE,
                   zerolinecolor="#222222"),
      plot_bgcolor  = "#333333",
      paper_bgcolor = "#333333",
      font = list(color = '#FFFFFF')
    ) %>% 
    
    add_markers(
      data=optimal_thresholds,
      type="scatter",
      mode="markers+lines",
      marker=list(size=15,
                  color = 'rgba(17, 157, 255,0)',
                  line=list(color="mediumspringgreen",
                            width=3.25)),
      x=~FalseAlarms,
      y=~Missed,
      line=list(color='rgba(17, 157, 255,0)'),
      name=paste("<b>Miss Tolerances:</b><br>",
                 paste0(names(missTolerance)," | ",missTolerance,collapse="\n")),
      showlegend=TRUE
    )) 
  
  return(list("plot"=trade_off,
              "optThresholds"=optimal_thresholds))}

evalPerformance <- function(
  numSamples=100,
  optimalThresholds) {
  
  # Aggregate point estimates over 3,4,5
  cat("\nPreparing data\n")
  aggData <- optimalThresholds[,.(aggMissed=sum(Missed),
                    aggFalseAlarms=sum(FalseAlarms),
                    aggWarnings=sum(WarningCount)),by=fragLabel]
  
  # Run monte carlo using the optimal thresholds ####
  cat("Sampling\n")
  aggSim <- purrr::map_dfr(.x=1:numSamples,
                           .f=function(x) {
    
      # Sample
      sampleInfo <- debrisInfo
      sampleInfo[["simData"]] <- debrisInfo$data[sample(.N,.N,replace=TRUE)]
      
      # Debris confusion function
      sim <- apply(
        X=optimalThresholds[,c("WarnThreshold","TCA_Bin","fragLabel")],
        MARGIN=1,
        FUN=function(x) {
          sampleInfo[["data"]] <- sampleInfo[["simData"]][fragLabel==x[3]]
          debrisConfusion(Pc_warn=as.numeric(x[1]),
                                        tca_of_interest=x[2],
                                        debrisInfo=sampleInfo)}) %>% rbindlist()
      
      # Total the FNs and FPs for the year
      sim[,.(aggMissed=sum(Missed),
             aggFalseAlarms=sum(FalseAlarms),
             aggWarnings=sum(WarningCount)),by=fragLabel]
  })

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
  
  finalPerformance <- aggSim_and_Data[,.(
    `fragLabel`=fragLabel,
    `Missed`=paste0(aggMissed," [",round(missedLow)," - ",round(missedHigh),"]"),
    `False Alarms`=paste0(aggFalseAlarms," [",round(falseAlarmLow)," - ",round(falseAlarmHigh),"]"),
    `Warnings`=paste0(aggWarnings," [",round(warningsLow)," - ",round(warningsHigh),"]")
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
      colors=colorRampPalette(c("dodgerblue","orange"))(5),
      marker=list(size=12),
      hoverinfo="text",
      text=~paste("<b>fragLabel:</b> ",fragLabel,"<br>",
                  "<b>Total Missed:</b> ",aggMissed,"<br>",
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
          fillcolor = "#00008B", line=list(color="#00008B"),opacity = 0.3,
          x0 = aggSim_and_Data[fragLabel==">= 1",falseAlarmLow], 
          x1 = aggSim_and_Data[fragLabel==">= 1",falseAlarmHigh], xref = "x",
          y0 = aggSim_and_Data[fragLabel==">= 1",missedLow], 
          y1 = aggSim_and_Data[fragLabel==">= 1",missedHigh], yref = "y"),
        list(
          type="rect",
          fillcolor = "#0F48C4",line=list(color="#0F48C4"),opacity = 0.3,
          x0 = aggSim_and_Data[fragLabel==">= 10",falseAlarmLow], 
          x1 = aggSim_and_Data[fragLabel==">= 10",falseAlarmHigh], xref = "x",
          y0 = aggSim_and_Data[fragLabel==">= 10",missedLow], 
          y1 = aggSim_and_Data[fragLabel==">= 10",missedHigh], yref = "y"),
        list(
          type="rect",
          fillcolor = "#1E90FF", line=list(color="#1E90FF"),opacity = 0.3,
          x0 = aggSim_and_Data[fragLabel==">= 100",falseAlarmLow], 
          x1 = aggSim_and_Data[fragLabel==">= 100",falseAlarmHigh], xref = "x",
          y0 = aggSim_and_Data[fragLabel==">= 100",missedLow], 
          y1 = aggSim_and_Data[fragLabel==">= 100",missedHigh], yref = "y"),
        list(
          type="rect",
          fillcolor = "#78B1E9", line=list(color="#78B1E9"),opacity = 0.3,
          x0 = aggSim_and_Data[fragLabel==">= 1000",falseAlarmLow], 
          x1 = aggSim_and_Data[fragLabel==">= 1000",falseAlarmHigh], xref = "x",
          y0 = aggSim_and_Data[fragLabel==">= 1000",missedLow], 
          y1 = aggSim_and_Data[fragLabel==">= 1000",missedHigh], yref = "y"),
        list(
          type="rect",
          fillcolor = "#D3D3D3", line=list(color="#D3D3D3"),opacity = 0.3,
          x0 = aggSim_and_Data[fragLabel==">= 10000",falseAlarmLow], 
          x1 = aggSim_and_Data[fragLabel==">= 10000",falseAlarmHigh], xref = "x",
          y0 = aggSim_and_Data[fragLabel==">= 10000",missedLow], 
          y1 = aggSim_and_Data[fragLabel==">= 10000",missedHigh], yref = "y")
      ),
      legend = list(title = list(text="<b>fragLabel</b>")),
      xaxis = list(title="<b>False Alarms (FP)</b>",
                   gridcolor="#222222"),
      yaxis = list(title="<b>Missed Concern Events (FN)</b>",
                   gridcolor="#222222",
                   zeroline = TRUE,
                   zerolinecolor="#222222"),
      plot_bgcolor  = "#333333",
      paper_bgcolor = "#333333",
      font = list(color = '#FFFFFF')
    ) %>% plotly_build()
  
  return(list("plot"=final,
              "finalPerformance"=finalPerformance))
}


# Testing
url <- "https://www.dropbox.com/s/4si129wa5kou67i/DebrisOnDebris3.csv?dl=0"
raw <- downloadFromDropBox(dropbox_csv_url=url) 
debrisInfo <- raw %>% processDebris()
trade_off_plot(debrisInfo=debrisInfo,
               missTolerance=list(">= 1"=15,
                                  ">= 10"=10,
                                  ">= 100"=5,
                                  ">= 1000"=1,
                                  ">= 10000"=0))

concernCount_5 <- trade_off_plot(debrisInfo=debrisInfo,
                                 missTolerance=list(">= 1"=15,
                                                    ">= 10"=10,
                                                    ">= 100"=5,
                                                    ">= 1000"=1,
                                                    ">= 10000"=0),
                                 tcaBin=5)
concernCount_4 <- trade_off_plot(debrisInfo=debrisInfo,
                                 missTolerance=list(">= 1"=15,
                                                    ">= 10"=10,
                                                    ">= 100"=5,
                                                    ">= 1000"=1,
                                                    ">= 10000"=0),
                                 tcaBin=4)
concernCount_3 <- trade_off_plot(debrisInfo=debrisInfo,
                                 missTolerance=list(">= 1"=15,
                                                    ">= 10"=10,
                                                    ">= 100"=5,
                                                    ">= 1000"=1,
                                                    ">= 10000"=0),
                                 tcaBin=3)

optimalThresholds <- list(
  concernCount_5$optThresholds,
  concernCount_4$optThresholds,
  concernCount_3$optThresholds
) %>% rbindlist()

evalPerformance(numSamples=10000,
                optimalThresholds=optimalThresholds)
