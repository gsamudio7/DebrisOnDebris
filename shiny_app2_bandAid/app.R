## app.R ##
library(shinydashboard)
library(data.table)
library(shiny)
library(plotly)
library(dplyr)
library(plotly)
library(RColorBrewer)
library(DT)
library(purrr)

# Load data
load("data/debrisData.RData")
events <- fread("data/events.csv")
scale <- 365/events[,uniqueN(round(time_of_screening))]
concernProbs <- debrisData[Pc_min >= quantile(Pc_min,.75),
                           .(`Fragment Size`=fragLabel,
                             `.75 Super Quantile`=mean(Pc_min) %>% 
                                 formatC(format="e",digits=2)),by=fragLabel][,!"fragLabel"]

# Source functions
# Debris Confusion Function ####
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
    }
    return(rbindlist(result))
}

# Trade off plot function
trade_off_plot <- function(tcaBins=c(5),
                           test_thresholds=seq(from=1e-7,to=1e-5,by=1e-7),
                           missTolerance) {
    
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
        function(x) {concernCount[fragLabel==x & Missed <= missTolerance[[x]]][
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
                colors=colorRampPalette(c("blue4","dodgerblue","lightgray"))(5),
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

evalPerformance <- function(
    numSamples=100,
    optimalThresholdList) {
    
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
                    Pc_warn=optimalThresholdList[["5"]]$WarnThreshold[j] %>% as.numeric(),
                    fragVector=optimalThresholdList[["5"]]$fragLabel[j])
            }) %>% rbindlist(),
            
            lapply(1:5,function(j) {
                debrisConfusion(
                    df=debrisData[TCA_Bin == 4][sample(.N,nrow(debrisData[TCA_Bin == 4]),replace=TRUE)],
                    Pc_warn=optimalThresholdList[["4"]]$WarnThreshold[j] %>% as.numeric(),
                    fragVector=optimalThresholdList[["4"]]$fragLabel[j])
            }) %>% rbindlist(),
            
            lapply(1:5,function(j) {
                debrisConfusion(
                    df=debrisData[TCA_Bin == 3][sample(.N,nrow(debrisData[TCA_Bin == 3]),replace=TRUE)],
                    Pc_warn=optimalThresholdList[["3"]]$WarnThreshold[j] %>% as.numeric(),
                    fragVector=optimalThresholdList[["3"]]$fragLabel[j])
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
            colors=colorRampPalette(c("blue4","dodgerblue","lightgray"))(5),
            marker=list(size=12),
            hoverinfo="text",
            text=~paste("<b>Fragment Size:</b> ",fragLabel,"<br>",
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

# UI ####

ui <- dashboardPage(

    dashboardHeader(title = "USSPACECOM Space Conjuction",
                    titleWidth = 350),
    dashboardSidebar(
        HTML('<center><img src="SPACECOM_logo.png" width="150"></center>'),
        h4("Miss Tolerances"),
        p( paste("Space conjunction analysis in support of USSPACECOM. 
        Choose the parameters then press update.")),
        
        # Concern inputs
        textInput("input_1", "Fragments >= 1", 15),
        textInput("input_10", "Fragments >= 10", 10),
        textInput("input_100", "Fragments >= 100", 5),
        textInput("input_1000", "Fragments >= 1,000", 1),
        textInput("input_10000", "Fragments >= 10,000", 0),

        
        actionButton("MakeUpdates", HTML("<b>Check Performance</b>")),
        actionButton("monte", HTML("<b>Check Variability</b>"))
    ),
    dashboardBody(
        # add css custom to the app
        tags$head(
            tags$link(rel = "stylesheet", type = "text/css", href = "custom.css")
        ),
        # Boxes need to be put in a row (or column)
        fluidRow(
            box(
                h4("5 Days to TCA"),
                plotlyOutput("plot1", height = 430)),
            box(
                h4("4 Days to TCA"),
                plotlyOutput("plot2", height = 430))
        
            
            ),
            

        fluidRow(
            box(
                h4("3 Days to TCA"),
                plotlyOutput("plot3", height = 430)),
            box(
                h4("Thresholds"),
                DTOutput('tbl_final_rec', height = 430))
            ),
        
        fluidRow(
            box(
            h4("Variability"),
            plotlyOutput("plot4", height = 430)),
            
            box(h4("Performance"),
                DTOutput('tbl_final_per', height = 430))
        )
        )
    ) # end UI


server <- function(input, output) {
    
    # Get input
    mt_1 <- eventReactive(input$input_1, {as.numeric(input$input_1)})
    mt_10 <- eventReactive(input$input_10, {as.numeric(input$input_10)})
    mt_100 <- eventReactive(input$input_100, {as.numeric(input$input_100)})
    mt_1000 <- eventReactive(input$input_1000, {as.numeric(input$input_1000)})
    mt_10000 <- eventReactive(input$input_10000, {as.numeric(input$input_10000)})
    
    # Make plots
    concernCount_5 <- eventReactive(input$MakeUpdates, {
        trade_off_plot(tcaBins=c(5),missTolerance = list(
            ">= 1"=mt_1(),
            ">= 10"=mt_10(),
            ">= 100"=mt_100(),
            ">= 1000"=mt_1000(),
            ">= 10000"=mt_10000())
        )})

    concernCount_4 <- eventReactive(input$MakeUpdates, {
        trade_off_plot(tcaBins=c(4),missTolerance = list(
            ">= 1"=mt_1(),
            ">= 10"=mt_10(),
            ">= 100"=mt_100(),
            ">= 1000"=mt_1000(),
            ">= 10000"=mt_10000())
        )})
    
    concernCount_3 <- eventReactive(input$MakeUpdates, {
        trade_off_plot(tcaBins=c(3),missTolerance = list(
            ">= 1"=mt_1(),
            ">= 10"=mt_10(),
            ">= 100"=mt_100(),
            ">= 1000"=mt_1000(),
            ">= 10000"=mt_10000())
        )})
    
    final <- eventReactive(input$monte, {
        evalPerformance(numSamples=100,
                        optimalThresholdList=list(
                            "5"=concernCount_5()[["optThresholds"]],
                            "4"=concernCount_4()[["optThresholds"]],
                            "3"=concernCount_3()[["optThresholds"]]
                        ))
    })
    
    output$plot1 <- renderPlotly({
        concernCount_5()[["plot"]]
    })
    
    output$plot2 <- renderPlotly({
        concernCount_4()[["plot"]]
    })
    
    output$plot3 <- renderPlotly({
        concernCount_3()[["plot"]]
    })
    
    output$plot4 <- renderPlotly({
        final()[["plot"]]
    })
    
    output$tbl_final_rec <- renderDT({
        rbindlist(list(concernCount_5()[["optThresholds"]],
                       concernCount_4()[["optThresholds"]],
                       concernCount_3()[["optThresholds"]]))[
                           ,c("TCA_Bin","fragLabel","WarnThreshold")] %>%
            datatable(style="bootstrap",
                      options = list(searching = FALSE))
    })

    output$tbl_final_per <- renderDT({
        datatable(final()[["finalPerformance"]], style = 'bootstrap',options = list(searching = FALSE))
    })

} # end server function

shinyApp(ui, server)
