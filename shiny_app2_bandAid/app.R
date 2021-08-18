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
library(shinycssloaders)

# Load data
load("data/debrisData.RData")
events <- fread("data/events.csv")
scale <- 365/events[,uniqueN(round(time_of_screening))]
concernProbs <- debrisData[Pc_min >= quantile(Pc_min,.75),
                           .(`Fragment Size`=fragLabel,
                             `.75 Super Quantile`=mean(Pc_min) %>% 
                                 formatC(format="e",digits=2)),by=fragLabel][,!"fragLabel"]

# Source functions
source("scripts/supportFunctions.R")

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
        downloadButton('downloadData', HTML("<b>Download Thresholds</b>")),
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
                plotlyOutput("plot1", height = 430) %>% withSpinner(color="#2359c4")),
            box(
                h4("4 Days to TCA"),
                plotlyOutput("plot2", height = 430) %>% withSpinner(color="#2359c4"))
        
            
            ),
            

        fluidRow(
            box(
                h4("3 Days to TCA"),
                plotlyOutput("plot3", height = 430) %>% withSpinner(color="#2359c4")),
            box(
                h4("Thresholds"),
                DTOutput('tbl_final_rec', height = 430) %>% withSpinner(color="#2359c4"))
            ),
        
        fluidRow(
            box(
            h4("Variability"),
            plotlyOutput("plot4", height = 430) %>% withSpinner(color="#2359c4")),
            
            box(h4("Performance"),
                DTOutput('tbl_final_per', height = 430) %>% withSpinner(color="#2359c4"))
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
    
    output$tbl_final_rec <- renderDT(server=FALSE, {
        rbindlist(list(concernCount_5()[["optThresholds"]],
                       concernCount_4()[["optThresholds"]],
                       concernCount_3()[["optThresholds"]]))[
                           ,c("TCA_Bin","fragLabel","WarnThreshold")] %>%
            datatable(style="bootstrap",
                      options = list(searching = FALSE))
    })

    output$tbl_final_per <- renderDataTable({
        datatable(final()[["finalPerformance"]], 
                  style = 'bootstrap',
                  options = list(searching = FALSE))
    })
    
    output$downloadData <- downloadHandler(
          filename = function() {
            paste('data-', Sys.Date(), '.csv', sep='')
          },
          content = function(con) {
            fwrite(rbindlist(list(concernCount_5()[["optThresholds"]],
                                  concernCount_4()[["optThresholds"]],
                                  concernCount_3()[["optThresholds"]]))[
                                      ,c("TCA_Bin","fragLabel","WarnThreshold")], con)
          }
        )

} # end server function

shinyApp(ui, server)
