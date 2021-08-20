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


# Source functions
source("scripts/supportFunctions.R")


# UI ####
ui <- dashboardPage(

    dashboardHeader(title = "USSPACECOM Space Debris Risk Analysis Tool",
                    titleWidth = 500),
    dashboardSidebar(
        HTML('<center><img src="SPACECOM_logo.png" width="150"></center>'),
        h4("Load data from dropbox"),
        textInput("url","Paste Dropbox data share url"),
        actionButton("loadData","Load data"),
        fluidRow(column(12, style = "margin-left: 13px; font-weight:bold",
                   textOutput("dataStatus") %>% withSpinner(color="orange"))),
        hr(),
        h4("Miss Tolerances"),
        p("Choose the minimum number of yearly missed events (False Negatives) required for each debris on debris collision type."),
        
        # Concern inputs
        textInput("input_1", ">= 1", 15),
        textInput("input_10", ">= 10", 10),
        textInput("input_100", ">= 100", 5),
        textInput("input_1000", ">= 1,000", 1),
        textInput("input_10000", ">= 10,000", 0),

        actionButton("MakeUpdates", HTML("<b>Check Performance</b>")),
        actionButton("monte", HTML("<b>Check Variability</b>")),
        
        fluidRow(column(3, style = "margin-left: 12px;",
                        downloadButton("downloadData", HTML("<b>Download Thresholds</b>"),
                           style="color: #fff; background-color: #0645ad; border-color: #fff;"))
        )
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
                plotlyOutput("plot1", height = 430) %>% withSpinner(color="orange")),
            box(
                h4("4 Days to TCA"),
                plotlyOutput("plot2", height = 430) %>% withSpinner(color="orange"))
        
            
            ),
            

        fluidRow(
            box(
                h4("3 Days to TCA"),
                plotlyOutput("plot3", height = 430) %>% withSpinner(color="orange")),
            box(
                h4("Recommended Thresholds"),
                DTOutput('tbl_final_rec', height = 430) %>% withSpinner(color="orange"))
            ),
        
        fluidRow(
            box(
            h4("Variability"),
            plotlyOutput("plot4", height = 430) %>% withSpinner(color="orange")),
            
            box(h4("Performance"),
                DTOutput('tbl_final_per', height = 430) %>% withSpinner(color="orange"))
        )
        )
    ) # end UI

# Server ####
server <- function(input, output) {
    
    # Pull data
    debrisInfo <- reactiveVal()
    observeEvent(input$loadData, {
        dataStatus("Uploading data from dropbox")
        debrisInfo(downloadFromDropBox(dropbox_csv_url=input$url) %>% processDebris())
        dataStatus("Data ready")
    })
    dataStatus <- reactiveVal("Ready to upload data")
    output$dataStatus <- renderText({dataStatus()})
    
    # Get input
    mt_1 <- eventReactive(input$input_1, {as.numeric(input$input_1)})
    mt_10 <- eventReactive(input$input_10, {as.numeric(input$input_10)})
    mt_100 <- eventReactive(input$input_100, {as.numeric(input$input_100)})
    mt_1000 <- eventReactive(input$input_1000, {as.numeric(input$input_1000)})
    mt_10000 <- eventReactive(input$input_10000, {as.numeric(input$input_10000)})
    
    # Make plots
    concernList <- eventReactive(input$MakeUpdates, {
        purrr::map(.x=c(5,4,3),
                   .f=trade_off_plot,
                   debrisInfo=debrisInfo(),
                   missTolerance = list(">= 1"=mt_1(),
                                        ">= 10"=mt_10(),
                                        ">= 100"=mt_100(),
                                        ">= 1000"=mt_1000(),
                                        ">= 10000"=mt_10000()))
    })

    final <- eventReactive(input$monte, {
        evalPerformance(numSamples=100,
                        debrisInfo=debrisInfo(),
                        optimalThresholds=purrr::map_dfr(concernList(), function(x) x[[2]]))
    })

    output$plot1 <- renderPlotly({
        concernList()[[1]][["plot"]]
    })

    output$plot2 <- renderPlotly({
        concernList()[[2]][["plot"]]
    })

    output$plot3 <- renderPlotly({
        concernList()[[3]][["plot"]]
    })

    output$plot4 <- renderPlotly({
        final()[["plot"]]
    })

    output$tbl_final_rec <- renderDT(server=FALSE, {
        purrr::map_dfr(concernList(), function(x) x[[2]])[
                           ,c("TCA_Bin","fragLabel","WarnThreshold")] %>%
            datatable(style="bootstrap",selection="none",
                      options = list(searching = FALSE))
    })

    output$tbl_final_per <- renderDataTable({
        datatable(final()[["finalPerformance"]],
                  style = 'bootstrap',selection="none",
                  options = list(searching = FALSE))
    })

    output$downloadData <- downloadHandler(
          filename = function() {
            paste('data-', Sys.Date(), '.csv', sep='')
          },
          content = function(con) {
            fwrite(purrr::map_dfr(concernList(), function(x) x[[2]])[
                ,c("TCA_Bin","fragLabel","WarnThreshold")], con)
          }
        )

} # end server function

shinyApp(ui, server)
