## app.R ##
library(shinydashboard)
library(data.table)
library(shiny)
library(plotly)
library(dplyr)
library(plotly)
library(RColorBrewer)


ui <- dashboardPage(

    dashboardHeader(title = "USSPACECOM Space Conjuction",
                    titleWidth = 350),
    dashboardSidebar(
        HTML('<center><img src="SPACECOM_logo.png" width="150"></center>'),
        p( paste("Space conjunction analysis in support of USSPACECOM. 
        Choose the parameters then press update.")),
        
        # Number of Fragments
        checkboxGroupInput("NumFrags", label = "Number of fragments:",
                            choices = c("PcBest","PcFrag10","PcFrag100"), 
                           selected = c("PcBest","PcFrag10","PcFrag100")),
        # Days to TCA 
        sliderInput("DaysTCA", label = "Days to TCA",
                    min = 1, max = 8, value = 5, step = 1),
        
        # Concern inputs
        textInput("inputA", "Concern Category A", 1e-5),
        textInput("inputB", "Concern Category B", 1e-6),
        textInput("inputC", "Concern Category C", 1e-7),
        
        actionButton("MakeUpdates", "Update")
    ),
    dashboardBody(
        # add css custom to the app
        tags$head(
            tags$link(rel = "stylesheet", type = "text/css", href = "custom.css")
        ),
        # Boxes need to be put in a row (or column)
        fluidRow(
            box(
                h4("Concern Category A"),
                plotlyOutput("plot1", height = 235), width = 250)),
            
        fluidRow(
            box(
                h4("Concern Category B"),
                plotlyOutput("plot2", height = 235), width = 250)),
        fluidRow(
            box(
                h4("Concern Category C"),
                plotlyOutput("plot3", height = 235), width = 250))
        )
    ) # end UI


server <- function(input, output) {
    
    load("data/concernData.RData",envir = .GlobalEnv)
    source("scripts/debrisConfusion.R")
    source("scripts/FN_FP.R")
    
    make_data <- function(){
        
        index = 1
        list_of_plots = list()
        
        for( i in c(input$inputA, input$inputB, input$inputC)){
            i <- as.numeric(i)
            
            list_of_plots[[index]] <-FN_FP_Plot(pc_concern=i,
                                                TCA_days =input$DaysTCA,
                                                frags=input$NumFrags)
            index = index + 1
  
        } # end for loop
        
        return(list_of_plots)}
    
    list_of_df = list()
    
    react <- eventReactive(input$MakeUpdates, {make_data()})
    
    output$plot1 <- renderPlotly({
        dfA <- react()[[1]]
        dfA
        })
    
    output$plot2 <- renderPlotly({
         dfB <- react()[[2]]
         dfB

    })
    
    output$plot3 <- renderPlotly({
        dfC <- react()[[3]]
        dfC
    })

} # end server function

shinyApp(ui, server)
