## app.R ##
library(shinydashboard)
library(data.table)
library(shiny)
library(plotly)
library(dplyr)
library(plotly)
library(RColorBrewer)
library(DT)


ui <- dashboardPage(

    dashboardHeader(title = "USSPACECOM Space Conjuction",
                    titleWidth = 350),
    dashboardSidebar(
        HTML('<center><img src="SPACECOM_logo.png" width="150"></center>'),
        # p( paste("Space conjunction analysis in support of USSPACECOM. 
        # Choose the parameters then press update.")),
        
        # Number of Fragments
        # checkboxGroupInput("NumFrags", label = "Number of fragments:",
        #                     choices = c("PcBest","PcFrag10","PcFrag100"), 
        #                    selected = c("PcBest","PcFrag10","PcFrag100")),
        # Days to TCA 
        # sliderInput("DaysTCA", label = "Days to TCA",
        #             min = 1, max = 8, value = 5, step = 1),
        h4("Miss Tolerances"),
        p( paste("Space conjunction analysis in support of USSPACECOM. 
        Choose the parameters then press update.")),
        
        # Concern inputs
        textInput("input_1", "Fragments >= 1", 15),
        textInput("input_10", "Fragments >= 10", 10),
        textInput("input_100", "Fragments >= 100", 5),
        textInput("input_1000", "Fragments >= 1,000", 1),
        textInput("input_10000", "Fragments >= 10,000", 0),

        
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
                h4("Time of Closest Approach 5 Days"),
                plotlyOutput("plot1", height = 430)),
            box(
                h4("Time of Closest Approach 4 Days"),
                plotlyOutput("plot2", height = 430))
        
            
            ),
            

        fluidRow(
            box(
                h4("Time of Closest Approach 3 Days"),
                plotlyOutput("plot3", height = 430)),
            box(
                h4("Final Plot"),
                plotlyOutput("plot4", height = 430))
            

            ),
        
        fluidRow(
            box(h4("Final Recommendations"),
                DTOutput('tbl_final_rec', height = 430)),
            
            box(h4("Final Performance"),
                DTOutput('tbl_final_per', height = 430))
        )
        )
    ) # end UI


server <- function(input, output) {
    
    #load("data/concernData.RData",envir = .GlobalEnv) # new data with new functions
    load("data/debrisData.RData", envir = .GlobalEnv)
    events <- fread("data/events.csv")
    scale <- 1.520833 # there is problem remove this tell Gabe
    source("scripts/derbrisConfusion_function.R")
    source("scripts/trade_off_plots_function.R")
    source("scripts/evalPerformance_function.R")
    #source("scripts/debrisConfusion.R")
    #source("scripts/FN_FP.R")

    
    
    make_data <- function(){
        
        index = 1
        list_of_plots = list()
        

        mt_list <- list(
            ">= 1"=input$input_1,
            ">= 10"=input$input_10,
            ">= 100"=input$input_100,
            ">= 1000"=input$input_1000,
            ">= 10000"=input$input_10000)
        
        tt <-c(seq(from=1e-7,to=1e-5,by=1e-6),1e-7) # by 1e-8 
        
        
        # assign(concernCount_5, trade_off_plot(tcaBins=c(5),test_thresholds=tt, missTolerance = mt_list),envir = .GlobalEnv)
        # assign(concernCount_4, trade_off_plot(tcaBins=c(4),test_thresholds=tt, missTolerance = mt_list),envir = .GlobalEnv)
        # assign(concernCount_3, trade_off_plot(tcaBins=c(3),test_thresholds=tt, missTolerance = mt_list),envir = .GlobalEnv)
        # 
        print(ls())
        concernCount_5 <- trade_off_plot(tcaBins=c(5),test_thresholds=tt, missTolerance = mt_list)
        concernCount_4 <- trade_off_plot(tcaBins=c(4),test_thresholds=tt, missTolerance = mt_list)
        concernCount_3 <- trade_off_plot(tcaBins=c(3),test_thresholds=tt, missTolerance = mt_list)
        
        # numsamples reduce to speed up app
        final <- evalPerformance(numSamples=100,
                                 optimalThresholdList=list(
                                     concernCount_5$optThresholds,
                                     concernCount_4$optThresholds,
                                     concernCount_3$optThresholds))
        
        list_of_plots[[1]] <- concernCount_5$plot %>% plotly_build()
        list_of_plots[[2]] <- concernCount_4$plot %>% plotly_build()
        list_of_plots[[3]] <- concernCount_3$plot %>% plotly_build()
        list_of_plots[[4]] <- final$finalRec
        list_of_plots[[5]] <- final$plot
        list_of_plots[[6]] <- final$finalPerformance
        
        return(list_of_plots)}
    
    list_of_df = list()
    
    react <- eventReactive(input$MakeUpdates, {make_data()})
    
    output$plot1 <- renderPlotly({
        #react()[[1]]
        TC_5 <- react()[[1]]
        TC_5
        })
    
    output$tbl_final_rec <- renderDT({
        datatable(react()[[4]], selection="multiple", escape=FALSE, 
        options = list(sDom  = '<"top">lrt<"bottom">ip'))
        })
    
    output$plot2 <- renderPlotly({
         TC_4 <- react()[[2]]
         TC_4

    })
    
    output$plot3 <- renderPlotly({
        TC_3 <- react()[[3]]
        TC_3
    })
    
    output$plot4 <- renderPlotly({
        TC_3 <- react()[[5]]
        TC_3
    })
    
    output$tbl_final_per<- renderDT({
        datatable(react()[[6]], selection="multiple", escape=FALSE, 
                  options = list(sDom  = '<"top">lrt<"bottom">ip'))
    })

} # end server function

shinyApp(ui, server)
