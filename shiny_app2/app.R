## app.R ##
library(shinydashboard)
library(cowplot)
library(remotes)
library(bslib)
library(tidyverse)
library(flexdashboard)
library(scales)
library(data.table)
library(shiny)
library(plotly)
library(data.table)
library(dplyr)
library(plotly)
library(RColorBrewer)

#setwd("shiny_app2")

#source("scripts/debrisConfusion.R")
#
# Read in debrisCounfuction function
#source("scripts/debrisConfusion.R")


ui <- dashboardPage(

    dashboardHeader(title = "USSPACECOM Space Conjuction",
                    titleWidth = 350),
    dashboardSidebar(
        HTML('<center><img src="SPACECOM_logo.png" width="150"></center>'),
        #img(src='SPACECOM_logo.png', width = 150, height = 150, align = "center"),
        p( paste("Space conjunction analysis in support of USSPACECOM. 
        Choose the parameters then press update.")),
        # Number of Fragments
        
        #checkboxGroupInput("NumFrags")
        
        checkboxGroupInput("NumFrags", label = "Number of fragments:",
                            choices = c("PcBest","PcFrag10","PcFrag100"), 
                           selected = c("PcBest","PcFrag10","PcFrag100")),
        # Days to TCA 
        sliderInput("DaysTCA", label = "Days to TCA",
                    min = 1, max = 8, value = 5, step = 1),
        
        #threshold_Pc
        
        # sliderInput("PcWarn", label = "Warning Threshold",
        #             min = 0.000000001, max = 0.001, value = 0.000001, step =0.00001),
        
        textInput("inputA", "Concern Category A", 1e-5),
        textInput("inputB", "Concern Category B", 1e-6),
        textInput("inputC", "Concern Category C", 1e-7),
        
        actionButton("MakeUpdates", "Update")
    ),
    dashboardBody(
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
    )


server <- function(input, output) {
    load("data/concernData.RData",envir = .GlobalEnv)
    conf_matrix_plot<-function(df = data_from_function, i){
        
        
        
        mat <-matrix(df[i,2:5],2,2) %>% data.frame() %>%  rename("No concern event"= X1, "concern event" = X2)
        
        data <-matrix(rep(0,12),4,3) %>% data.frame()%>%  rename("Warning"= X1, "Reference" = X2, "Freq" = X3)
        
        data$Warning <- c("concern event", "No concern event", "concern event","No concern event")
        #data$Warning <- c("No concern event", "concern event", "No concern event","concern event")
        data$Reference <- c("concern event", "concern event", "No concern event", "No concern event")
        data$discript <-c("Warning given and concern event","No warning and concern event", "Warning given and no concern event" ,"No warning and no concern event")
        
        
        #data$Reference <- c("No concern event", "No concern event", "concern event", "concern event")
        data[1,3] <- mat[2,2]
        data[3,3] <- mat[1,2]
        data[2,3] <- mat[2,1]
        data[4,3] <- mat[1,1]
        data <-data %>% 
            mutate(Warning = factor(Warning, levels = c( "No concern event","concern event")),  
                   Reference= factor(Reference, levels = c("concern event", "No concern event"))) 
        
        num_obs <- sum(data$Freq)
        num_ben <- sum(data[3,3]+ data[4,3])
        num_mal <- sum(data[1,3] + data[2,3])
        
        data <- data  %>%group_by(Reference) %>% 
            mutate(
                total = sum(Freq),
                frac_fill = if_else(Warning == Reference, Freq / total, 0),
                frac = Freq / num_obs) %>% 
            ungroup() %>% 
            mutate(dist = ifelse(Reference == "concern event", Freq/num_mal, Freq/num_ben )) 
        
        
        ggplot(data,aes(Warning,Reference,fill = frac)) +
            geom_tile() +
            scale_fill_gradient(low = "#FFFFFF", high = "#0d4077",guide = FALSE) +
            scale_x_discrete(position = "top") +
            geom_text(aes(color = frac>0.1,label = str_c(discript,"\n",comma(as.integer(Freq)),"\n ", round(frac * 100,2), "% of all Observations \n", round(dist *100,2), "% of True ", Reference  )), size = 6)+
            scale_color_manual(values=c("#000000", "#FFFFFF"),guide = FALSE)+
            geom_tile(color = "black", fill = "black", alpha = 0)+
            theme_cowplot()+
            theme(axis.text=element_text(size=16),
                  axis.title=element_text(size=16,face="bold"))
    }
    
    make_data <- function(){
        
        index = 1
        list_of_plots = list()
        
        for( i in c(input$inputA, input$inputB, input$inputC)){
            i <- as.numeric(i)
            
            list_of_plots[[index]] <-FN_FP_Plot(pc_concern=i,
                                                TCA_days =input$DaysTCA,
                                                frags=input$NumFrags)
            index = index + 1

            
            # Temp_df <- debrisConfusion(concernData_at_TOI=concernSummary_at_TOI,
            #                            concernData=concernSummary,
            #                            Pc_concern=i,
            #                            days_to_TCA=input$DaysTCA,
            #                            Pc_warn=input$PcWarn,
            #                            frag_number_Pc=input$NumFrags,
            #                            verbose=FALSE,
            #                            rate = FALSE)
            # 
            # Temp_df <- Temp_df %>% 
            #     pivot_wider(names_from = `Type`, values_from = Count)
            # 
            # names(Temp_df)<-  c("name", "false_neg", "false_pos","true_pos" ,"true_neg")
            # 
            # Temp_df <- Temp_df %>% select(name, true_neg,false_neg ,false_pos ,true_pos )
            # 
            # index = index + 1
            # list_of_df[[index]] <- Temp_df
        }
        
        return(list_of_plots)}
    
    
    
    source("scripts/debrisConfusion.R")
    source("scripts/FN_FP.R")
    list_of_df = list()
    
    react <- eventReactive(input$MakeUpdates, {make_data()})
    
    output$plot1 <- renderPlotly({
        dfA <- react()[[1]]
        dfA
        })
    
    output$plot2 <- renderPlotly({
         dfB <- react()[[2]]
         dfB
        # FN_FP_Plot(pc_concern =input$inputB,
        #            TCA_days =input$DaysTCA )
    })
    
    output$plot3 <- renderPlotly({
        dfC <- react()[[3]]
        dfC
    })

    
}

shinyApp(ui, server)
