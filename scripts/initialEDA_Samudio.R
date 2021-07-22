################################
# Exploratory Data Analysis ####
################################


# Read in trimmed data
events <- fread("events.csv")

# Insert PcNom value when there is no PcBest ####

# Verify the "when there is no Pc Best"
events[is.na(`Pc Best`),.N] # 148950
events[is.nan(`Pc Best`),.N] # 0
events[!is.na(`Pc Best`),.N] # 261939
nrow(events)==events[is.na(`Pc Best`),.N] + events[!is.na(`Pc Best`),.N]

# Generate "Pc" column
events[,"Pc" := ifelse(is.na(`Pc Best`),`Pc Nom`,`Pc Best`)]

# Drop redundant Pc columns so we have only one Pc column
events <- events[,!c("Pc Best","Pc Nom")]

# How many NAs for Pc still left?
events[is.na(Pc),.N] # 13

# Inspect
events[is.na(Pc)]

# NAs across ProbCatIfColl / NumFragIfCatColl / NumFragIfNonCatColl / Pc
# Drop them!
events <- events[!is.na(Pc)]

# Write and push to Git hub
fwrite(events,"events.csv")


# Divide entries up into bins by time to TCA ####
library(data.table)
library(dplyr)
library(fst)

# Read in trimmed data
events <- fread("events.csv")
events[,`Days_to_TCA`] %>% summary()

# Crude histogram of TCA distribution (will pretty this up for the final product)
events[,`Days_to_TCA`] %>% hist() # Uniform distribution

# Create TCA Bin column
events[,"TCA_Bin" := as.factor(round(`Days_to_TCA`))]

# Plot PDF of the Pc values for each bin ####
events[,unique(`TCA_Bin`)]

# How many NA values for TCA?
events[is.na(`TCA_Bin`),.N] # Just one
events[is.na(`TCA_Bin`)] # Inspect

# How many other events with the same event number (118.47)
events[`EventNumber`==118.47,.N] # Just one

# Drop it!
events <- events[`EventNumber`!=118.47]

# Configure data set to show Pc Values grouped by TCA bin

# How many Pc at 0?
events[Pc==0,.N] # 231677

# How many Pc above 1e-4?
events[Pc > 1e-4,.N]

# Percentage of events with Pc above 1e-4?
events[Pc > 1e-4,.N]/events[,.N]*100 # 0.0954


# Plot
library(plotly)
library(RColorBrewer)

# Configure data set to show Pc Values grouped by TCA bin (fix me - this looks so ugly)
events[,Pc,by=`TCA_Bin`] %>% 
  
  # Plot
  plot_ly(type="histogram",
          x=~Pc, 
          color =~`TCA_Bin`,
          nbinsx=35,
          colors=colorRampPalette(brewer.pal(10,"Spectral"))(10)) 

# Lots of 0 - plot without Pc == 0

events[Pc > 0,Pc,by=`TCA_Bin`] %>% 
  
  # Plot
  plot_ly(type="histogram",
          x=~log(Pc), 
          color =~`TCA_Bin`,
          nbinsx=35,
          colors=colorRampPalette(brewer.pal(10,"Spectral"))(10)) %>%
  
  layout(
    title = '<b>Histogram of Pc by TCA Bin</b>',
    xaxis = list(title="<b>Probability of Collision</b>",
                 tickvals = ~Pc %>% log() %>% quantile(0:3/3),
                 ticktext = ~Pc %>% quantile(0:3/3)),
    yaxis = list(gridcolor = "#D3D3D3"),
    shapes = list(type ="line",
                  line = list(color="white"),
                  x0 = log(1e-4), x1 = log(1e-4),
                  y0 = 0, y1 = 6000),
    annotations = list(text = paste("<b>1e-4 Pc</b>"),
                       x = log(1e-35), y = 5800, showarrow=FALSE),
    legend=list(title=list(text='<b> TCA Bin </b>')),
    plot_bgcolor  = "#3F3F3F",
    paper_bgcolor = "#3F3F3F",
    font = list(color = '#D3D3D3'))

# How many warnings would be issued at each TCA bin?
events[Pc > 1e-5,.("Warning Threshold"=ifelse(Pc > 1e-4,"1e-4","1e-5"),
                   "Warnings Issued"=.N),by=`TCA_Bin`] %>% 
  plot_ly(type="bar",
          x=~`TCA_Bin`,
          y=~`Warnings Issued`,
          color=~`Warning Threshold`,
          colors=colorRampPalette(brewer.pal(3,"Set2"))(2))




# Introduce the # of fragments as a second consideration ####

# Label catastrophic events!
events[,"Catastrophic" := as.factor(ifelse(ProbCatIfColl > 0.5,"Catastrophic","Not catastrophic"))]

# Label num of fragments
events[,"NumFrag" := ifelse(Catastrophic==TRUE,NumFragIfCatColl,NumFragIfNonCatColl)]

# Drop two redundant columns
events <- events[,!c("NumFragIfCatColl","NumFragIfNonCatColl")]

# Write and push to Git
fwrite(events,"events.csv")

# Plot histogram of # of fragments for each catastrophic event
events[!is.na(NumFrag),NumFrag,by=c("EventNumber","Catastrophic")] %>%
  
  plot_ly(type="histogram",
          x=~log(NumFrag),
          color=~Catastrophic,
          nbinsx=35) %>%
  layout(
    title = '<b>Histogram of Number of Fragments by Event</b>',
    xaxis = list(title="<b>Fragment Count</b>",
                 tickvals = ~NumFrag %>% log() %>% quantile(0:7/7),
                 ticktext = ~NumFrag %>% quantile(0:7/7) %>% round(2)),
    yaxis = list(title="<b>Number of Events</b>",
                 gridcolor = "#D3D3D3"),
    plot_bgcolor  = "#3F3F3F",
    paper_bgcolor = "#3F3F3F",
    font = list(color = '#D3D3D3'))


# Time series!!!
# Use 0.5 < Days_to_TCA < 1 as proxy for an actual collision
events[,"Collision" := as.factor(Days_to_TCA > 0.5 & Days_to_TCA < 1.0)]




# If the General wants a five day warning time, 
# what threshold would have to be applied at five days to TCA to warn
# for every conjunction with a Pc > 1e-4 or 1e-5 at 1 day to TCA?

events[
  # Find the events that have a TCA of at least 5 days out
  Days_to_TCA > 5 & Collision==TRUE,Pc] %>% 
  
  plot_ly(type="histogram",
          x=~Pc) %>%
  layout(
    title = "<b>Histogram of Probability of Collision at 5 days\nfor Events with a Pc > 1e-5 at 1 day out"
    yaxis = list(title="<b>Number of Events</b>",
                 gridcolor = "#D3D3D3"),
    plot_bgcolor  = "#3F3F3F",
    paper_bgcolor = "#3F3F3F",
    font = list(color = '#D3D3D3'))


# 
