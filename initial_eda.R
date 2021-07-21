# 
# ###################
# # Process Data ####
# ###################
# 
# 
# Imports
# library(data.table)
# library(dplyr)
# library(fst)

# # Read in the raw data 
# debrisData <- fread("DebrisOnDebris.csv", fill=TRUE)
# dim(debrisData) # 410889 rows and 49 columns
# 
# # Write and compress to store in Git repo
# write_fst(debrisData, "DebrisOnDebris.fst", compress = 100)
# 
# # Remove raw large data set 
# unlink("DebrisOnDebris.csv")
# 
# # Now entire data set is compressed and can be stored in Git repo
# 
# 
# # Read in .fst file and convert to data.table for fast compute
# events <- read_fst("DebrisOnDebris.fst") %>% data.table()
# dim(events)
# 
# # Transform data to only the relevant columns (based on Mr. Hejduk's guidance)
# events <- events[,c(48,49,9,6,12,14,15)]
# events <- events %>%
#   setnames(c("V48","V49","V9","V6","V12","V14","V15"),
#            c("Event Number","Days to TCA","Pc Best","Pc Nom","ProbCatIfColl",
#              "NumFragIfCatColl","NumFragIfNonCatColl"))
# 
# # Write
# fwrite(events,"events.csv")



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
events[,`Days to TCA`] %>% summary()

# Crude histogram of TCA distribution (will pretty this up for the final product)
events[,`Days to TCA`] %>% hist() # Uniform distribution

# Create TCA Bin column
events[,"TCA Bin" := as.factor(round(`Days to TCA`))]

# Plot PDF of the Pc values for each bin ####
events[,unique(`TCA Bin`)]

# How many NA values for TCA?
events[is.na(`TCA Bin`),.N] # Just one
events[is.na(`TCA Bin`)] # Inspect

# How many other events with the same event number (118.47)
events[`Event Number`==118.47,.N] # Just one

# Drop it!
events <- events[`Event Number`!=118.47]

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
events[,Pc,by=`TCA Bin`] %>% 

  # Plot
  plot_ly(type="histogram",
          x=~Pc, 
          color =~`TCA Bin`,
          nbinsx=35,
          colors=colorRampPalette(brewer.pal(10,"Spectral"))(10)) 
  
# Lots of 0 - plot without Pc == 0

events[Pc > 0,Pc,by=`TCA Bin`] %>% 
  
  # Plot
  plot_ly(type="histogram",
          x=~log(Pc), 
          color =~`TCA Bin`,
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
                   "Warnings Issued"=.N),by=`TCA Bin`] %>% 
  plot_ly(type="bar",
          x=~`TCA Bin`,
          y=~`Warnings Issued`,
          color=~`Warning Threshold`,
          colors=colorRampPalette(brewer.pal(3,"Spectral"))(2))




# Introduce the # of fragments as a second consideration ####



