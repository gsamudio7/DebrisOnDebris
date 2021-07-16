
# Imports
library(data.table)
library(dplyr)


# Read in the data
observations <- fread("DebrisOnDebris.csv", fill=TRUE)


# Transform data to only the relavant columns (based on Mr. Hejduk's guidance)

observations <- observations[,c(49,48,9,6,12,14,15)]
observations <- observations %>% 
  setnames(c("V49","V48","V9","V6","V12","V14","V15"),
           c("Event Number","Days to TCA","Pc Best","Pc Nom","ProbCatIfColl",
             "NumFragIfCatColl","NumFragIfNonCatColl"))

# Write
fwrite(observations,"trimmedData.csv")