
# Imports
library(data.table)
library(dplyr)


## Read in the data (executed before initial clean)
# observations <- fread("DebrisOnDebris.csv", fill=TRUE)


# # Transform data to only the relevant columns (based on Mr. Hejduk's guidance)
# observations <- observations[,c(48,49,9,6,12,14,15)]
# observations <- observations %>% 
#   setnames(c("V48","V49","V9","V6","V12","V14","V15"),
#            c("Event Number","Days to TCA","Pc Best","Pc Nom","ProbCatIfColl",
#              "NumFragIfCatColl","NumFragIfNonCatColl"))
# 
# # Write
# fwrite(observations,"trimmedData.csv")


# Read in trimmed data
observations <- fread("trimmedData.csv")

# Insert PcNom value when there is no PcBest
observations[,"Pc" := ifelse(is.na(`Pc Best`),`Pc Nom`,`Pc Best`)]

# Drop redundant Pc columns so we have only one Pc column
observations <- observations[,!c("Pc Best","Pc Nom")]

# Write
fwrite(observations,"trimmedData.csv")
