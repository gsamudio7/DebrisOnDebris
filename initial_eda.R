
###################
# Process Data ####
###################


# Imports
library(data.table)
library(dplyr)
library(fst)

# Read in the raw data 
debrisData <- fread("DebrisOnDebris.csv", fill=TRUE)
dim(debrisData) # 410889 rows and 49 columns

# Write and compress to store in Git repo
write_fst(debrisData, "DebrisOnDebris.fst", compress = 100)

# Remove raw large data set 
unlink("DebrisOnDebris.csv")

# Now entire data set is compressed and can be stored in Git repo


# Read in .fst file and convert to data.table for fast compute
events <- read_fst("DebrisOnDebris.fst") %>% data.table()

# Transform data to only the relevant columns (based on Mr. Hejduk's guidance)
events <- events[,c(48,49,9,6,12,14,15)]
events <- events %>%
  setnames(c("V48","V49","V9","V6","V12","V14","V15"),
           c("Event Number","Days to TCA","Pc Best","Pc Nom","ProbCatIfColl",
             "NumFragIfCatColl","NumFragIfNonCatColl"))

# Write
fwrite(events,"trimmedData.csv")


# Read in trimmed data
events <- fread("trimmedData.csv")

# Insert PcNom value when there is no PcBest ####
events[,"Pc" := ifelse(is.na(`Pc Best`),`Pc Nom`,`Pc Best`)]

# Drop redundant Pc columns so we have only one Pc column
events <- events[,!c("Pc Best","Pc Nom")]

# Write
fwrite(events,"trimmedData.csv")


# Divide entries up into bins by time to TCA ####
# Read in trimmed data
events <- fread("trimmedData.csv")
events[,`Days to TCA`] %>% summary()

# Crude histogram of TCA distribution
events[,`Days to TCA`] %>% hist()


