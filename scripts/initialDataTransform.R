# 
# ###################
# # Process Data ####
# ###################
# 
# 
# Imports
library(data.table)
library(dplyr)
library(fst)

# Read in the raw data
debrisData <- fread("data/DebrisOnDebris3.txt")
dim(debrisData) # 707188 rows and 58 columns

# Write and compress to store in Git repo
write_fst(debrisData, "data/DebrisOnDebris3.fst", compress = 100)

# Now entire data set is compressed and can be stored in Git repo

# Read in .fst file and convert to data.table for fast compute
events <- read_fst("data/DebrisOnDebris3.fst") %>% data.table()
dim(events) # 707188/58 -> no loss of data


# Transform data to only the relevant columns (based on Mr. Hejduk's guidance)
events <- events[,c(1,2,3,11,14,15,16,17,45,46,58)]

# Name columns
events <- events %>%
  setnames(c("V1","V2","V3",
             "V11","V14","V15",
             "V16","V17","V45",
             "V46","V58"),
           c("Primary","Secondary","MissDistance",
             "PcBest","PcFrag10","PcFrag100",
             "PcFrag1000","PcFrag10000",
             "time_of_screening","TCA",
             "eventNumber"))

# Generate time to TCA column and remove redundant columns
events[,"time2TCA" := TCA - time_of_screening]

# Write and push to Git
fwrite(events,"data/events.csv")

# Git commands in the terminal:
# git add .
# git commit -m "uploading transform script and events.csv"
# git push


