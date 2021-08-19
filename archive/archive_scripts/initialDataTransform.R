# 
# ###################
# # Process Data ####
# ###################

# Imports
rm(list=ls())
library(data.table)
library(dplyr)
library(fst)
library(readxl)

########################################
# This part completed - do not re run 
# DebrisOnDebris3.xlsx not in repo 

# Read in the raw data
# debrisData <- read_excel("data/DebrisOnDebris3.xlsx")
# dim(debrisData) # 707188 rows and 58 columns

# Write and compress to store in Git repo and delete original .xlsx file
# write_fst(debrisData, "data/DebrisOnDebris3.fst", compress = 100)
# unlink("data/DebrisOnDebris3.xlsx")

# Now entire data set is compressed and can be stored in Git repo
###########################################

#############################################
# Data transform begins here from .fst file 
#############################################

# Read in .fst file and convert to data.table for fast compute
read_fst("data/DebrisOnDebris3.fst") %>% data.table() %>% fwrite("data/debrisData.csv")
events <- fread("data/debrisData.csv",header=FALSE)
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

# Remove redundant data in data file
unlink(c("data/DebrisOnDebris_v2.fst",
         "data/DebrisOnDebris.fst",
         "data/debrisData.csv"))

# Git commands in the terminal:
# git add .
# git commit -m "uploading transform script and events.csv"
# git push


