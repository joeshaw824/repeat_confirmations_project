################################################################################
## DNA Aliquoting by Mark Gaskin
################################################################################

##############################
# Load libraries
##############################

library(tidyverse)
library(readxl)
library(janitor)

setwd(dir = "W:/MolecularGenetics/Neurogenetics/Research/Joe Shaw Translational Post 2022/RFC1 R code/RFC1_analysis/")

source("functions/rfc1_functions.R")

##############################
# Load pullsheets
##############################

all_pullsheets <- list.files("W:/MolecularGenetics/Neurogenetics/Research/DNA_aliquots_for_research/Pull sheets/")

pullsheet_merge <- data.frame()

for (i in all_pullsheets) {
  
  tmp_sheet <- read_pullsheet(i)
  
  pullsheet_merge <- rbind(pullsheet_merge, tmp_sheet)
  
  rm(tmp_sheet)
  
}

##############################