################################################################################
## Frozen Blood Samples for Neurology Research
################################################################################

##############################
# Load libraries
##############################

library(tidyverse)
library(readxl)
library(janitor)
library(lubridate)

##############################

bioresource_2020_files <- list.files("Bioresources/2020")

read_bioresource <- function(file, year) {
  
  stopifnot(typeof(file) == "character")
  
  stopifnot(typeof(year) == "double")
  
  stopifnot(year %in% c(2020, 2021, 2022))
  
  bioresource_excel <- read_excel(paste0("I:/Genetics/DNA Lab/Prep Lab/Research samples from NHNN/Bioresources/", year, "/", file),
                          sheet = "Bioresource") %>%
    janitor::clean_names()
  
  return(bioresource_excel)
  
}

bioresource_collated <- data.frame()

for (file in bioresource_2020_files) {
  
  temp_file <- read_bioresource(file, 2020)
  
  bioresource_collated <- rbind(bioresource_collated, temp_file)
  
  rm(temp_file)
  
  return(bioresource_collated)
  
}


##############################