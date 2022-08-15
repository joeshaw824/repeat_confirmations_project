################################################################################
# DNA Sample Locations for Research Aliquots
# August 2022
# Joseph.Shaw@gosh.nhs.uk
################################################################################

# This is a script for quickly identifying tray locations for DNA samples required for research aliquots.
# Because there is no consistent test set/feature of these samples, an Epic MyReport can't be used.
# Instead, the approach I took was to pull all DNA sample locations out of Epic, and then query these
# with the sample numbers of interest.

##############################
# Load libraries
##############################

library(tidyverse)
library(readxl)
library(janitor)

setwd(dir = "W:/MolecularGenetics/Neurogenetics/Research/Joe Shaw Translational Post 2022/RFC1 R code/RFC1_analysis")

##############################
# Getting information out of Epic
##############################

# Epic MyReports has an export limit of 10000 records.
# So I pulled out sample locations in 10000 sample chunks, using the "collection date" of the 
# most recent sample in each chunk as a new threshold for the next query.

# The sample selection logic was:
# Collection Date: Greater than or equal to [modify]
# Test group: TestingLaboratory: GOSH REGIONAL GENETICS LAB 

##############################
# Load DNA locations from Epic outputs
##############################

dna_location_files <- list.files("data/dna_locations")

dna_locations <- data.frame()

for (file in dna_location_files) {
  
  tmp_file <- readxl::read_excel(path = paste0("data/dna_locations/", file)) %>%
    janitor::clean_names()
  
  dna_locations <- rbind(dna_locations, tmp_file)
  
  rm(tmp_file)
  
}

dna_locations_mod <- dna_locations %>%
  filter(!base::duplicated(specimen_id) &
           # Ensure only regional genetics (RG) samples included
           substr(dna_locations$specimen_id, 3, 4) == "RG")

##############################
# Select samples for Mark
##############################

sample_request <- read_excel("resources/mark_gaskin_sample_location_request.xlsx")

sample_request_with_locations <- sample_request %>%
  left_join(dna_locations_mod, by = "specimen_id")

write.csv(x = sample_request_with_locations, 
          file = paste0("outputs/mark_gaskin_samples_with_locations_", format(Sys.time(), "%Y%m%d"), ".csv"),
          row.names = FALSE)

##############################
