################################################################################
## Selecting RFC1 samples for confirmations
################################################################################

##############################
# Load libraries
##############################

library(tidyverse)
library(readxl)
library(janitor)
library(lubridate)

setwd(dir = "W:/MolecularGenetics/Neurogenetics/Research/Joe Shaw Translational Post 2022/RFC1 R code/RFC1_analysis")

source("scripts/rfc1_database.R")

##############################
# Read in files
##############################

rfc1_epic_export <- read_excel("data/CANVAS_Referrals_with_Clinicians_20220621_1608.xlsx") %>%
  janitor::clean_names() %>%
  mutate(full_name = paste(toupper(pt_first_nm), toupper(patient_surname), sep = " ")) 

automatic_reports <- read_excel( path = "data/AC_CANVAS_Screeninglist_2021_GOSH.xlsx") %>%
  janitor::clean_names() %>%
  mutate(full_name = paste(toupper(first_name), toupper(surname), sep = " ")) %>%
  filter(!base::duplicated(full_name))

##############################
# Select samples for next worksheet
##############################

# Identify samples from previous runs
previous_worksheets <- c("22-2268", "22-2325")

samples_tested_already <- rfc1_db %>%
  filter(worksheet %in% previous_worksheets)

# We only want to test samples with research results that we haven't already tested
samples_with_locations <- automatic_reports %>%
  left_join(rfc1_epic_export, by = "full_name") %>%
  select(test_specimen_id, pt_first_nm, patient_surname, full_name,
         storage_location, final_dna_conc_ng_ml, interpretation) %>%
  filter(!is.na(test_specimen_id) &
           # Check that the patient name hasn't been tested before (the same patient
           # may have multiple sample numbers)
           !full_name %in% samples_tested_already$full_name &
           interpretation != "pending") %>%
  # Remove duplicates
  filter(!base::duplicated(full_name))

samples_to_test <- dplyr::slice_sample(samples_with_locations, n = 40) %>%
  select(-full_name)

write.csv(samples_to_test, "W:/MolecularGenetics/Neurogenetics/Research/Joe Shaw Translational Post 2022/RFC1 worksheets/22-2543/22_2543_samples.csv",
          row.names = FALSE)

##############################