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

rfc1_epic_export <- read_excel("data/Checking_CANVAS_referrals_20220902.xlsx",
                               sheet = "epic_output") %>%
  janitor::clean_names() %>%
  mutate(full_name = paste(toupper(pt_first_nm), toupper(patient_surname), sep = " ")) 

# Updated spreadsheet from Riccardo Curro
updated_research_results <- read_excel(path = "data/CANVAS_Screeninglist_update_Aug2022_GOSH.xlsx") %>%
  janitor::clean_names() %>%
  mutate(full_name = paste0(first_name, " ", surname))

##############################
# Samples with a positive research result
##############################

positives <- c("biallelic RFC1 expansion", "likely biallelic RFC1 expansion", 
               "likely biallelic RFC1 expansion (rechecked by Joe)", "confirmed",
               "already confirmed", "patient already confirmed")

positives_for_confirmation <- updated_research_results %>%
  filter(interpretation %in% positives) %>%
  mutate(full_name = paste(toupper(first_name), toupper(surname), sep = " "))

# 112 positive results.
# We have had 29 positive results so far.

##############################
# Samples already tested
##############################

# Identify samples from previous runs
previous_worksheets <- c("22-2268", "22-2325", "22-2543", "22-2649")

samples_tested_already <- rfc1_db %>%
  filter(worksheet %in% previous_worksheets)

# We only want to test samples with research results that we haven't already tested
samples_for_confirmation <- positives_for_confirmation %>%
  left_join(rfc1_epic_export, by = "full_name") %>%
  select(test_specimen_id, pt_first_nm, patient_surname, full_name,
         storage_location, final_dna_conc_ng_ml, interpretation) %>%
  filter(!is.na(test_specimen_id) &
           # Check that the patient name hasn't been tested before (the same patient
           # may have multiple sample numbers)
           !full_name %in% samples_tested_already$full_name) %>%
  # Remove duplicates
  filter(!base::duplicated(full_name)) %>%
  select(-full_name)

write.csv(samples_for_confirmation, "W:/MolecularGenetics/Neurogenetics/Research/Joe Shaw Translational Post 2022/RFC1 worksheets/22-3206/22_3206_samples.csv",
          row.names = FALSE)

##############################
# Checking research results
##############################

# Cases where duplicate samples tested
duplicate_tests <- updated_research_results %>%
  filter(duplicated(full_name, fromLast = FALSE) | duplicated(full_name, fromLast = TRUE)) %>%
  arrange(full_name)

# Cases with no PCR products on flanking PCR

no_product <- c("no PCR products", "no PCR product", "no PCR products (checked twice)", "No PCR products",
                "No PCR product", "no product", "no pcr product", "no pcr prod",
                "no PCR product (empty?)", "no pcr Product" )

# Did they retest the discordant sample?

##############################