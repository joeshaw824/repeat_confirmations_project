################################################################################
## RFC1 Epic exports
################################################################################

# This script is to check whether there are discrepancies between the list of samples for RFC1 testing
# extracted from Epic and the list of samples emailed to Andrea Cortese.

##############################
# Load libraries
##############################

library(tidyverse)
library(readxl)
library(janitor)
library(lubridate)

setwd(dir = "W:/MolecularGenetics/Neurogenetics/Research/Joe Shaw Translational Post 2022/RFC1 R code/RFC1_analysis")

##############################
# Load spreadsheets
##############################

rfc1_epic_export <- read_excel("data/CANVAS_Referrals_with_Clinicians_20220621_1608.xlsx") %>%
  janitor::clean_names()

email_list <- read_excel("W:/MolecularGenetics/Neurogenetics/Research/CANVAS/Canvas list sent to Andrea/Combined_Beaker-EPIC_CANVAS_list_sent.xlsx")

##############################
# Check for differences
##############################

colnames(rfc1_epic_export)

# Samples on Epic but not in the email list
missed_samples <- setdiff(rfc1_epic_export$test_specimen_id, email_list$`Specimen ID`)

missed_sample_info <- rfc1_epic_export %>%
  filter(test_specimen_id %in% missed_samples) %>%
  select(test_specimen_id, patient_surname, pt_first_nm, batch)

##############################