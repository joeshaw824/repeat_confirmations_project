################################################################################
## RFC1 DNA Aliquots for Research Testing
################################################################################

##############################
# Load libraries
##############################

library(tidyverse)
library(readxl)
library(janitor)

setwd(dir = "W:/MolecularGenetics/Neurogenetics/Research/Joe Shaw Translational Post 2022/RFC1 R code/RFC1_analysis")

##############################
# Load resources
##############################

# This is an Excel that I manually typed out from the paper forms in the prep lab research
# aliquots folder
prep_pab_samples <- read_excel("data/samples_requested_by_cortese_group.xlsx",
                               sheet = "aliquots") %>%
  dplyr::rename(specimen_id = episode)

samples_with_ids <- prep_pab_samples %>%
  filter(!is.na(specimen_id))

epic_details <- read_excel("data/CANVAS_Referrals_with_Clinicians_20220621_1608.xlsx") %>%
  janitor::clean_names() %>%
  dplyr::rename(specimen_id = test_specimen_id)

# Check to see if manual entry of numbers was incorrect

samples_with_details <- samples_with_ids %>%
  left_join(epic_details, by = "specimen_id") %>%
  filter(is.na(patient_surname))

##############################
