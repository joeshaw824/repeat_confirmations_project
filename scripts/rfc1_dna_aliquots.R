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
# Which samples had RFC1 testing requested?
##############################

# "CANVAS Referrals with Clinicians" Epic MyReport
# Extracts all samples with an RFC1 test set added

rfc1_epic_requests <- read_excel(path = "data/Checking_CANVAS_referrals_20220902.xlsx",
                                 sheet = "epic_output") %>%
  janitor::clean_names() %>%
  dplyr::rename(specimen_id = test_specimen_id)

##############################
# Which samples have had DNA sent to the Cortese group?
##############################

# This is an Excel that I manually typed out from the paper forms in the prep lab research
# aliquots folder.
# These samples were aliquotted by the GOSH prep lab team and sent to the Cortese group.
prep_lab_samples <- read_excel("data/samples_requested_by_cortese_group.xlsx",
                               sheet = "aliquots") %>%
  dplyr::rename(specimen_id = episode) %>%
  mutate(aliquotted_by = "GOSH Prep Lab") %>%
  left_join(rfc1_epic_requests %>%
              select(specimen_id, mrn),
            by = "specimen_id")

# This is Mark Gaskin's pull sheet for DNA samples that he aliquotted on 11/08/2022
gaskin_samples <- read_excel(path = "W:/MolecularGenetics/Neurogenetics/Research/Joe Shaw Translational Post 2022/DNA_aliquots_for_research/T1611-Pull sheet.xlsx",
                             sheet = "Pull sheet",
                             skip = 2) %>%
  janitor::clean_names() %>%
  select(-c("x13", "x14", "x15")) %>%
  filter(!is.na(original_sample_id)) %>%
  rename(specimen_id = original_sample_id,
         volume_ul = amount_taken_ul) %>%
  mutate(date_aliquotted = as.Date("2022-08-11 UTC"),
         aliquotted_by = "Mark Gaskin") %>%
  left_join(rfc1_epic_requests %>%
              select(specimen_id, pt_first_nm, patient_surname, dob, mrn),
            by = "specimen_id") %>%
  mutate(patient_name = paste0(pt_first_nm, " ", patient_surname))


# Join dataframes together to give a summary table of all aliquotting since
# April 2021
all_aliquots <- rbind(gaskin_samples %>%
                        select(specimen_id, patient_name, dob, mrn, 
                               date_aliquotted, aliquotted_by, volume_ul),
                      prep_lab_samples %>%
                        select(specimen_id, patient_name, dob, mrn, 
                               date_aliquotted, aliquotted_by, volume_ul))

##############################
# 02/09/2022 - Manually checking for missed samples
##############################

# I manually checked through all the RFC1 requests from Epic to find any that had not been aliquotted
# either by the GOSH prep lab or by Mark Gaskin.
# I chose this method because there can be multiple episode numbers for the same patient, only one
# of which may have the RFC1 test set added.

##############################
# How many samples are from the same patient?
##############################

duplicated_requests <- rfc1_epic_requests %>%
  filter(duplicated(mrn, fromLast = FALSE) |
           duplicated(mrn, fromLast = TRUE)) %>%
  arrange(mrn) %>%
  select(specimen_id, pt_first_nm, patient_surname, dob, mrn, collection_instant)

# Number of patients with duplicated samples
length(unique(duplicated_requests$mrn))

##############################
# Which samples need to be aliquotted?
##############################

checked_epic_samples <- read_excel(path = "data/Checking_CANVAS_referrals_20220902.xlsx",
                                 sheet = "epic_output_annotated") %>%
  janitor::clean_names() %>%
  dplyr::rename(specimen_id = test_specimen_id)

additional_samples <- read_excel(path = "data/Checking_CANVAS_referrals_20220902.xlsx",
                                 sheet = "additional_samples") %>%
  janitor::clean_names() %>%
  filter(notes != "Doesn't appear to be for CANVAS") %>%
  mutate(storage_location = "-80 freezers in basement") %>%
  dplyr::rename(specimen_id = sample)

epic_samples_to_aliquot <- checked_epic_samples %>%
  filter(dna_required == "Yes" | is.na(dna_required))

extra_samples <- data.frame(
  specimen_id = c("21RG-304G0081", "21RG-294G0055"),
  storage_location = c("RG DNA Tray 233 slot 19-E", "RG DNA Tray 230 slot 21-F"))

samples_for_aliqotting <- additional_samples %>%
  select(specimen_id, storage_location) %>%
  rbind(extra_samples) %>%
  rbind(epic_samples_to_aliquot %>%
          select(specimen_id, storage_location)) %>%
  # "Empty" on Mark's list, other sample is 21RG-304G0081
  filter(specimen_id != "21RG-070G0100") %>%
  # Blood sample with no location, DNA for this patient is 21RG-294G0055
  filter(specimen_id != "21RG-294G0054")

write.csv(x = samples_for_aliqotting, 
          file = paste0("outputs/rfc1_samples_for_aliqotting_patient_identifiers_removed_", format(Sys.time(), "%Y_%m_%d_%H_%M_%S"), ".csv"),
          row.names = FALSE)

##############################
