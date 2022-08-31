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

rfc1_epic_requests <- read_excel(path = "data/CANVAS_Referrals_with_Clinicians_20220831_1507.xlsx") %>%
  janitor::clean_names()

july_referrals <- grep("2022-07", rfc1_epic_requests$collection_instant, value = TRUE)

rfc1_epic_requests %>% 
  filter(substr(collection_instant, 1, 7) == "2022-07")

# 27 samples per month
# By December 2022 there will be 455 + (27*4) = 563 samples for testing


# Roy's spreadsheet of sample details he sent to Andrea
sent_samples <- read_excel("W:/MolecularGenetics/Neurogenetics/Research/CANVAS/Canvas list sent to Andrea/Combined_Beaker-EPIC_CANVAS_list_sent.xlsx") %>%
  janitor::clean_names()

##############################
# Which samples have had DNA sent to the Cortese group?
##############################

# This is an Excel that I manually typed out from the paper forms in the prep lab research
# aliquots folder.
# These samples were aliquotted by the GOSH prep lab team and sent to the Cortese group.
prep_pab_samples <- read_excel("data/samples_requested_by_cortese_group.xlsx",
                               sheet = "aliquots") %>%
  dplyr::rename(specimen_id = episode)

aliquotted_prep_lab_samples <- prep_pab_samples %>%
  # Remove samples where no sample was sent
  filter(volume_ul != 0)

# This is Mark Gaskin's pull sheet for DNA samples that he aliquotted on 11/08/2022
gaskin_samples <- read_excel(path = "W:/MolecularGenetics/Neurogenetics/Research/Joe Shaw Translational Post 2022/DNA_aliquots_for_research/T1611-Pull sheet.xlsx",
                             sheet = "Pull sheet",
                             skip = 2) %>%
  janitor::clean_names() %>%
  select(-c("x13", "x14", "x15")) %>%
  filter(!is.na(original_sample_id)) %>%
  rename(specimen_id = original_sample_id)

aliquotted_gaskin_samples <- gaskin_samples %>%
  # Remove samples which were "missing" or "empty"
  filter(sample_status == "Present")

##############################
# Which samples haven't been aliquotted?
##############################

aliquotted_samples <- c(aliquotted_gaskin_samples$specimen_id, aliquotted_prep_lab_samples$specimen_id)

not_aliquotted <- setdiff(rfc1_epic_requests$test_specimen_id, aliquotted_samples)

not_aliquotted_details <- rfc1_epic_requests %>%
  filter(test_specimen_id %in% not_aliquotted) %>%
  select(-c("submitter_ordering_dept", "external_requisition_comment", "storage_location",
            "final_dna_conc_ng_ml", "batch"))

not_on_sent_spreadsheet <- setdiff(rfc1_epic_requests$test_specimen_id, sent_samples$specimen_id)

not_on_spreadsheet_samples <- rfc1_epic_requests %>%
  filter(test_specimen_id %in% not_on_sent_spreadsheet) %>%
  select(-c("submitter_ordering_dept", "external_requisition_comment", "storage_location",
            "final_dna_conc_ng_ml", "batch"))


##############################
