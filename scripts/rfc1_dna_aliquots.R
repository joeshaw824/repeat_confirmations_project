################################################################################
## RFC1 DNA Aliquots for Research Testing
################################################################################

# This script is used for generating a list of samples which require RFC1
# testing on a research basis by the Cortese group at the Institute of 
# Neurology.

# This whole process is complicated by:
# 2 labs having recently merged
# 4 different information management systems (GOSH-Epic, UCLH-Epic, Winpath and PDMS)
# Covid disruption
# The researchers don't issue any reports
# The researchers don't consistently use a patient-specific identifier
# Mark Gaskin does not want to see patient-specific identifiers
# Clinicians sent multiple samples for the same patient

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

most_recent_export <- "CANVAS_Referrals_with_Clinicians_20221205_0842.xlsx"

rfc1_epic_requests <- read_excel(path = paste0("data/", most_recent_export)) %>%
  janitor::clean_names() %>%
  dplyr::rename(specimen_id = test_specimen_id)

##############################
# Samples sent to the Cortese group by the GOSH prep lab
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

##############################
# Samples aliquotted by Mark Gaskin
##############################

# This script reads all the pull sheets sent by Mark Gaskin and compiles as a single
# dataframe called "pullsheet_merge".
source("scripts/dna_aliquotting_mark_gaskin.R")

# Join dataframes together to give a summary table of all aliquotting since
# April 2021
all_aliquots <- rbind(pullsheet_merge %>%
                        mutate(aliquotted_by = "Mark Gaskin") %>%
                        dplyr::rename(specimen_id = original_sample_id) %>%
                        select(specimen_id, amount_taken_ul, aliquotted_by, sheet),
                      prep_lab_samples %>%
                        dplyr::rename(amount_taken_ul = volume_ul) %>%
                        mutate(sheet = "") %>% 
                        select(specimen_id, amount_taken_ul, aliquotted_by, sheet))

all_successful_aliquots <- all_aliquots %>% filter(amount_taken_ul > 0)

##############################
# 05/12/2022 - Which samples need to be aliquotted?
##############################

# Find samples on the Epic list which aren't on the list of successfully aliquotted samples.

# Because a consistent patient-specific identifier isn't used and because the clinicians
# keep sending multiple samples, I will rely on my previous checks on 02/09/22 to identify
# samples which are duplicates for the same patient.

previous_check <- read_excel("data/Checking_CANVAS_referrals_20220902.xlsx",
                                     sheet = "epic_output_annotated") %>%
  janitor::clean_names()

samples_already_tested <- previous_check %>% filter(dna_required == "No")


samples_to_aliquot <- rfc1_epic_requests %>%
  filter(!specimen_id %in% all_successful_aliquots$specimen_id &
           !specimen_id %in% samples_already_tested$test_specimen_id)
  
suspicious_samples <- samples_to_aliquot %>%
  filter(collection_instant < "2022-09-03")


to_check <- previous_check %>%
  filter(test_specimen_id %in% suspicious_samples$specimen_id)


# 22RG-333G0126
# 22RG-332G0070 - sample from Dr Anna Latorre

##############################
# 02/09/2022 - Manually checking for missed samples
##############################

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
          file = paste0("outputs/rfc1_samples_for_aliqotting_patient_identifiers_removed_", 
                        format(Sys.time(), "%Y_%m_%d_%H_%M_%S"), ".csv"),
          row.names = FALSE)

##############################
# The summary
##############################

# On Epic, there are:
nrow(rfc1_epic_requests)
# requests for RFC1 testing.

# Since 2019, we have sent
nrow(all_successful_aliquots)
# aliquots of DNA to Andrea Cortese's group


