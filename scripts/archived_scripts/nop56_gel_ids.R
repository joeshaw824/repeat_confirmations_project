################################################################################
## NOP56 Samples for WGS
################################################################################

# Script for adding GEL identifiers to list of samples for Dalia Kasperaviciute

library(tidyverse)
library(janitor)
library(readxl)

source("functions/rfc1_functions.R")

setwd(dir = "W:/MolecularGenetics/Neurogenetics/Research/Joe Shaw Translational Post 2022/RFC1 R code/RFC1_analysis")

gel_ids <- read_csv("data/glh_samples_filtered-2022-10-07.csv") %>%
  janitor::clean_names()

ataxia_cases <- read_csv("data/North Thames R54.3 Cases to 06.10.22.csv") %>%
  dplyr::rename(gel1001_primary_sample_id_glh_lims = specimen_id)

ataxia_cases_with_ids <- ataxia_cases %>%
  left_join(gel_ids, by = "gel1001_primary_sample_id_glh_lims") %>%
  select(gel1001_primary_sample_id_glh_lims, gel1001_patient_ngis_id, 
         gel1001_referral_id,
         gel1001_dispatched_sample_lsid, gel1001_gmc_rack_id, 
         gel1001_gmc_rack_well, )

write.csv(ataxia_cases_with_ids, "outputs/north_thames_ataxia_cases_with_ids.csv",
          row.names = FALSE)


#################################
# 21/11/2022: get patient names and DNA concentrations
#################################

ataxia_names <- read_excel("data/Ataxia_panel_referrals_20221121_1636.xlsx") %>%
  janitor::clean_names()

ataxia_names_split_locations <- split_dna_location(ataxia_names)

ataxia_cases_with_names <- ataxia_cases %>%
  dplyr::rename(test_specimen_id = gel1001_primary_sample_id_glh_lims) %>%
  left_join(ataxia_names_split_locations, by = "test_specimen_id" ) %>%
  mutate(storage_location2 = gsub("RG DNA ", "", storage_location)) %>%
  arrange(tray) %>%
  select(test_specimen_id, pt_first_nm, patient_surname, storage_location2, final_dna_conc_ng_ml) %>%
  filter(!is.na(final_dna_conc_ng_ml))

write.csv(ataxia_cases_with_names, "outputs/ataxia_cases_with_names.csv",
          row.names = FALSE)

################################################################################