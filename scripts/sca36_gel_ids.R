################################################################################
## SCA36 Samples for WGS
################################################################################

# Script for adding GEL identifiers to list of samples for Dalia Kasperaviciute

library(tidyverse)
library(janitor)

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

################################################################################