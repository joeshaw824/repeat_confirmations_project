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

##############################
# Select samples
##############################

rfc1_epic_export <- read_excel("data/CANVAS_Referrals_with_Clinicians_20220621_1608.xlsx") %>%
  janitor::clean_names() %>%
  mutate(full_name = paste(toupper(pt_first_nm), toupper(patient_surname), sep = " ")) 

automatic_reports <- read_excel( path = "data/AC_CANVAS_Screeninglist_2021_GOSH.xlsx") %>%
  janitor::clean_names() %>%
  mutate(full_name = paste(toupper(first_name), toupper(surname), sep = " ")) %>%
  filter(!base::duplicated(full_name))

ws_22_2268 <- read_excel("W:/MolecularGenetics/Neurogenetics/Research/Joe Shaw Translational Post 2022/RFC1 worksheets/22-2268/22-2268.xlsx",
                         skip = 3,
                         sheet = "Data Sheet")

samples_with_locations <- automatic_reports %>%
  left_join(rfc1_epic_export, by = "full_name") %>%
  select(test_specimen_id, pt_first_nm, patient_surname,
         storage_location, final_dna_conc_ng_ml, interpretation) %>%
  filter(!is.na(test_specimen_id) &
           !test_specimen_id %in% ws_22_2268$Episode &
           interpretation != "pending")

samples_22_2325 <- dplyr::slice_sample(samples_with_locations, n = 35)

write.csv(samples_22_2325, "W:/MolecularGenetics/Neurogenetics/Research/Joe Shaw Translational Post 2022/RFC1 worksheets/22-2325/22_2325_samples.csv",
          row.names = FALSE)

##############################
