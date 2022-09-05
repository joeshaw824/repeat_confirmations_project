################################################################################
## RFC1 DNA Aliquoting by Mark Gaskin
################################################################################

##############################
# Load libraries
##############################

library(tidyverse)
library(readxl)
library(janitor)

setwd(dir = "W:/MolecularGenetics/Neurogenetics/Research/Joe Shaw Translational Post 2022/DNA_aliquots_for_research")

##############################
# Load pullsheets
##############################

T1611 <- read_excel(path = "T1611-Pull sheet.xlsx",
                             sheet = "Pull sheet",
                             skip = 2) %>%
  janitor::clean_names() %>%
  filter(!is.na(original_sample_id)) %>%
  mutate(epic_comment = paste0(amount_taken_ul, "ul aliquotted for Dr Andrea Cortese at the Institute of Neurology by Mark Gaskin on 11/08/2022 (batch T1611)"))

T1630 <- read_excel(path = "T1630-Pull sheet.xlsx",
                    sheet = "Pull sheet",
                    skip = 2) %>%
  janitor::clean_names() %>%
  filter(!is.na(original_sample_id)) %>%
  mutate(epic_comment = paste0(amount_taken_ul, "ul aliquotted for Dr Pietro Fratta at the Institute of Neurology by Mark Gaskin on 25/08/2022 (batch T1630)"))

##############################
# Add comments
##############################

epic_comments_for_update <- T1630 %>%
  select(original_sample_id, epic_comment) %>%
  rbind(T1611 %>%
          select(original_sample_id, epic_comment)) %>%
  filter(grepl("RG", original_sample_id)) %>%
  arrange(original_sample_id)

##############################
# Export
##############################

write.csv(epic_comments_for_update,
          file = "epic_comments_for_update.csv",
          row.names = FALSE)

##############################