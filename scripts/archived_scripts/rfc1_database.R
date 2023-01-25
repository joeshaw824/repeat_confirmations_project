################################################################################
## RFC1 Sample Testing Database
################################################################################

# Source this function to load the database and avoid code duplication.

library(tidyverse)
library(readxl)
library(janitor)

rfc1_db <-readxl::read_excel("W:/MolecularGenetics/Neurogenetics/Research/Joe Shaw Translational Post 2022/RFC1 R code/RFC1_analysis/data/RFC1_testing_database.xlsx",
                              sheet = "rfc1_db") %>%
    janitor::clean_names() %>%
  mutate(full_name = paste(toupper(forename), toupper(surname), sep = " "))

################################################################################