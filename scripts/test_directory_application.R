################################################################################
## RFC1 Test Directory Application
################################################################################

##############################
# Load libraries
##############################

library(tidyverse)
library(readxl)
library(janitor)
library(lubridate)

setwd(dir = "W:/MolecularGenetics/Neurogenetics/Research/Joe Shaw Translational Post 2022/RFC1 R code/RFC1_analysis")

referral_centres <- read.csv("resources/referring_centres.csv")

non_english_centres <- referral_centres %>%
  filter(location == "not england")

consistent_columns <- c("dna_number", "episode_number", "forename", "surname", 
                        "dob", "referring_centre", "collection_date")

##############################
# Winpath Referrals
##############################

# Winpath spreadsheet is split into 4 sections which are formatted differently.

winpath_spreadsheet_path <- "W:/MolecularGenetics/Neurogenetics/Research/CANVAS/Canvas list sent to Andrea/Winpath_samples/Combined Canvas List to Andrea from 2019 to 2021.xlsx"

winpath1 <- read_excel(path = winpath_spreadsheet_path,
                       range = "B1:F99") %>%
  janitor::clean_names() %>%
  dplyr::rename(
    "forename" = "first_name",
    "dob" = "do_b") %>%
  mutate(referring_centre = "",
         collection_date = format(as.POSIXct("2019-12-20"))) %>%
  select(all_of(consistent_columns))

winpath2 <- read_excel(path = winpath_spreadsheet_path,
                       range = "B101:K109") %>%
  janitor::clean_names() %>%
  mutate(collection_date = format(as.POSIXct(request_date), '%Y-%m-%d')) %>%
  dplyr::rename(
    "referring_centre" = "source") %>%
  select(all_of(consistent_columns))

winpath3 <- read_excel(path = winpath_spreadsheet_path,
                       range = "B111:K218") %>%
  janitor::clean_names() %>%
  mutate(collection_date = format(as.POSIXct(date_received), '%Y-%m-%d')) %>%
  dplyr::rename(
    "forename" = "first_name",
    "dob" = "do_b",
    "referring_centre" = "source") %>%
  select(all_of(consistent_columns))

winpath4 <- read_excel(path = winpath_spreadsheet_path,
                       range = "A220:K280",
                       col_types = c("text","text", "text", 
                                     "text", "text", "date", 
                                     "text", "text", "text", "text", "date")) %>%
  janitor::clean_names() %>%
  mutate(collection_date = format(as.POSIXct(booked_date), '%Y-%m-%d')) %>%
  filter(dna_number != "DNA Number") %>%
  dplyr::rename(
    "referring_centre" = "source") %>%
  select(all_of(consistent_columns))

##############################
# Epic Referrals
##############################

epic_referrals <- read_excel(path = "data/CANVAS_referrals_20220324.xlsx") %>%
  janitor::clean_names() %>%
  mutate(collection_date = format(as.POSIXct(collection_instant), '%Y-%m-%d')) %>%
  dplyr::rename(
    "episode_number" = "test_specimen_id",
    "surname" = "patient_surname",
    "forename" = "pt_first_nm",
    "referring_centre" = "submitter_ordering_dept") %>%
  mutate(dna_number = "") %>%
  select(all_of(consistent_columns))

##############################
# All RFC1 Referrals
##############################

rfc1_referrals <- rbind(
  epic_referrals, winpath1, winpath2, winpath3, winpath4) %>%
  mutate(
    forename = toupper(forename),
    surname = toupper(surname),
    name_dob = paste0(forename, " ", surname, " ", dob)) %>%
  filter(!base::duplicated(name_dob)) %>%
  filter(!referring_centre %in% c(non_english_centres$centre)) %>%
  mutate(collection_month = format(as.POSIXct(collection_date), '%Y-%m'),
         disease = "CANVAS")

# Application text
paste0("From December 2019 to March 2022 (28 months), there have been ",
        nrow(rfc1_referrals),
       " patient samples sent for this testing from neurology departments in England.")

##############################
# Tests requested
##############################

non_method_tests <- c("", 
                      "(RG) Specimen Reception", 
                      "NGEN_DNA to be Stored (R346.1)",
                      "DNA Extraction - Kurabo Large (Method)",
                      "DNA Extraction - Fuji (Method)",
                      "DNA Extraction - Fuji Large (Method)",
                      "DNA Extraction - Chemagic 360D (Method)",
                      "(Genetics) Inappropriate Request",
                      "(Genetics) Inappropriate Sample",
                      "(RG) Specimen Reception - WGS",
                      "DNA Extraction - WGS",
                      "DNA to be stored",
                      "GOSH DNA Extraction - Chemagic Star (Method)",
                      "(RG) Specimen Reception - WGS")

rfc1_tests <- read_excel("data/RFC1_tests_20220328.xlsx") %>%
  janitor::clean_names()

test_matrix <- stringr::str_split_fixed(rfc1_tests$test, "\r\n", n = 19)

test_df <- as.data.frame(test_matrix)

rfc1_long <- rfc1_tests %>%
  select(specimen_id) %>%
  cbind(test_df) %>%
  pivot_longer(cols = -specimen_id, 
               names_to = "column",
               values_to = "test_set") %>%
  filter(!test_set %in% non_method_tests)

rfc1_summary <- rfc1_long %>%
  group_by(specimen_id) %>%
  summarise(n = n()) %>%
  arrange(n)

# Text for application

paste0("For the ",
       length(unique(rfc1_long$specimen_id)),
       " samples for which such data is readily available, ",
       nrow(rfc1_summary %>%
              filter(n == 1)),
       " were referred with RFC1 as the only requested test.") 

##############################
# Comparison to other STR referrals
##############################

read_epic_list <- function(file_input, disease_input) {
  
  epic_list <- read_excel(path = paste0("data/", file_input, ".xlsx")) %>%
    janitor::clean_names() %>%
    mutate(disease = disease_input,
           collection_date = format(as.POSIXct(collection_instant), '%Y-%m-%d'),
           collection_month = format(as.POSIXct(collection_instant), '%Y-%m')) %>%
    dplyr::rename(
      "episode_number" = "test_specimen_id",
      "surname" = "patient_surname",
      "forename" = "pt_first_nm")
  
  return(epic_list)
}

HD_referrals <- read_epic_list("HD_20220330", "HD")

DRPLA_referrals <- read_epic_list("DRPLA_20220330", "DRPLA")

FRDA_referrals <- read_epic_list("FRDA_20220330", "FRDA")

C9orf72_referrals <- read_epic_list("C9orf72_20220330", "C9orf72")

SCA1_7_referrals <- read_epic_list("SCA1_7_20220330", "SCA1_7")

SBMA_referrals <- read_epic_list("SBMA_20220330", "SBMA")

JPH3_referrals <- read_epic_list("JPH3_20220330", "JPH3")

summarise_epic_list <- function(epic_list) {
  
  new_list <- epic_list %>%
    filter(collection_date >= "2019-12-31") %>%
    group_by(collection_month, disease) %>%
    summarise(monthly_total = n())
  
  return(new_list)
  
}

monthly_referrals <- rbind(summarise_epic_list(HD_referrals),
                summarise_epic_list(DRPLA_referrals),
                summarise_epic_list(FRDA_referrals),
                summarise_epic_list(rfc1_referrals),
                summarise_epic_list(C9orf72_referrals),
                summarise_epic_list(SBMA_referrals),
                summarise_epic_list(SCA1_7_referrals),
                summarise_epic_list(JPH3_referrals))


ggplot(monthly_referrals, aes(x = collection_month, y = monthly_total, group = disease, 
                              colour = disease)) +
  geom_point() +
  geom_line() +
  theme_bw() +
  theme(panel.grid = element_blank(),
        axis.text.x = element_text(angle = 45, vjust = 0.5)) +
  labs(x = "Month", y = "Samples received", 
       title = "UCLH Neurogenetics STR Referrals from January 2020")

##############################
# Have samples been requested?
##############################

# Requests (non-English centres included)
rfc1_requests <- rbind(
  epic_referrals, winpath1, winpath2, winpath3, winpath4) %>%
  mutate(
    forename = toupper(forename),
    surname = toupper(surname),
    name_dob = paste0(forename, " ", surname, " ", dob),
    full_name = paste(forename, surname, sep = " ")) %>%
  filter(!base::duplicated(full_name))

# Results provided by Cortese group
research_results <- read_excel(path = "data/RFC1 summary_AC Dec2021.xlsx",
                               sheet = "RFC1 all tested") %>%
  janitor::clean_names() %>%
  mutate(full_name = paste(toupper(forename), toupper(surname), sep = " ")) %>%
  filter(!base::duplicated(full_name))

# Check if results are available for the requested samples

requests_and_results <-  left_join(x = rfc1_requests, y = research_results, 
            by = "full_name",
            na_matches = "never") %>%
  select(dna_number, episode_number, full_name, flanking_pcr, aaggg, aaagg, aaaag,
         sb_result, interpretation, notes)

# Samples without results

nrow(requests_and_results %>%
       filter(is.na(flanking_pcr)))

##############################