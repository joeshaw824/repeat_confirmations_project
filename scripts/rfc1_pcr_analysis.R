################################################################################
## RFC1 PCR Analysis
################################################################################

##############################
# Load libraries
##############################

library(tidyverse)
library(readxl)
library(janitor)
library(lubridate)
library(ggpubr)

setwd(dir = "W:/MolecularGenetics/Neurogenetics/Research/Joe Shaw Translational Post 2022/RFC1 R code/RFC1_analysis")

source("functions/rfc1_functions.R")

##############################
# Read and collate worksheets
##############################

ws_22_2268 <- read_rfc1_ws("22-2268")

ws_22_2325 <- read_rfc1_ws("22-2325")

ws_22_2543 <- read_rfc1_ws("22-2543")

ws_22_2649 <- read_rfc1_ws("22-2649")

ws_22_3206 <- read_rfc1_ws("22-3206")

ws_22_3382 <- read_rfc1_ws("22-3382")

collated_diagnostic_results <- rbind(ws_22_2268, ws_22_2325, ws_22_2543, 
                          ws_22_2649, ws_22_3206, ws_22_3382) %>%
  filter(!report_type %in% c(NA, "Water control") ) %>%
  mutate(full_name = paste(toupper(first_name), toupper(surname), sep = " "))

positive_diagnostic_results <- collated_diagnostic_results %>%
  filter(coded_result %in% c("AAGGG expansion presumed homozygous", "Other")  &
           !base::duplicated(dna_no))

##############################
# Total tested and positives
##############################

# Total number of DNA samples tested
length(unique(collated_diagnostic_results$dna_no))

# Total with a positive or unusual result
length(unique(positive_diagnostic_results$dna_no))

##############################
# How many do we still need to confirm?
##############################

# Updated spreadsheet from Riccardo Curro
updated_research_results <- read_excel(path = "data/CANVAS_Screeninglist_update_Aug2022_GOSH.xlsx") %>%
  janitor::clean_names() %>%
  mutate(full_name = paste(toupper(first_name), toupper(surname), sep = " "))

positives <- c("biallelic RFC1 expansion", "likely biallelic RFC1 expansion", 
               "likely biallelic RFC1 expansion (rechecked by Joe)", "confirmed",
               "already confirmed", "patient already confirmed")

positive_research_results <- updated_research_results %>%
  filter(interpretation %in% positives) %>%
  select(full_name, interpretation, consultant)

length(unique(positive_research_results$full_name))

# Read in Epic information for DNA locations
rfc1_epic_export <- read_excel("data/Checking_CANVAS_referrals_20220902.xlsx",
                               sheet = "epic_output") %>%
  janitor::clean_names() %>%
  mutate(full_name = paste(toupper(pt_first_nm), toupper(patient_surname), sep = " ")) 

# Collated winpath referrals
winpath_rfc1_referrals <- read_csv("outputs/winpath_rfc1_referrals.csv") %>%
  janitor::clean_names() %>%
  mutate(full_name = paste0(forename, " ", surname))

# Check which samples have already been confirmed
to_manally_check <- split_dna_location(positive_research_results %>%
  # Add on DNA locations
  left_join(rfc1_epic_export %>%
              select(full_name, storage_location, 
                     test_specimen_id, dob), by = "full_name") %>%
  left_join(collated_diagnostic_results, by = "full_name")) %>%
  arrange(tray) %>%
  filter(!base::duplicated(full_name))%>%
  filter(is.na(coded_result)) %>%
  # Add Winpath DNA numbers
  left_join(winpath_rfc1_referrals %>%
              select(dna_number, full_name), by = "full_name") %>%
  mutate(query_already_tested = "",
         worksheet_if_tested = "",
         dna_volume = "",
         query_other_sample_available = "",
         other_sample_id = "") %>%
  select(full_name, dob, test_specimen_id, dna_number, tray, ycoord, xcoord, query_already_tested,
         worksheet_if_tested, dna_volume, query_other_sample_available, other_sample_id)

write.csv(to_manally_check, "outputs/manually_checking_dnas.csv", row.names = FALSE)

##############################
# Comparing results with research results - out of date
##############################

automatic_reports <- read_excel( path = "data/AC_CANVAS_Screeninglist_2021_GOSH.xlsx") %>%
  janitor::clean_names() %>%
  mutate(full_name = paste(toupper(first_name), toupper(surname), sep = " ")) %>%
  filter(!base::duplicated(full_name))

control_results <- data.frame(
  # Removed patient names as these should not go on Github.
  full_name = c(),
  interpretation_2 = c("positive", "carrier", "negative"))

result_comparison1 <- collated_results %>%
  filter(!base::duplicated(full_name)) %>%
  inner_join(automatic_reports, by = "full_name") %>%
  select(full_name, report_type, interpretation_2) 

result_comparison2 <- collated_results %>%
  filter(!base::duplicated(full_name)) %>%
  inner_join(control_results, by = "full_name") %>%
  select(full_name, report_type, interpretation_2)

neg_results <- c("negative", "carrier")
pos_results <- c("positive", "VUS (please contact me)")

result_comparison <- rbind(result_comparison1, result_comparison2) %>%
  mutate(category = case_when(
    report_type == "RFC1 disorder not confirmed" & interpretation_2 %in% neg_results ~"TN",
    report_type == "RFC1 disorder not confirmed" & interpretation_2 %in% pos_results ~"FN",
    report_type == "Consistent with RFC1 disorder" & interpretation_2 %in% pos_results ~"TP",
    report_type == "Consistent with RFC1 disorder" & interpretation_2 %in% neg_results ~"FP")) 

result_table <- result_comparison %>%
  group_by(category) %>%
  summarise(total = n())

##############################
# Plotting RFC1 amplicon sizes
##############################

flanking <- collated_results %>%
  # Need to filter repeated samples
  filter(marker == "RFC1_FL") %>%
  select(dna_no, size_1, size_2) %>%
  pivot_longer(cols = -dna_no,
               names_to = "allele",
               values_to = "size_bp") %>%
  mutate(allele_size = round(size_bp, 0)) %>%
  filter(!is.na(allele_size))

num_alleles <- nrow(flanking)

large_alleles <- nrow(flanking %>% 
                        filter(size_bp > 400))

(large_alleles / num_alleles) * 100
# 50% of alleles are greater than 400bp (~21 pentanucleotide repeats)
# (400-293) / 5 = 21.4 repeats

flanking_amplicon_abi_plot <- ggplot(flanking, aes(x = size_bp)) +
  geom_bar(stat = "bin", binwidth = 5) +
  labs(x = "", y = "Frequency",
       title = "RFC1 flanking PCR amplicon sizes- ABI-3730") +
  theme_bw() +
  xlim(250, 2100) +
  geom_vline(xintercept = 860, linetype = "dashed")

gel_sizes <- read_excel("data/rfc1_flanking_pcr_gel_sizes.xlsx",
                        col_types = c("text", "text", "text", "numeric", "numeric"))

flanking_amplicon_gel_plot <- gel_sizes %>%
  select(-c(worksheet, tube)) %>%
  pivot_longer(cols = -sample,
               names_to = "allele",
               values_to = "size_bp") %>%
  filter(!is.na(size_bp)) %>%
  ggplot(aes(x = size_bp)) +
  geom_bar() +
  theme_bw() +
  xlim(250, 2100) +
  labs(x = "Flanking PCR amplion (bp)", y = "Frequency",
       title = "RFC1 flanking PCR amplicon sizes - agarose gel") +
  geom_vline(xintercept = 860, linetype = "dashed")

sizes_plot <- ggpubr::ggarrange(flanking_amplicon_abi_plot,
                  flanking_amplicon_gel_plot,
                  ncol = 1, nrow = 2, align = "v")

##############################
# Tables for STR working group presentation
##############################

rfc1_table1 <- collated_results %>%
  filter(!is.na(coded_result)) %>%
  filter(!base::duplicated(dna_no)) %>%
  group_by(coded_result) %>%
  summarise(N = n()) %>%
  mutate("%" = round((N/sum(N)*100))) %>%
  arrange(desc(N))

rfc1_table2 <- collated_results %>%
  filter(!is.na(report_type)) %>%
  filter(!base::duplicated(dna_no)) %>%
  group_by(report_type) %>%
  summarise(N = n()) %>%
  mutate("%" = round((N/sum(N)*100))) %>%
  arrange(desc(N))

write.csv(rfc1_table1, "outputs/rfc1_table1.csv", row.names = FALSE)

write.csv(rfc1_table2, "outputs/rfc1_table2.csv", row.names = FALSE)

collated_results %>%
  filter(coded_result == "Other") %>%
  select(dna_no, result)

##############################