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

setwd(dir = "W:/MolecularGenetics/Neurogenetics/Research/Joe Shaw Translational Post 2022/RFC1 R code/RFC1_analysis")

source("scripts/rfc1_database.R")

##############################
# Read in worksheets
##############################

read_rfc1_ws <- function(worksheet_number) {
  
  output <- readxl::read_excel(paste0("W:/MolecularGenetics/Neurogenetics/Research/Joe Shaw Translational Post 2022/RFC1 worksheets/", 
  worksheet_number, "/", worksheet_number, ".xlsx"),
                     sheet = "results_sheet",
                     skip = 2) %>%
    janitor::clean_names() %>%
    select(sample, dna_no, first_name, surname, sample_file, marker, allele_1, size_1, height_1, 
           allele_2, size_2, height_2, result, report_type, coded_result)
}

ws_22_2268 <- read_rfc1_ws("22-2268")

ws_22_2325 <- read_rfc1_ws("22-2325")

ws_22_2543 <- read_rfc1_ws("22-2543")

collated_results <- rbind(ws_22_2268, ws_22_2325, ws_22_2543) %>%
  filter(!is.na(report_type)) %>%
  mutate(full_name = paste(toupper(first_name), toupper(surname), sep = " "))

##############################
# Comparing results with research results
##############################

automatic_reports <- read_excel( path = "data/AC_CANVAS_Screeninglist_2021_GOSH.xlsx") %>%
  janitor::clean_names() %>%
  mutate(full_name = paste(toupper(first_name), toupper(surname), sep = " ")) %>%
  filter(!base::duplicated(full_name))

result_comparison <- collated_results %>%
  left_join(automatic_reports, by = "full_name") %>%
  select(full_name, result, coded_result, report_type, interpretation, interpretation_2, flanking, aaaag, aaagg,
         aaggg) %>%
  filter(report_type == "Consistent with RFC1 disorder")

##############################
# Plotting RFC1 amplicon sizes
##############################

flanking <- collated_results %>%
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
# 38% of alleles are greater than 400bp

flanking_amplicon_plot <- ggplot(flanking, aes(x = size_bp)) +
  geom_bar(stat = "bin", binwidth = 5) +
  labs(x = "Flanking PCR amplion (bp) - measured on ABI-3730", y = "Frequency",
       title = "RFC1 flanking PCR amplicon sizes") +
  theme_bw()

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

##############################