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

#source("scripts/rfc1_database.R")

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

ws_22_2649 <- read_rfc1_ws("22-2649")

collated_results <- rbind(ws_22_2268, ws_22_2325, ws_22_2543, ws_22_2649) %>%
  filter(!is.na(report_type)) %>%
  mutate(full_name = paste(toupper(first_name), toupper(surname), sep = " "))

##############################
# Positive results for sequencing
##############################

positive_results <- collated_results %>%
  filter(coded_result == "AAGGG expansion presumed homozygous" &
           !base::duplicated(dna_no)) %>%
  select(dna_no, first_name, surname)

write.csv(positive_results, 
          "W:/MolecularGenetics/Neurogenetics/Research/Joe Shaw Translational Post 2022/RFC1 worksheets/22-2998/positive_results.csv",
          row.names = FALSE)
        
##############################
# Comparing results with research results
##############################

automatic_reports <- read_excel( path = "data/AC_CANVAS_Screeninglist_2021_GOSH.xlsx") %>%
  janitor::clean_names() %>%
  mutate(full_name = paste(toupper(first_name), toupper(surname), sep = " ")) %>%
  filter(!base::duplicated(full_name))

control_results <- data.frame(
  full_name = c("ALISON CLARKE", "GERALDINE JONES", "IRMA BALDENWEG"),
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
# Uncertainty of measurement
##############################

collated_results %>%
  filter(dna_no == "109437")



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
# 38% of alleles are greater than 400bp

flanking_amplicon_abi_plot <- ggplot(flanking, aes(x = size_bp)) +
  geom_bar(stat = "bin", binwidth = 5) +
  labs(x = "", y = "Frequency",
       title = "RFC1 flanking PCR amplicon sizes- ABI-3730") +
  theme_bw() +
  xlim(250, 2100)

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
       title = "RFC1 flanking PCR amplicon sizes - agarose gel")

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