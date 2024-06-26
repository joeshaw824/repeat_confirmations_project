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

# Worksheets 22-3206 and 22-3382 were runs including only samples with positive research results.

ws_22_2268 <- read_rfc1_ws("22-2268")

ws_22_2325 <- read_rfc1_ws("22-2325")

ws_22_2543 <- read_rfc1_ws("22-2543")

ws_22_2649 <- read_rfc1_ws("22-2649")

ws_22_3206 <- read_rfc1_ws("22-3206")

ws_22_3382 <- read_rfc1_ws("22-3382")

ws_22_4370 <- read_rfc1_ws("22-4370")

ws_23_0124 <- read_rfc1_ws("23-0124")

collated_diagnostic_results <- rbind(ws_22_2268, ws_22_2325, ws_22_2543, 
                          ws_22_2649, ws_22_3206, ws_22_3382, ws_22_4370,
                          ws_23_0124) %>%
  mutate(full_name = paste(toupper(first_name), toupper(surname), sep = " ")) %>%
  # Gets rid of empty lines, but not samples with "NA" as report_type
  filter(!is.na(report_type))

positive_diagnostic_results <- collated_diagnostic_results %>%
  filter(coded_result %in% c("AAGGG expansion presumed homozygous", "Other")  &
           !base::duplicated(dna_no))

## Total tested

water_controls <- unique(grep("Water", collated_diagnostic_results$dna_no, ignore.case = TRUE,
                       value = TRUE))

results_for_summary <- collated_diagnostic_results %>%
  filter(!dna_no %in% water_controls) %>%
  mutate(repeat_sample = ifelse(base::duplicated(dna_no, fromLast = TRUE) |
                                  base::duplicated(dna_no, fromLast = FALSE),
                                "repeated",
                                "not repeated")) %>%
  filter(!(repeat_sample == "repeated" & coded_result %in% c("Fail", "Further analysis")))

summary_table <- results_for_summary %>%
  filter(!base::duplicated(dna_no)) %>%
  group_by(coded_result) %>%
  summarise(N = n()) %>%
  mutate("%" = round((N/sum(N)*100))) %>%
  arrange(desc(N))

sum(summary_table$N)

##############################
# AAGGG peak heights
##############################

rp3_peaks <- read_csv("data/RP3_files_highest_peak_calls.csv",
                      skip = 12) %>%
  janitor::clean_names() %>%
  mutate(specimen_id = str_extract(sample, pattern = "..RG-...G...."),
         worksheet = str_extract(sample, pattern = "22-...."),
         worksheet_number = as.numeric(substr(worksheet, 4, 7)))
  

rp3_peaks_results <- rp3_peaks %>%
  filter(!is.na(allele_number_1) & !is.na(specimen_id)) %>%
  left_join(collated_diagnostic_results %>%
              select(dna_no, worksheet, coded_result) %>%
              dplyr::rename(specimen_id = dna_no),
            by = c("specimen_id", "worksheet"))

# Colour by coded result
ggplot(rp3_peaks_results, aes(x = reorder(sample, desc(height_number_1)),
             y = height_number_1, colour = coded_result)) +
  geom_point(size = 3) +
  theme_bw() +
  theme(axis.text.x = element_blank(),
        panel.grid = element_blank()) +
  labs(x = "", y = "Heighest peak (RFU)",
       title = "Highest AAGGG peak heights in samples with expansions") +
  ylim(0, 20000) +
  geom_hline(yintercept = 1000, linetype ="dashed") +
  # This is the sample with an unusual repeat motif on the second allele
  geom_text(aes(label=ifelse(specimen_id == "20RG-323G0124", 
                             "20RG-323G0124",'')),
            hjust=2,vjust=0, size=3, colour = "black")
  
# Colour by worksheet
ggplot(rp3_peaks %>%
         filter(!is.na(allele_number_1)), aes(x = worksheet,
                              y = height_number_1, colour = worksheet)) +
  geom_jitter(size = 3) +
  theme_bw() +
  theme(axis.text.x = element_blank(),
        panel.grid = element_blank()) +
  labs(x = "", y = "Heighest peak (RFU)",
       title = "Fluorescence values on RFC1 worksheets are increasing over time",
       caption = "Data for AAGGG highest peak in samples with AAGGG expansions") +
  ylim(0, 20000)

pos_control_wells <- grep("108808", rp3_peaks$sample, value = TRUE)

rp3_peaks %>%
  filter(sample %in% pos_control_wells & worksheet != "22-2268") %>%
  ggplot(aes(x = reorder(sample, worksheet_number), y = height_number_1,
             colour = worksheet)) +
  geom_point(size = 3) +
  ylim(0, 9000) +
  theme_bw() +
  theme(panel.grid = element_blank()) +
  labs(x = "", y = "Heighest peak (RFU)",
       title = "Same sample dilution tested on 3 different worksheets")


##############################
# For update meeting: total tested and positives
##############################

# Total number of DNA samples tested
length(unique(collated_diagnostic_results$dna_no))

# Total with a positive or unusual result
length(unique(positive_diagnostic_results$dna_no))

cohort_table_1 <- collated_diagnostic_results %>%
  filter(!base::duplicated(dna_no) & coded_result != "Water blank clear") %>%
  group_by(coded_result) %>%
  summarise(N = n()) %>%
  mutate("%" = round((N/sum(N)*100))) %>%
  arrange(desc(N))

# 175 samples tested on 7 worksheets - 25 samples per worksheet

84/8

write.csv(cohort_table_1, "outputs/cohort_table_1.csv", row.names = FALSE)

# How many samples with AAGGG expansions have the strange dip pattern at every 4th repeat?
length(grep("Lower peaks every 4th peak later in the trace.", collated_diagnostic_results$result))

collated_diagnostic_results %>%
  filter(coded_result %in% c("Further analysis")) %>%
  select(worksheet, sample, dna_no, result) %>%
  filter(!base::duplicated(dna_no))

##############################
# How many samples could be reported currently?
##############################

# Action from update meeting on 21/10/2022

random_sample_set <- collated_diagnostic_results %>%
  filter(!base::duplicated(dna_no)) %>%
  filter(worksheet %in% c("22-2268", "22-2325", "22-2543", "22-2649")) %>%
  filter(report_type != "NA")

random_total <- nrow(random_sample_set)

random_sample_set_summary <- random_sample_set %>%
  group_by(report_type) %>%
  summarise(total = n()) %>%
  mutate(percent = round((total/random_total)*100,0)) %>%
  arrange(desc(total))

write.csv(random_sample_set_summary, "outputs/random_sample_set_summary.csv",
          row.names = FALSE)

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
  select(full_name, interpretation, consultant) %>%
  filter(!base::duplicated(full_name))

# Read in Epic information for DNA locations
rfc1_epic_export <- read_excel("data/Checking_CANVAS_referrals_20220902.xlsx",
                               sheet = "epic_output") %>%
  janitor::clean_names() %>%
  mutate(full_name = paste(toupper(pt_first_nm), toupper(patient_surname), sep = " ")) 

# Collated winpath referrals
winpath_rfc1_referrals <- read_csv("outputs/winpath_rfc1_referrals.csv") %>%
  janitor::clean_names() %>%
  mutate(full_name = paste0(forename, " ", surname))

# Research positives compared to diagnostic testing
research_vs_diagnostic <- positive_research_results %>%
  # Add on DNA locations
  left_join(rfc1_epic_export %>%
              filter(!base::duplicated(full_name)) %>%
              select(full_name, storage_location, 
                     test_specimen_id, dob), by = "full_name") %>%
  left_join(collated_diagnostic_results %>%
              filter(!base::duplicated(full_name)), by = "full_name") %>%
  # Add Winpath DNA numbers
  left_join(winpath_rfc1_referrals %>%
              select(dna_number, full_name), by = "full_name") %>%
  select(full_name, test_specimen_id, storage_location, interpretation, coded_result) %>%
  arrange(coded_result)

write.csv(research_vs_diagnostic, "outputs/research_vs_diagnostic.csv",
          row.names = FALSE)


research_vs_diagnostic %>%
  group_by(coded_result) %>%
  summarise(total = n())

# 29 still to be tested

##############################
# Sequencing normal controls
##############################

ws_22_2268_info <- get_ws_sample_info("22-2268")

ws_22_2325_info <- get_ws_sample_info("22-2325")

ws_22_2543_info <- get_ws_sample_info("22-2543")

ws_22_2649_info <- get_ws_sample_info("22-2649")

ws_22_3206_info <- get_ws_sample_info("22-3206")

ws_22_3382_info <- get_ws_sample_info("22-3382")

collated_info <- rbind(ws_22_2268_info, ws_22_2325_info, ws_22_2543_info,
                       ws_22_2649_info, ws_22_3206_info, ws_22_3382_info)


to_sequence <- c("20RG-325G0154", "20RG-286G0123", "20RG-318G0122", "20RG-293G0108", "20RG-273G0073", 
                 "20RG-218G0174", "21RG-183G0092", "20RG-321G0068", "21RG-212G0054", "20RG-307G0077", 
                 "21RG-219G0090", "22RG-075G0096", "20RG-276G0109", "21RG-166G0059", "21RG-204G0118",
                 "20RG-314G0053", "21RG-168G0006")

to_export <- collated_info %>%
  filter(episode %in% to_sequence) %>%
  dplyr::rename(dna_no = episode) %>%
  left_join(collated_diagnostic_results %>%
              select(dna_no, result), by = "dna_no", keep = FALSE)

write.csv(to_export, "outputs/for_sequencing.csv")

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

flanking <- collated_diagnostic_results %>%
  # Don't include runs which were full of expansion confirmations
  filter(worksheet %in% c("22-2268", "22-2325", "22-2543", "22-2649") &
           marker == "RFC1_FL") %>%
  select(dna_no, size_1, size_2, allele_1, allele_2) %>%
  pivot_longer(cols = -c(dna_no,allele_1, allele_2),
               names_to = "allele",
               values_to = "size_bp") %>%
  mutate(allele_size = round(size_bp, 0)) %>%
  filter(!is.na(allele_size) & !base::duplicated(dna_no))

num_alleles <- nrow(flanking)

# 293bp non-repeat sequence
# 150bp is the read length of whole genome sequencing
# 293+150 = 443bp

large_alleles <- nrow(flanking %>% 
                        filter(size_bp > 443))

large_allele_percentage <- (large_alleles / num_alleles) * 100
# 64% of alleles are greater than 443bp (~30 pentanucleotide repeats)

flanking_amplicon_abi_plot <- ggplot(flanking, aes(x = size_bp)) +
  geom_bar(stat = "bin", binwidth = 5) +
  labs(x = "Size (bp)", y = "Frequency (number of amplicons)",
       title = "RFC1 flanking PCR amplicon sizes: ABI-3730",
       subtitle = paste0(large_allele_percentage, "% of sized alleles are too large to be sequenced by a single WGS read (150bp)"),
       caption = "Data for randomly selected RFC1 referrals") +
  theme_bw() +
  theme(panel.grid = element_blank()) +
  xlim(250, 900) +
  geom_vline(xintercept = 443, linetype = "dashed")

ggsave(plot = flanking_amplicon_abi_plot, 
       filename = paste0(format(Sys.time(), "%Y%m%d_%H%M%S"),
                         "_",
                         "flanking_amplicons",
                         ".jpeg"),
       path = "plots/", 
       device='jpeg', 
       dpi=600,
       units = "cm",
       width = 15,
       height = 15)



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
# Tables for STR working group presentation - out of date
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
# Andrea's large spreadsheet

large_spreadsheet <- read_excel("data/RFC1 summary_AC Dec2021.xlsx") %>%
  janitor::clean_names()

interpretation_terms <- c("likely biallelic RFC1 expansion", "confirmed", "likley biallelic RFC1 expansion",
                          "confirmed 27/9/21", "confirmed, ad long read" )

english_centres <- c("Sheffield", "NHNN", "Salford Royal Hosptial", "Royal Victoria Infirmary Newcastle",
                     "South Tyneside Hospital", "Southmead Hospital Britstol", "James Cook University Hospital",
                     "Nuffield Health Hereford Hospital", "Salford Royal Hospital",
                     "Queen Elizabeth Hospital","Stepping Hill Hospital", "Norwich University Hospital",
                     "York Hospital", "Wessex Clinical Genetics", "Royal Hallamshire Hospital", NA)

test <- large_spreadsheet %>%
  filter(interpretation %in% interpretation_terms & centre %in% english_centres & !is.na(forename)) %>%
  mutate(full_name = paste(toupper(forename), toupper(surname), sep = " "))

test2 <- test %>%
  left_join(collated_diagnostic_results, by = "full_name") %>%
  # Samples without a result and with a Winpath/Epic DNA number
  filter(is.na(coded_result) & !is.na(internal_id)) %>%
  select(forename, surname.x, internal_id, external_id, hospital_no, interpretation)

# Do we have samples logged onto Epic already?
write.csv(test2, "outputs/test2.csv",
          row.names = FALSE)
