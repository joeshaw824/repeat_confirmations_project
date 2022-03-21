################################################################################
## Loading and Cleaning RFC1 research results
################################################################################

##############################
# Load libraries
##############################

library(tidyverse)
library(readxl)
library(janitor)
library(ggplot2)

setwd(dir = "W:/MolecularGenetics/Neurogenetics/Research/Joe Shaw Translational Post 2022/RFC1 R code/RFC1_analysis")

##############################
# Load research results
##############################

research_results <- read_excel(path = "data/AC_CANVAS_Screeninglist_2021_GOSH.xlsx") %>%
  janitor::clean_names()

##############################
# Standardise results
##############################

flanking_noamp_variants <- c("No PCR product", "no PCR product", "no PCR products", "no bands",
                             "no PCR products (checked twice)", "no pCR products", "Negative", "no product")

flanking_amp_variants <- c("reference", "Intermediate", "Intermediate - faint", "intermediate",
                           "2 bands", "reference (rechecked, first time no PCR products)",
                           "reference (faint band)", "reference (re-checked)", "Positive",
                           "product")

aaggg_positive_variants <- c("positive", "pos", "Positive")

aaggg_negative_variants <- c("negative", "neg", "Negative")

##############################
# Clean research results
##############################

cleaned_results <- research_results %>%
  
  # Remove duplicates
  filter(!base::duplicated(name_string)) %>%
  
  # Remove "patient names" containing numbers
  filter(!name_string %in% grep("1|2|3|4|5|6|7|8|9|0", 
                                research_results$name_string, value = TRUE)) %>%
  
  mutate(
    flanking_clean = case_when(
      flanking_pcr %in% flanking_noamp_variants ~"no amplification",
      flanking_pcr %in% flanking_amp_variants ~"amplification",
      TRUE ~"other"),
    
    aaggg_clean = case_when(
      aaggg %in% aaggg_positive_variants ~"positive",
      aaggg %in% aaggg_negative_variants ~"negative",
      TRUE ~"other"),
    
    result_clean = case_when(
      flanking_clean == "no amplification" & aaggg_clean == "positive" ~"biallelic AAGGG repeat expansion detected",
      flanking_clean == "amplification" & aaggg_clean == "negative" ~"biallelic AAGGG repeat expansion not detected",
      TRUE ~"other"),
    
    external_id_new = case_when(
      is.na(external_id) ~ " ",
      TRUE ~external_id),
    
    internal_id_new = case_when(
      is.na(internal_id) ~ " ",
      TRUE ~internal_id),
    
    hospital_no_new = case_when(
      is.na(hospital_no) ~ " ",
      TRUE ~hospital_no),
    
    identifiers = paste(external_id_new, internal_id_new, hospital_no_new, sep = " ")
    )

##############################
# Southern blotting alleles
##############################

patient_alleles <- cleaned_results %>%
  select(name_string, large, small) %>%
  pivot_longer(
    cols = c("large", "small"),
    names_to = "allele",
    values_to = "size") %>%
  filter(!is.na(size))

ggplot(patient_alleles, aes(x = reorder(name_string, size), y = size)) +
  geom_point() +
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank()) +
  labs(y = "Allele pentanucleotide repeats", x = "Sample", 
       title = "Southern blotting allele sizes for RFC1")

