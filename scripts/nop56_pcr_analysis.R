################################################################################
## SCA36 PCR Analysis
################################################################################

##############################
# Load libraries
##############################

library(tidyverse)
library(readxl)
library(janitor)

setwd(dir = "W:/MolecularGenetics/Neurogenetics/Research/Joe Shaw Translational Post 2022/RFC1 R code/RFC1_analysis")

source("functions/rfc1_functions.R")

##############################
# Plot NOP56 normal allele size frequencies
##############################


ws_22_4357 <- readxl::read_excel(path = "W:/MolecularGenetics/Neurogenetics/Research/Joe Shaw Translational Post 2022/NOP56 worksheets/22-4357/22-4357.xlsx",
                                 sheet = "results_sheet",
                                 skip = 5) %>%
  janitor::clean_names()

allele_sizes <- ws_22_4357 %>%
  filter(!sample %in% c("02", "03", "05", "06")) %>%
  filter(marker != "NOP56_RP") %>%
  select(sample, dna_no, marker, allele_1, allele_2) %>%
  pivot_longer(cols = -c(sample, dna_no, marker),
               names_to = "allele",
               values_to = "size") %>%
  filter(!is.na(size))

ggplot(allele_sizes, aes(x = size, y = )) +
  geom_bar() +
  facet_wrap(~marker) +
  theme(axis.text.x = element_text(angle = 90)) +
  labs(x = "", y = "Frequency",
       title = "NOP56 GGCCTG repeats in normal alleles")

##############################