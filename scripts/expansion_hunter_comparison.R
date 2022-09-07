################################################################################
## Comparison of ExpansionHunter and Orthogonal Testing
################################################################################

##############################
# Load libraries
##############################

library(tidyverse)
library(readxl)
library(janitor)

setwd(dir = "W:/MolecularGenetics/Neurogenetics/Research/Joe Shaw Translational Post 2022/RFC1 R code/RFC1_analysis")

##############################
# Load data
##############################

comparison_data <- read_excel(path = "data/combined_othogonal_comparison.xlsx",
                              sheet = "combined_othogonal_output") %>%
  janitor::clean_names()


##############################
# Plot data
##############################

allele_1_results <- comparison_data %>%
  select(-c(pcr_allele_size_2, eh_min_ci_allele_2, eh_allele_2, eh_max_ci_allele_2, 
            pcr_eh_max_ci_allele_2, pcr_eh_max_ci_allele_1, investigation)) %>%
  dplyr::rename(pcr_allele_size = pcr_allele_size_1,
                eh_min_ci = eh_min_ci_allele_1,
                eh_allele = eh_allele_1,
                eh_max_ci = eh_max_ci_allele_1)

allele_2_results <- comparison_data %>%
  select(-c(pcr_allele_size_1, eh_min_ci_allele_1, eh_allele_1, eh_max_ci_allele_1, 
            pcr_eh_max_ci_allele_2, pcr_eh_max_ci_allele_1, investigation)) %>%
  dplyr::rename(pcr_allele_size = pcr_allele_size_2,
                eh_min_ci = eh_min_ci_allele_2,
                eh_allele = eh_allele_2,
                eh_max_ci = eh_max_ci_allele_2)

longer_results <- rbind(allele_1_results, allele_2_results)

longer_results %>%
  ggplot(aes(x = pcr_allele_size, y = eh_allele)) +
  geom_point(alpha = 0.5) +
  geom_errorbar(aes(ymin = eh_min_ci, ymax = eh_max_ci)) +
  geom_abline(linetype="dashed") +
  facet_wrap(~repeat_id)

##############################
