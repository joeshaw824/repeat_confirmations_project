################################################################################
## RFC1 Allele Sizes
################################################################################

##############################
# Load libraries
##############################

library(tidyverse)
library(readxl)
library(janitor)

setwd(dir = "W:/MolecularGenetics/Neurogenetics/Research/Joe Shaw Translational Post 2022/RFC1 R code/RFC1_analysis")

##############################
# Load Southern blot data
##############################

southern_blot_data <- read_excel(path = "data/RFC1 summary_AC Dec2021.xlsx",
                               sheet = "RFC1 all tested") %>%
  janitor::clean_names() %>%
  mutate(name_string = paste0(forename, "_", surname),
         name_string = toupper(name_string)) %>%
  filter(name_string != "NA_NA")

non_uk <- c("Paris", "Israel", "Genova", "Torino", "Sydney", "Berciano", "Waldemann",
            "Tubingen", "Pavia", "Padova", "Swisse", "Napoli", "Auckland City Hospital",
            "Hospital Roger Salengro", "Finland", "Halle", "France")

##############################
# Southern blotting alleles
##############################

patient_alleles <- southern_blot_data %>%
  filter(!centre %in% non_uk) %>%
  mutate(small_bp = round((small*5) + 293, 0),
         large_bp = round((large*5) + 293, 0))
  
alleles_longer <- patient_alleles %>%
  select(name_string, large_bp, small_bp) %>%
  pivot_longer(
    cols = c("large_bp", "small_bp"),
    names_to = "allele",
    values_to = "size") %>%
  filter(!is.na(size))

mean(alleles_longer$size)

# Scatter plot
ggplot(alleles_longer, aes(x = reorder(name_string, size), y = size)) +
  geom_point() +
  theme_bw() +
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank()) +
  labs(y = "Amplicon size (bp)", x = "Patient samples", 
       title = "Southern blotting allele sizes in NHS patients")

# Boxplot
ggplot(alleles_longer, aes(y = size)) +
  geom_boxplot() +
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank()) +
  labs(y = "Amplicon size (bp)", x = "Sample", 
       title = "Southern blotting allele sizes in NHS patients") +
  ylim(0, 30000)

# Violin plot
ggplot(alleles_longer, aes(x = "sample", y = size)) +
  geom_violin() +
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank()) +
  labs(y = "Amplicon size (bp)", x = "Sample", 
       title = "Southern blotting allele sizes in NHS patients") +
  ylim(0, 30000)

##############################
