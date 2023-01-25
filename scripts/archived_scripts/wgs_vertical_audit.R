################################################################################
## Whole Genome Sequencing: Vertical Audit
################################################################################

##############################
# Load libraries
##############################

library(tidyverse)
library(readxl)
library(janitor)
library(lubridate)

##############################
# Load data and analysis
##############################

wgs_folder <- "I:/Genetics/DNA Lab/databases/Specialist_Services/WGS/"

wgs_analysis <- read_excel(path = paste0(wgs_folder, "WGS Analysis.xlsx"),
                          sheet = "Sample_list") %>%
  janitor::clean_names() %>%
  mutate(tat = difftime(reported_date, activation_date, units = "days"),
         tat_numeric = as.numeric(tat))

ngen_cases <- wgs_analysis %>%
  filter(analysis_lab == "NGEN" & !is.na(tat_numeric))

cases_to_inspect <- sample_n(ngen_cases, 3)

# Plot
wgs_analysis %>%
  filter(analysis_lab == "NGEN" & !is.na(tat_numeric)) %>%
  ggplot(aes(x = reorder(lab_number, tat_numeric), y = tat_numeric)) +
  geom_point() +
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        panel.grid = element_blank()) +
  geom_hline(yintercept = 42, linetype = "dashed") +
  labs(x = "Specimens", y = "Turnaround time (days)",
       title = "Neurogenetics WGS turnaround times")
  
##############################