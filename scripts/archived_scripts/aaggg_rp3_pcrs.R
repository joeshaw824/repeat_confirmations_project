################################################################################
## RFC1 AAGGG Repeat Primed PCR Analysis
################################################################################

# 23/09/2022: I have a noticed a strange pattern in the AAGGG RP-PCR (RP3) for
# some samples which appear to have expansions.
# There is a decrease in peak signal with 4 repeat periodicity:
# i.e. peak heights follow the pattern:
# normal, normal, normal, low, normal, normal, normal, low etc

library(tidyverse)
library(janitor)

setwd(dir = "W:/MolecularGenetics/Neurogenetics/Research/Joe Shaw Translational Post 2022/RFC1 R code/RFC1_analysis")

rp3_results <- read_csv("data/RP3_files_meta_analysis_AlleleReport.csv",
                        skip = 12) %>%
  janitor::clean_names() %>%
  select(-x1) %>%
  mutate(specimen_id = str_extract(sample, pattern = "..RG-...G...."))

# Add this step in to exclude samples
samples_to_exclude <- c(
  # Weird stutter pattern, like a bite taken out of an apple.
  "21RG-202G0062", 
  # Seems to have a smaller expansion and a larger one.
  "20RG-318G0185")

rp3_expansions <- rp3_results %>%
  filter(!is.na(marker)) %>%
  mutate(pattern_repeat_start = as.numeric(gsub("rpts", "", allele_number_1)),
         pattern_repeat_end = as.numeric(gsub("rpts", "", allele_number_2)),
         category = ifelse(allele_number_1 == "Pileup peak", "No interruptions", 
                           "Interruption pattern")) %>%
  filter(!specimen_id %in% samples_to_exclude)

# How many samples have this pattern versus don't?

rp3_expansions %>%
  group_by(category) %>%
  summarise(total = n())

# Where does this pattern start and stop?

rp3_expansions %>%
  filter(category == "Interruption pattern") %>%
  select(sample, pattern_repeat_start, pattern_repeat_end) %>%
  pivot_longer(cols = -sample,
               names_to = "allele",
               values_to = "repeat_number") %>%
  ggplot(aes(x = repeat_number, y = , fill = allele)) +
  geom_bar() +
  xlim(1, 160) +
  labs(title = "Unusual pattern position in RP3 expansion data")
  

# Is the pattern seen in the same samples when repeated?

duplicates <- rp3_expansions %>%
  filter(base::duplicated(specimen_id, fromLast = TRUE) |
           base::duplicated(specimen_id, fromLast = FALSE))

################################################################################