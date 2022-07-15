################################################################################
## RFC1 Sanger sequencing
################################################################################

##############################
# Load libraries
##############################

library(tidyverse)
library(janitor)
# spgs for complement function
library(spgs)

# Reference file – negative strand
# Traces exported from Mutation Surveyor

setwd(dir = "W:/MolecularGenetics/Neurogenetics/Research/Joe Shaw Translational Post 2022/RFC1 R code/RFC1_analysis")

sequencing_ids <- read.csv("data/sequencing_positions.csv")

##############################
# Read sequencing traces function
##############################

get_sequence <- function(seq_file, phred_threshold = 30,
                         filepath = "data/sanger_sequencing/") {
  
  seq_data <- read.table(paste0(filepath, seq_file),
                         skip=1,
                         sep = "\t",
                         header = TRUE) %>%
    janitor::clean_names() %>%
    filter(phred_score > phred_threshold)
  
  sequence <- paste(seq_data$base_call, collapse = "")
  
  output <- data.frame("file" = seq_file, "sequence" = sequence)
  
  return(output)
  
}

##############################
# Collate sequence data
##############################

seq_files <- list.files("data/sanger_sequencing")

seq_df <- data.frame()

for (i in seq_files) {
  
  temp_data <- get_sequence(i, phred_threshold = 20)
  
  seq_df <- rbind(seq_df, temp_data)
  
  rm(temp_data)
  
}

##############################
# Modify sequence data and create complement
##############################

seq_df_mod <- seq_df %>%
  mutate(direction = case_when(
    str_detect(file, "R") ~"reverse",     
    TRUE ~"forward"),
    worksheet = paste0("22-", substr(file, 1, 4)),
    position = as.numeric(substr(file, 6,7)),
    sequence_complement = spgs::complement(sequence, case="upper")) %>%
  left_join(sequencing_ids, by = c("worksheet", "position"))

##############################
# Export complement sequences
##############################

complement_sequences <- seq_df_mod %>%
  select(worksheet, position, direction, file, episode, patient, sequence_complement)

write.csv(complement_sequences, "outputs/complement_sequences.csv",
          row.names = FALSE)

##############################