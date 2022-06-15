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

setwd(dir = "W:/MolecularGenetics/Neurogenetics/Research/Joe Shaw Translational Post 2022/RFC1 R code/RFC1_analysis")

##############################
# Read sequencing traces
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
  
  return(sequence)
  
}

##############################
# Read forward sequencing on the negative strand
##############################

forward_seq_files <- list.files("data/rfc1_forward_sequencing")

get_repeat_seq <- function(seq_file){
  
  seq_output <- spgs::complement(get_sequence(seq_file, phred_threshold = 20), 
                                 case = "upper")
  
  final_output <- data.frame(file = seq_file, 
                    sequence = seq_output)
  
  return(final_output)
}

# Collate all forward data together

seq_df <- data.frame()

for (i in forward_seq_files) {
  
  temp_data <- get_repeat_seq(i)
  
  seq_df <- rbind(seq_df, temp_data)
  
  rm(temp_data)
  
}

mod_seq_df <- seq_df %>%
  mutate(repeat_sequence = sub(".*TACG", "", sequence))


write.csv(mod_seq_df, "outputs/forward_sequencing.csv",
          row.names = FALSE)

##############################


list.files("data/sanger_sequencing")

get_sequence("1277-06-upper-R_G08.txt", phred_threshold = 0)
