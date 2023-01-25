################################################################################
## Sanger Quality Control Analysis
################################################################################

##############################
# Load libraries
##############################

library(tidyverse)
library(readxl)
library(janitor)

##############################
# Functions
##############################

# Function to get Excel filepaths for each folder
get_excel_filepaths <- function(input_filepath) {
  
  excel_files <- list.files(input_filepath, pattern = ".xlsx")
  
  folder_excel_paths <- data.frame("path"= c())
  
  for (j in excel_files) {
    
    temp_excel_filepath <- data.frame("path" = paste0(input_filepath, j))
    
    folder_excel_paths <- rbind(folder_excel_paths, temp_excel_filepath)
    
    rm(temp_excel_filepath)
  }
  
  return(folder_excel_paths)
  
}

# Function for reading sequencing results sheet
read_seq_excel <- function(input_excel_path) {
  
  excel_tbl <- read_excel(path = input_excel_path, 
                          skip = 5) %>%
    janitor::clean_names()
  
  #stopifnot(colnames(excel_tbl) == seq_table_cols)
  
  return(excel_tbl)
  
}

##############################
# Get filepaths
##############################

ngen_seq_filepath <- "W:/MolecularGenetics/Mutation Surveyor/Neurogenetics/Runs/2022/"

seq_folder_names <- list.files(ngen_seq_filepath)

# For files in these filepaths, list the ones with .xlsx filetypes

seq_folder_paths <- data.frame("path" = c())

# Get all the folder filepaths
for (i in seq_folder_names) {
  
  tmp_filepath <- data.frame("path" = paste0(ngen_seq_filepath, i, "/"))
  
  seq_folder_paths <- rbind(seq_folder_paths, tmp_filepath)
  
  rm(tmp_filepath)
  
}

# Get all Excel filepaths and collate

all_excel_filepaths <- data.frame("path"= c())

for (i in seq_folder_paths$path) {
  
  temp_df <- get_excel_filepaths(i)
  
  all_excel_filepaths <- rbind(all_excel_filepaths, temp_df)
  
  rm(temp_df)
  
}

seq_excels <- all_excel_filepaths %>%
  # Remove the QC files
  filter(!path %in% grep("qc", path, ignore.case = TRUE,
               value = TRUE) &
           
           !path %in% grep("~", path, ignore.case = TRUE,
                           value = TRUE))

##############################
# Read and collate sequencing Excel files
##############################

all_seq_results <- data.frame()

for (l in seq_excels$path) {
  
  seq_table_cols <- c("w_s_no", "well_no", "patient_name", "gene_exon", "raw_data",
                      "results_1st_read", "results_2nd_read")
  
  tmp_excel <- read_seq_excel(l)
  
  if (length(setdiff(colnames(tmp_excel), seq_table_cols)) == 0) {
    
    all_seq_results <- rbind(all_seq_results, tmp_excel)
  }

  rm(tmp_excel)

}

first_read_fail <- unique(grep("fail", all_seq_results$results_1st_read, ignore.case = TRUE, value = TRUE))

second_read_fail <- unique(grep("fail", all_seq_results$results_2nd_read, ignore.case = TRUE, value = TRUE))

all_seq_results_edit <- all_seq_results %>%
  # Remove empty rows
  filter(!w_s_no %in% c("Do not use", "0", NA)) %>%
  mutate(status = case_when(
    results_1st_read %in% first_read_fail ~"fail",
    results_2nd_read %in% second_read_fail ~"fail",
    TRUE ~"pass"))

# Stacked bar plot
fail_bar_plot <- all_seq_results_edit %>%
  group_by(w_s_no, status) %>%
  summarise(total = n()) %>%
  ggplot(aes(x = w_s_no, y = total, fill = status))+
  scale_fill_manual(values = c("#666666", "#FFFFFF")) +
    geom_col(colour = "black") +
    theme_bw() +
    theme(axis.text.x = element_blank(),
          panel.grid = element_blank()) +
    labs(x = "Worksheets", y = "Number of wells",
         title = "Neurogenetics robot Sanger worksheets in 2022")

ggsave(plot = fail_bar_plot, 
       filename = paste0(format(Sys.time(), "%Y%m%d"),
                         "_",
                         "fail_bar",
                         ".jpeg"),
       path = "plots/", 
       device='jpeg', 
       dpi=600,
       units = "cm",
       width = 24,
       height = 10)


# Which exons fail most often?
all_seq_results_edit %>%
  filter(status == "fail") %>%
  group_by(gene_exon) %>%
  summarise(total = n()) %>%
  arrange(desc(total))

# Which samples fail most often?
all_seq_results_edit %>%
  filter(status == "fail") %>%
  group_by(patient_name) %>%
  summarise(total = n()) %>%
  arrange(desc(total))

fail_position <- all_seq_results_edit %>%
  filter(status == "fail") %>%
  group_by(well_no) %>%
  summarise(total = n()) %>%
  ggplot(aes(x = well_no, y = total))+
  #scale_fill_manual(values = c("#666666", "#FFFFFF")) +
  geom_col(colour = "black") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90),
        panel.grid = element_blank()) +
  labs(x = "Well position", y = "Number of failed wells",
       title = "Frequency of fails by ABI well position")

ggsave(plot = fail_position, 
       filename = paste0(format(Sys.time(), "%Y%m%d"),
                         "_",
                         "fail_position",
                         ".jpeg"),
       path = "plots/", 
       device='jpeg', 
       dpi=600,
       units = "cm",
       width = 24,
       height = 10)

#################################################################################