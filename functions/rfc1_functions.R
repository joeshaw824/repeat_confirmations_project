################################################################################
## Useful functions
################################################################################

##############################
# Epic functions
##############################

# Useful function for splitting DNA storage location string used in Epic

# Regex notes:
# \\s = space
# \\d = digits
# \\d{3} = 3 digits
# \\d{1,2} = at least 1 digit and at most 2 digits
# () = grouping, i.e. select this group

split_dna_location <- function(input_df) {
    
  stopifnot("storage_location" %in% colnames(input_df))
  
  new_df <- input_df %>%
    
    mutate(
      # Remove line breaks
      location_stripped = gsub("[\r\n]", "", storage_location),
    
      # Get tray coordinate
      tray = as.numeric(sub(pattern = ".*?RG DNA Tray\\s(.*?)\\sslot.*", 
                 replacement = "\\1", 
                 x = location_stripped,
                 # ignore.case required as locations can have "Tray" or "tray"
                 ignore.case = TRUE)),
      
      # Get y coordinate
      ycoord = as.numeric(sub(pattern = ".*?RG DNA Tray\\s\\d{3}\\sslot\\s(.*?)-.", 
                   replacement = "\\1", 
                   x = location_stripped,
                   ignore.case = TRUE)),
      
      # Get x coordinate
      xcoord = sub(pattern = ".*?RG DNA Tray\\s\\d{3}\\sslot\\s\\d{1,2}-(.*?)", 
                   replacement = "\\1", 
                   x = location_stripped,
                   ignore.case = TRUE))
  
  return(new_df)
}

##############################
# Worksheet functions
##############################

read_rfc1_ws <- function(worksheet_number) {
  
  output <- readxl::read_excel(paste0("W:/MolecularGenetics/Neurogenetics/Research/Joe Shaw Translational Post 2022/RFC1 worksheets/", 
                                      worksheet_number, "/", worksheet_number, ".xlsx"),
                               sheet = "results_sheet",
                               skip = 2) %>%
    janitor::clean_names() %>%
    mutate(worksheet = worksheet_number) %>%
    select(worksheet, sample, dna_no, first_name, surname, sample_file, marker, allele_1, size_1, height_1, 
           allele_2, size_2, height_2, result, report_type, coded_result)
}


get_ws_sample_info <- function(worksheet_number) {
  
  output <- readxl::read_excel(paste0("W:/MolecularGenetics/Neurogenetics/Research/Joe Shaw Translational Post 2022/RFC1 worksheets/", 
                          worksheet_number, "/", worksheet_number, ".xlsx"),
                   sheet = "data_sheet",
                   skip = 2,
                   n_max = 33)%>%
  janitor::clean_names() %>%
  select(episode, first_name, surname, location, dna_ng_ul) %>%
  filter(!episode %in% c("Water1", "Water2", "Water3", "Water4"))

  return(output)
  
}

##############################
# Pullsheet functions
##############################

read_pullsheet <- function(pullsheet) {
  
  pull_sheet_filepath <- "W:/MolecularGenetics/Neurogenetics/Research/DNA_aliquots_for_research/Pull sheets/"
  
  output <- read_excel(paste0(pull_sheet_filepath, pullsheet),
                       skip = 2) %>%
    janitor::clean_names() %>%
    mutate(sheet = pullsheet) %>%
    select(sheet, original_sample_id, amount_taken_ul) %>%
    filter(!is.na(original_sample_id))
  
  return(output)
}

##############################