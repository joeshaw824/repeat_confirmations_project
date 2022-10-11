################################################################################
## Useful functions
################################################################################

##############################
# Epic functions
##############################

# Useful function for splitting DNA storage location string used in Epic
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

##############################