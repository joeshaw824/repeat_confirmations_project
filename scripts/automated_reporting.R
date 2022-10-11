################################################################################
## Automated reporting of RFC1 research results
## Joseph.Shaw@gosh.nhs.uk
################################################################################

##############################
# Load libraries
##############################

library(rmarkdown)
library(tidyverse)
library(readxl)
library(janitor)

setwd(dir = "W:/MolecularGenetics/Neurogenetics/Research/Joe Shaw Translational Post 2022/RFC1 R code/RFC1_analysis")

##############################
# User input
##############################

results_file <- "CANVAS_Screeninglist_update_Aug2022_GOSH"
positive_text <- "RFC1 CANVAS Spectrum Disorder confirmed"
negative_text <- "RFC1 CANVAS Spectrum Disorder NOT confirmed"
result_texts <- c(positive_text, negative_text)

##############################
# Load research results
##############################

# Additional formatting steps required for new results Excel.

negative_strings <- c("negative", "carrier")

positive_strings <-c("biallelic RFC1 expansion", "likely biallelic RFC1 expansion",
                     "likely biallelic RFC1 expansion (rechecked by Joe)",
                     "confirmed",
                     "already confirmed",
                     "patient already confirmed")

research_results <- read_excel(path = paste0("data/",results_file,".xlsx")) %>%
  janitor::clean_names() %>%
  mutate(
    name_string = toupper(paste0(first_name," ",surname)),
    result = case_when(
      interpretation %in% positive_strings ~positive_text,
      interpretation %in% negative_strings ~negative_text,
      TRUE ~"other")) %>%
  filter(!base::duplicated(name_string) &
           result %in% result_texts)

##############################
# Automated reporting function
##############################

generate_report <- function(name_input) {
  
  stopifnot(c("name_string", "consultant",
              "result") %in% colnames(research_results))
  
  # Isolate patient details
  patient_details <- research_results %>%
    dplyr::filter(name_string == name_input)

  stopifnot(nrow(patient_details) == 1)
  
  # Patient details for report
  patient_name  <- as.character(patient_details$name_string)
  
  clinician     <- as.character(patient_details$consultant)
  
  result        <- as.character(patient_details$result)
  
  file_name     <- paste0(patient_name, " RFC1 research report ", 
                          format(Sys.time(), '%d_%m_%Y'))

  stopifnot(result %in% c(positive_text, negative_text))

  # Automation Code
  rmarkdown::render(
    input         = "scripts/rfc1_report_template.Rmd",
    output_file   = file_name,
    output_dir    = "automated_reports/",
    output_format = "word_document",
    params        = list(
      patient_name = patient_name,
      clinician = clinician,
      result = result,
      show_code      = FALSE)
  )

}

##############################
# Generating reports
##############################

lapply(research_results$name_string, generate_report)

##############################
