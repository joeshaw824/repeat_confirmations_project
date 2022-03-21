################################################################################
## Automated reporting of RFC1 research results
################################################################################

##############################
# Load libraries
##############################

library(rmarkdown)

source("scripts/load_research_results.R")

##############################
# Automated reporting
##############################

generate_report <- function(name_input) {
  
  stopifnot(c("name_string", "identifiers", "consultant",
              "result_clean") %in% colnames(cleaned_results))
  
  # Isolate patient details
  patient_details <- cleaned_results %>%
    dplyr::filter(name_string == name_input)

  stopifnot(nrow(patient_details) == 1)
  
  # Patient details for report
  patient_name  <- as.character(patient_details$name_string)
  
  identifiers   <- as.character(patient_details$identifiers)
  
  clinician     <- as.character(patient_details$consultant)
  
  result        <- as.character(patient_details$result_clean)
  
  file_name     <- paste0(patient_name, "  RFC1 research report")

  stopifnot(result %in% c("biallelic AAGGG repeat expansion detected",
                          "biallelic AAGGG repeat expansion not detected"))

  # Automation Code
  rmarkdown::render(
    input         = "scripts/rfc1_report_template.Rmd",
    output_format = "word_document",
    output_file   = file_name,
    output_dir    = "automated_reports/",
    params        = list(
      patient_name = patient_name,
      identifiers = identifiers,
      clinician = clinician,
      result = result,
      show_code      = FALSE)
  )

}

##############################
# Testing the reporting
##############################

test_patients <- c("21RG-228G0027", "20RG-252G0140", "21RG-236G0106", "21RG-222G0087")

test_results <- cleaned_results %>%
  filter(internal_id %in% test_patients)

lapply(test_results$name_string, generate_report)

##############################
