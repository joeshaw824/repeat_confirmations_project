################################################################################
## Unexplained Cases for Bionano Optical Genome Mapping
## Joseph.Shaw@gosh.nhs.uk / joseph.shaw3@nhs.net
################################################################################

# This script is for selecting unexplained ataxia cases for optical genome mapping 
# with Bionano technology. Inputs are a variant table exported from GOSHG2P and
# patient details exported from Epic.

##############################
# Load libraries
##############################

library(tidyverse)
library(readxl)
library(janitor)

setwd(dir = "W:/MolecularGenetics/Neurogenetics/Research/Joe Shaw Translational Post 2022/RFC1 R code/RFC1_analysis")

##############################
# Read data
##############################

goshg2p_cases <- read_csv("data/goshg2p_export_notepad_20220614.csv") %>%
  janitor::clean_names()

# "Ataxia panel referrals" [6303714] report on Epic. Results exported from Epic for these test codes:

# NGEN NGS - Adult onset hereditary spastic paraplegia
# NGEN_Hereditary neuropathy or pain disorder – NOT PMP22 copy number (Panel) (R78.1)
# NGEN_Hereditary ataxia with onset in adulthood (Panel) (R54.1)

case_info <- read_excel("data/Ataxia_panel_referrals_20220722_1606.xlsx") %>%
  janitor::clean_names() %>%
  dplyr::rename(sample_name = test_specimen_id)

##############################
# Initial filtering
##############################

filtered_cases <- goshg2p_cases %>%
  filter(inheritance %in% c("AR", "AD/AR", "NULL") &
           acmg_classification %in% c("4", "5") &
           # Remove cases with homozygous variants
           alt_freq < 0.8) %>%
  arrange(sample_name) %>%
  # Only "NULL" values in these columns
  select(-c(gene_constraint, inheritance_pattern_id, assigned_reviewer_id, 
            aaf_exac_all, clinvar_sig, impact)) %>%
  # Not interested in these columns
  select(-c(date_reported, user_reported_id, id_transcripts_and_inheritance,
            acmg_assigned, acmg_checked, acmg_date_checked, results_id, assigned, checked, date_checked,
            gt_types, gene_roi, gene_2, transcript, acmg_variant_id,
            hgvs_nomenclature, protein_change, classification,
            chrom, start, end, gosh_id, ref, alt, total_alt,
            panel_id, pipeline_run_id))

##############################
# Cases to exclude from manual inspection on Epic
##############################

# Diagnosis confirmed/consistent with condition
diag_confirmed <- c("20RG-331G0039", "19RG-350G0090", "20RG-014G0059", "20RG-227G0100", 
                    "20RG-344G0069", "21RG-124G0137", "20RG-310G0020", "21RG-035G0046", 
                    "21RG-183G0164", "21RG-057G0015", "21RG-057G0059", "21RG-077G0260", 
                    "21RG-098G0002", "21RG-102G0018", "21RG-266G0115", "21RG-275G0112", 
                    "22RG-054G0062", "22RG-056G0023")

# Class 5 and class 3 variants detected in same gene (unlikely to have structural variant as underlying cause)
snvs_same_gene <- c("20RG-328G0168", "21RG-104G0080", "21RG-118G0161", "21RG-104G0080", 
                    "20RG-289G0066",
                    # 2 class 3 variants in same gene
                    "22RG-026G0054")

# Class 3 with dominant inheritance noted on report
dom_inheritance <- c("21RG-117G0050", "19RG-357G0098")

##############################
# Further filtering and add case identifiers
##############################

potential_bionano_cases <- filtered_cases %>%
  filter(!sample_name %in% c(diag_confirmed, snvs_same_gene, dom_inheritance)) %>%
  # Remove duplicate calls of the same variant for the same sample due to multiple pipeline runs
  filter(!(base::duplicated(sample_name) & base::duplicated(acmg_cdna))) %>%
  # Add identifiers
  left_join(case_info, by = "sample_name") %>%
  filter(!is.na(patient_surname)) %>%
  select(sample_name, pt_first_nm, patient_surname, dob, nhs_number, mrn, storage_location,
       gene, acmg_cdna, acmg_protein_change, alt_freq, acmg_classification, inheritance)

# How many patients?
length(unique(potential_bionano_cases$sample_name))

write.csv(potential_bionano_cases, 
          paste0("outputs/potential_cases_for_bionano_", format(Sys.time(), "%Y%m%d"), ".csv"),
          row.names= FALSE)

##############################



