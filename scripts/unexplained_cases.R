################################################################################
## Unexplained Cases
## Selecting samples for optical genome mapping with Bionano technology
################################################################################

##############################
# Load libraries
##############################

library(tidyverse)

setwd(dir = "W:/MolecularGenetics/Neurogenetics/Research/Joe Shaw Translational Post 2022/RFC1 R code/RFC1_analysis")

##############################
# Read data
##############################

goshg2p_cases <- read_csv("data/goshg2p_export_notepad_20220614.csv")

##############################
# Cases to exclude
##############################

# Diagnosis confirmed/consistent
diag_confirmed <- c("20RG-331G0039", "19RG-350G0090", "20RG-014G0059", "20RG-227G0100", 
                    "20RG-344G0069", "21RG-124G0137", "20RG-310G0020", "21RG-035G0046", 
                    "21RG-183G0164", "21RG-057G0015", "21RG-057G0059", "21RG-077G0260", 
                    "21RG-098G0002", "21RG-102G0018", "21RG-266G0115", "21RG-275G0112", 
                    "22RG-054G0062", "22RG-056G0023")

# Class 5 and class 3 variants detected in same gene
snvs_same_gene <- c("20RG-328G0168", "21RG-104G0080", "21RG-118G0161", "21RG-104G0080", 
                    "20RG-289G0066",
                    # 2 class 3 variants in same gene
                    "22RG-026G0054")

# Dominant inheritance
dom_inheritance <- c("21RG-117G0050", "19RG-357G0098")

##############################
# Potential cases
##############################

filtered_cases <- goshg2p_cases %>%
  filter(Inheritance %in% c("AR", "AD/AR", "NULL") &
           acmg_classification %in% c("3 hot", "4", "5") &
           alt_freq < 0.8) %>%
  arrange(Sample_name) %>%
  # Multiple pipeline runs lead to duplicated rows, so remove with distinct
  dplyr::distinct(.keep_all = TRUE) %>%
  # Only "NULL" values in these columns
  select(-c(gene_constraint, inheritance_pattern_id, assigned_reviewer_id, 
            aaf_exac_all, clinvar_sig, impact)) %>%
  # Not interested in these columns
  select(-c(date_reported, user_reported_id, idTranscripts_and_inheritance,
            acmg_assigned, acmg_checked, acmg_date_checked, results_id, assigned, checked, date_checked,
            gt_types, gene_roi, Gene, Transcript, acmg_variant_id,
            hgvs_nomenclature, protein_change, classification,
            chrom, start, end, gosh_id, ref, alt, total_alt,
            PanelID, Pipeline_Run_ID)) %>%
  filter(!Sample_name %in% c(diag_confirmed, snvs_same_gene, dom_inheritance)) %>%
  filter(reported == 1)

colnames(filtered_cases)

# How many patients
length(unique(filtered_cases$Sample_name))

candidate_genes <- (filtered_cases$gene)

write.csv(candidate_genes, "candidate_genes.csv", row.names = FALSE)

paste0(unique(filtered_cases$gene), collapse = " OR ")

length(unique(filtered_cases$gene))

##############################
# Potential cases
##############################

# Cases with a heterozygous hot 3, 4 or 5 in a recessive gene, with no other genetic explanation of phenotype.

"20RG-016G0016"
"20RG-245G0075"
"20RG-262G0026"
"20RG-322G0082"
"20RG-325G0004"
"20RG-332G0044"
"20RG-351G0085"
"21RG-049G0043"
"21RG-063G0217"
"21RG-063G0219"
"21RG-098G0005"
"21RG-125G0227"
"21RG-159G0125"
"21RG-176G0047"
"21RG-182G0188"
"21RG-183G0172"
"21RG-189G0130"
"21RG-214G0059"
"21RG-305G0092"
"21RG-308G0021"
"21RG-348G0152"
"22RG-038G0053"
"22RG-054G0090"

##############################
# Cases where I'm unsure
##############################

"21RG-057G0018"
# Variant not on report - may not fit clinical phenotype

"21RG-172G0023"
"21RG-231G0073"
"21RG-258G0022"
"21RG-277G0023"
"22RG-044G0092"

##############################