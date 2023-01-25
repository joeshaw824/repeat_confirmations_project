################################################################################
## NOP56 PCR Validation
## joseph.shaw3@nhs.net
################################################################################

##############################
# Load libraries
##############################

library(tidyverse)
library(readxl)
library(janitor)
library(ggpubr)

setwd(dir = "W:/MolecularGenetics/Neurogenetics/Research/Joe Shaw Translational Post 2022/RFC1 R code/RFC1_analysis")

##############################
# Load data
##############################

read_nop56_allele_report <- function(report_file) {
  
  report <- read_excel(paste0("data/", report_file), skip = 9) %>%
    select(-c(`...1`)) %>%
    janitor::clean_names() %>%
    mutate(specimen_id = str_extract(sample, pattern = "..RG-...G...."),
           control_id = str_extract(sample, pattern = "C275"),
           allele_1_repeat = parse_number(allele_number_1),
           allele_2_repeat = parse_number(allele_number_2),
           worksheet = str_extract(report_file, "2.-...."))
}

allele_report_22_4357 <- read_nop56_allele_report("Diseases_Genes22-4357_FL_PCR_1_MergeProj_AlleleReport.xlsx")

allele_report_23_0086 <- read_nop56_allele_report("Diseases_Genes23-0086_FL1_MergeProj_AlleleReport.xlsx")

allele_report_22_4159 <- read_nop56_allele_report("Diseases_Genes22-1459_FL1_MergeProj_AlleleReport.xlsx")

allele_report_22_4170 <- read_nop56_allele_report("Diseases_Genes22-4170_FL1_MergeProj_AlleleReport.xlsx")

sequenced_ids <- c("22RG-273G0151", "22RG-273G0153", "22RG-273G0145", "22RG-283G0052",
                   "22RG-304G0056", "22RG-304G0056", "21RG-203G0002", "C275")

##############################
# PCR Amplicon Sizing
##############################
#######################
# Manipulate data
#######################

allele_repeats_longer <- allele_report_22_4357 %>%
  select(c(sample, specimen_id, marker, allele_1_repeat, allele_2_repeat)) %>%
  pivot_longer(cols = c(allele_1_repeat, allele_2_repeat),
             names_to = "category",
             values_to = "repeats") %>%
  filter(!is.na(marker)) %>%
  mutate(allele = parse_number(category),
         predicted_size_bp = case_when(
           marker == "NOP56_FL_PCR1" ~round((repeats*6)+113, 0),
           marker == "NOP56_flanking_PCR_2" ~round((repeats*6)+66,0),
           marker == "NOP56_RP" ~round((repeats*6)+49, 0)))

allele_sizes_longer <- allele_report_22_4357 %>%
  select(c(sample, marker, size_number_1, size_number_2)) %>%
  pivot_longer(cols = c(size_number_1, size_number_2),
               names_to = "category",
               values_to = "measured_size_bp") %>%
  mutate(allele = parse_number(category))

# Combine the two tables together
combined <- allele_repeats_longer %>%
  left_join(allele_sizes_longer, by = c("sample", "marker", "allele")) %>%
  filter(!is.na(repeats) & repeats < 40 & specimen_id %in% sequenced_ids)

#######################
# Plot results
#######################

x_label <- "ABI-3730xl amplicon (bp)"
y_label <- "Amplicon size (bp) from sequence"

fl1_plot <- ggplot(combined %>% 
                    filter(marker == "NOP56_FL_PCR1"), aes(x = measured_size_bp, y = predicted_size_bp)) +
  geom_point(size = 3, pch = 21) +
  geom_abline(linetype = "dashed") +
  theme_bw() +
  labs(x = "", y = y_label, title = "Flanking PCR 1") +
  ylim(140, 180) +
  xlim(140,180)

fl2_plot <- ggplot(combined %>% 
                     filter(marker == "NOP56_flanking_PCR_2"), 
                   aes(x = measured_size_bp, y = predicted_size_bp)) +
  geom_point(size = 3, pch = 21) +
  geom_abline(linetype = "dashed") +
  theme_bw() +
  labs(x = x_label, y = "", title = "Flanking PCR 2") +
  xlim(90, 140) +
  ylim(90, 140)

rp_plot <- ggplot(combined %>% 
                    filter(marker == "NOP56_RP"), aes(x = measured_size_bp, y = predicted_size_bp)) +
  geom_point(size = 3, pch = 21) +
  geom_abline(linetype = "dashed") +
  theme_bw() +
  labs(x = "", y = "", title = "Repeat primed PCR") +
  xlim(60, 110) +
  ylim(60, 110)

all_plots <- ggpubr::ggarrange(fl1_plot, fl2_plot, rp_plot,
                               ncol = 3, nrow = 1)

ggsave(plot = all_plots, 
       filename = "nop56_pcr_sizing_plot.png",
       path = "plots/",
       device = "png",
       dpi = 600,
       units = "cm",
       width = 18.5,
       height = 8)

##############################
# Normal allele distribution
##############################

affected_samples <- c("21RG-203G0002", "21RG-026G0059", "22RG-273G0143", "22RG-273G0135",
                      "22RG-273G0151", "22RG-273G0153", "22RG-273G0145", "22RG-273G0154",
                      "22RG-273G0156", "22RG-283G0052", "22RG-304G0056")

combined_reports <- rbind(allele_report_22_4357, allele_report_23_0086)

allele_sizes <- combined_reports %>%
  filter(marker == "NOP56_flanking_PCR_2" &
           !specimen_id %in% affected_samples &
           !is.na(specimen_id)) %>%
  select(sample, specimen_id, allele_1_repeat, allele_2_repeat) %>%
  pivot_longer(cols = -c(sample, specimen_id),
               names_to = "allele",
               values_to = "size") %>%
  filter(!is.na(size))

# Filter normal controls
num_controls <- length(unique(allele_sizes$specimen_id))

nop56_normal_alleles_plot <- ggplot(allele_sizes, aes(x = size, y = )) +
  geom_bar() +
  theme(axis.text.x = element_text(angle = 90)) +
  labs(x = "NOP56 hexanucleotide repeats", y = "Frequency",
       title = paste0("Size distribution for NOP56 repeats in ", num_controls, " normal controls"),
       subtitle = "Data from flanking PCR 2") +
  theme_bw() +
  scale_x_continuous(breaks = seq(15:15, 1)) +
  ylim(0, 30)

ggsave("nop56_normal_alleles_plot.png",
       nop56_normal_alleles_plot, 
       path = "plots/",
       units = "cm",
       width = 15,
       height = 10,
       dpi = 600)

##############################
# Uncertainty of Measurement
##############################

combined_data <- rbind(allele_report_22_4170, allele_report_22_4357, 
                       allele_report_22_4159,  allele_report_23_0086) %>%
  filter((specimen_id %in% c("21RG-203G0002", "22RG-304G0056") | control_id == "C275") 
         & marker %in% c("NOP56_FL_PCR1", "NOP56_flanking_PCR_2")) %>%
  mutate(identifier = case_when(
    is.na(specimen_id) ~control_id,
    is.na(control_id) ~specimen_id)) %>%
  select(sample, identifier, marker, worksheet, size_number_1, size_number_2) %>%
  arrange(identifier, marker, worksheet)

pcr_1_uom <- combined_data %>%
  filter(marker == "NOP56_FL_PCR1" & !sample %in% c(
    # Not 0.2uM primers
    "02_FL_1_21RG-203G0002_0.4_22-4170_B01_013.fsa", 
    "06_FL_1_21RG-203G0002_0.3_22-4170_F01_005.fsa",
    "01_FL_1_22RG-304G0056_0.4_22-4170_A01_015.fsa",
    "05_FL_1_22RG-304G0056_0.3_22-4170_E01_007.fsa",
    "03_FL_1_C275_0.4_22-4170_C01_011.fsa",
    "07_FL_1_C275_0.3_22-4170_G01_003.fsa",
    # Not 50ng input
    "02_FL1_PCR_21RG-203G0002_B01_22-4357_B01_013.fsa",
    "05_FL1_PCR_C275_E01_22-4357_E01_007.fsa"))


pcr_2_uom <- combined_data %>%
  filter(marker == "NOP56_flanking_PCR_2" & !sample %in% 
           c("04_22-4159_21RG-203G0002_FL2_D01_009.fsa",
            "18_FL_2_21RG-203G0002_0.4_22-4170_B03_029.fsa",
            "22_FL_2_21RG-203G0002_0.3_22-4170_F03_021.fsa",
            "02_FL2_PCR_21RG-203G0002_B05_22-4357_B05_045.fsa",
            "17_FL_2_22RG-304G0056_0.4_22-4170_A03_031.fsa",
            "21_FL_2_22RG-304G0056_0.3_22-4170_E03_023.fsa",
            "05_22-4159_C275_FL2_E01_007.fsa",
            "19_FL_2_C275_0.4_22-4170_C03_027.fsa",
            "23_FL_2_C275_0.3_22-4170_G03_019.fsa",
            "05_FL2_PCR_C275_E05_22-4357_E05_039.fsa"))

generate_uom_table <- function(input_table) {
  
  table_short <- input_table %>%
    select(-c(sample, marker, size_number_2)) %>%
    pivot_wider(id_cols = c(identifier),
                names_from = worksheet, 
                values_from = c(size_number_1)) %>%
    mutate(allele = "Short") %>%
    select(identifier, allele, `22-1459`, `22-4170`, `22-4357`, `23-0086`)
  
  table_long <- input_table %>%
    select(-c(sample, marker, size_number_1)) %>%
    pivot_wider(id_cols = c(identifier),
                names_from = worksheet, 
                values_from = c(size_number_2)) %>%
    mutate(allele = "Long") %>%
    select(identifier, allele, `22-1459`, `22-4170`, `22-4357`, `23-0086`)
  
  table_combined <- rbind(table_short, table_long) %>%
    arrange(identifier, desc(allele))
  
  return(table_combined)
  
}

pcr1_table <- generate_uom_table(pcr_1_uom)
pcr2_table <- generate_uom_table(pcr_2_uom)

write.csv(pcr1_table, "outputs/pcr1_table.csv", row.names = FALSE)
write.csv(pcr2_table, "outputs/pcr2_table.csv", row.names = FALSE)

##############################