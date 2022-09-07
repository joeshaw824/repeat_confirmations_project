################################################################################
## Frozen Blood Samples for Neurology Research
################################################################################

##############################
# Load libraries
##############################

library(tidyverse)
library(readxl)
library(janitor)
library(lubridate)

##############################

setwd("I:/Genetics/DNA Lab/Prep Lab/Research samples from NHNN")

bioresource_files <- c(list.files("Bioresources/2020"), list.files("Bioresources/2021"), list.files("Bioresources/2022"))

read_excel("Bioresources/2020/Bioresources List 1- August 2020_Collected by MG on 23.10.2020.xlsx",
           sheet = "Bioresource")

##############################