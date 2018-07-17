# This script will read in data with corresponding labels in last column and analyze the cluster groups using the patient data
library(SNFtool)
library(tidyverse)
library(tibble)
library(readr)

# read in results
dat_normal <- readRDS('../data/normal_distance_results.rda')
dat_daisy <- readRDS('../data/daisy_distance_results.rda')


normal_results <- summarise_cluster_groups(dat_normal)
daisy_results <- summarise_cluster_groups(dat_normal)

