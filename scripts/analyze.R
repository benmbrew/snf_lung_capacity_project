# This script will read in data with corresponding labels in last column and analyze the cluster groups using the patient data
library(SNFtool)
library(tidyverse)
library(tibble)
library(readr)

# read in results
dat_normal <- readRDS('../data/normal_distance_results.rda')
dat_daisy <- readRDS('../data/daisy_distance_results.rda')

# group by cluster label and summarise mean results, do same thing again for factors 
summarise_cluster_groups <- function(temp_dat){
  
  # summarise numeric and factor variables and then join by clustering_label
  temp_num <- temp_dat %>% group_by(clustering_label) %>% summarise_all(funs(mean))
  temp_num$childgender_group_1_binary <- temp_num$lbf_3m_status_group_2_factor <- temp_num$lbf_6m_status_group_2_factor <- NULL
  temp_fac <- temp_dat %>% 
    group_by(clustering_label) %>% 
    summarise(counts = n(),
              percent_male = round((sum(childgender_group_1_binary == 'M')/counts)*100, 2),
              percent_lbf_3m_exclusive = round((sum(lbf_3m_status_group_2_factor == 'Exclusive')/counts)*100, 2),
              percentlbf_3m_exclusive_after_hospital = round((sum(lbf_3m_status_group_2_factor == 'Exclusive After Hospital')/counts)*100, 2),
              percentlbf_3m_partial = round((sum(lbf_3m_status_group_2_factor == 'Partial')/counts)*100, 2),
              percentlbf_3m_zero = round((sum(lbf_3m_status_group_2_factor == 'Zero')/counts)*100, 2),
              percentlbf_6m_exclusive = round((sum(lbf_6m_status_group_2_factor == 'Exclusive')/counts)*100, 2),
              percentlbf_6m_exclusive_after_hospital = round((sum(lbf_6m_status_group_2_factor == 'Exclusive After Hospital')/counts)*100, 2),
              percentlbf_6m_partial = round((sum(lbf_6m_status_group_2_factor == 'Partial')/counts)*100, 2),
              percentlbf_6m_zero = round((sum(lbf_6m_status_group_2_factor == 'Zero')/counts)*100, 2))
  # join by clustering label
  temp_final <- inner_join(temp_fac, temp_num, by = 'clustering_label')
    
  return(temp_final)
}

normal_results <- summarise_cluster_groups(dat_normal)
daisy_results <- summarise_cluster_groups(dat_normal)

