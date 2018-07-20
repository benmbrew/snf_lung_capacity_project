# This script will read in data with corresponding labels in last column and analyze the cluster groups using the patient data

# NOTES
###########
# ideas for analyzing results 
# use the summriase_cluster_groups function (needs editing), linear regression 
# regressing each variable on groups to find significance
# MANOVA can take multiple y - 
# https://www.r-bloggers.com/multiple-analysis-of-variance-manova/
# kruskal.test() non parametic anova (numeric)
# chi square and fisher exact test (categorical)
# longitudinal data - variable = random effect is individual and time, fixed effect is group.
# all variables for each group, kruskal.test() is mean fev higher in one group compared to another
###########
library(SNFtool)
library(tidyverse)
library(tibble)
library(readr)
library(kml3d)
library(longitudinalData)
library(broom)
library(reshape2)
library(ggplot2)
library(ggthemes)

# source functions
source('functions.R')

# read in results
dat_normal <- readRDS('../data/normal_distance_results.rda')
dat_daisy <- readRDS('../data/daisy_distance_results.rda')

# make any_adclinic_g6 variable numeric
dat_normal$any_adclinic_g6 <- as.numeric(as.character(dat_normal$any_adclinic_g6))
dat_daisy$any_adclinic_g6 <- as.numeric(as.character(dat_daisy$any_adclinic_g6))

# convert all integers into factors
dat_normal <- dat_normal %>% mutate_if(sapply(as.data.frame(dat_normal), is.integer), as.factor)
dat_daisy <- dat_daisy %>% mutate_if(sapply(as.data.frame(dat_daisy), is.integer), as.factor)

# regress each variable on our groups and find signficant relationships
lm_normal <- get_group_stats(dat_normal, method = 'lm')
lm_daisy <- get_group_stats(dat_daisy, method = 'lm')

# use kruskal test, non parametri anova, on numeric variables 
kruskal_normal <- get_group_stats(dat_normal, method = 'kruskal')
kruskal_daisy <- get_group_stats(dat_daisy, method = 'kruskal')

# use chi squared of fisher exact test on categorical variables
chi_normal <- get_group_stats(dat_normal, method = 'chi_squared')
chi_daisy <- get_group_stats(dat_daisy, method = 'chi_squared')

# find significant variables 
lm_normal$sig <- ifelse(lm_normal$p.value <=0.05, TRUE, FALSE)
lm_daisy$sig <- ifelse(lm_daisy$p.value <=0.05, TRUE, FALSE)

kruskal_normal$sig <- ifelse(kruskal_normal$p.value <=0.05, TRUE, FALSE)
kruskal_daisy$sig <- ifelse(kruskal_daisy$p.value <=0.05, TRUE, FALSE)

chi_normal$sig <- ifelse(chi_normal$p.value <=0.05, TRUE, FALSE)
chi_daisy$sig <- ifelse(chi_daisy$p.value <=0.05, TRUE, FALSE)

# subset by only significant variables 
lm_sig_norm <- lm_normal$var_name[lm_normal$sig == TRUE]
lm_sig_norm <- lm_sig_norm[!duplicated(lm_sig_norm)]
lm_sig_daisy <- lm_daisy$var_name[lm_daisy$sig == TRUE]
lm_sig_daisy <- lm_sig_daisy[!duplicated(lm_sig_daisy)]

kruskal_sig_norm <- kruskal_normal$var_name[kruskal_normal$sig == TRUE]
kruskal_sig_norm <- kruskal_sig_norm[!duplicated(kruskal_sig_norm)]
kruskal_sig_daisy <- kruskal_daisy$var_name[kruskal_daisy$sig == TRUE]
kruskal_sig_daisy <- kruskal_sig_daisy[!duplicated(kruskal_sig_daisy)]


chi_sig_norm <- chi_normal$var_name[chi_normal$sig == TRUE]
chi_sig_norm <- chi_sig_norm[!duplicated(chi_sig_norm)]
chi_sig_daisy <- chi_daisy$var_name[chi_daisy$sig == TRUE]
chi_sig_daisy <- chi_sig_daisy[!duplicated(chi_sig_daisy)]

# use the summarise_cluster_groups to get mean difference and plots
temp_normal <- summarise_cluster_groups(dat_normal)
temp_daisy <- summarise_cluster_groups(dat_daisy)

# make box plot for each type of analysis: lm, kruskal, chi

# plot_vars <- c('ebf_duration_g4', "lbf_3m_status_g4", "bronchodilator_2_3y_g6", "wheeze3m_g6", "percent_lbf_3m_partial" )
temp_dat <- dat_normal
# subset dat aby sig_Vars
temp_dat <- temp_dat[, c("percent_lbf_3m_exclusive", "percent_lbf_3m_exclusive_after_hospital", "percent_lbf_3m_partial",
                         "percent_lbf_3m_zero", "ebf_duration_g4", "percent_bronchodilator_2_3y_g6", "percent_wheeze3m_g6",
                         "fev_zscore_18m_g3","fev_zscore_1y_g3","fev_zscore_3m_g3","fev_zscore_3y_g3", 'clustering_label')]

temp_male <- temp_dat[, grepl('male|clustering_label', names(temp_dat))]
temp_lbf_3m <- temp_dat[, grepl('lbf_3m|clustering_label', names(temp_dat))]
temp_broncho <- temp_dat[, grepl('broncho|clustering_label', names(temp_dat))]
temp_wheeze_3m <- temp_dat[, grepl('wheeze|clustering_label', names(temp_dat))]
fev_18m <- temp_dat[, grepl('fev_zscore_18m|clustering_label', names(temp_dat))]
fev_1y <- temp_dat[, grepl('fev_zscore_1y|clustering_label', names(temp_dat))]
fev_3m <- temp_dat[, grepl('fev_zscore_3m|clustering_label', names(temp_dat))]
fev_3y <- temp_dat[, grepl('fev_zscore_3y|clustering_label', names(temp_dat))]


# boxplot on numeric data - using full data to get points and true avg, median
ggplot(dat_daisy, aes(clustering_label, fev_zscore_3y_g3)) + 
  geom_boxplot(col = '#5182D9', size = 1) + geom_jitter(size = 2, col = 'white', width = 0.1, alpha = 0.5) +
  theme_dark(base_size = 16) +  
  labs(x = 'Clustering group', 
       y = 'Months', 
       title = '')


# barplot for group variable
temp_melt <- melt(temp_lbf_3m, id.vars = 'clustering_label')

ggplot(temp_melt, aes(clustering_label, value, group = variable, fill = variable)) + 
  geom_bar(stat = 'identity',position = 'dodge',  alpha = 0.7) +
  theme_solarized_2(base_size = 16) + 
  scale_fill_manual(values=c("#999999", "#E69F00", "#56B4E9", 'black'), 
                    name="LBF at 3m",
                    breaks=c("percent_lbf_3m_exclusive", "percent_lbf_3m_exclusive_after_hospital", "percent_lbf_3m_partial", "percent_lbf_3m_zero"),
                    labels=c("% exclusive", "% exclusive after hospital", "% partial", "% zero")) +
  labs(x = 'Clustering group', 
       y = 'Percent',
       title = '', 
       subtitle = '')

# barplot for singe variable
ggplot(temp_wheeze_3m, aes(clustering_label, percent_wheeze3m_g6)) + 
  geom_bar(stat = 'identity',fill = 'blue', col = 'black', alpha = 0.4) +
  theme_solarized_2(base_size = 16) + 
  labs(x = 'Clustering group', 
       y = 'Percent',
       title = '',
       subtitle = '')

