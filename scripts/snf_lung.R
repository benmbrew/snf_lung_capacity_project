###### this script will run snf over lung function data
library(SNFtool)
library(tidyverse)
library(tibble)
library(readr)
library(cluster)
library(analogue)

# source functions
source('functions.R')

# read in data 
temp_dat <- read.csv('../data/lung_data.csv', 
                     na.strings = c("777", NA), 
                     stringsAsFactors = FALSE)

# determine which variables we can throw out - keep only the Z-score FEVs
# get the variables that have a z score transformation 
cols_temp <- colnames(temp_dat)[grepl('LCI|FEV|FVC', colnames(temp_dat))]

# only keep zscores
cols_remove <- cols_temp[!grepl('zscore', cols_temp)]

# remove those variables from the data set
temp_dat <- temp_dat[, !colnames(temp_dat) %in% cols_remove]

# remove StudyCenter column
temp_dat$StudyCenter <- NULL

# make column names lowerase 
colnames(temp_dat) <- tolower(colnames(temp_dat))

# recode factor variables to deal with blank space and NA
temp_dat$lbf_3m_status <- ifelse(temp_dat$lbf_3m_status == '', NA,
                                 ifelse(temp_dat$lbf_3m_status == 'unknown', NA, temp_dat$lbf_3m_status))

temp_dat$lbf_6m_status <- ifelse(temp_dat$lbf_6m_status == '', NA,
                                 ifelse(temp_dat$lbf_6m_status == 'unknown', NA, temp_dat$lbf_6m_status))


# implement funciton that combines hospital and er visits into one column
# then remove the variables used

# combine lbf_3m and lbf_6m
temp_dat$lbf_6m_status <-NULL
temp_dat$sf6 <- temp_dat$n2 <- NULL

# combine adclinic1y and adclinic3y into one variable 
temp_dat$any_adclinic <- combine_with_sum(temp_dat$adclinic1y, temp_dat$adclinic3y)
temp_dat$adclinic1y <- temp_dat$adclinic3y <- NULL

# use on recurrentwheeze variables 
temp_dat$any_recurrent_wheeze <- as.numeric(combine_with_sum(temp_dat$recwheeze1y, temp_dat$recwheeze3y))
temp_dat$recwheeze1y <- temp_dat$recwheeze3y <- NULL

# for 3m 
temp_dat$any_hos_3m <- as.numeric(combine_with_sum(temp_dat$numer_3m, temp_dat$numhos_3m))
temp_dat$numer_3m <- temp_dat$numhos_3m <- NULL

# for 6m 
temp_dat$any_hos_6m <- as.numeric(combine_with_sum(temp_dat$numer_6m, temp_dat$numhos_6m))
temp_dat$numer_6m <- temp_dat$numhos_6m <- NULL

# for 12m 
temp_dat$any_hos_12m <- as.numeric(combine_with_sum(temp_dat$numer_12m, temp_dat$numhos_12m))
temp_dat$numer_12m <- temp_dat$numhos_12m <- NULL

# for 18m 
temp_dat$any_hos_18m <- as.numeric(combine_with_sum(temp_dat$numer_18m, temp_dat$numhos_18m))
temp_dat$numer_18m <- temp_dat$numhos_18m <- NULL

# for 24m 
temp_dat$any_hos_24m <- as.numeric(combine_with_sum(temp_dat$numer_24m, temp_dat$numhos_24m))
temp_dat$numer_24m <- temp_dat$numhos_24m <- NULL

# for 30m 
temp_dat$any_hos_30m <- as.numeric(combine_with_sum(temp_dat$numer_30m, temp_dat$numhos_30m))
temp_dat$numer_30m <- temp_dat$numhos_30m <- NULL

# for 36m 
temp_dat$any_hos_36m <- as.numeric(combine_with_sum(temp_dat$numer_36m, temp_dat$numhos_36m))
temp_dat$numer_36m <- temp_dat$numhos_36m <- NULL

#remove subject type and move subject number (ID) to rownames, then remove
rownames(temp_dat) <- temp_dat$subjectnumber
temp_dat$subjectnumber <- temp_dat$subjecttype <- NULL

# remove maternal variables with too much NAs
temp_dat$bmi_mom <- temp_dat$egwg_mom <- temp_dat$gwg_mom <- NULL

# remove cateogrical variables with too much missingness
temp_dat$smk_inside_6m <- temp_dat$smk_inside_12m <- temp_dat$smk_inside_18m <- 
  temp_dat$smk_inside_24m <- temp_dat$smk_inside_30m <- NULL

temp_dat$wheeze6m <- temp_dat$wheeze1y <- temp_dat$wheeze18m <- 
  temp_dat$wheeze2hy <- temp_dat$wheeze2y <- NULL

temp_dat$any_hos_6m <- temp_dat$any_hos_12m <-  temp_dat$any_hos_18m <- temp_dat$any_hos_24m <-
  temp_dat$any_hos_30m <-NULL

# remove ocs_2_3y and ics_2_3y
temp_dat$ocs_2_3y <- temp_dat$ics_2_3y <- NULL

# impute on time series trajectories using mean or splines
# column trajectories: 
# 1) age_3m, age_1y, age_18m, age_3y, weight_3m, weight_1y, weight_18m, weight_3y, height_3m, height_1y, height_18m, height_3y - ok
# 2) smk_inside_3m, smk_inside_6m, smk_inside_12m, smk_inside_18m, smk_inside_24m, smk_inside_30m, smk_inside_36m - bad
# 3)'no2_t1', 'no2_t2', 'no2_t3', 'no2_t4', 'no2_t5', 'no2_t6', 'no2_t7' - ok 
# 4) wheeze3m, wheeze6m, wheeze1y, wheeze18m, wheeze2y, wheeze2hy, wheeze3y - categorical
# 5) recwheeze3y, recwheeze1y - combined into one
# 6) any_hos_3m, any_hos_6m, any_hos_12m, any_hos_18m, any_hos_24m, any_hos_30m, any_hos_36m - categorical
# 7) lci_zscore_3m, lci_zscore_1y, lci_zscore_18m, lci_zscore_3y ok 679
# 8) fev_zscore_3m, fev_zscore_1y, fev_zscore_18m, fev_zscore_3y
# 9) fevfvc_zscore_3m, fevfvc_zscore_1y, fevfvc_zscore_18m, fevfvc_zscore_3y


# loop through trajectories and impute
# first get NA index for age, weight, and height, so we can impute those back to NA
# age_all_na <-temp_dat[, grepl('age', temp_dat)]
# which(rowSums(is.na(temp_dat))==ncol(temp_dat))

imputation_list <- list()
unique_traj <- c('age', 'height', 'weight', 'no2', 'lci_zscore', 'fev_zscore', 'fevfvc_zscore')
for(traj in 1:length(unique_traj)){
  this_traj <- unique_traj[traj]
  traj_vector <- names(temp_dat)[grepl(this_traj, names(temp_dat))]
  imputation_list[[traj]] <- impute_trajectory(traj_vector)
}

# get data frame 
temp_imputed <- do.call('cbind', imputation_list)

# joint temp_imputed with with temp_dat
temp_imputed$ids <- rownames(temp_imputed)
temp_dat$ids <- rownames(temp_dat)
temp_dat <- inner_join(temp_imputed, temp_dat, by = 'ids')
rm(temp_imputed)
rownames(temp_dat) <- temp_dat$ids 
temp_dat$ids <- NULL

# remove columns .y and then remove the .x from the remaining columns (.x is our new imputed columns)
temp_dat <- temp_dat[, !grepl('.y', names(temp_dat), fixed = TRUE)]
names(temp_dat) <- gsub('.x', '', names(temp_dat), fixed = TRUE)

# keep only fev_zcore
temp_dat <- temp_dat[,!grepl('lci|fvc', names(temp_dat))]

# save data and manually enter data groups
# write_csv(temp_cols, '../data/variable_percent_missing_imputed.csv')

# read in variable - group dictionary to map each variable to its data group
var_group_dict <- read.csv('../data/variable_group_dictionary.csv', stringsAsFactors = FALSE)

# save data for exploration 
# write.csv(temp_dat, '../data/data_exploration.csv')

# sort each by variable name 
var_group_dict <- var_group_dict[order(var_group_dict$variable),]
temp_dat <- temp_dat[, order(names(temp_dat))]

# paste variable and groupvector together and then overlay the temp_dat variable names 
new_var_names <- paste0(var_group_dict$variable, '_', var_group_dict$group)
names(temp_dat) <- new_var_names

# determine the amount of missingness for each column
temp_cols <-as.data.frame(unlist(apply(temp_dat, 2, function(x) {
  round(length(which(is.na(x)))/nrow(temp_dat), 2)
})))

temp_cols$variable <- rownames(temp_cols)
colnames(temp_cols)[1] <- 'percent_missing'
rownames(temp_cols) <- NULL

# this removes 571 sampeles
temp_full <- temp_dat[complete.cases(temp_dat),]

# create a data list with each group consituting an element in the list
data_list <- list()
unique_groups <- c('g1', 'g2', 'g3', 'g4', 'g5', 'g6', 'g7')

for(i in 1:length(unique_groups)){
  temp_group <- unique_groups[i]
  temp_sub <- temp_full[, grepl(temp_group, names(temp_full))]
  temp_ids <- rownames(temp_full)
  temp_sub <- temp_sub %>% mutate_if(sapply(as.data.frame(temp_sub), is.character), as.factor)
  rownames(temp_sub) <- temp_ids
  data_list[[i]] <- temp_sub
}


# make heatmap for each data type 
i = 1
for(i in 1:length(data_list)){
  sub_dat <- data_list[[i]]
  sub_dat <- sub_dat[,!grepl('childgen', names(sub_dat))]
  sub_dat <- as.matrix(dist(sub_dat))
  sub_aff <- affinityMatrix(sub_dat)
  displayClustersWithHeatmap(as.matrix(sub_aff))

}

# # look for time series variables and impute only if an NA is sandwhiched between two time 
# # series variables with same value
# # possible variables: (adclinic1y, adclinicc3y, lbf_3m_status, lbf_6m_status, smk_inside_3m, smk_inside_6m,
# # 'smk_inside_12m, smk_inside_18m, smk_inside_24m, smk_inside_30m, smk_inside_36m, recwheeze3y, recwheeze1y, wheeze1y, wheeze3y, wheeze2hy)
# temp <- as.data.frame(cbind(n2 = temp_dat$n2 , sf6 = temp_dat$sf6, three = temp_dat$smk_inside_3m, 
#                             six = temp_dat$smk_inside_6m, twelve=temp_dat$smk_inside_12m,
#                             eighteen = temp_dat$smk_inside_18, twenty_four = temp_dat$smk_inside_24m, 
#                             thirty = temp_dat$smk_inside_30m, thirty_six =temp_dat$smk_inside_36,
#                             rec_wheeze1y = temp_dat$recwheeze1y_group_6_binary, rec_wheeze3y = temp_dat$recwheeze3y_group_6_binary,
#                             wheeze_1y = temp_dat$wheeze1y_group_6_binary, wheeze_2hy = temp_dat$wheeze2hy_group_6_binary,
#                             wheez_3y = temp_dat$wheeze3y_group_6_binary))





# 
# # start by removing variables with over 50%  NA
# temp_complete <- temp_dat[, -which(colMeans(is.na(temp_dat)) > 0.5)]
# 
# temp_full <- temp_complete[complete.cases(temp_complete), ]

###### 
# implement SNF
######
# run SNF using normal distance function
final_data <- snf_clustering(data = data_list, 
                            numClus = 4, 
                            sampleRows = TRUE, 
                            dist_type = 'normal')

# run SNF using daisy (interval scaled) distance function
final_data_daisy <- snf_clustering(data = data_list, 
                                  numClus = 4, 
                                  sampleRows = TRUE, 
                                  dist_type = 'daisy')


# save data for analysis
saveRDS(final_data, '../data/normal_distance_results.rda')
saveRDS(final_data_daisy, '../data/daisy_distance_results.rda')
