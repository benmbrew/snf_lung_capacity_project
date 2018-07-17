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

#remove subject type and move subject number (ID) to rownames, then remove
rownames(temp_dat) <- temp_dat$subjectnumber
temp_dat$subjectnumber <- temp_dat$subjecttype <- NULL

# implement funciton that combines hospital and er visits into one column
# then remove the variables used

# for 3m 
temp_dat$any_hos_3m <- combine_with_sum(temp_dat$numer_3m, temp_dat$numhos_3m)
temp_dat$numer_3m <- temp_dat$numhos_3m <- NULL

# for 6m 
temp_dat$any_hos_6m <- combine_with_sum(temp_dat$numer_6m, temp_dat$numhos_6m)
temp_dat$numer_6m <- temp_dat$numhos_6m <- NULL

# for 12m 
temp_dat$any_hos_12m <- combine_with_sum(temp_dat$numer_12m, temp_dat$numhos_12m)
temp_dat$numer_12m <- temp_dat$numhos_12m <- NULL

# for 18m 
temp_dat$any_hos_18m <- combine_with_sum(temp_dat$numer_18m, temp_dat$numhos_18m)
temp_dat$numer_18m <- temp_dat$numhos_18m <- NULL

# for 24m 
temp_dat$any_hos_24m <- combine_with_sum(temp_dat$numer_24m, temp_dat$numhos_24m)
temp_dat$numer_24m <- temp_dat$numhos_24m <- NULL

# for 30m 
temp_dat$any_hos_30m <- combine_with_sum(temp_dat$numer_30m, temp_dat$numhos_30m)
temp_dat$numer_30m <- temp_dat$numhos_30m <- NULL

# for 36m 
temp_dat$any_hos_36m <- combine_with_sum(temp_dat$numer_36m, temp_dat$numhos_36m)
temp_dat$numer_36m <- temp_dat$numhos_36m <- NULL


# impute on time series trajectories using mean or splines
temp_data <- temp_dat
# column trajectories: num_
impute_trajectory <- function(temp_data){
  
}

# create group vector to merge to column names
group_vector <- c('group_1_binary', 'group_1', 'group_1', 'group_1', 'group_1', 'group_1', 'group_1', 
                  'group_1', 'group_1', 'group_1', 'group_1', 'group_1', 'group_1', 'group_1_binary', 'group_1_binary', 
                  'group_2', 'group_2', 'group_2', 'group_2', 'group_2', 'group_2_binary', 'group_2', 'group_2', 
                  'group_2_factor', 'group_2_factor', 
                  'group_3', 'group_3', 'group_3', 'group_3', 'group_3', 'group_3', 'group_3', 
                  
                  'group_4_binary', 'group_4_binary', 'group_4_binary', 'group_4_binary', 'group_4_binary', 
                  'group_4_binary', 'group_4_binary', 'group_4_binary', 
                  'group_5', 'group_5', 'group_5', 'group_5', 'group_5', 'group_5', 'group_5', 
                  'group_6_binary','group_6_binary', 'group_6_binary', 'group_6_binary', 'group_6_binary',  
                  'group_6_binary', 'group_6_binary', 'group_6_binary', 'group_6_binary', 'group_6_binary',  
                  'group_6_binary', 'group_6_binary', 'group_7_binary', 'group_7_binary', 'group_7',  
                  'group_7', 'group_7', 'group_7', 'group_7', 'group_7','group_7',  'group_7', 'group_7', 
                  'group_7', 'group_7', 'group_7')

# add groupvector to names
names(temp_dat) <- paste0(names(temp_dat), '_', group_vector)

# determine the amount of missingness for each column
temp_cols <-as.data.frame(unlist(apply(temp_dat, 2, function(x) {
  round(length(which(is.na(x)))/nrow(temp_dat), 2)
})))

temp_cols$variable <- rownames(temp_cols)
colnames(temp_cols)[1] <- 'percent_missing'
rownames(temp_cols) <- NULL

# save data and manually enter data groups
write_csv(temp_cols, '../data/variable_percent_missing.csv')

# remove columns with over 30% missing
temp_30 <- temp_cols$variable[temp_cols$percent_missing < .41]
# removes 46 variables
temp_full <- temp_dat[, temp_30]
# this removes 571 sampeles
temp_full <- temp_full[complete.cases(temp_full),]

# create a data list with each group consituting an element in the list
data_list <- list()
unique_groups <- c('group_1', 'group_2', 'group_3', 'group_4', 'group_5', 'group_6')
i = 1

for(i in 1:length(unique_groups)){
  temp_group <- unique_groups[i]
  temp_sub <- temp_full[, grepl(temp_group, names(temp_full))]
  temp_ids <- rownames(temp_full)
  temp_sub <- temp_sub %>% mutate_if(sapply(as.data.frame(temp_sub), is.character), as.factor)
  rownames(temp_sub) <- temp_ids
  data_list[[i]] <- temp_sub
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




# plot the distribution of NA
# ggplot(temp_cols, 
#        aes(reorder(variable, 
#                    -percent_missing), 
#            percent_missing)) +
#   geom_bar(stat = 'identity', 
#            fill = 'red', 
#            alpha = 0.4) +
#   labs(x = 'Variables', 
#        y = 'Percent missing') +
#   theme_gray(base_size = 10) +
#   theme(axis.text.x = element_text(angle=90)) 
#   
# 
# hist(temp_cols$percent_missing, 
#      xlab = 'Percent NA', 
#      main = 'Distribution of NAs', 
#      col = 'lightblue',
#      breaks = 15)
# 
# # start by removing variables with over 50%  NA
# temp_complete <- temp_dat[, -which(colMeans(is.na(temp_dat)) > 0.5)]
# 
# temp_full <- temp_complete[complete.cases(temp_complete), ]

###### 
# implement SNF
######
# run SNF using normal distance function
final_data <- SNFClustering(data = data_list, 
                            numClus = 4, 
                            sampleRows = TRUE, 
                            dist_type = 'normal')

# run SNF using daisy (interval scaled) distance function
final_data_daisy <- SNFClustering(data = data_list, 
                                  numClus = 4, 
                                  sampleRows = TRUE, 
                                  dist_type = 'daisy')


# save data for analysis
saveRDS(final_data, '../data/normal_distance_results.rda')
saveRDS(final_data_daisy, '../data/daisy_distance_results.rda')
