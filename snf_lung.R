###### this script will run snf over lung function data
library(SNFtool)
library(tidyverse)
library(tibble)
library(readr)
library(cluster)
library(analogue)

# read in data 
temp_dat <- read.csv('lung_data.csv', 
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

# create group vector to merge to column names
group_vector <- c('group_1_binary', 'group_1', 'group_1', 'group_1', 'group_1', 'group_1', 'group_1', 
                  'group_1', 'group_1', 'group_1', 'group_1', 'group_1', 'group_1', 'group_1_binary', 'group_1_binary', 
                  'group_2', 'group_2', 'group_2', 'group_2', 'group_2', 'group_2_binary', 'group_2', 'group_2', 
                  'group_2_factor', 'group_2_factor', 
                  'group_3', 'group_3', 'group_3', 'group_3', 'group_3', 'group_3', 'group_3', 'group_3', 
                  'group_3', 'group_3',  'group_3', 'group_3', 'group_3', 'group_3', 
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
write_csv(temp_cols, 'variable_percent_missing.csv')

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

# HERE- remember if the factor is only one level (like gender) then maybe you can just code it as 0 and 1 numeric - the only non binary facotrs are lbf_3m, lbf_6m
# HERE- anytime there is an ER visit, does that emply hospitalization? 
# HERE - how many layers can SNF handle
# HERE - look at recurrent wheeze and its relationship to normal wheeze
# HERE - Normalize across datatype? Gender, survey, etc. Test for variance.


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
data <- data_list
# create a list of data frames based on the types of data 
# maybe observable clinical data and numerical data
SNFClustering <- function(data, 
                          numClus = 4, 
                          sampleRows = TRUE,
                          dist_type) {
  # Transpose data for the dist function
  # The dist function takes distances between rows
  sampleRowsRequired <- TRUE
  transposeData <- sampleRows != sampleRowsRequired
  if (transposeData) {
    data <- lapply(data, t)
  }
  
  # Calculate the distance between samples
  if(dist_type == 'normal') {
    distances <- lapply(data, function(x) as.matrix(dist(x)))
  } 
  if (dist_type == 'daisy') {
    distances <- lapply(data, function(x) as.matrix(daisy(x, metric = 'gower')))
  } 
  
  # Convert the distances to affinities
  affinities <- lapply(distances, affinityMatrix)
  
  # Fuse the affinity matrices
  fusedMatrix <- SNF(affinities)
  
  # Cluster the fused matrix
  labels <- spectralClustering(fusedMatrix, numClus)
  
  # cbind data to create a final data frame with corresponding cluster labels
  final_data <- do.call('cbind', data)

  # add a column to data for the labels
  final_data$clustering_label <- labels
  
  return(final_data)
}

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


