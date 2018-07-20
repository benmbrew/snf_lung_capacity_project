# This script will get summary stats for the dataset

###### this script will run snf over lung function data
library(SNFtool)
library(tidyverse)
library(tibble)
library(readr)
library(cluster)
library(analogue)

# source functions
source('functions.R')

# explore original data
# read in data 
temp_dat <- read.csv('../data/lung_data.csv', 
                     na.strings = c("777", NA), 
                     stringsAsFactors = FALSE)


# keep only zscore
cols_temp <- colnames(temp_dat)[grepl('LCI|FEV|FVC', colnames(temp_dat))]

# only keep zscores
cols_remove <- cols_temp[!grepl('zscore', cols_temp)]

# remove those variables from the data set
temp_dat <- temp_dat[, !colnames(temp_dat) %in% cols_remove]

# determine the amount of missingness for each column
temp_cols <-as.data.frame(unlist(apply(temp_dat, 2, function(x) {
  round(length(which(is.na(x)))/nrow(temp_dat), 2)
})))

temp_cols$variable <- rownames(temp_cols)
names(temp_cols)[1] <- 'percent_missing'

# plot the distribution of NA
ggplot(temp_cols,
       aes(reorder(variable,
                   -percent_missing),
           percent_missing)) +
  geom_bar(stat = 'identity',
           fill = 'red',
           alpha = 0.4) +
  labs(x = 'Variables',
       y = 'Percent missing') +
  theme_gray(base_size = 14) +
  theme(axis.text.x = element_text(angle=90))


hist(temp_cols$percent_missing,
     xlab = 'Percent NA',
     main = 'Distribution of NAs',
     col = 'lightblue',
     breaks = 15)


# plot that shows sample size on y axis and number of variables removed (from most to least NAs)
# on the x axis
# FEVFCV_zscore_18m, FEV_zscore_18m, FEVFVC_zscore_3m, FEV_zscore_3m, LCI_zscore_18m, LCI_zscore_3m, FEVFVC_zscore_1y,
# sort data by percent missing 
first <- temp_dat[complete.cases(temp_dat),]
temp_cols <- temp_cols[order(temp_cols$percent_missing, decreasing = T),]
result_list <- list()
sub_var_list <- list()
vars <- temp_cols$variable
i = 1
for(i in 1:80) {
  sub_var_list[[i]] <- vars[i]
  sub_vars <- unlist(sub_var_list)
  temp <- temp_dat[, !names(temp_dat) %in% sub_vars]
  n_rows <- nrow(temp[complete.cases(temp),])
  result_list[[i]] <- n_rows
  print(i)
}

temp_final <- as.data.frame(do.call('rbind', result_list))
temp_final$var_removed <- rownames(temp_final)
names(temp_final) <- c('sample_size', 'var_removed')

temp_final <- as.data.frame(temp_final)
temp_final$sample_size <- as.numeric(temp_final$sample_size)
temp_final$var_removed <- as.numeric(temp_final$var_removed)
# plot vars remved against sample size
# plot the distribution of NA
ggplot(temp_final,
       aes(var_removed,
           sample_size)) +
  geom_bar(stat = 'identity',
           fill = 'blue',
           alpha = 0.4) +
  labs(x = 'Varibles removed',
       y = 'Sample size',
       title = 'Complete data',
       subtitle = 'Removing one variable at a time') +
  theme_gray(base_size = 18) 


write.csv(as.data.frame(temp_cols),'~/Desktop/missing.csv')

#------------------------------------------------------------------------------
# - explore curated data

# read in data 
temp_dat <- read.csv('../data/data_exploration.csv', stringsAsFactors = FALSE)
var_group_dict <- read.csv('../data/variable_group_dictionary.csv', stringsAsFactors = FALSE)

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

# remove columns with over 30% missing
temp_50 <- temp_cols$variable[temp_cols$percent_missing <= .50]

# removes 46 variables
temp_full <- temp_dat[, temp_50]

# this removes 571 sampeles
temp_full1 <- temp_full[complete.cases(temp_full),]
