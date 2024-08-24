## title: "Analyzing Elisa Assays"
## author: "Merai Dandouch"

# load packages
library('tidyverse')
library('ggplot2')
library('dplyr')

#' Load a csv located at specific location `filename` into a tibble
#' and merge PlateDay and Read columns into col called platedayNO
#'
#'
#' @param filename (str): the path to a specific file (ie 'file/path/to/file.tsv')
#'
#' @return tibble: a (n x n) tibble with a PlateDay, Read, Description columns followed by 
#' Concentration <dbl>, Signal <dbl> column
#' 
#' @note PlateDay and Read are two separate columns and it's better to conjoin the two columns
#' into one 
#' 
#' 
load_data <- function(filename){
  data <- read_csv(filename)
  data <- data %>% mutate(platedayNO = paste(PlateDay, Read), .after = Read) %>% 
    select(-c(PlateDay, Read))
  return(data)
}

#' Blank correction by subtracting blank from every other signal measured for standards,
#' test subjects, and quality control, not entirely sure if this needs to be done 
#'
#'
#' @param data (tibble): a (n x n) tibble with a platedayNO, Description, Signal, Concentrations columns 
#' 
#' @return tibble: a (n x n) tibble with a platedayNO, Description, Signal, Concentrations columns
#' followed by adj_signal columns containing the true absorbency signal 
#' 
#' 
blank_correction <- function(data){
  data <- data %>% 
    group_by( platedayNO_A) %>%
    mutate(mean_Signal = abs(mean_Signal[Description_A == 'BLANK'] - mean_Signal)) %>% 
    mutate(Signal_A = abs(Signal_A[Description_A == 'BLANK'] - Signal_A)) %>% 
    mutate(Signal_B = abs(Signal_B[Description_A == 'BLANK'] - Signal_B))
  return(data)
}


#' Because data is formatted row-wise, this makes it difficult to analyze and compute means 
#' for duplicates. This function converts the data so that it pivots wider and assigns 
#' each duplicate their own column
#'
#' @param data (tibble): a (n x n) tibble with a platedayNO, Description, Signal, Concentrations columns 
#' @param grp_dup (vector): a grouping scheme to assign duplicates their own row 
#' @param grp_AB (vector): a grouping scheme to assign names for columns 
#' @param colnames (vector): a vector containing a list of columns from which to take values from 
#' 
#' @return tibble: a (n x n) wider tibble with a platedayNO, Description, Signal, Concentrations columns for groups
#' A and B. 
#' 
#' 
data_wrangle <- function(data, gr_dup, gr_AB, colnames) {
  data$grp_dup <- gr_dup
  data$grp_AB <- gr_AB
  data <- data %>% 
    pivot_wider(id_cols = grp_dup, 
                names_from = grp_AB, 
                values_from = colnames) %>% 
    select (-c(Description_B, platedayNO_B))
  return(data)
}

#' This function aims to solve for the concentration of duplicates by solving 
#' for x using the slope and b and plugging it into the equation of the line.
#' 
#' @param data (tibble): a (n x n) tibble with a platedayNO, Description, Signal, Concentrations columns 
#' @param lin_mod (tibble): a (n x n) tibble with linear regression model describing the linear relationship between
#' signal and concentration. 
#' 
#' @return tibble: a (n x n) tibble with Concentration A and Concentration B of the technical replicates. 
#' 
#' 
solve_for_x <- function(data, lin_mod){
  data <- data %>%
    group_by(grp_dup, platedayNO_A) %>%
    mutate(Concentration_A = ((Signal_A - lin_mod[lin_mod$plateday_lm == platedayNO_A,]$b) / lin_mod[lin_mod$plateday_lm == platedayNO_A,]$slope)) %>% 
    mutate(Concentration_B = ((Signal_B - lin_mod[lin_mod$plateday_lm == platedayNO_A,]$b) / lin_mod[lin_mod$plateday_lm == platedayNO_A,]$slope))
  return(data)
}

#' CV is a way to measure the spread of a group in relation to the mean. To do this the std dev of technical 
#' replicates is divided by the mean. 
#' 
#' @param data (tibble): a (n x n) tibble with a platedayNO, Description, Signal, Concentrations columns 
#' 
#' @return tibble: a (n x n) tibble with std dev and cv of the technical replicates. 
#' 
find_cv <- function(data){
  data <- data %>%  
    group_by(grp_dup) %>% 
    mutate(mean_Concentration = mean(c(Concentration_A, Concentration_B))) %>% 
    mutate(sd_c = sd(c(Concentration_A, Concentration_B))) %>% 
    mutate(cv = (sd_c / mean_Concentration) * 100 )
  return(data)
}

#'  Filter data to extract only QCs for each plate type
#' 
#' @param data (tibble): a (n x n) tibble with a platedayNO, Description, Signal, Concentrations columns 
#' 
#' @return tibble: a (n x n) tibble with std dev and cv of the technical replicates. 
#' 
extract_qc <- function(data, pattern){
  data <- data %>% 
    group_by(platedayNO_A) %>% 
    filter(Description_A == 'Quality Control Samples', str_detect(platedayNO_A, pattern))
  return(data)
}

#' This function aims to solve for the upper and lower control limits of using the
#' corrected std dev value called sigma.
#' 
#' @param data (tibble): a (n x n) tibble with containing only QC values  
#' @param N (int): an integer with the total number of plates ran 
#'  
#' 
#' @return tibble: a (n x n) with low, mean, and high ranges for each QC level 
#' 
#' 
calculate_LUCL <- function(data, N){
  x_grandmean <- data %>% group_by(qc_lvl) %>% summarise(qc_mean = mean(mean_Concentration))
  s_bar <- data %>% group_by(qc_lvl) %>% summarise(sd_sk = mean(sd_c))
  num.an <- sqrt(2) * gamma(N/2)
  den.an <- sqrt(N-1) * gamma((N-1)/2)
  an <- num.an / den.an
  LCL <- x_grandmean - (3 * (s_bar / (an * sqrt(N))))
  UCL <- x_grandmean + (3 * (s_bar / (an * sqrt(N))))
  LU_CL <- cbind(low = LCL$qc_mean, grand_mean = x_grandmean$qc_mean,high = UCL$qc_mean)
  LU_CL <- as_tibble(LU_CL)
  LU_CL$qc_lvl <- rep(1:(nrow(LU_CL)))
  return(LU_CL)
}
