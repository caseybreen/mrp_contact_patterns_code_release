# This R.script contains helpful functions for working with contact matrices.
#
#

###########################################
# Library all Packages for analysis      #
###########################################
library(tidyverse)
library(ipumsr)
library(here)
library(haven)
library(janitor)
library(ggthemes)
library(cowplot)
library(broom)
library(brms)
library(tidybayes)
library(tictoc)
library(lubridate)
library(loo)
library(cmdstanr) # install.packages("cmdstanr", repos = c('https://stan-dev.r-universe.dev', getOption("repos")))
library(cowplot)
library(geofacet)
library(tidyr)
library(Matrix)
library(RSpectra)
library(socialmixr)
library(mice)
library(factoextra)
library(gghighlight)
library(glue)


###########################################
# global options                          #
###########################################

options(dplyr.summarise.inform = FALSE)

###########################################
# calculate leading eigenvalue function   #
###########################################

calculate_leading_eigenvalue <- function(input_matrix) {
  
  if (class(input_matrix)[[1]] != "matrix")
  {input_matrix = convert_df_to_matrix(input_matrix)}

  eigen = Re(RSpectra::eigs(input_matrix, 2, which = "LM", opts = list(retvec = FALSE))$values[[1]]) # faster than base::eigen()
  
  return(eigen)
}

###########################################
# calculate second eigenvalue function   #
###########################################

calculate_second_eigenvalue <- function(input_matrix) {
  
  eigen = Re(RSpectra::eigs(input_matrix, 2, which = "LM", opts = list(retvec = FALSE))$values[[2]]) # faster than base::eigen()
  
  return(eigen)
}

###########################################
# calculate left eigenvectors           #
###########################################

calculate_left_eigenvector <- function(input_matrix) {
  
  eigen = Re(RSpectra::eigs(input_matrix, 2, which = "LM")$vectors[,0:1]) # faster than base::eigen()]
  
  eigen <- data.frame(agecat =rownames(input_matrix), eigenvector = eigen)
  
  return(eigen)
  
}

###########################################
# calculate right eigenvectors           #
###########################################

calculate_right_eigenvector <- function(input_matrix) {
  
  input_matrix <- t(input_matrix)
  
  eigen = Re(RSpectra::eigs(input_matrix, 2, which = "LM")$vectors[,0:1]) # faster than base::eigen()]
  
  eigen <- data.frame(agecat =rownames(input_matrix), eigenvector = eigen)
  
  return(eigen)
  
}

###########################################
# enforce reciprocity function            # 
###########################################

enforce_reciprocity <- function(df) {
  
  ## rename variables 
  df_sym <- df %>%
    mutate(.alter_age = case_when(
      .category == "alterage018" ~ "[0,18)",
      .category == "alterage1825" ~ "[18,25)",
      .category == "alterage2535" ~ "[25,35)",
      .category == "alterage3545" ~ "[35,45)",
      .category == "alterage4555" ~ "[45,55)",
      .category == "alterage5565" ~ "[55,65)",
      .category == "alterage65100" ~ "[65,100]"
    )) %>%
    mutate(.ego_age = as.character(agecat)) %>%
    mutate(unadj_avg_per_ego = mean)
  
  ## add on youngest age category 
  df_sym <- df_sym %>%
    bind_rows(df_sym %>% mutate(.ego_age = "[0,18)") %>%
      tidyr::expand(.ego_age, .alter_age, state, day_std, .draw))

  ## reciprocal matrix
  reciprocal_matrix <- df_sym %>%
    select(state, day_std, .draw,
           .alter_age = .ego_age,
           .ego_age = .alter_age,
           other_unadj_avg_per_ego = unadj_avg_per_ego
    )
  
  df_sym <- df_sym %>%
    left_join(acs_state_symmetry %>%
                select(ego_acs_N = count, .ego_age = agecat, state), by = c(".ego_age", "state")) %>%
    left_join(acs_state_symmetry %>%
                select(alter_acs_N = count, .alter_age = agecat, state), by = c(".alter_age", "state"))
  
  matrix_sym <- df_sym %>%
    left_join(reciprocal_matrix,
              by = c(".ego_age", ".alter_age", "state", "day_std", ".draw")
    ) %>%
    mutate(sym_avg_per_ego = case_when(
      .alter_age %in% "[0,18)" ~ unadj_avg_per_ego,
      .ego_age %in% "[0,18)" ~ (other_unadj_avg_per_ego * alter_acs_N) / ego_acs_N,
      TRUE ~ (1 / (2 * ego_acs_N)) *
        ((unadj_avg_per_ego * ego_acs_N) + (other_unadj_avg_per_ego * alter_acs_N))
    ))
  
  leading_eigenvalues <- matrix_sym %>%
    filter(.alter_age != "[0,18)" & .ego_age != "[0,18)") %>%
    group_by(state, day_std, .draw) %>%
    nest() %>%
    mutate(leading_eigen = map(data, ~ calculate_leading_eigenvalue(input_matrix = .))) %>%
    unnest(leading_eigen) %>%
    select(state, day_std, leading_eigen, .draw)
  
  ## polymod data
  polymod_country <- "United Kingdom"
  polymod_age_limits <- c(0, 18, 25, 35, 45, 55, 65)
  
  ## get polynmod data 
  ## supress warnings so it doesn't clog up pipeline 
  polymod_mat <- suppressMessages(suppressWarnings(getPolymodMatrixNoSchool(
    polymod_country = polymod_country,
    polymod_age_limits = polymod_age_limits
  )))
  
  polymod_mat_modified <- polymod_mat[-1, -1]
  
  youngest_age_data <- leading_eigenvalues %>%
    mutate(sym_avg_per_ego = polymod_mat[1, 1] * leading_eigen / eigen(polymod_mat_modified)$values[[1]]) %>%
    mutate(
      .alter_age = "[0,18)",
      .ego_age = "[0,18)"
    )
  
  matrix_sym <- youngest_age_data %>%
    bind_rows(matrix_sym %>%
                filter(!(.alter_age == "[0,18)" & .ego_age == "[0,18)"))) %>%
    arrange(day_std, state) %>%
    ungroup()
  
  return(matrix_sym)
}


###########################################
# Calculate eigenvalues                   #
###########################################

eigenvalue_ratio <- function(df) {
  
  ## calculate eigenvalues
  eigenvalue_df <- df %>% 
    group_by(state, day_std, .draw) %>%
    nest() %>%
    mutate(summary_stats = map(data, ~summarize_contact(input_df = .))) %>%
    unnest(summary_stats) %>%
    select(state,
           day_std,
           leading_eigenvalue,
           second_eigenvalue, 
           right_eigenvector,
           left_eigenvector,
           R0_agespecific,
           .draw
    )
  
  ##
  eigenvalue_df <- eigenvalue_df %>%
    mutate(q_index = second_eigenvalue/leading_eigenvalue)
  
  return(eigenvalue_df)
}


###################################
# Total Contact                   # 
###################################

average_contact <- function(input_df ) {
  
  input_df <- input_df %>% 
    group_by(state, day_std, .ego_age, .alter_age) %>% 
    summarize(contacts = mean(sym_avg_per_ego),
              lower = quantile(sym_avg_per_ego, 0.1),
              upper = quantile(sym_avg_per_ego, 0.9))
  
  return(input_df)
} 


###################################
# summarize_contact               # 
###################################


summarize_contact <- function(input_df){
  
  input_matrix <- convert_df_to_matrix(input_df)
  
  ## calculate summary stats
  summary_stats <- tibble(
    leading_eigenvalue = calculate_leading_eigenvalue(input_matrix = input_matrix),
    second_eigenvalue = calculate_second_eigenvalue(input_matrix = input_matrix),
    right_eigenvector = list(calculate_right_eigenvector(input_matrix = input_matrix)),
    left_eigenvector = list(calculate_left_eigenvector(input_matrix = input_matrix)),
    R0_agespecific = compute_R0_relative(C = input_matrix)
  )
  
  return(summary_stats)
  
}

## convert data.frame to matrix 

## faster 
convert_df_to_matrix <- function(input) {
  
 # m_old <- Matrix(nrow = sqrt(length(input$sym_avg_per_ego)), ncol = sqrt(length(input$sym_avg_per_ego)), data = input$sym_avg_per_ego)
  matrix <- input %>% select(sym_avg_per_ego, .ego_age, .alter_age) %>% 
    pivot_wider(names_from = .ego_age, values_from = sym_avg_per_ego) %>%
    column_to_rownames(var = ".alter_age") %>% 
    as.matrix()
  
  return(matrix)
  
}

# convert_df_to_matrix2 <- function(input) {
# 


## Function for creating dummy weekday / weekend variable 
is_weekday = function(timestamp){
  lubridate::wday(timestamp, week_start = 1) < 6
}

## Source .Rmd 

source_rmd = function(file, ...) {
  tmp_file = tempfile(fileext=".R")
  on.exit(unlink(tmp_file), add = TRUE)
  knitr::purl(file, output=tmp_file, quiet = T)
  source(file = tmp_file, ...)
}

## function to get polymod contact matrix , with school contacts dropped, based on country and age groupings

getPolymodMatrixNoSchool <- function(polymod_country = "United Kingdom",
                                     polymod_age_limits = c(0, 18, 25, 35, 45, 55,  65)){
  #need the socialmixr package loaded
  data(polymod)
  data_part <- polymod$participants
  data_cnt <- polymod$contacts %>% filter(cnt_school == 0)
  
  poly_no_school <- socialmixr::contact_matrix(survey(data_part, data_cnt), 
                                               countries = polymod_country, 
                                               age.limits = polymod_age_limits, symmetric = TRUE)
  polymod_no_school_mat <- poly_no_school$matrix
  dimnames(polymod_no_school_mat)[[1]] <- dimnames(polymod_no_school_mat)[[2]]
  return(polymod_no_school_mat)
}

## Compute R0 Relative 


compute_R0_relative = function(C){
  #set fixed parameters:
  
  gamma <- 1/6 # recovery period (I -> R), ref: Davies
  rel.symp <- c(0.26, 0.24, 0.3, 0.36, 0.44, 0.48, 0.69) # rho 
  u = c(0.39, 0.62, 0.81, 0.83, 0.81, 0.74)/5 # susceptibility for each age class
  
  alpha <- 0.5
  
  # Davies NGM
  Du <- diag(u, 7)
  Dy <- diag(1/gamma, 7)
  Dx <- diag((1-rel.symp)*alpha+rel.symp, 7)
  NGM <- Du %*% C %*% Dy %*% Dx 
  R0  <- abs(eigen(NGM)$values[1])
  return(R0)
}


