
# Setup -------------------------------------------------------------------

# load packages
library(tidyverse)

# load data
input_map <- read_tsv('data/raw_data/My_SNPS.map', 
                      col_names = c('chromosome', 'snp_id', 'genetic_distance',
                                    'chromosome_position'))
input_ped <- read_delim('data/raw_data/My_SNPS.ped', delim = ' ',
                        col_names = c('family_id', 'individual_id', 'father_id',
                                      'mother_id', 'sex', 'affected_status'))


# Data manipulation -----------------------------------------------------

## basic mutations first

processed_ped <- input_ped %>%
  # add prefix to family ID
  mutate(family_id = paste0('INCH_', family_id), 
  # make a column stating if a sample ID is even or odd, for clarity in the
  # later steps
  is_even = individual_id %% 2 == 0,
  # change phenotype status, making all even-numbered samples OTHER THAN
  # sample #27 have a status of '2' (affected). Keep even sample_ids unaltered
  affected_status = ifelse(is_even == T & individual_id !=27, 2, affected_status),
  # now change the phenotype status, making all odd-numbered samples OTHER THAN 
  # sample #27 have a status of 1 (affected). Keep odd sample_ids unaltered
  affected_status = ifelse(is_even == F & individual_id !=27, 1, affected_status))


