
# Setup -------------------------------------------------------------------

# load packages
library(tidyverse)

# load data
input_map <- read_tsv('data/raw_data/My_SNPS.map', 
                      col_names = c('chromosome', 'snp_id', 'genetic_distance',
                                    'chromosome_position')) %>%
  # add a variable for easier combining later
  mutate(SNP_name = paste0('SNP', seq(1,5)))


input_ped <- read_delim('data/raw_data/My_SNPS.ped', delim = ' ',
                        col_names = c('family_id', 'individual_id', 'father_id',
                                      'mother_id', 'sex', 'affected_status', 
                                      # unique identifier for the SNPs columns
                                      paste0('SNP', rep(1:5, each=2), '_', c(1,2))
                                      ))


# Data manipulation -----------------------------------------------------

## basic mutations for objectives 1 and 2 first

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

## General data munging for objective 3

# currently we have a load of redundancy as we have columns such as SNP1_1, SNP1_2
# etc, and if I understand correctly we want them to be just one column for each,
# containing both values
for_transposing <- processed_ped %>%
  unite('SNP1', SNP1_1, SNP1_2, sep = ' ') %>%
  unite('SNP2', SNP2_1, SNP2_2, sep = ' ') %>%
  unite('SNP3', SNP3_1, SNP3_2, sep = ' ') %>%
  unite('SNP4', SNP4_1, SNP4_2, sep = ' ') %>%
  unite('SNP5', SNP5_1, SNP5_2, sep = ' ')

trimmed_ped <- for_transposing %>%
  # select only SNP columns
  select(starts_with('SNP')) %>%
  # transpose it
  t() %>%
  # convert it to a tibble, for easier manipulation. Store SNP_name in a column
  as_tibble(rownames = 'SNP_name') %>%
  # combine all of the DNA columns into one
  unite('dna_string', starts_with('V'), sep = ' ')

# create the output trans tped file by combining the above object with the 
# initial map file
trans_tped <- input_map %>%
  left_join(trimmed_ped) %>%
  # remove redundant 'SNP_name' column
  select(-SNP_name)

# save output
write_delim(trans_tped,
            file = 'data/processed_data/trans.tped',
            delim = ' ',
            col_names = F,
            quote = 'none')



## now create trans.tfam file (just the first six columns of the ped file)
processed_ped %>%
  select(!starts_with('SNP')) %>%
  select(!c(is_even, family_id)) %>%
  write_tsv('data/processed_data/trans.tfam',
            col_names = F)
