library(dplyr)
library(magrittr)
library(stringr)
library(tidyr)
source('~/miRNA_targeting/NP_Pipeline_MiR.R')
source('~/miRNA_targeting/NP_Funcs_MiR.R')

#Function to load and prep the St. Jude data for entry to our pipeline
#going to rename the columns to match the columns for our other data - this will create some confusing column names but I'll rename them back
#after it goes through the pipeline
prep_sj <- function() {
  load("~/miRNA_targeting/data_files/stjude_counts.Rda")
  df <- filter(jude_df, sj_diseases == "EWS") %>% 
    select(-scaled_log_count, -sj_diseases, -attr_oncotree_disease_code,
           -STDIN) %>% 
    rename(gene_name = Geneid, 
           experiment_id = sample_name, 
           cell_line = subject_name, 
           expression = log_count) %>% 
    mutate(condition = "patient", 
           experiment_num = 1)
}

run_sj <- function() {
  df <- prep_sj()
  df_np <- gibbs_pipeline(df, ncores = 20, filename = "sj_ews")
}

run_sj()

