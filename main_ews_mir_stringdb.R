library(disruptr)

##Define Globals
cache = "~/miRNA_targeting/data_files/"
experiment_name = "EWS_NP_string2"
ppi = "stringdb"
ncores = 20
n_genes = 50
n = 1000

#Function to load and prep the St. Jude data for entry to our pipeline
#going to rename the columns to match the columns for our other data - this will create some confusing column names but I'll rename them back
#after it goes through the pipeline
df <- read.csv(system.file("test_data/rld_Counts.csv", package = "disruptr"))

df <- compute_np(cache = cache, exp_mat = df,
                 mir_paper = TRUE, ncores = ncores, 
                 experiment_name = experiment_name, 
                 ppi = ppi, min_score = 700)

df <- compute_dnp(cache = cache, df = df,
                  experiment_name = experiment_name, ppi = ppi,
                  ncores = ncores)

df_null <- compute_null(cache = cache, experiment_name = experiment_name,
                        df = df, ppi = ppi, ncores = ncores,
                        n_genes = n_genes, n = n)
