# Program to reproduce doi:10.1038/nature14468

# Analyse antibody-associated vasculitis AAV
# Systemic lupus erythematosus SLE
# Inflammatory bowel disease IBD
# Murine lymphocytic choriomeningitis virus LCMV
# malaria vaccination MV
# HCV
# influenza vaccination IV
# yellow fever vaccination YFV
# dengue fever DF
# Idiopathic pulmonary fibrosis IPF
# type 1 diabetes T1D
# NOD # A type of mice
# rheumatoid arthritis RA
# in vitro CD8 stimulation 

# Provide with objects with where is the data stored
source("download_datasets.R")

sdrf <- list.files(data.folder, pattern = "(sdrf.txt)", full.names = T)
idf <- list.files(data.folder, pattern = "(idf.txt)", full.names = T)
raw <- list.files(data.folder, pattern = "(.raw.*.zip)", full.names = T)
adf <- list.files(data.folder, pattern = "(adf)", full.names = T)
txt <- list.files(data.folder, pattern = ".txt", full.names = T)
aav <- list.files(data.folder, pattern = "aav", full.names = T)
processed <- list.files(data.folder, pattern = "(.processed.*.zip)", 
                        full.names = T)
names(aav) <-basename(aav)

# Relays on the fact that there is only one file with CDX
cd8 <- read.delim(grep("CD8", aav, value = T), row.names = 1)
cd4 <- read.delim(grep("CD4", aav, value = T), row.names = 1)

