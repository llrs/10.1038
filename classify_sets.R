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

########################### AAV ###########################
names(aav) <-basename(aav)

# Relays on the fact that there is only one file with CDX
cd8 <- read.delim(grep("CD8", aav, value = T), row.names = 1)
cd4 <- read.delim(grep("CD4", aav, value = T), row.names = 1)
cd8 <- cd8[-1, ]
cd4 <- cd4[-1, ]

cd8.names <- strsplit(colnames(cd8), ".", fixed = T)
cd4.names <- strsplit(colnames(cd4), ".", fixed = T)
  
cd8.patients <- as.character(lapply(cd8.names, `[[`, 1))
cd4.patients <- as.character(lapply(cd4.names, `[[`, 1))

# Intersect of the samples in both datasets
samples_a <- intersect(cd4.patients, cd8.patients)
# Selecting those that start with V
# They are the n= 44 showed on the letter
samples <- samples[grepl("^V", samples.a)]
samples <- paste0(samples, ".")
samples.a <- paste0(samples.a, ".")

# Filtering for just this 44 samples
cd8.44 <- cd8[, as.numeric(lapply(samples, grep, colnames(cd8), fixed = T))]
cd4.44 <- cd4[, as.numeric(lapply(samples, grep, colnames(cd4), fixed = T))]

# Filtering for just this 58 samples
cd8.58 <- cd8[, as.numeric(lapply(samples.a, grep, colnames(cd8), fixed = T))]
cd4.58 <- cd4[, as.numeric(lapply(samples.a, grep, colnames(cd4), fixed = T))]
########################### SLE ###########################


