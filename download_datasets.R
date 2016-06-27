# Program to reproduce doi:10.1038/nature14468

data.folder <- "~/Documents/clinic_correlation"

AAV.SLE <- c("E-MTAB-2452", "E-MTAB-157", "E-MTAB-145")
IBD <- "E-MTAB-331"
T1D <- "E-TABM-666"
CD8s <- "E-MTAB-3470"

arr.exp <- c(AAV.SLE, IBD, T1D, CD8s)
paths.ae <- file.path(data.folder, arr.exp)

library("ArrayExpress", quietly = T, verbose = F)

# Download from ArrayExpress
# lapply(arr.exp, getAE, type = "full", extract = T, path = data.folder)

library("GEOquery", quietly = T, verbose = F)

NOD <- "GSE21897"
arthritis <- c("GSE15258", "GSE33377")
LCMV <- "GSE9650"
HCV <- "GSE7123"
malaria <- "GSE18323"
influenza <- "GSE29619"
YFW <- "GSE13486"
dengue <- "GSE25001"
IPF <- "GSE28221"

gse.list <- c(NOD, arthritis, LCMV, HCV, malaria, influenza, YFW, dengue, IPF)

# lapply(gse.list, getGEOSuppFiles, baseDir = data.folder)

# List of files and data:
names(paths.ae) <- names(arr.exp)
paths.gse <- file.path(data.folder, gse.list)
names(paths.gse) <- names(gse.list)
paths.data <- c(paths.gse, paths.ae)

all.files <- list.files(path = data.folder, full.names = T)
