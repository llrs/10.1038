# WGCNA of the AAV data

source(classify_sets.R)

library("WGCNA")
options(stringsAsFactors = FALSE)
enableWGCNAThreads()

cd4.t <- t(cd4.44)
cd8.t <- t(cd8.44)

exp2 <- function(x){
  # Exponential of 2/ reverse of log2
  return(2 ^ x)
}

cd4.exp <- apply(cd4.t, 1:2, exp2)
cd8.exp <- apply(cd8.t, 1:2, exp2)

gsg.4 <-  goodSamplesGenes(cd4.exp, verbose = 3)
gsg.8 <-  goodSamplesGenes(cd8.exp, verbose = 3)
if (!gsg.4$allOK) {
  stop("Should be checked sample of cd4")
} else if (!gsg.8$allOK) {
  stop("Should be checked sample of cd8")
}