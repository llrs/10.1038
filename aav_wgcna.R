# WGCNA of the AAV data

source("classify_sets.R")

library("WGCNA")
options(stringsAsFactors = FALSE)
enableWGCNAThreads()

# Transpose to obtain the desired input for WGCNA
cd4.t <- t(cd4.44)
cd8.t <- t(cd8.44)

# Transform to numeric values
cd4.t <- apply(cd4.t, 2, as.numeric)
cd8.t <- apply(cd8.t, 2, as.numeric)

# Mantain the rownames
rownames(cd4.t) <- colnames(cd4.44)
rownames(cd8.t) <- colnames(cd8.44)


exp2 <- function(x) {
  # Exponential of 2/ reverse of log2
  return(2 ^ (x))
}

# Transform the log expression to raw expression
cd4.exp <- apply(cd4.t, 1:2, exp2)
cd8.exp <- apply(cd8.t, 1:2, exp2)

exp_conditions <- list("CD4" = cd4.exp, "CD8" = cd8.exp)
n <- 0
for (exp in exp_conditions) {
  n <- n + 1
  cat("Working with new data\n")
  # Checking if genes expression is ok
  gsg <- goodSamplesGenes(exp, verbose = 1)
  if (!gsg$allOK) {
    stop("Should check the sample.")
  }
  sampleTree <- hclust(dist(exp), method = "average")
  plot(sampleTree)
  clust <- cutreeStatic(sampleTree)
  table(clust)
  
  # Choose a set of soft-thresholding powers
  powers = c(c(1:10), seq(from = 12, to = 20, by = 2))
  # Call the network topology analysis function
  sft <- pickSoftThreshold(exp, powerVector = powers, verbose = 5)
  
  cat(paste0("Value of power used is ", sft$powerEstimate, "\n"))
  
  # Do an automatic blocking
  net <- blockwiseModules(exp, 
            power = sft$powerEstimate, # The default suggested value
            minModuleSize = 30, 
            mergeCutHeight = 0.01,
            saveTOMs = TRUE,
            saveTOMFileBase = paste0("sample_", n, "_TOM"),
            verbose = 1)
  
  table(net$colors)
  
  
  # Make the table of correlation between modules and clinical data
  traits <- read.delim("fileinexistente.com")
  MEs0 <- moduleEigengenes(exp, net$colors)
  MEs <- orderMEs(MEs0$eigengenes)
  moduleTraitCor <- cor(MEs, traits, use = "p")
  moduleTraitPvalue <- corPvalueStudent(moduleTraitCor, ncol(exp))
  
  # plotMat()
  # plotModuleSignificance
  barplot(as.matrix(MEs0$varExplained), 
          names = colnames(MEs0$eigengenes), las = 2)
  boxplot(t(as.matrix(MEs0$varExplained)), main = "Explanation of the Eigengenes",
          xlab = "Eigengenes", ylab = "% Explained", outline = FALSE)
  stripchart(t(as.matrix(MEs0$varExplained)), add = TRUE, vertical = T,
             pch = 21, bg = substring(colnames(MEs0$eigengenes), 3), method = "jitter")
  sizeGrWindow(10,6)
  # Will display correlations and their p-values
  textMatrix <-  paste(ifelse(
    moduleTraitPvalue <= 0.1 & abs(moduleTraitCor) >= 0.3, "* ", ""),
    signif(moduleTraitCor, 2), 
    " (",
    signif(moduleTraitPvalue, 1), ")", 
    sep = "")
  # Force to get the structure of teh matrix
  dim(textMatrix) <-  dim(moduleTraitCor)
  par(mar = c(3, 10, 4, 2));
  
  # Display the correlation values within a heatmap plot
  labeledHeatmap(Matrix = moduleTraitCor,
                 xLabels = "Shedding",
                 yLabels = names(MEs),
                 ySymbols = names(MEs),
                 colorLabels = F,
                 colors = blueWhiteRed(50),
                 textMatrix = textMatrix,
                 setStdMargins = F,
                 cex.text = 0.5,
                 zlim = c(-1,1),
                 main = "Module-trait relationships",
                 width = 1)
}
