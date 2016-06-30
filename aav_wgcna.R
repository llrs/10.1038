# WGCNA of the AAV data

source("classify_sets.R")

library("WGCNA")
options(stringsAsFactors = FALSE)
enableWGCNAThreads()

# Transpose to obtain the desired input for WGCNA
cd4.t <- t(cd4.44)
cd8.t <- t(cd8.44)
cd4.t <- as.matrix(apply(cd4.t, 2, as.numeric), nrow(cd4.t), ncol(cd4.t))
cd8.t <- as.matrix(apply(cd8.t, 2, as.numeric), nrow(cd8.t), ncol(cd8.t))

# Mantain the rownames
rownames(cd4.t) <- colnames(cd4.44)
rownames(cd8.t) <- colnames(cd8.44)


# exp2 <- function(x) {
#   # Exponential of 2/ reverse of log2
#   return(2 ^ (x))
# }
#
# # Transform the log expression to raw expression
# cd4.exp <- apply(cd4.t, 1:2, exp2)
# cd8.exp <- apply(cd8.t, 1:2, exp2)

# Calculate the median absolute deviation values for all proves
# cd4.mad <- apply(cd4.t, 2, mad, na.rm = T)
# cd4.mad <- cd4.mad[!is.na(cd4.mad)]
# cd8.mad <- apply(cd8.t, 2, mad, na.rm = T)
# cd8.mad <- cd8.mad[!is.na(cd8.mad)]
# plot(cd4.mad, 1 - rank(cd4.mad)/length(cd4.mad),
# ylab = "% of points over", main = "CD4)
# plot(cd8.mad, 1 - rank(cd8.mad)/length(cd8.mad),
#      ylab = "% of points over", main = "CD8")

# Lacks of the filtering using inflexion point of the ranked list of mad.

exp_conditions <- list("CD4" = cd4.t, "CD8" = cd8.t)

n <- 0
for (exp in exp_conditions) {
  nam <- names(exp_conditions)[n]
  cat(paste("Working with new data", nam, "\n"))
  n <- n + 1
  # Checking if genes expression is ok
  gsg <- goodSamplesGenes(exp, verbose = 3)
  if (!gsg$allOK) {
    if (sum(!gsg$goodGenes) > 0) {
      cat(paste(sum(!gsg$goodGenes), "genes don't pass the filter.\n"))
      # genes <- paste(names(exp)[!gsg$goodGenes], collapse = ", ")
      # cat(paste("Removing genes:", genes))
    }
    if (sum(!gsg$goodSamples) > 0){
      cat(paste(sum(!gsg$goodSamples)), "samples don't pass the filter.\n")
      # samples <- paste(rownames(exp)[!gsg$goodSamples], collapse = ", ")
      # cat(paste("Removing samples:", samples))
    }
    cat("Removing genes and samples which don't pass the filter.\n")
    # Remove the offending genes and samples from the data:
    exp <- exp[gsg$goodSamples, gsg$goodGenes]

  }
  sampleTree <- hclust(dist(exp), method = "average")
  pdf(paste0("sampletree_", nam, ".pdf"))
  plot(sampleTree)
  dev.off()
  clust <- cutreeStatic(sampleTree)
  print(table(clust))

  # Call the network topology analysis function
  sft <- pickSoftThreshold(exp, verbose = 5, RsquaredCut = 0.6)
  print(str(sft))

  # Plot them
  sizeGrWindow(9, 5)
  pars_org <- par(mfrow = c(1,2));
  cex1 = 0.9
  # Choose a set of soft-thresholding powers
  powers = c(c(1:10), seq(from = 12, to = 20, by = 2))
  pdf(paste0("threshold_", nam, ".pdf"))
  # Scale-free topology fit index as a function of the soft-thresholding power
  plot(sft$fitIndices[, "Power"],
       -sign(sft$fitIndices[,"slope"])*sft$fitIndices[, "SFT.R.sq"],
       xlab = "Soft Threshold (power)",
       ylab = "Scale Free Topology Model Fit, signed R^2", type = "n",
       main = "Scale independence")
  text(sft$fitIndices[, "Power"],
       -sign(sft$fitIndices[,"slope"])*sft$fitIndices[, "SFT.R.sq"],
       labels = powers, cex = cex1, col = "red")
  # this line corresponds to using an R^2 cut-off of h
  abline(h = 0.90, col = "red")
  # Mean connectivity as a function of the soft-thresholding power
  plot(sft$fitIndices[, "Power"], sft$fitIndices[, "mean.k."],
       xlab = "Soft Threshold (power)", ylab = "Mean Connectivity", type = "n",
       main = "Mean connectivity")
  text(sft$fitIndices[, "Power"], sft$fitIndices[, "mean.k."], labels = powers,
       cex = cex1, col = "red")
  dev.off()

  if (is.null(sft$powerEstimate)) {
    stop("RsquaredCut option of the pickSoftThreshold to hight.\n")
  }
  cat(paste0("Value of power used is ", sft$powerEstimate, "\n"))

  # Do an automatic blocking
  net <- blockwiseModules(exp,
            power = sft$powerEstimate, # The default suggested value
            minModuleSize = 30,
            mergeCutHeight = 0.01,
            saveTOMs = TRUE,
            saveTOMFileBase = paste0("samples_", nam, "_TOM"),
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
