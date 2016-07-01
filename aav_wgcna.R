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



cd4.mad <- apply(cd4.t, 2, mad, na.rm = T)
cd4.mad <- cd4.mad[!is.na(cd4.mad)]
cd4.mad <- sort(cd4.mad, decreasing = T)
cd8.mad <- apply(cd8.t, 2, mad, na.rm = T)
cd8.mad <- cd8.mad[!is.na(cd8.mad)]
cd8.mad <- sort(cd8.mad, decreasing = T)
# Calculate the inflection point
#library("RootsExtremaInflections")
v2 <- 1 - rank(cd4.mad, ties.method = "random")/length(cd4.mad)
rv2 <- sort(v2, decreasing = TRUE)

#infl4 <- inflexi(cd4.mad[rv2], v2[rv2], 3000, 13000, 2, 5, "x", "y", 3,3)
#v4 <- cd4.mad[infl4$finfl][1]
#v2 <- 1 - rank(cd8.mad, ties.method = "random")/length(cd8.mad)
#rv2 <- rank(v2)
#infl8 <- inflexi(cd8.mad[rv2], v2[rv2], 3000, 13000, 2, 5, "x", "y", 3,3)
#v8 <- cd8.mad[infl8$finfl][1]



count.p <- function(data, per){
  # Calculate the percentatge above per of data
  sum(data >= per)/length(data)
}

perc.cd4 <- unlist(lapply(cd4.mad, count.p, data = cd4.mad))
perc.cd8 <- unlist(lapply(cd8.mad, count.p, data = cd8.mad))

# Calculating the inflection points
library("inflection")
infl.4 <- ede(cd4.mad, perc.cd4, 0)
infl.8 <- ede(cd8.mad, perc.cd8, 0)

print(infl.4)
print(infl.8)

png("MAD_filtering.png", width = 1200, height = 1200)
plot(cd4.mad, perc.cd4,
  ylab = "Proportion of points over", main = "MAD score", col = "blue")
points(cd8.mad, perc.cd8, col = "red")
abline(v = c(infl.8[1, 3], infl.4[1, 3]), col = c("blue","red"))
legend("topright", legend = c("CD8", "CD4"),
       fill = c("red", "blue"))
# Lacks of the filtering using inflexion point of the ranked list of mad.
dev.off()
#
# fitl <- loess(perc.cd4~cd4.mad)
# xl <- seq(min(cd4.mad), max(cd4.mad), (max(cd4.mad) - min(cd4.mad))/1000)
# out <- predict(fitl, xl)
# infl <- c(FALSE, diff(diff(out) > 0) != 0)
# plot(fitl)
#
# secant_method <- function(x, fun, x0, x1, iteration = 500, precision = 0.05, ...) {
#   for (i in 1:iteration ) {
#     x2 <- x[x1] - fun(x1, ...)*(x1 - x0)/(fun(x1, ...) - fun(x0, ...))
#     if (abs(x2 - x1) < precision) {
#       return(x2)
#     }
#     x0 <- x1
#     x1 <- x2
#   }
#   stop("Exceeded allowed number of interactions")
# }

exp_conditions <- list("CD4" = cd4.t, "CD8" = cd8.t)
n <- 0
for (exp in exp_conditions) {
  n <- n + 1
  nam <- names(exp_conditions)[n]
  cat(paste("Working with new data", nam, "\n"))
  # Checking if genes expression is ok
  gsg <- goodSamplesGenes(exp, verbose = 3)
  if (!gsg$allOK) {
    if (sum(!gsg$goodGenes) > 0) {
      cat(paste(sum(!gsg$goodGenes), "genes don't pass the filter.\n"))
      # genes <- paste(names(exp)[!gsg$goodGenes], collapse = ", ")
      # cat(paste("Removing genes:", genes))
    }
    if (sum(!gsg$goodSamples) > 0) {
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
  plot(sampleTree, main = nam)
  dev.off()
  clust <- cutreeStatic(sampleTree)
  print(table(clust))

  # Call the network topology analysis function
  sft <- pickSoftThreshold(exp, verbose = 5, RsquaredCut = 0.6)

  # Plot them
  png(paste0("threshold_", nam, ".png"))
  pars_org <- par(mfrow = c(1,2));
  cex1 = 0.9
  # Choose a set of soft-thresholding powers
  powers = c(c(1:10), seq(from = 12, to = 20, by = 2))
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
    stop("RsquaredCut option of the pickSoftThreshold too height.\n")
  }
  cat(paste0("Value of power used is ", sft$powerEstimate, ".\n"))

  # Do an automatic blocking
  net <- blockwiseModules(exp,
            power = sft$powerEstimate, # The default suggested value
            minModuleSize = 30,
            mergeCutHeight = 0.01,
            saveTOMs = TRUE,
            saveTOMFileBase = paste0("samples_", nam, "_TOM"),
            verbose = 1, loadTOM = TRUE)
  print(table(net$colors))

  pdf(paste0("modules_", nam, ".pdf"))

  # Plot the dendrogram and the module colors underneath
  plotDendroAndColors(net$dendrograms[[1]],
                      net$colors[net$blockGenes[[1]]],
                      "Module colors",
                      dendroLabels = FALSE, hang = 0.03,
                      addGuide = TRUE, guideHang = 0.05)
  dev.off()

  # Make the table of correlation between modules and clinical data
  # traits <- read.delim("fileinexistente.com")
  # MEs0 <- moduleEigengenes(exp, net$colors)
  # MEs <- orderMEs(MEs0$eigengenes)
  # moduleTraitCor <- cor(MEs, traits, use = "p")
  # moduleTraitPvalue <- corPvalueStudent(moduleTraitCor, ncol(exp))
  #
  # # plotMat()
  # # plotModuleSignificance
  # barplot(as.matrix(MEs0$varExplained),
  #         names = colnames(MEs0$eigengenes), las = 2)
  # boxplot(t(as.matrix(MEs0$varExplained)), main = "Explanation of the Eigengenes",
  #         xlab = "Eigengenes", ylab = "% Explained", outline = FALSE)
  # stripchart(t(as.matrix(MEs0$varExplained)), add = TRUE, vertical = T,
  #            pch = 21, bg = substring(colnames(MEs0$eigengenes), 3), method = "jitter")
  # # Will display correlations and their p-values
  # textMatrix <-  paste(ifelse(
  #   moduleTraitPvalue <= 0.1 & abs(moduleTraitCor) >= 0.3, "* ", ""),
  #   signif(moduleTraitCor, 2),
  #   " (",
  #   signif(moduleTraitPvalue, 1), ")",
  #   sep = "")
  # # Force to get the structure of teh matrix
  # dim(textMatrix) <-  dim(moduleTraitCor)
  # par(mar = c(3, 10, 4, 2));
  #
  # # Display the correlation values within a heatmap plot
  # labeledHeatmap(Matrix = moduleTraitCor,
  #                xLabels = "Shedding",
  #                yLabels = names(MEs),
  #                ySymbols = names(MEs),
  #                colorLabels = F,
  #                colors = blueWhiteRed(50),
  #                textMatrix = textMatrix,
  #                setStdMargins = F,
  #                cex.text = 0.5,
  #                zlim = c(-1,1),
  #                main = "Module-trait relationships",
  #                width = 1)
}
