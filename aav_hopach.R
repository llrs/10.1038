# hopach of the AAV data

source("classify_sets.R")

library("hopach")

cd4.44 <- apply(cd4.44, 1:2, as.numeric)
gene.dist <- distancematrix(cd4.44, "cosangle", verbose = TRUE)
gene.hobj <- hopach(cd4.44, dmat = gene.dist)

print(table(gene.hobj$clust$labels))

png("cluster_hopach_CD4.png", width = height = 1500)
dplot(gene.dist,gene.hobj, ord = "final",
      main = "Cluster hopach CD4", showclusters = F)
dev.off()
