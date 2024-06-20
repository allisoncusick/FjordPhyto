total_clust <- data.frame(total_clust[,-1], row.names= total_clust$Sample)
#Cluster LMG1510 con distancia Euclidea
matrixd <- dist(total_clust, method = "euclidean")
hd <- hclust(matrixd, method = "average")
library(dendextend)
dend10 <- as.dendrogram(hd)
# Color the branches based on the clusters:
dend10 <- color_branches(dend10, k=4) 
# We hang the dendrogram a bit:
dend10 <- hang.dendrogram(dend10,hang_height=0.1)
# reduce the size of the labels:
#dend <- assign_values_to_leaves_nodePar(dend, 0.5, "lab.cex")
dend10 <- set(dend10, "labels_cex", 0.7)
dend10 <- highlight_branches_lwd(dend10, values =4)
dend10 <- rotate(dend10, 1)
#Plot
plot(dend10, 
     main = "Spring and Autumn", 
     horiz =  TRUE,  nodePar = list(cex = .007), xlab = "Euclidean distance", axes = FALSE)
lines(x = c(0,0), y = c(0,100), type = "n") 
axis(side = 1, at = seq(0,100,10), labels = seq(100,0,-10))

require(vegan)
total_clust2 <- data.frame(total_clust2[,-1], row.names= total_clust2$Sample)
ano1<- anosim(matrixd, grouping = total_clust2$Clust, permutations = 999, distance = "euclidean")
plot(ano1)
