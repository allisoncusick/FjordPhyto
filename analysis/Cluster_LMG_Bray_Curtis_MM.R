clust2 <- data.frame(clust1_LMG1510[,-1], row.names= clust1_LMG1510$Sample)
library(vegan)
d = (1 - vegdist(clust2, method="bray")) * 100
h = hclust(d, method = "average")
plot(h, main = "LMG1510 Cluster using Bray Curtis method", sub = "", xlab="", axes = FALSE, hang = -1)
lines(x = c(0,0), y = c(0,100), type = "n") # force extension of y axis
axis(side = 2, at = seq(0,100,10), labels = seq(100,0,-10))

d = (1 - vegdist(clust, method="euclidean")) * 100
h = hclust(d, method = "average")
plot(h, main = "LMG1510 Cluster using Bray Curtis method", sub = "", xlab="", axes = FALSE, hang = -1)
lines(x = c(0,0), y = c(0,100), type = "n") # force extension of y axis
axis(side = 2, at = seq(0,100,10), labels = seq(100,0,-10))

#Cluster LMG1510 con distancia Euclidea
matrixd <- dist(clust, method = "euclidean")
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
#Plot
plot(dend10, 
     main = "LMG1510", 
     horiz =  TRUE,  nodePar = list(cex = .007), xlab = "Euclidean distance", axes = FALSE)
lines(x = c(0,0), y = c(0,100), type = "n") 
axis(side = 1, at = seq(0,100,10), labels = seq(100,0,-10))
