#Cluster Completo
#Para poner como nombre la primer columna
Complete <- data.frame(Complete[,-1], row.names= Complete$Sample)
d <- dist(Complete)
cluster_comp <- hclust(d, method = "average")
plot(cluster_comp)

library(dendextend)
dend <- as.dendrogram(cluster_comp)
# order it the closest we can to the order of the observations:
dend <- rotate(dend, 1:45)

# Color the branches based on the clusters:
dend <- color_branches(dend, k=4) 
# We hang the dendrogram a bit:
dend <- hang.dendrogram(dend,hang_height=0.1)
# reduce the size of the labels:
#dend <- assign_values_to_leaves_nodePar(dend, 0.5, "lab.cex")
dend <- set(dend, "labels_cex", 0.4)

#Plot
plot(dend, 
     main = "LMG1510", 
     horiz =  TRUE,  nodePar = list(cex = .007))

#Cluster Promedio
#Para poner como nombre la primer columna
samp.with.rownames <- data.frame(Promedio[,-1], row.names= Promedio$Sample)
d2 <- dist(samp.with.rownames)
cluster_prom <- hclust(d2, method = "complete")
plot(cluster_prom)

library(dendextend)
dend2 <- as.dendrogram(cluster_prom)
# order it the closest we can to the order of the observations:
dend2 <- rotate(dend2, 1:15)

# Color the branches based on the clusters:
dend2 <- color_branches(dend2, k=4) 
# We hang the dendrogram a bit:
dend2 <- hang.dendrogram(dend2,hang_height=0.1)
# reduce the size of the labels:
#dend2 <- assign_values_to_leaves_nodePar(dend2, 0.5, "lab.cex")
dend2 <- set(dend2, "labels_cex", 0.6)

#Plot
plot(dend2, main = "LMG1510", horiz =  TRUE,  nodePar = list(cex = .007))


#Cluster Superficial
#Para poner como nombre la primer columna
sup <- data.frame(Superficial[,-1], row.names= Superficial$Sample)
d3 <- dist(sup)
cluster_sup <- hclust(d3, method = "complete")
plot(cluster_sup)

library(dendextend)
dend3 <- as.dendrogram(cluster_sup)
# order it the closest we can to the order of the observations:
dend3 <- rotate(dend3, 1:15)

# Color the branches based on the clusters:
dend3 <- color_branches(dend3, k=4) 
# We hang the dendrogram a bit:
dend3 <- hang.dendrogram(dend3,hang_height=0.1)
# reduce the size of the labels:
#dend2 <- assign_values_to_leaves_nodePar(dend2, 0.5, "lab.cex")
dend3 <- set(dend3, "labels_cex", 0.6)

#Plot
plot(dend3, main = "LMG1510", horiz =  TRUE,  nodePar = list(cex = .007))

#Para comparar 2 dendrogramas
tanglegram(dend2, dend3, main_left = "Average", main_right = "50% Light")

#Cluster Bio superficial
#Para poner como nombre la primer columna
bio_sup <- data.frame(bio_sup [,-1], row.names= bio_sup$Sample)
d3 <- dist(bio_sup)
cluster_bio <- hclust(d3, method = "complete")
plot(cluster_bio)

library(dendextend)
dend3 <- as.dendrogram(cluster_bio)
# order it the closest we can to the order of the observations:
dend3 <- rotate(dend3, 1:45)

# Color the branches based on the clusters:
dend3 <- color_branches(dend3, k=4) 
# We hang the dendrogram a bit:
dend3 <- hang.dendrogram(dend3,hang_height=0.1)
# reduce the size of the labels:
#dend <- assign_values_to_leaves_nodePar(dend, 0.5, "lab.cex")
dend3 <- set(dend3, "labels_cex", 0.5)

#Plot
plot(dend3, 
     main = "LMG1510", 
     horiz =  TRUE,  nodePar = list(cex = .007))

#Cluster Promedio superficial grupos
#Para poner como nombre la primer columna
abun_sup <- data.frame(abun_sup [,-1], row.names= abun_sup$Sample)
d4 <- dist(abun_sup)
cluster_abun <- hclust(d4, method = "complete")
plot(cluster_abun)

library(dendextend)
dend4 <- as.dendrogram(cluster_abun)
# order it the closest we can to the order of the observations:
dend4 <- rotate(dend4, 1:45)

# Color the branches based on the clusters:
dend4 <- color_branches(dend4, k=4) 
# We hang the dendrogram a bit:
dend4 <- hang.dendrogram(dend4,hang_height=0.1)
# reduce the size of the labels:
#dend <- assign_values_to_leaves_nodePar(dend, 0.5, "lab.cex")
dend4 <- set(dend4, "labels_cex", 0.5)

#Plot
plot(dend4, 
     main = "LMG1510", 
     horiz =  TRUE,  nodePar = list(cex = .007))

#Cluster Bio superficial
#Para poner como nombre la primer columna
bio_sup_pro <- data.frame(bio_sup_pro [,-1], row.names= bio_sup_pro$Sample)
d5 <- dist(bio_sup_pro)
cluster_bio2 <- hclust(d5, method = "complete")
plot(cluster_bio2)

library(dendextend)
dend5 <- as.dendrogram(cluster_bio2)
# order it the closest we can to the order of the observations:
dend3 <- rotate(dend3, 1:45)

# Color the branches based on the clusters:
dend5 <- color_branches(dend5, k=4) 
# We hang the dendrogram a bit:
dend5 <- hang.dendrogram(dend5,hang_height=0.1)
# reduce the size of the labels:
#dend <- assign_values_to_leaves_nodePar(dend, 0.5, "lab.cex")
dend5 <- set(dend5, "labels_cex", 0.5)

#Plot
plot(dend5, 
     main = "LMG1510", 
     horiz =  TRUE,  nodePar = list(cex = .007))

#Para comparar 2 dendrogramas
tanglegram(dend4, dend5, main_left = "Abundancia", main_right = "Biomasa")

#Cluster Abundancia relativa superficial grupos
#Para poner como nombre la primer columna
abun_rel_sup <- data.frame(abun_rel_sup [,-1], row.names= abun_rel_sup$Sample)
d6 <- dist(abun_rel_sup)
cluster_abun_rel <- hclust(d6, method = "complete")
plot(cluster_abun_rel)

library(dendextend)
dend6 <- as.dendrogram(cluster_abun_rel)
# order it the closest we can to the order of the observations:
dend6 <- rotate(dend6, 1:45)

# Color the branches based on the clusters:
dend6 <- color_branches(dend6, k=4) 
# We hang the dendrogram a bit:
dend6 <- hang.dendrogram(dend6,hang_height=0.1)
# reduce the size of the labels:
#dend <- assign_values_to_leaves_nodePar(dend, 0.5, "lab.cex")
dend6 <- set(dend6, "labels_cex", 0.5)

#Plot
plot(dend6, 
     main = "LMG1510", 
     horiz =  TRUE,  nodePar = list(cex = .007))

#Cluster Abundancia relativa superficial fitoplancton
#Para poner como nombre la primer columna
abun_rel_fito <- data.frame(abun_rel_fito [,-1], row.names= abun_rel_fito$Sample)
d7 <- dist(abun_rel_fito)
cluster_abun_rel_fito <- hclust(d7, method = "complete")
plot(cluster_abun_rel_fito)

library(dendextend)
dend7 <- as.dendrogram(cluster_abun_rel_fito)
# Color the branches based on the clusters:
dend7 <- color_branches(dend7, k=4) 
# We hang the dendrogram a bit:
dend7 <- hang.dendrogram(dend7,hang_height=0.1)
# reduce the size of the labels:
#dend <- assign_values_to_leaves_nodePar(dend, 0.5, "lab.cex")
dend7 <- set(dend7, "labels_cex", 0.7)

#Plot
plot(dend7, 
     main = "LMG1510 abundancias relativas", 
     horiz =  TRUE,  nodePar = list(cex = .007))
#Para comparar 2 dendrogramas
tanglegram(dend6, dend7, main_left = "Total", main_right = "Fito")

#Cluster Bio relativa superficial (solo fitoplancton)
#Para poner como nombre la primer columna
bio_rel_sup <- data.frame(bio_rel_sup [,-1], row.names= bio_rel_sup$Sample)
d8 <- dist(bio_rel_sup)
cluster_bio3 <- hclust(d8, method = "average")#usar method = "average" es UPGMA
plot(cluster_bio3)

library(dendextend)
dend8 <- as.dendrogram(cluster_bio3)
# Color the branches based on the clusters:
dend8 <- color_branches(dend8, k=4) 
# We hang the dendrogram a bit:
dend8 <- hang.dendrogram(dend8,hang_height=0.1)
# reduce the size of the labels:
#dend <- assign_values_to_leaves_nodePar(dend, 0.5, "lab.cex")
dend8 <- set(dend8, "labels_cex", 0.7)
dend8 <- highlight_branches_lwd(dend8, values =4)

#Plot
plot(dend8, 
     main = "LMG1510", 
     horiz =  TRUE,  nodePar = list(cex = .007))
summary(dend8)
dend8

#Cluster Bio superficial (solo fitoplancton)
#Para poner como nombre la primer columna
bio_sup <- data.frame(bio_sup [,-1], row.names= bio_sup$Sample)
d10 <- dist(bio_sup)
cluster_bio5 <- hclust(d10, method = "average")#usar method = "average" es UPGMA
plot(cluster_bio5)

library(dendextend)
dend10 <- as.dendrogram(cluster_bio5)
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
     horiz =  TRUE,  nodePar = list(cex = .007))

#Para comparar 2 dendrogramas
tanglegram(dend8, dend10, main_left = "Relativa", main_right = "Total")
#Cluster Bio total superficial quitando spp "raras" que solo aparecian en 1 st
#Para poner como nombre la primer columna
bio_complete <- data.frame(bio_com [,-1], row.names= bio_com$Sample)
d9 <- dist(bio_complete)
cluster_bio4 <- hclust(d9, method = "average")#usar method = "average" es UPGMA
plot(cluster_bio4)

library(dendextend)
dend8 <- as.dendrogram(cluster_bio4)
# Color the branches based on the clusters:
dend8 <- color_branches(dend8, k=4) 
# We hang the dendrogram a bit:
dend8 <- hang.dendrogram(dend8,hang_height=0.1)
# reduce the size of the labels:
#dend <- assign_values_to_leaves_nodePar(dend, 0.5, "lab.cex")
dend8 <- set(dend8, "labels_cex", 0.5)
#Plot
plot(dend8, 
     main = "LMG1510", 
     horiz =  TRUE,  nodePar = list(cex = .007))

#Cluster Bio superficial subgrupos
#Para poner como nombre la primer columna
bio_subg_rel <- data.frame(bio_subg_rel [,-1], row.names= bio_subg_rel$Muestra)
d11 <- dist(bio_subg_rel)
cluster_bio6 <- hclust(d11, method = "average")#usar method = "average" es UPGMA
plot(cluster_bio6)

library(dendextend)
dend10 <- as.dendrogram(cluster_bio6)
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
     horiz =  TRUE,  nodePar = list(cex = .007))

#Cluster LMG1510 con distancia Euclidea
matrixd <- dist(Complete, method = "euclidean")
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
dend10 <- rotate(dend10, c(26:33, 1))
#Plot
plot(dend10, 
     main = "Primavera y Otoño", 
     horiz =  TRUE,  nodePar = list(cex = .007), xlab = "Euclidean distance", axes = FALSE)
lines(x = c(0,0), y = c(0,100), type = "n") 
axis(side = 1, at = seq(0,100,10), labels = seq(100,0,-10))
plot(rotate(dend10, c(26:33, 1)))
