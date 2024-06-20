#NMDS English version
#Se necesitan 2 set de datos, uno de spp y otro de "grupos"
#Para poner la primer columna como titulo
bio <- data.frame(bio_table[,-1], row.names= bio_table$Sample)
#bio1 <- data.frame(data_table1[,-1], row.names= data_table1$Sample)
ID <- data.frame(ID_table[,-1], row.names= ID_table$Sample)

library(vegan) # for vegetation and community analysis
library(dplyr) # for data manipulation

#Conviene usar datos relativos!! Ver Rpubs
#Para transformar las variables (no siempre necesario)
#bio1 <- decostand(bio, "normalize")
# Calcular la matriz de distancia
bio_distmat1 <- vegdist(bio, binary=FALSE, method = "bray")
#bio_distmat1 <- vegdist(bio, binary=FALSE, method = "jaccard")
# Correr el NMDS con vegan (metaMDS)
bio_NMS1 <-  metaMDS(bio_distmat1,
                     distance = "bray",
                     k = 2,
                     maxit = 999, 
                     trymax = 500,
                     wascores= TRUE, expand = TRUE, autotransform = FALSE)

#bio_NMS1 <-  monoMDS(bio_distmat1, model= "hybrid")
plot(bio_NMS1)
bio_NMS1

#Graficos basicos
ordiplot(bio_NMS1, type = "n", main = "ellipses")
orditorp(bio_NMS1, display = "sites", labels = F, pch = c(16, 8, 1, 2, 4, 6, 7, 9, 10, 11, 12, 13, 14, 15, 17, 18, 19, 20) [as.factor(ID_table$`Mes+A?o+Zona`)],
         col = c("green", "blue") [as.factor(ID_table$Group)], cex = 1)
#Con Elipses
ordiellipse(bio_NMS1, groups = ID_table$Group, draw = "polygon", lty = 3, col = "grey90")
#Con poligonos
ordihull(bio_NMS1, groups = ID_table$Group, draw = "polygon", lty = 1, col = "grey90")

# Shepards test/goodness of fit (para comprobar que tan bien las graficas de ordenacion representan los datos)
goodness(bio_NMS1) # Produces a results of test statistics for goodness of fit for each point
stressplot(bio_NMS1) # Produces a Shepards diagram

#Para calcular los spp scores
bio.spp_scrs1 <-   sppscores(bio_NMS1) <- bio

#Para agregar las especies significativas al dataframe
bio.spp.fit <- envfit(bio_NMS1, bio, permutations = 999) # this fits species vectors
spp.scrs <- as.data.frame(scores(bio.spp.fit, display = "vectors", arrow.mul=2.5)) #save species intrinsic values into dataframe #arrow.mul es un vector que multiplica el largo de las flechas
spp.scrs <- cbind(spp.scrs, Species = rownames(spp.scrs)) #add species names to dataframe
spp.scrs <- cbind(spp.scrs, pval = bio.spp.fit$vectors$pvals) #add pvalues to dataframe so you can select species which are significant
#spp.scrs<- cbind(spp.scrs, abrev = abbreviate(spp.scrs$Species, minlength = 6)) #abbreviate species names
sig.spp.scrs <- subset(spp.scrs, pval<=0.05) #subset data to show species significant at 0.05

#extract NMDS scores (x and y coordinates)
data.scores1 <- as.data.frame(scores(bio_NMS1$points))

#Factorizar los anos
Ano_mes<-as.factor(ID_table$`Mes+Año`)
Ano_mes_zona<-as.factor(ID_table$`Mes+Año+Zona`)
Group<-as.factor(ID_table$Group)
Season<-as.factor(ID_table$Subgroup)


#Hacemos un nuevo data frame con las coordenadas x e y
data.scores1$Sample = ID_table$Sample
data.scores1$Year = Ano_mes
data.scores1$Month = ID_table$Month
data.scores1$Area = ID_table$Area
data.scores1$'Year+Z' = Ano_mes_zona
data.scores1$Group = Group
data.scores1$Season = Season

#Para graficar!
library(ggplot2)


xx = ggplot(data.scores1, aes(x = MDS1, y = MDS2)) + 
  geom_point(size = 5, aes( shape = Season, colour = Group))+ 
  theme(axis.text.y = element_text(colour = "black", size = 12, face = "bold"), 
        axis.text.x = element_text(colour = "black", face = "bold", size = 12), 
        legend.text = element_text(size = 12, face ="bold", colour ="black"), 
        legend.position = "right", axis.title.y = element_text(face = "bold", size = 14), 
        axis.title.x = element_text(face = "bold", size = 14, colour = "black"), 
        legend.title = element_text(size = 14, colour = "black", face = "bold"), 
        panel.background = element_blank(), panel.border = element_rect(colour = "black", fill = NA, size = 1.2),
        legend.key=element_blank()) + 
  labs(x = "NMDS1", colour = "Group", y = "NMDS2", shape = "Time frame")  + 
  scale_colour_manual(values = c("lightseagreen", "indianred3")) +
  scale_shape_manual(values=seq(0,9))

xx


#Para agregar las especies significativas 
#options(ggrepel.max.overlaps = Inf) para quitar el overlap
xx+
  geom_segment(data =sig.spp.scrs, aes(x = 0, xend=NMDS1, y=0, yend=NMDS2), arrow = arrow(length = unit(0.25, "cm")), colour = "grey10", lwd=0.3) + #add vector arrows of significant species
  ggrepel::geom_text_repel(data = sig.spp.scrs, aes(x=NMDS1, y=NMDS2, label = Species), cex = 3, direction = "both")+ #add labels for species, use ggrepel::geom_text_repel so that labels do not overlap
  labs(title = "Ordination with species vectors")

#Para ver las spp que son significativas en el analisis
bio.spp.fit1 <- envfit(bio_NMS1, bio_table, permutations = 999)
head(bio.spp.fit1)

#ANOSIM (A Significance value less than 0.05 is generally considered to be statistically significant)
ano <- anosim(bio_distmat1, data.scores1$Group, distance = "bray", permutations = 9999)
ano
plot(ano)

ano1 <- anosim(bio_distmat1, data.scores1$Season, distance = "bray", permutations = 9999)
ano1
plot(ano1)
