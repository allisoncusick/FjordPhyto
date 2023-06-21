#NMDS Martina Mascioni English version
#to edit final plots use INKSCAPE (downloaded to my mac already)

#Se necesitan 2 set de datos, uno de spp y otro de "grupos"
#first column in each file has to be ID sample and has to be exactly the same order sequence for each file
#at some point I want to categorize sites into categories "northern" "middle" southern" 
#this code is for THE ENTIRE dataset, so make smaller subset dataframes to do this again
#this code worked for my entire dataset so thats great! 

#Para poner la primer columna como titulo
bio <- data.frame(sumdf_wider[,-1], row.names= sumdf_wider$samples)
#create the ID file from my working_df
#I want columns: samples (which exactly have to match the bio dataframe), date, month, Site_Name
colnames(working_df)

ID <- working_df %>% 
  distinct(samples, Date, Site_Name, month,season)
    
ID2 <- data.frame(ID[,-1], row.names= ID$samples)

library(vegan) # for vegetation and community analysis
library(dplyr) # for data manipulation

#Conviene usar datos relativos!! Ver Rpubs
#Para transformar las variables (no siempre necesario)
#bio1 <- decostand(bio, "normalize")
# Calcular la matriz de distancia
bio_distmat1 <- vegdist(bio, binary=FALSE, method = "bray")
bio_distmatjaccard <- vegdist(bio, binary=FALSE, method = "jaccard")
# Correr el NMDS con vegan (metaMDS)
bio_NMS1 <-  metaMDS(bio_distmat1,
                     distance = "bray",
                     k = 2,
                     maxit = 999, 
                     trymax = 500,
                     wascores= TRUE, expand = TRUE, autotransform = FALSE)

#bio_NMS1 <-  monoMDS(bio_distmat1, model= "hybrid")
plot(bio_NMS1)
#error said it cant create species scores, dont worry about that now, will calculate later.
bio_NMS1

#Graficos basicos
ordiplot(bio_NMS1, type = "n", main = "ellipses")
ordi
orditorp(bio_NMS1, display = "sites", labels = F, pch = c(16, 8, 1, 2, 3, 4, 5, 6, 7, 9, 10, 11, 12, 13, 14, 15, 17, 18, 19, 20, 21, 22, 23, 24,25) [as.factor(ID2$Site_Name)],
         col = c("green", "blue", "red") [as.factor(ID2$season)], cex = 1)

#Con Elipses
ordiellipse(bio_NMS1, groups = ID2$Site_Name, draw = "polygon", lty = 3, col = "grey90")

#Con poligonos
ordihull(bio_NMS1, groups = ID2$Site_Name, draw = "polygon", lty = 1, col = "grey90")

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
spp.scrs<- cbind(spp.scrs, abrev = abbreviate(spp.scrs$Species, minlength = 6)) #abbreviate species names
sig.spp.scrs <- subset(spp.scrs, pval<=0.05) #subset data to show species significant at 0.05

#extract NMDS scores (x and y coordinates)
data.scores1 <- as.data.frame(scores(bio_NMS1$points))

#Factorizar los anos
month<-as.factor(ID2$month)
site<-as.factor(ID2$Site_Name)
season<-as.factor(ID2$season)


#Hacemos un nuevo data frame con las coordenadas x e y
data.scores1$Sample = ID2$samples
data.scores1$Month = ID2$month
data.scores1$Site = ID2$Site_Name
#data.scores1$'Year+Z' = Ano_mes_zona
data.scores1$Season = ID2$season

#Para graficar!
library(ggplot2)


colorsMDS <- c("#519935",
               "#bc4fac",
               "#5acd7c",
               "#543586",
               "#acb839",
               "#6d71d8",
               "#ce9c2d",
               "#5e8bd5",
               "#c76d27",
               "#43c8ac",
               "#d14c84",
               "#59b275",
               "#862a63",
               "#9ac666",
               "#bc81d5",
               "#98a03c",
               "#d97cb8",
               "#447228",
               "#d64c55",
               "#cfa14c",
               "#b94b65",
               "#caa662",
               "#a4473a",
               "#856521",
               "#c75e3c")

#plot by site and month
xx = ggplot(data.scores1, aes(x = MDS1, y = MDS2)) + 
  geom_point(size = 5, aes( shape = Site, colour = month))+ 
  theme(axis.text.y = element_text(colour = "black", size = 12, face = "bold"), 
        axis.text.x = element_text(colour = "black", face = "bold", size = 12), 
        legend.text = element_text(size = 12, face ="bold", colour ="black"), 
        legend.position = "right", axis.title.y = element_text(face = "bold", size = 14), 
        axis.title.x = element_text(face = "bold", size = 14, colour = "black"), 
        legend.title = element_text(size = 14, colour = "black", face = "bold"), 
        panel.background = element_blank(), panel.border = element_rect(colour = "black", fill = NA, size = 1.2),
        legend.key=element_blank()) + 
  labs(x = "NMDS1", colour = "month", y = "NMDS2", shape = "Site")  + 
  scale_colour_manual(values = colorsMDS) +
  scale_shape_manual(values=seq(0,24))

xx

#Para agregar las especies significativas 
#options(ggrepel.max.overlaps = Inf) para quitar el overlap
xx +
  geom_segment(data = sig.spp.scrs, aes(x = 0, xend=NMDS1, y=0, yend=NMDS2), arrow = arrow(length = unit(0.25, "cm")), colour = "grey10", lwd=0.3) + #add vector arrows of significant species
  ggrepel::geom_text_repel(data = sig.spp.scrs, aes(x=NMDS1, y=NMDS2, label = Species), cex = 3, direction = "both")+ #add labels for species, use ggrepel::geom_text_repel so that labels do not overlap
  labs(title = "Ordination with phytogroup vectors")

#plot by site and season
yy = ggplot(data.scores1, aes(x = MDS1, y = MDS2)) + 
  geom_point(size = 5, aes( shape = Site, colour = Season))+ 
  theme(axis.text.y = element_text(colour = "black", size = 12, face = "bold"), 
        axis.text.x = element_text(colour = "black", face = "bold", size = 12), 
        legend.text = element_text(size = 12, face ="bold", colour ="black"), 
        legend.position = "right", axis.title.y = element_text(face = "bold", size = 14), 
        axis.title.x = element_text(face = "bold", size = 14, colour = "black"), 
        legend.title = element_text(size = 14, colour = "black", face = "bold"), 
        panel.background = element_blank(), panel.border = element_rect(colour = "black", fill = NA, size = 1.2),
        legend.key=element_blank()) + 
  labs(x = "NMDS1", colour = "Season", y = "NMDS2", shape = "Site")  + 
  scale_colour_manual(values = colorsMDS) +
  scale_shape_manual(values=seq(0,24))

yy

#Para agregar las especies significativas 
#options(ggrepel.max.overlaps = Inf) para quitar el overlap
yy +
  geom_segment(data = sig.spp.scrs, aes(x = 0, xend=NMDS1, y=0, yend=NMDS2), arrow = arrow(length = unit(0.25, "cm")), colour = "grey10", lwd=0.3) + #add vector arrows of significant species
  ggrepel::geom_text_repel(data = sig.spp.scrs, aes(x=NMDS1, y=NMDS2, label = Species), cex = 3, direction = "both")+ #add labels for species, use ggrepel::geom_text_repel so that labels do not overlap
  labs(title = "Ordination with phytogroup vectors")

#plot by month and season
zz = ggplot(data.scores1, aes(x = MDS1, y = MDS2)) + 
  geom_point(size = 5, aes( shape = month, colour = Season))+ 
  theme(axis.text.y = element_text(colour = "black", size = 12, face = "bold"), 
        axis.text.x = element_text(colour = "black", face = "bold", size = 12), 
        legend.text = element_text(size = 12, face ="bold", colour ="black"), 
        legend.position = "right", axis.title.y = element_text(face = "bold", size = 14), 
        axis.title.x = element_text(face = "bold", size = 14, colour = "black"), 
        legend.title = element_text(size = 14, colour = "black", face = "bold"), 
        panel.background = element_blank(), panel.border = element_rect(colour = "black", fill = NA, size = 1.2),
        legend.key=element_blank()) + 
  labs(x = "NMDS1", colour = "Season", y = "NMDS2", shape = "month")  + 
  scale_colour_manual(values = colorsMDS) +
  scale_shape_manual(values=seq(0,24))

zz

#Para agregar las especies significativas 
#options(ggrepel.max.overlaps = Inf) para quitar el overlap
zz +
  geom_segment(data = sig.spp.scrs, aes(x = 0, xend=NMDS1, y=0, yend=NMDS2), arrow = arrow(length = unit(0.25, "cm")), colour = "grey10", lwd=0.3) + #add vector arrows of significant species
  ggrepel::geom_text_repel(data = sig.spp.scrs, aes(x=NMDS1, y=NMDS2, label = Species), cex = 3, direction = "both")+ #add labels for species, use ggrepel::geom_text_repel so that labels do not overlap
  labs(title = "Ordination with phytogroup vectors")

#Para ver las spp que son significativas en el analisis
bio.spp.fit1 <- envfit(bio_NMS1, bio, permutations = 999)
head(bio.spp.fit1)

#ANOSIM (A Significance value less than 0.05 is generally considered to be statistically significant)
ano <- anosim(bio_distmat1, data.scores1$Month, distance = "bray", permutations = 9999)
ano
plot(ano) #takes a while to plot, patience. 

ano1 <- anosim(bio_distmat1, data.scores1$Season, distance = "bray", permutations = 9999)
ano1
plot(ano1)

ano2 <- anosim(bio_distmat1, data.scores1$Site, distance = "bray", permutations = 9999)
ano2
plot(ano2)

#not clustering at all, so split dataframes into subsets and try all code again

