library(sf)
library(dplyr)
library(ggplot2)
library(reshape2)
library(gplots)
library(RColorBrewer)

load("C:/Users/davidle.WISMAIN//Box/lab folder/hotspots/maps raw data/biomes/wwf_simple_b2.RData")
load("C:/Users/davidle.WISMAIN//Box/lab folder/hotspots/paper/raw datat/4. quadrant freq/counts.coord.inland.norm.RData")
load("C:/Users/davidle.WISMAIN//Box/lab folder/hotspots/paper/raw datat/4. quadrant freq/local_G II.RData")
load("C:/Users/davidle.WISMAIN//Box/lab folder/hotspots/paper/raw datat/4. quadrant freq/global.inland.RData")
load("C:/Users/davidle.WISMAIN//Box/lab folder/hotspots/paper/raw datat/4. quadrant freq/counts.coord.inland.norm.RData")
load("C:/Users/davidle.WISMAIN/Box/lab folder/hotspots/maps raw data/biomes/names_biome.RData")
load("C:/Users/davidle.WISMAIN/Box/lab folder/hotspots/paper/raw datat/biome.int/interserctions.b2.RData")


int <- which(lengths(intersections.b$origins)!=1)
biom <- which(lengths(intersections.b$origins)==1)

internal <- intersections.b[int,]
outernal <- intersections.b[biom,]

#### preparing cells for analysis
Gi.cells <- Gi.prep.fun(local_G,global.inland)
cont.cells <- counts.prep.fun(counts.coord.inland2,global.inland)

#### intersection df and matrix
intrsct.df <- intrsct.fun(Gi.cells,internal = internal)

intrsct.mat <- intrsct.mat.fun(intrsct.df,internal)

#### calculating the number of edges at edge of biomes
perc.df <- calc.int.fun(intrsct.df,Gi.cells)

btsrp.mat <- btsrp.edges(intrsct.df,RE.df = cont.cells,reps = 100)

####calculate bootstrab for each biome intersect 

df2 <- data.frame()
for (i in 1:1005) {
  df <- data.frame()
  df <- fct.btstrp.biome.cont.a(cont.cells,internal)
  df2 <- rbind(df2,df)
}

df3 <- df2[complete.cases(df2),]
names(df3) <- c("ID1","ID2","sum")

#####create normalized matricies using the mean of the bootstrap

btsrp.mat2 <- fct.btstrp.biome.cont.b(intrsct.df,df3)
comb.mat <- intrsct.mat$sum.mat2
pval.sum <- btsrp.mat2$pval.sum
pval <- matrix(ncol = 14, nrow = 14)
pval[(which(pval.sum<0.01))] = "*"
diag(pval) <- NA

at <- gplots::heatmap.2(comb.mat, trace = 'none',
                        density.info = 'none', 
                        col = my_palette, distfun = dist_no_na,
                        scale = "none",symm = T,na.color = "grey",
                        cellnote = pval, notecol="black", notecex=2,
                        cexRow = 1, cexCol = 1)



################
mat.sum.norm <- intrsct.mat$sum.mat/btsrp.mat2$sum.mat.btsrp

cont.lst <- total.cont.fun(mat.sum.norm,0.05)
cont.df <- cont.lst$comb.df
tk <- cont.lst$tk
reord <- c(1,  7,  8,  9, 10, 11, 12, 13, 14,  2,  3,  4,  5,  6)
myPalette <- colorRampPalette(rev(brewer.pal(11, "Spectral")))
sd <- myPalette(14)[reord]
ggplot(cont.df, aes(reorder(biome,value),value)) + 
  geom_boxplot(aes(fill = biome),alpha = 0.5) + 
  #scale_fill_manual(values=c("bisque2","bisque2","bisque2","darkgoldenrod1","darkgoldenrod1","darkgoldenrod1","bisque2","orange3","cornsilk","darkgoldenrod1","orange4","bisque2","orange3","cornsilk"))+
  scale_fill_manual(values=sd) +
  theme_test()+
  geom_text(data = tk,  aes(label = cld, x = biome, y = quant),
            vjust = -0.35, hjust = -0.2, size = 3) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        axis.title.x = element_blank(),
        legend.position = "none") +
  labs( y = "log(normalized contribution)") 


  
