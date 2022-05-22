WwfLoad <-
  function(x){
    if ( x == ""){x <- getwd()}
    download.file("http://assets.worldwildlife.org/publications/15/files/original/official_teow.zip",
                  destfile = file.path(x, "wwf_ecoregions.zip"))#?1349272619")
    unzip(file.path(x, "wwf_ecoregions.zip"), exdir = file.path(x, "WWF_ecoregions"))
    file.remove(file.path(x, "wwf_ecoregions.zip"))
    #wwf <- sf::st_read(file.path(x, "WWF_ecoregions", "official", "wwf_terr_ecos.shp"))
    wwf2 <- sf::st_read(file.path(x, "WWF_ecoregions", "official", "wwf_terr_ecos.shp"))
    return(wwf2)
  }


wwf_eco2 = WwfLoad("")
wwf_eco2$BIOME = as.factor(wwf_eco2$BIOME)
wwf_eco2 = wwf_eco2%>% mutate(biome_num = as.numeric(BIOME))

a = ggplot()+
geom_sf(data = wwf_eco, aes(fill = BIOME), color=NA, alpha = 0.6 )
coords = c("N","S","E","W")
colors = c("forestgreen", "cornflowerblue","tomato1","darkorchid1")


names_biome = (c("Tropical and subtropical moist broadleaf forests","Tropical and subtropical moist broadleaf forests","	Tropical and subtropical coniferous forests","Temperate broadleaf and mixed forests","Temperate coniferous forests","Boreal forests/taiga","Tropical and subtropical grasslands, savannas, and shrublands","Temperate grasslands, savannas, and shrublands","Flooded grasslands and savannas","Montane grasslands and shrublands","Tundra","Mediterranean forests, woodlands, and shrub","Deserts and xeric shrublands","Mangrove")) %>% as.data.frame()
names_biome[,2] = seq.int(nrow(names_biome))
wwf_eco3 = merge(wwf_eco2,names_biome, by.x = "BIOME", by.y = "V2") 
wwf_eco3 = wwf_eco3 %>% arrange(desc(biome_num))
wwf_eco4 = wwf_eco3[22:24] %>% group_by(biome_num) %>% summarize(geometry = st_union(st_make_valid(geometry)))
wwf_eco4$BIOME = names_biome

for (i in 1:4){
  
  sp<- ggplot() + 
    geom_sf(data = wwf_eco3, aes(fill = .), alpha = 0.4, color = NA) +
  #(data=grid,      aes(x=long, y=lat, group=group, fill = count), alpha=1)    +
    #geom_path   (data=grid,      aes(x=long, y=lat, group=group), alpha=0.4, color="white") +
    
    #geom_sf(data = noncoastline_intersect, aes(size = count, colour = coord)) +
    geom_point(data = grid_point_f, aes(x=lon_deg, y=lat_deg,size = count, colour = coord ), ) +
    
    #geom_sf(data = noncoastline_intersect[noncoastline_intersect$coord == coords[i],], aes(size = count), colour = colors[i]) +
    #geom_point(data = grid_point_f[grid_point_f$coord ==coords[i],], aes(x=lon_deg, y=lat_deg,size = count), fill = colors[i], color = "black" ) +
    
    scale_size_continuous(range=c(1, 6)) +
    coord_sf(xlim = c(-25,160), ylim = c(-45,80)) 
  #coord_sf(xlim = c(-150,-10), ylim = c(-50,75)) 
  sp  <- sp + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                 panel.background = element_blank(), axis.line = element_line(colour = "black"))
  print(sp)
  ggsave(sp, file=paste0("Biomes",".png"))
  
}

a + scale_fill_manual(values = getPalette(colorCount))

##################
library(rmapshaper)
library(ggplot2)
library(sf)
library(dplyr)

wwf_eco4 = wwf_eco3[22:24] %>% group_by(biome_num) %>% summarize(geometry = st_combine(geometry))
wwf_simpl <- ms_simplify(wwf_eco4,keep = 0.0005, keep_shapes = T) %>% st_make_valid(st_union())

wwf_ecoreg <- wwf_eco3[c(5,7,24)] %>% group_by(ECO_NAME)%>% summarize(geometry = st_combine(geometry))
wwf_ecoreg.simple <- ms_simplify(wwf_ecoreg,keep = 0.0005, keep_shapes = T) %>% st_make_valid(st_union())

ggplot() + 
  geom_sf(wwf_ecoreg.simple, mapping = aes(fill = ECO_NAME, geometry = geometry),)


ggplot() + 
  geom_sf(desert, mapping = aes( geometry = geometry),) + 
  theme(legend.position = "none")
  ggsave(filename = "ecoregions_k=0.001.jpg",width = 12.6, height = 11.8)
  
tundra <- ms_simplify(wwf_eco4[11,], keep = 0.001, keep_shapes = T)
wwf_simpl$geometry[11] <- tundra$geometry

#when doing desert independently, we get a higher result of interaction with 
#the deserts and the other biomes... it definately acts weird
desert <- (wwf_eco4[13,])
wwf_simpl$geometry[13] <- desert$geometry

wwf_simpl <- ms_simplify(wwf_simpl)
