#Description: Download the data for biomes from the WWF webpage and simplify the maps to reduced complexity of edges using ms_simply function

library(rmapshaper)
library(ggplot2)
library(sf)
library(dplyr)

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

names_biome = (c("Tropical and subtropical moist broadleaf forests","Tropical and subtropical moist broadleaf forests","	Tropical and subtropical coniferous forests","Temperate broadleaf and mixed forests","Temperate coniferous forests","Boreal forests/taiga","Tropical and subtropical grasslands, savannas, and shrublands","Temperate grasslands, savannas, and shrublands","Flooded grasslands and savannas","Montane grasslands and shrublands","Tundra","Mediterranean forests, woodlands, and shrub","Deserts and xeric shrublands","Mangrove")) %>% as.data.frame()
names_biome[,2] = seq.int(nrow(names_biome))
wwf_eco3 = merge(wwf_eco2,names_biome, by.x = "BIOME", by.y = "V2") 
wwf_eco3 = wwf_eco3 %>% arrange(desc(biome_num))
wwf_eco4 = wwf_eco3[22:24] %>% group_by(biome_num) %>% summarize(geometry = st_union(st_make_valid(geometry)))
wwf_eco4$BIOME = names_biome

wwf_eco4 = wwf_eco3[22:24] %>% group_by(biome_num) %>% summarize(geometry = st_combine(geometry))
wwf_simpl <- ms_simplify(wwf_eco4,keep = 0.0005, keep_shapes = T) %>% st_make_valid(st_union())

wwf_ecoreg <- wwf_eco3[c(5,7,24)] %>% group_by(ECO_NAME)%>% summarize(geometry = st_combine(geometry))
wwf_ecoreg.simple <- ms_simplify(wwf_ecoreg,keep = 0.0005, keep_shapes = T) %>% st_make_valid(st_union())

