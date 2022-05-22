library(dggridR)
library(ggplot2)
library(sf)
#Get polygons for each country of the world
#Construct a global grid with cells approximately 1000 miles across

dggs          <- dgconstruct(aperture = 3, res = 7, metric=FALSE, show_info = TRUE, resround = "down")
#Old world
#a <- dgrectgrid(dggs, minlat = -50,minlon = -20,maxlat = 75,maxlon = 175, frame = FALSE)
#New world

#global
a <- dgrectgrid(dggs, minlat = -60,minlon = -180,maxlat = 75,maxlon = 179.9, frame = FALSE)

global.cell <- data.frame(cell=getSpPPolygonsIDSlots(a), row.names=getSpPPolygonsIDSlots(a))
global <- SpatialPolygonsDataFrame(a, global.cell)
 
 for(i in 1:length(global@polygons))
   {if(max(global@polygons[[i]]@Polygons[[1]]@coords[,1]) -  min(global@polygons[[i]]@Polygons[[1]]@coords[,1]) > 180) {
         global@polygons[[i]]@Polygons[[1]]@coords[,1] <- (global@polygons[[i]]@Polygons[[1]]@coords[,1] +360) %% 360}
      }
  
global <- st_as_sf(global)

global$cell <- 1:nrow(global)

#Get the corresponding grid cells for all the x,y coordinates in the tree (lat-long pair)
marked_cells_in.all <- data.frame()
specie_in <- data.frame()
marked_cells_peri.all <- data.frame()
specie_peri <- data.frame()

for (i in 1:nrow(GBIFpolygon_groups))
{
  GBIFchose = GBIFpolygon_groups[i,] %>% st_buffer(0.01)
  
  marked_cells_in <- as.data.frame(st_intersects(GBIFchose,global))
  marked_cells_in <- as.data.frame(marked_cells_in$col.id)
  
  marked_cells_in.all <- rbind(marked_cells_in.all,marked_cells_in) 
  
  marked_cells_in <- list(as.numeric(marked_cells_in$`marked_cells_in$col.id`))  
  specie_in[i,1] <- paste(GBIFchose$specie[1],GBIFchose$group[1])
  specie_in[i,2] <-  list(marked_cells_in)  
  
  
  GBIFperi = GBIFchose %>% st_cast("POLYGON") %>% st_cast("LINESTRING")
  
  marked_cells_peri  <- as.data.frame(st_intersects(GBIFperi,global))
  marked_cells_peri <- as.data.frame(marked_cells_peri$col.id)
  
  
  marked_cells_peri.all <- rbind(marked_cells_peri.all,marked_cells_peri) 
  
  marked_cells_peri <- list(as.numeric(marked_cells_peri$`marked_cells_peri$col.id`))  
  specie_peri[i,1] <- paste(GBIFchose$specie[1],GBIFchose$group[1])
  specie_peri[i,2] <-  list(marked_cells_peri)  
  
}

specie_in2 <- data.frame()
marked_cells_in2.all <- data.frame()

#identify those grids that are internal and not intersecting between internal and peripheral
for (i in 7622:nrow(specie_in)){
  ins <- as.data.frame(specie_in$V2[[i]])
  peri<- as.data.frame(specie_peri$V2[[i]])
  
  names(ins) <- "V1"
  names(peri) <- "V1"
  
  marked_cells_in.2 <- subset(ins, !(V1 %in% peri$V1))
  
  marked_cells_in2.all <- rbind(marked_cells_in2.all,marked_cells_in.2) 
  
  marked_cells_in.2 <- list(as.numeric(marked_cells_in.2$V1))  
  specie_in2[i,1] <- specie_in[i,1]
  specie_in2[i,2] <- list(marked_cells_in.2) 
}


counts.in <- as.data.frame(1:nrow(global))
counts.peri <- as.data.frame(1:nrow(global))
counts.in2 <- as.data.frame(1:nrow(global))

for (i in 1:nrow(global)){
  counts.in[i,2] <- length(which(marked_cells_in.all == counts.in[i,1]))
  counts.peri[i,2] <- length(which(marked_cells_peri.all == counts.peri[i,1]))
  counts.in2[i,2] <- length(which(marked_cells_in2.all== counts.in2[i,1]))
}

#counts.in <- counts.in[counts.in$V2>0,]
#counts.peri <- counts.peri[counts.peri$V2>0,]
#countss.in2 <- counts.in2[counts.in2$V2>0,]
save(specie_in, file = "specie_in.RData")
save(specie_in2, file = "specie_in2.RData")
save(specie_peri, file = "specie_peri.RData")
save(marked_cells_in.all, file = "marked_cells_in.all.RData")
save(marked_cells_in2.all, file = "marked_cells_in2.all.RData")
save(marked_cells_peri.all, file = "marked_cells_peri.all.RData")

save(counts.in, file = "counts.in.RData")
save(counts.in2, file = "counts.in2.RData")
save(counts.peri, file = "counts.peri.RData")

margin_allN$coord = "N"
margin_allS$coord = "S"
margin_allE$coord = "E"
margin_allW$coord = "W"
margin_all <- rbind(margin_allN,margin_allS,margin_allE,margin_allW)
coords = c("N","S","E","W")

marked_cells_coord.all <- data.frame()
specie_coord <- data.frame()

for (i in (1:4))
{
  #use one of the margins (cardinal_)
  margin_chose = as.data.frame(margin_all[margin_all$coord == coords[i],]) #chose the specific margin
  margin_chose <- st_as_sf(margin_chose, coords = c("V1","V2"))
  st_crs(margin_chose) <- st_crs(global)
  
  #which polygons intersect
  marked_cells_coord  <- as.data.frame(st_intersects(margin_chose,global))
  #which polygons intersect
  marked_cells_coord2 <- as.data.frame(marked_cells_coord$col.id)
  marked_cells_coord2[,2] <- coords[i]
  
  marked_cells_coord.all <- rbind(marked_cells_coord.all,marked_cells_coord2)
  
  #which specie intersects where
  specie_c <- as.data.frame(margin_chose$species[marked_cells_coord$row.id])
  specie_c[,2] <- margin_chose$ID[marked_cells_coord$row.id]
  specie_c[,3] <- coords[i]
  specie_c[,4] <- marked_cells_coord$col.id
  

  specie_coord <- rbind(specie_coord,specie_c)  
}



counts.coord<- as.data.frame(1:nrow(global))
for (j in 1:4){
  for (i in 1:nrow(global)){
    
    counts.coord[i,j+1] <- length(which(marked_cells_coord.all[marked_cells_coord.all$V2 == coords[j],] == counts.coord[i,1]))
    
  }
}

save(specie_coord, file = "specie_coord.RData")
save(marked_cells_coord.all, file = "marked_cells_coord.all.RData")
save(counts.coord, file = "counts.coord.RData")



away = c(which(margin_all$coord == "N" & margin_all$V2>0),which(margin_all$coord == "S" & margin_all$V2<0))
towards = c(which(margin_all$coord == "N" & margin_all$V2<0),which(margin_all$coord == "S" & margin_all$V2>0))
margin_all$direction = "t"
margin_all$direction[away] = "a"
margin_all$direction[towards] = "t"
direc <- c("a","t")

marked_cells_dir.all <- data.frame()
specie_dir <- data.frame()

for (i in (1:2))
{

  margin_chose = as.data.frame(margin_all[margin_all$direction == direc[i],])
  margin_chose <- margin_chose[complete.cases(margin_chose),]
  margin_chose <- st_as_sf(margin_chose, coords = c("V1","V2"))
  st_crs(margin_chose) <- st_crs(global)
  
  #which polygons intersect
  marked_cells_dir  <- as.data.frame(st_intersects(margin_chose,global))
  #which polygons intersect
  marked_cells_dir2 <- as.data.frame(marked_cells_dir$col.id)
  marked_cells_dir2[,2] <- direc[i]
  
  marked_cells_dir.all <- rbind(marked_cells_dir.all,marked_cells_dir2)
  
  #which specie intersects where
  specie_d <- as.data.frame(margin_chose$species[marked_cells_dir$row.id])
  specie_d[,2] <- margin_chose$ID[marked_cells_dir$row.id]
  specie_d[,3] <- direc[i]
  specie_d[,4] <- marked_cells_dir$col.id
  
  
  specie_dir <- rbind(specie_dir,specie_d)  
}


counts.dir<- as.data.frame(1:nrow(global))
for (j in 1:2){
  for (i in 1:nrow(global)){
    
    counts.dir[i,j+1] <- length(which(marked_cells_dir.all[marked_cells_dir.all$V2 == direc[j],] == counts.dir[i,1]))
    
  }
}

names(counts.dir) <- c("count","away","towards")


###plot 

names(counts.coord)[1] <- "cell"
global1 <- merge(global, counts.coord, by= "cell")

ggplot() + geom_polygon(data=countries, aes(long, lat, group=group), fill = NA, col = "black") + geom_sf(data = global1[global1$countss>0,], aes(fill = log(countss)))

                                                                                                         