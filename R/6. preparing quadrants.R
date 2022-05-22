library(dggridR)
library(ggplot2)
library(sf)

#############Get polygons for each country of the world
#Construct a global grid 

dggs          <- dgconstruct(aperture = 3, res = 7, metric=FALSE, show_info = TRUE, resround = "down")

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

#############Get the corresponding grid cells for the different cardinal coordinates

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
  #use one of the margins (cardinal)
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

