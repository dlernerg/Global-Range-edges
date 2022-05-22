#find clusters using dbscan
library(concaveman)
library(dbscan)
library(parallel)
library(sf)

#find what epsilon to use
i=8
spec_to_pull = as.character(count_species[i,1])
GBIFtree_dat = data.frame(GBIFtree_f5) %>% subset(species == spec_to_pull)
GBIFtree_dat = GBIFtree_dat[,2:3]

kNNdistplot(centr_matrix, k=3)
abline(h=4, col = "black", lty=2)

#a = kNN(GBIFtree_dat, k=5)
#plot(a$dist)

#carry out dbscan
#i=12

numCores <- detectCores()
cl <- makeCluster(numCores-1)

DBSCAN_clust = list()
for (i in 15:nrow(count_species_f5))
{
  spec_to_pull = as.character(count_species_f5[i,1])
  GBIFtree_dat = data.frame(GBIFtree_f5) %>% subset(species == spec_to_pull)
  GBIFtree_dat = GBIFtree_dat[,22:23]
  
  if (nrow(GBIFtree_dat)>200000){
    e = 2  
  } else if (nrow(GBIFtree_dat) <=200000 & nrow(GBIFtree_dat) >7200){
    e = 3.8
  } else if (nrow(GBIFtree_dat) <=7200){
  e = 8.41
  } 
  
  res = dbscan(GBIFtree_dat, eps = e, minPts = 2)
  
  DBSCAN_clust[[i]] = res
  
  }
stopCluster(cl)


###get the INliers using the function remove_outliers
numCores <- detectCores()
cl <- makeCluster(numCores-1)
a = data.frame()
IN = list()
for (i in 15:nrow(count_species_f5))
{
  GBIFtree_dat = GBIFpull(GBIFtree_f5, count_species_f5, i)
  clus = DBSCAN_clust[[i]]$cluster
  
  #determine the outliers from DBSCAN that will return to the function
  IN[[i]] = remove_outliers(GBIFtree_dat,clus)
  #GBIFclean = rbind(all,IN)
  
}
stopCluster(cl)


#####concave hulls
#polygons = list()

cl <- makeCluster(numCores)
polygon_dat = data.frame()
polygons = data.frame()
for (i in 15:nrow(count_species_f5))
{

    GBIFtree_dat = GBIFpull(GBIFtree_f5, count_species_f5, i)
    clus = (DBSCAN_clust[[i]]$cluster)
    maxn = max(clus)
    noutlier = GBIFtree_dat[which(clus!=0),]
    clus_clean = clus[clus!=0]
    all = cbind(noutlier,clus_clean)
    
    if (is.na(IN[[i]])){
      GBIFclean = all
    }else{
      colnames(IN[[i]]) = colnames(all)
      GBIFclean = rbind(all,IN[[i]])
      
    }
  
    polygon_dat = data.frame()
    
    for (j in 1:(maxn))
    {
      pts <- st_as_sf(GBIFclean[(GBIFclean$clus_clean==j),1:2], coords=c('decimalLongitude','decimalLatitude'), crs= "WGS84" )
      #polygon_dat = (concaveman(pts, concavity = 2.5))
      polygon_dat[j,1] = (concaveman(pts, concavity = 2.5))
      polygon_dat[j,2] = i
    } 
    polygons = rbind(polygons,polygon_dat)
}
stopCluster(cl)

polygons = st_sf(polygons, crs="WGS84")
polygons$specie = count_species_f5[polygons$V2,1]

#remove the polygons that are from an unspecified specie 
a = which(polygons$specie == "")
polygons = polygons[-a,]

save(polygons, file = "polygons.RData")


######plot polygon
i=20
rand = sample(1:nrow(count_species_f5), 20, replace = F)
rand = round(rand)


for (i in 1:length(rand))
{
  
  j = rand[i]  
  
#polygon_dat = polygons[(polygons$V2==i),1]
#a = which(GBIFpolygon_clipped2 == flag_by_specie.dpl[30,1])
  polygon_dat = GBIFpolygon_clipped$tree2u[j]
  polygon_dat = st_sf(st_sfc(polygon_dat))
  
#b = which(count_species_f5 == flag_by_specie.dpl[91,1])
  
 # GBIFtree_dat = GBIFpull(GBIFtree_f5, count_species_f5, b)
#  clus = (DBSCAN_clust[[b]]$cluster)
#  maxn = max(clus)
 # noutlier = GBIFtree_dat[which(clus!=0),]
  #clus_clean = clus[clus!=0]
#  all = cbind(noutlier,clus_clean)
  
  
nlrg = 15
  lims = st_bbox(polygon_dat)
  sp.map = ggplot() +
    geom_polygon(data=countries, aes(x=long, y=lat, group=group), fill=NA, color="black") +
    geom_sf(data = polygon_dat,fill = "blue") +
    
    coord_sf(xlim = c(lims[1]-nlrg,lims[3]+nlrg), ylim = c(lims[2]-nlrg,lims[4]+nlrg)) +
    #scale_fill_gradient(low="yellow", high="blue") +
   ggtitle(paste(GBIFpolygon_clipped2[j,1]))
  
  plot(sp.map)
}

########plot DBSCAN and RO function
rand = sample(1:nrow(count_species_f5), 20, replace = F)
rand = round(rand)

  
  for (i in 1:length(rand))
{
    
  j = rand[i]
GBIFtree_dat = GBIFpull(GBIFtree_f5, count_species_f5, a)
clus = DBSCAN_clust[[a]]$cluster
INa = remove_outliers(GBIFtree_dat,clus)
all = cbind(GBIFtree_dat,clus)
all$clus = as.factor(all$clus)

if (is.na(INa))
{
  sp.map = ggplot() +
    geom_polygon(data=countries, aes(x=long, y=lat, group=group), fill=NA, color="black") +
    #geom_sf(data = polygon_dat) +
    geom_point(data = all, aes(y=decimalLatitude,x=decimalLongitude, color = clus)) +
    # geom_point(data = IN, aes(y=decimalLatitude,x=decimalLongitude), color = "black") +
    
    #geom_point(data = GBIFtree_dat, aes(y=decimalLatitude,x=decimalLongitude), color = "red") +
    #geom_point(data = gbif_filtered, aes(y=decimalLatitude,x=decimalLongitude), color = "red") +
    
    #geom_raster(data = kde, aes(x = x, y=y, alpha = layer)) +
    #  coord_sf(xlim = c(lims[1]-nlrg,lims[3]+nlrg), ylim = c(lims[2]-nlrg,lims[4]+nlrg)) +
    #scale_fill_gradient(low="yellow", high="blue") +
    ggtitle(paste(count_species_f5[a,1]))
  
}else{
  
  sp.map = ggplot() +
    geom_polygon(data=countries, aes(x=long, y=lat, group=group), fill=NA, color="black") +
    #geom_sf(data = polygon_dat) +
    geom_point(data = all, aes(y=decimalLatitude,x=decimalLongitude, color = clus)) +
     geom_point(data = INa, aes(y=decimalLatitude,x=decimalLongitude), color = "black") +
    
    #geom_point(data = gbif_filtered, aes(y=decimalLatitude,x=decimalLongitude), color = "red") +
    
    #geom_raster(data = kde, aes(x = x, y=y, alpha = layer)) +
    #  coord_sf(xlim = c(lims[1]-nlrg,lims[3]+nlrg), ylim = c(lims[2]-nlrg,lims[4]+nlrg)) +
    #scale_fill_gradient(low="yellow", high="blue") +
    ggtitle(paste(count_species_f5[j,1]))
  

}
plot(sp.map)



  }

