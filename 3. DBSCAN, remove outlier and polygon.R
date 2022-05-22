#find clusters using dbscan
library(concaveman)
library(dbscan)
library(parallel)
library(sf)


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


}

