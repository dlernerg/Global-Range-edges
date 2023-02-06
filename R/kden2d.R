kden.highr <- list()
for (i in 448:nrow(count_species_f5)){

  spec_to_pull = as.character(count_species_f5[i,1])
  GBIFtree_dat = data.frame(GBIFtree_all) %>% subset(species == spec_to_pull)
  GBIFtree_dat = GBIFtree_dat[,22:23]
  
  kden1 <- MASS::kde2d(GBIFtree_dat$decimalLongitude, GBIFtree_dat$decimalLatitude, n = 800, lims = c(range(GBIFtree_dat$decimalLongitude), range(GBIFtree_dat$decimalLatitude)))
  kden.highr[i] <- (list(kden1))
  
}



# -------------------------------------------------------------------------

kden.z <- reshape2::melt(kden[[i]]$z)  

kden.all.f <- kden.z %>%
  filter(!is.na(value)) %>%
  filter(value> quantile(.$value,0.96))%>%
  tbl_df() 

ggplot(kden.all.f, aes(x = Var1, y = Var2, z = value, fill = value)) + 
  geom_tile() + 
  geom_contour(color = "white", alpha = 0.05) + 
  scale_fill_distiller(palette="Spectral", na.value="white") + 
  theme_bw()



# get contour -------------------------------------------------------------

getContour <- function(kden, prob){
 # kde.chose <- kden
  #kde <- raster::raster(kde.chose)
  
  kde[kde == 0]<-NA
  kde_values <- raster::getValues(kde)
  
  sortedValues <- sort(kde_values[!is.na(kde_values)],decreasing = TRUE)
  # find cumulative sum up to ith location
  sums <- cumsum(as.numeric(sortedValues))
  # binary response is value in the probabily zone or not
  p <- sum(sums <= prob * sums[length(sums)])
  if(p == 0){
    kdeprob <- NA
  }else{
    # Set values in raster to 1 or 0
    kdeprob <- raster::setValues(kde, kde_values >= sortedValues[p])
    # return new kde
    return(kdeprob)  
  }
}


#GBIF.cont <- lapply(kden,FUN = getContour,prob = 0.90)


GBIF.poly <- list()
GBIF.cont <- list()
for (i in 3519:length(kden.highr)){
  
  kde.chose <- kden.highr[[i]]
  kde <- raster::raster(kde.chose)
  kdeprob <- list(getContour(kde,prob = 0.99))
  GBIF.cont[i] <- kdeprob
}

#GBIF.cont.samp <- GBIF.cont[125]
for (i in 2184:length(GBIF.cont)){
  GBIF.c <- GBIF.cont[[i]]
  GBIF.c[GBIF.c == 0] <- NA
  GBIF.poly[i] <- rasterToPolygons(GBIF.c, dissolve = T)
}

GBIF.poly <- lapply(GBIF.cont, 
                      FUN = function(x){
                        x[x==0]<-NA
                        y <- rasterToPolygons(x, dissolve = TRUE)
                        return(y)
                      })

names(GBIF.poly) <- count_species_f5$species[1:length(GBIF.poly)]

a <- sapply(GBIF.poly, class)
b <- which(a != a[1])
GBIF.poly.new <- GBIF.poly[-b]
GBIF.poly.sf <- lapply(GBIF.poly.new, st_as_sf) 

GBIF.poly.sf <- do.call(rbind,GBIF.poly.sf) %>% 
  mutate(species = rownames(.))



# Separate into clusters --------------------------------------------------

kden.p_groups_pre_ii = data.frame()


for (i in 941:nrow(GBIF.poly.sf)){
  
  polygon_dat = GBIF.poly.sf[i,] %>% st_cast("POLYGON") 
  
  if (nrow(polygon_dat) == 0){
    next
  }else{
    
    polygon_dat$ID = as.character(1:nrow(polygon_dat))
    
    #filter for anythign that is just a line...
    #area_p = st_area(polygon_dat) %>% set_units(.,km^2)
    #a = set_units(5.5e+1, km^2)
    #polygon_dat = polygon_dat[area_p>a,]
    
    
    if (nrow(polygon_dat) ==2){
      polygon_groups = data.frame()
      polyg_dist = st_distance(x = polygon_dat[1,1], y = polygon_dat[2,1])
      
      d=set_units(5e+05, m)
      if(polyg_dist<d){
        polygon_groups[1,1] = st_sf(st_union(polygon_dat$geometry))
        polygon_groups[1,2] = "1"
        polygon_groups[1,3] = polygon_dat$species[1]
      }else{
        polygon_groups[1,1] = st_sf(st_union(polygon_dat$geometry[1])) 
        polygon_groups[1,2] = "1"
        polygon_groups[2,1] = st_sf(st_union(polygon_dat$geometry[2]))
        polygon_groups[2,2] = "2"
        polygon_groups[,3] = polygon_dat$species[1]
      }
    }else{
      #pairwise distances in relation to the centers
      polyg_dist= as.matrix(st_distance(x =polygon_dat[,1], y = polygon_dat[,1], by_element = FALSE))
      polyg_dist[lower.tri(polyg_dist)] <- 0
      
      zer = set_units(0,m)
      med = median(polyg_dist[polyg_dist!=zer])
      dist = 500000 %>% set_units(m)
      
      #identify if there is any polygon that is at a distance larger than 2xmedian distance (outliers)
      outl = (which(polyg_dist>dist,arr.ind = TRUE))
      
      outllin = as.vector(outl)
      
      polygon_groups = data.frame()
      if (length(outllin) ==0){
        polygon_groups[1,1] = st_sf(st_union(polygon_dat$geometry))
        polygon_groups[1,2] = "1"
        polygon_groups[1,3] = polygon_dat$species[1]
      }else{
        
        outl_true =  data.frame(table(outllin))
        a = outl_true[outl_true$Freq > 1,]
        a = as.numeric(as.vector(a$outllin))
        b = 1:nrow(polygon_dat)
        b = b[b!=a]
        
        if (length(a) ==1){
          polygon_groups[1,1] = st_sf(st_union(polygon_dat$geometry[a])) 
          polygon_groups[1,2] = "1"
          polygon_groups[2,1] = st_sf(st_union(polygon_dat$geometry[b]))
          polygon_groups[2,2] = "2"
          polygon_groups[,3] = polygon_dat$species[1]
          
          
          
        } else {
          
          polyg_dist= as.matrix(st_distance(x =polygon_dat[,1], y = polygon_dat[,1], by_element = FALSE))
          hc <- hclust(as.dist(polyg_dist), method="ward.D2")
          dist = 500000
          clust <- cutree(hc, h=dist)
          maxn = max(clust)
          
          for (j in 1:maxn){
            polygon_groups[j,1] = st_sf(st_union(polygon_dat$geometry[which(clust ==j)])) 
            polygon_groups[j,2] = as.character(j)
            
          }
          polygon_groups[,3] = polygon_dat$species[1]
          
        }
      }
    }
    colnames(polygon_groups) = c("polygon","group","species")
    polygon_groups = st_sf(polygon_groups)
    kden.p_groups_pre_ii = rbind(kden.p_groups_pre_ii,polygon_groups)
    
    
  }

}



# test trials ------------------------------------------------------------
st_crs(polygon_groups) <- st_crs(polygon_dat)

polyg_dist = as.matrix(st_distance(x =polygon_groups[,1], y = polygon_groups[,1], by_element = FALSE))
hc <- hclust(as.dist(polyg_dist), method="ward.D2")
dist = 1000000
clust <- cutree(hc, h=dist)
maxn = max(clust)

polygon_groups1 <- data.frame()
for (j in 1:maxn){
  polygon_groups1[j,1] = st_sf(st_union(polygon_groups$polygon[which(clust ==j)])) 
  polygon_groups1[j,2] = as.character(j)
  
}
polygon_groups1[,3] = polygon_groups$specie[1]
colnames(polygon_groups1) = c("polygon","group","species")
polygon_groups1 = st_sf(polygon_groups1)

polygon_groups1.spat <- as_Spatial(polygon_groups)
proj4string(polygon_groups1.spat) <-  CRS("+proj=longlat +datum=WGS84")
sp.dist <- spDists(polygon_groups1.spat)
hc <- hclust(as.dist(sp.dist), method="ward.D2")
dist = 1000
clust <- cutree(hc, h=dist)
maxn = max(clust)

polygon_groups1 <- data.frame()
for (j in 1:maxn){
  polygon_groups1[j,1] = st_sf(st_union(polygon_groups$polygon[which(clust ==j)])) 
  polygon_groups1[j,2] = as.character(j)
  
}
polygon_groups1[,3] = polygon_groups$specie[1]
colnames(polygon_groups1) = c("polygon","group","species")
polygon_groups1 = st_sf(polygon_groups1)

#############

polygon_groups1.SL <- as(polygon_groups1.spat, 'SpatialLinesDataFrame')

p.pts = sapply(polygon_groups1.SL@lines, FUN = function(x){
  spsample(x, n = 1000, type = "regular") 
})

proj4string(p.pts[[1]]) <-  CRS("+proj=longlat +datum=WGS84")

find.min <- function (a){
  b <- a[[1]]
  min(apply(coordinates((b)), 1, function(p) 
    spDistsN1(coordinates((b)), p )))}

crossmap::xwalk(p.pts,find.min)



# filter out --------------------------------------------------------------
st_crs(kden.p_groups_pre_ii) <- st_crs(polygon_dat)

kden.p_groups_II_ii <- data.frame()

for (i in 741:length(unique(kden.p_groups_pre_ii$specie)))
  {
  
  spec <- unique(kden.p_groups_pre_ii$specie)[i]
  polygon_groups <- kden.p_groups_pre_ii[which(kden.p_groups_pre_ii$specie == spec),]
  if (nrow(polygon_groups) == 0){
    next
  }else if (nrow(polygon_groups) == 1){
    polygon_groups1 <- polygon_groups[,c(3,1,2)]
    st_crs(polygon_groups1) <- st_crs(kden.p_groups_II)
  }else{
    
    polygon_groups1.spat <- as_Spatial(polygon_groups)
    proj4string(polygon_groups1.spat) <-  CRS("+proj=longlat +datum=WGS84")
    polygon_groups1.SL <- as(polygon_groups1.spat, 'SpatialLinesDataFrame')
    
    p.pts = sapply(polygon_groups1.SL@lines, FUN = function(x){
      spsample(x, n = 1000, type = "regular") })
    #proj4string(p.pts[[1]]) <-  CRS("+proj=longlat +datum=WGS84")
    
    min.dist <- matrix(ncol = length(p.pts), nrow = length(p.pts))
    for (ii in 1:length(p.pts)){
      a <- p.pts[[ii]]
      for (jj in 1:length(p.pts)){
        b <- p.pts[[jj]]
        
        min.dist[ii,jj] <- min(apply(coordinates(a), 1, function(p) 
          spDistsN1(coordinates(b), p )))
        
      }
      
    }
    
    hc <- hclust(as.dist(min.dist), method="ward.D2")
    dist = 50
    clust <- cutree(hc, h=dist)
    maxn = max(clust)
    
    
    polygon_groups1 <- data.frame()
    for (k in 1:maxn){
      polygon_groups1[k,1] = st_sf(st_union(polygon_groups$polygon[which(clust ==k)])) 
      polygon_groups1[k,2] = as.character(k)
      
    }
    polygon_groups1[,3] = polygon_groups$specie[1]
    
    
  }
  colnames(polygon_groups1) = c("polygon","group","species")
  polygon_groups1 = st_sf(polygon_groups1)
  kden.p_groups_II = rbind(kden.p_groups_II,polygon_groups1)
  
}


# plot --------------------------------------------------------------------

for (i in 1:30){
  num <- sample(1:nrow(kden.p_groups_trim),1)
  #dataset.chose <- kden.p_groups_II[which(kden.p_groups_II$specie == kden.p_groups_II$specie[num]),]
  #dataset.chose <- kden.p_groups_pre_ii[which(kden.p_groups_pre_ii$specie == kden.p_groups_pre_ii$specie[num]),]
  
  dataset.chose <- kden.p_groups_trim[which(kden.p_groups_trim$species == "Prunus persica"),]
  
  
  lims <- st_bbox(dataset.chose)
  nlrg <- 2
  p1 <- ggplot(dataset.chose, aes())+
    geom_polygon(data=countries, aes(long, lat, group=group), fill = "ivory2", col = "grey") +
    geom_sf(aes(fill = group), alpha = 0.9) +
    coord_sf(xlim = c(lims[1]-nlrg,lims[3]+nlrg), ylim = c(lims[2]-nlrg,lims[4]+nlrg))+
    ggtitle(dataset.chose$species[1])
  
  print(p1)
  
  
}
ggsave(filename = paste("example_kden2d",dataset.chose$species[1],".jpg"), width = 15, height = 10)

