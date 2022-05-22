library(units)

GBIFpolygon_groups_pre = data.frame()

for (i in 1:max(polygons$V2)){

  polygon_dat = polygons[(polygons$V2 == i),]
  polygon_dat = st_sf(polygon_dat)# %>%   st_cast("MULTIPOLYGON") 
  
  if (nrow(polygon_dat) == 0){
    next
  }else{
    
    polygon_dat$ID = as.character(1:nrow(polygon_dat))
    
    #filter for anythign that is just a line...
    area_p = st_area(polygon_dat) %>% set_units(.,km^2)
    a = set_units(5.5e+2, km^2)
    polygon_dat = polygon_dat[area_p>a,]
    
    
    if (nrow(polygon_dat) ==2){
      polygon_groups = data.frame()
      polyg_dist = st_distance(x = polygon_dat[1,1], y = polygon_dat[2,1])
      
      d=set_units(5e+05, m)
      if(polyg_dist<d){
        polygon_groups[1,1] = st_sf(st_union(polygon_dat$polygons))
        polygon_groups[1,2] = "1"
        polygon_groups[1,3] = polygon_dat$specie[1]
      }else{
        polygon_groups[1,1] = st_sf(st_union(polygon_dat$polygons[1])) 
        polygon_groups[1,2] = "1"
        polygon_groups[2,1] = st_sf(st_union(polygon_dat$polygons[2]))
        polygon_groups[2,2] = "2"
        polygon_groups[,3] = polygon_dat$specie[1]
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
          polygon_groups[1,1] = st_sf(st_union(polygon_dat$polygons))
          polygon_groups[1,2] = "1"
          polygon_groups[1,3] = polygon_dat$specie[1]
        }else{
          
          outl_true =  data.frame(table(outllin))
          a = outl_true[outl_true$Freq > 1,]
          a = as.numeric(as.vector(a$outllin))
          b = 1:nrow(polygon_dat)
          b = b[b!=a]
          
          if (length(a) ==1){
            polygon_groups[1,1] = st_sf(st_union(polygon_dat$polygons[a])) 
            polygon_groups[1,2] = "1"
            polygon_groups[2,1] = st_sf(st_union(polygon_dat$polygons[b]))
            polygon_groups[2,2] = "2"
            polygon_groups[,3] = polygon_dat$specie[1]
            
            
            
          } else {
            
            polyg_dist= as.matrix(st_distance(x =polygon_dat[,1], y = polygon_dat[,1], by_element = FALSE))
            hc <- hclust(as.dist(polyg_dist), method="complete")
            dist = 500000
            clust <- cutree(hc, h=dist)
            maxn = max(clust)
            
            for (j in 1:maxn){
              polygon_groups[j,1] = st_sf(st_union(polygon_dat$polygons[which(clust ==j)])) 
              polygon_groups[j,2] = as.character(j)
              
            }
            polygon_groups[,3] = polygon_dat$specie[1]
            
           }
         }
       }
  
    colnames(polygon_groups) = c("polygon","group","specie")
    polygon_groups = st_sf(polygon_groups)
    GBIFpolygon_groups_pre = rbind(GBIFpolygon_groups_pre,polygon_groups)
    
  }
}

