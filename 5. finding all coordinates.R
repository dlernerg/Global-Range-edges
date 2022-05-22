library(units)
library(lwgeom)
library(ggplot2)
library(sf)


#1:
#create the four quartile polygons for each of the distributions (only the largest area)

rownames(GBIFpolygon_groups) <- c(1:nrow(GBIFpolygon_groups))

polygon_dat = GBIFpolygon_groups
polygon_dat2 = GBIFpolygon_groups


SEcoord_polygon = data.frame()
SWcoord_polygon = data.frame()
NEcoord_polygon = data.frame()
NWcoord_polygon = data.frame()

for (i in 1:nrow(GBIFpolygon_groups))
{
  
  tree = data.frame()
  tree = GBIFpolygon_groups[i,]
  tree_centroid = st_centroid(tree)
  
  a = st_bbox(tree)
  b = st_bbox(tree_centroid)
  
  SEpolygon = rbind(c(b[1],b[2]),c(b[1],a[2]),c(a[3],a[2]),c(a[3],b[2]),c(b[1],b[2]))
  SEpolygon_true = st_polygon(list(SEpolygon))
  SEpolygon_true = st_sfc(SEpolygon_true)
  
  SWpolygon = rbind(c(b[1],b[2]),c(b[1],a[2]),c(a[1],a[2]),c(a[1],b[2]),c(b[1],b[2]))
  SWpolygon_true = st_polygon(list(SWpolygon))
  SWpolygon_true = st_sfc(SWpolygon_true)
  
  NEpolygon = rbind(c(b[1],b[2]),c(b[1],a[4]),c(a[3],a[4]),c(a[3],b[2]),c(b[1],b[2]))
  NEpolygon_true = st_polygon(list(NEpolygon))
  NEpolygon_true = st_sfc(NEpolygon_true)
  
  NWpolygon = rbind(c(b[1],b[2]),c(b[1],a[4]),c(a[1],a[4]),c(a[1],b[2]),c(b[1],b[2]))
  NWpolygon_true = st_polygon(list(NWpolygon))
  NWpolygon_true = st_sfc(NWpolygon_true)
  
  SEcoord_polygon[i,1] = st_sf(SEpolygon_true, crs = "WGS84") 
  SWcoord_polygon[i,1] = st_sf(SWpolygon_true, crs = "WGS84") 
  NEcoord_polygon[i,1] = st_sf(NEpolygon_true, crs = "WGS84") 
  NWcoord_polygon[i,1] = st_sf(NWpolygon_true, crs = "WGS84")
  SEcoord_polygon[i,2] = polygon_dat$specie[i]
  SWcoord_polygon[i,2] = polygon_dat$specie[i]
  NEcoord_polygon[i,2] = polygon_dat$specie[i]
  NWcoord_polygon[i,2] = polygon_dat$specie[i]
  
}
SEcoord_polygon = st_sf(SEcoord_polygon,crs="WGS84")
SWcoord_polygon = st_sf(SWcoord_polygon,crs="WGS84")
NEcoord_polygon = st_sf(NEcoord_polygon,crs="WGS84")
NWcoord_polygon = st_sf(NWcoord_polygon,crs="WGS84")


#2:
#making N S E $ W data points
{
  marginS = data.frame()
  marginN = data.frame()
  marginE = data.frame()
  marginW = data.frame()
  
  for (i in 1:nrow(GBIFpolygon_groups))
  {
    
    tree = data.frame()
    tree =  GBIFpolygon_groups[i,]
    bbox_tree = st_bbox(tree)
    
    mat_tree= st_coordinates(tree)
    a = which(mat_tree[,1]==bbox_tree[3])
    if (is.null(nrow(mat_tree[a,1:2])) == TRUE) {
      marginE[i,1:2] = mat_tree[a,1:2]
    } else{
      marginE[i,1:2] = mat_tree[a,1:2][1,]
    }
    
    b = which(mat_tree[,1]==bbox_tree[1])
    if (is.null(nrow(mat_tree[b,1:2])) == TRUE){
      marginW[i,1:2] = mat_tree[b,1:2]
    } else{
      marginW[i,1:2] = mat_tree[b,1:2][1,]  
    }
    
    c =   which(mat_tree[,2]==bbox_tree[4])
    if (is.null(nrow(mat_tree[c,1:2])) == TRUE){
      marginN[i,1:2] = mat_tree[c,1:2]
    } else{
      marginN[i,1:2] = mat_tree[c,1:2][1,]
    }
    
    
    d =   which(mat_tree[,2]==bbox_tree[2])
    if (is.null(nrow(mat_tree[d,1:2])) == TRUE){
      marginS[i,1:2] = mat_tree[d,1:2]
    } else{
      marginS[i,1:2] = mat_tree[d,1:2][1,]
    }
  }
  
}

#3
#break up each polygon into their quartiles and find each's two most extreme

SEbbox = data.frame()
SWbbox = data.frame()
NEbbox = data.frame()
polygonSE = data.frame()
polygonSW = data.frame()
marginSE_1 = data.frame()
marginSE_2 = data.frame()
marginSW_1 = data.frame()
marginSW_2 = data.frame()
marginNE_1 = data.frame()
marginNE_2 = data.frame()
marginNW_1 = data.frame()
marginNW_2 = data.frame()


for (i in 1:nrow(GBIFpolygon_groups))
{
  
  #function to find the largest area of i
  tree = data.frame()
  tree =  GBIFpolygon_groups[i,]
  
  #break up each polygon into their quartiles
  croppedSE = st_crop(tree,SEcoord_polygon[i,])
  croppedSW = st_crop(tree,SWcoord_polygon[i,])
  croppedNE = st_crop(tree,NEcoord_polygon[i,])
  croppedNW = st_crop(tree,NWcoord_polygon[i,])
  
  bboxSE = st_bbox(croppedSE)
  bboxSW = st_bbox(croppedSW)
  bboxNE = st_bbox(croppedNE)
  bboxNW = st_bbox(croppedNW)
  
  
  matSE = sf::st_coordinates(croppedSE)
  a_1 = which(matSE[,1]==bboxSE[3])
  a_2 = which(matSE[,2]==bboxSE[2])
  if (is.integer0(a_1)==TRUE) {  #in case there is no overlap between the polygon and the S or E. 
    marginSE_1[i,1:2] = NA
    
  } else {
    if (is.null(nrow(matSE[a_1,1:2])) == TRUE ){
      marginSE_1[i,1:2] = matSE[a_1,1:2]
    } else {
      marginSE_1[i,1:2] = colMeans(matSE[a_1,1:2])
    }}
  
  if (is.integer0(a_2)==TRUE) {  #in case there is no overlap between the polygon and the S or E. 
    marginSE_2[i,1:2] = NA
  } else {
    if (is.null(nrow(matSE[a_2,1:2])) == TRUE ) {
      marginSE_2[i,1:2] = matSE[a_2,1:2]
      
    }else{
      marginSE_2[i,1:2] = colMeans(matSE[a_2,1:2])
      
    }
    
  }
  
  matSW = sf::st_coordinates(croppedSW)
  b_1 = which(matSW[,1]==bboxSW[1])
  b_2 = which(matSW[,2]==bboxSW[2])
  if (is.integer0(b_1)==TRUE) {
    marginSW_1[i,1:2] = NA
    marginSW_2[i,1:2] = NA
    
    
  } else {
    if (is.null(nrow(matSW[b_1,1:2])) == TRUE ){
      marginSW_1[i,1:2] = matSW[b_1,1:2]
    } else{
      marginSW_1[i,1:2] = colMeans(matSW[b_1,1:2])
      
    }
    if (is.null(nrow(matSW[b_2,1:2])) == TRUE) {
      marginSW_2[i,1:2] = matSW[b_2,1:2]
      
    }else{
      marginSW_2[i,1:2] = colMeans(matSW[b_2,1:2])
      
    }
  }
  
  matNE = sf::st_coordinates(croppedNE)
  c_1 = which(matNE[,1]==bboxNE[3])
  c_2 = which(matNE[,2]==bboxNE[4])
  if (is.integer0(c_1)==TRUE) {
    marginNE_1[i,1:2] = NA
    
  } else {
    if (is.null(nrow(matNE[c_1,1:2])) == TRUE ){
      marginNE_1[i,1:2] = matNE[c_1,1:2]
    } else{
      marginNE_1[i,1:2] = colMeans(matNE[c_1,1:2])
      
    }}
  
  if (is.integer0(c_2)==TRUE) {  #in case there is no overlap between the polygon and the S or E. 
    marginNE_2[i,1:2] = NA
  } else {
    if (is.null(nrow(matNE[c_2,1:2])) == TRUE ) {
      marginNE_2[i,1:2] = matNE[c_2,1:2]
      
    }else{
      marginNE_2[i,1:2] = colMeans(matNE[c_2,1:2])
      
    }
    
  }
  
  matNW = sf::st_coordinates(croppedNW)
  d_1 = which(matNW[,1]==bboxNW[1])
  d_2 = which(matNW[,2]==bboxNW[4])
  if (is.integer0(d_1)==TRUE) {
    marginNW_1[i,1:2] = NA
    
  } else {
    if (is.null(nrow(matNW[d_1,1:2])) == TRUE ){
      marginNW_1[i,1:2] = matNW[d_1,1:2]
    } else{
      marginNW_1[i,1:2] = colMeans(matNW[d_1,1:2])
      
    }}
  
  if (is.integer0(d_2)==TRUE) {  #in case there is no overlap between the polygon and the S or E. 
    marginNW_2[i,1:2] = NA
  } else {
    
    if (is.null(nrow(matNW[d_2,1:2])) == TRUE) {
      marginNW_2[i,1:2] = matNW[d_2,1:2]
      
    }else{
      marginNW_2[i,1:2] = colMeans(matNW[d_2,1:2])
      
    }
  }
}

#4
#Identify all of the edge coordinates. The two RE from the same cardinal distance are accounted for if d>15 arc-deg, otherwise only one   
  
  margin_ref_1= data.frame()
  margin_ref_2 = data.frame()
  margin_ref_3 = data.frame()
  
  for (i in 1:nrow(GBIFpolygon_groups)){
    m2 = data.frame()
    m1 = data.frame()
    m3 = data.frame() 
    
    #alter this number when analysing the different cardinal coordinates 
    #1 = E/W #2 = N/S
    a = marginSW_1[i,1:2] 
    b = marginNW_1[i,1:2]
    
    #will define if we account the long or lat axis: 2 when dealing with E/W margins and 1 when dealing with N/S
    c = 2 
    d = marginW[i,1:2]
    
    
    
    if (is.na(abs(b[c]-a[c])) == TRUE || abs(b[c]-a[c])<15) { #add the is.na because there are many situations where there is no marginSE and margin SW (tiny distributions)
      
      m1[1,1:2] = d 
      
    } else {
      m2[1,1:2] = a 
      m3[1,1:2] = b 
      
      
    }
    margin_ref_1 = rbind(margin_ref_1,m1)
    margin_ref_2 = rbind(margin_ref_2,m2)
    margin_ref_3 = rbind(margin_ref_3,m3)
    
    
  }
  
  #add the name of the tree to each margin ID. Use the fact that row name was saved for each. 
  #GBIFpolygon_groups$specie = paste(GBIFpolygon_groups$specie, GBIFpolygon_groups$group)
  V1 = paste(GBIFpolygon_groups$specie[as.numeric(rownames(margin_ref_1))],GBIFpolygon_groups$group[as.numeric(rownames(margin_ref_1))],GBIFpolygon_groups$world[as.numeric(rownames(margin_ref_1))])
  V2 = paste(GBIFpolygon_groups$specie[as.numeric(rownames(margin_ref_2))],GBIFpolygon_groups$group[as.numeric(rownames(margin_ref_2))],GBIFpolygon_groups$world[as.numeric(rownames(margin_ref_2))])
  length(c(V1,V2))
  
  margin_ref_1$ID = "all"
  margin_ref_1$species = V1
  margin_ref_2$ID = "side a"
  margin_ref_2$species = V2
  margin_ref_3$ID = "side b"
  margin_ref_3$species = V2
  
  margin_allW = rbind(margin_ref_1,margin_ref_2,margin_ref_3)
 
###########identify INLAND coordinates

library(dplyr)
library(rnaturalearth)
library(sf)

margin_allE$coord = "E"
margin_allN$coord = "N" 
margin_allS$coord = "S" 
margin_allW$coord = "W"
margin_all = rbind(margin_allE,margin_allN,margin_allS,margin_allW)
margin_all$species = paste(margin_all$species, margin_all$coord)
margin_all$species2 <- paste(margin_all$species,margin_all$ID)


options("download.file.method" = "libcurl")
oceans <- ne_download(scale = 110, type = 'ocean',category = 'physical', returnclass =  "sf")

coastline <- ne_coastline(scale = 110, returnclass =  "sf") %>% st_buffer(.,3)
margin_all_sf = st_as_sf(margin_all,coords = c("V1","V2"), crs = "WGS84")

margin_all = margin_all %>% mutate(coasts = apply(st_intersects(margin_all_sf,coastline), 1, any))

margin_all_inland = margin_all[margin_all$coasts==FALSE,]
margin_all_coast = margin_all[margin_all$coasts==TRUE,]
margin_all_coast_sf = st_as_sf(margin_all_coast,coords = c("V1","V2"), crs = "WGS84")

fract_out = data.frame()

for (i in 1:nrow(margin_all_coast)){
 
  n = margin_all_coast_sf[i,]
  n_b = st_buffer(n,3)
  lims = st_bbox(n_b)
  xy = st_coordinates(n)
  
  if (n$coord == "N"){
    cardinal = rbind(c(lims[1],xy[2]),c(lims[1],lims[4]),c(lims[3],lims[4]),c(lims[3],xy[2]),c(lims[1],xy[2]))
    cardinal = st_sfc(st_polygon(list(cardinal)), crs = "WGS84")
  }else if (n$coord == "S"){
    cardinal = rbind(c(lims[1],xy[2]),c(lims[1],lims[2]),c(lims[3],lims[2]),c(lims[3],xy[2]),c(lims[1],xy[2]))
    cardinal = st_sfc(st_polygon(list(cardinal)), crs = "WGS84")
  } else if (n$coord == "E"){
    cardinal = rbind(c(xy[1],lims[4]),c(lims[3],lims[4]),c(lims[3],lims[2]),c(xy[1],lims[2]),c(xy[1],lims[4]))
    cardinal = st_sfc(st_polygon(list(cardinal)), crs = "WGS84")
  }else if (n$coord == "W"){
    cardinal = rbind(c(xy[1],lims[2]),c(lims[1],lims[2]),c(lims[1],lims[4]),c(xy[1],lims[4]),c(xy[1],lims[2]))
    cardinal = st_sfc(st_polygon(list(cardinal)), crs = "WGS84")
  }
    
  cardinal = st_sf(cardinal)
  

  
  #get the fraction of out vs in  
  n_int = st_intersection(n_b,cardinal)
  coast = st_intersection(n_int,oceans)
 
    if (nrow(coast) == 0) {fract_out[i,1] = 0} else{fract_out[i,1] = st_area(coast)/st_area(n_int)}
}

  
  
  margin_all_coast$fraction = fract_out
  
  margin_coast = margin_all_coast[margin_all_coast$fraction>=0.5,]
  margin_inland = rbind(margin_all_coast[margin_all_coast$fraction<0.5,1:7],margin_all_inland)


