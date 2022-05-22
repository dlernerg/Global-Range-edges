library(spdep)
library(sf)
library(dplyr)
counts.coord.inland2 <- counts.coord.inland[,1:5]
counts.coord.inland2 <- counts.coord.inland2[,2:5]/((counts.in.inland$V2))
counts.coord.inland2 <- cbind(counts.coord.inland$`1:nrow(global)`,counts.coord.inland2)
names(counts.coord.inland2) <- c("cell","N","S","E","W")

counts.coord.inland2[which(is.na(counts.coord.inland2$S)),2:5] = 0
global.centr <- as.data.frame(st_centroid(global.inland$geometry)) 

counts.coord.inland.all <- as.data.frame(rowSums(counts.coord.inland[,2:5]))
counts.coord.inland.all2 <- as.data.frame(rowSums(counts.coord.inland2[,2:5]))
counts.coord.inland.all$cell <- counts.coord.inland2$cell
counts.coord.inland.all2$cell <- counts.coord.inland2$cell

names(counts.coord.inland.all) <- c("all","cell")
names(counts.coord.inland.all2) <- c("all","cell")

counts.coord.inland2 <- counts.coord.inland[,2:5]
counts.coord.inland2 <- counts.coord.inland2/((counts.in.inland$V2))
names(counts.coord.inland2) <- c("N","S","E","W")
count.which <- counts.coord.inland2[complete.cases(counts.coord.inland2),]

global.inland2 <- global.inland[complete.cases(counts.coord.inland2),]
st_geometry(count.which) <- st_geometry(global.inland2)

neighbours <- poly2nb(count.which, queen = F)
listw <- nb2listw(include.self(neighbours), style = "W", zero.policy = TRUE)
coords <- c("N","S","E","W")

local_G <- data.frame()
for (i in 1:4){
  a <- count.which
  #globalMoran <- moran.test(a$N,listw)
  #moran <- moran.plot(a$all, listw = nb2listw(neighbours, style = "W", zero.policy = TRUE))
  #local <- localmoran(x = a$E, listw,zero.policy = TRUE) %>% as.data.frame()
  #local = local %>% 
  #  mutate(quantilegroup = as.numeric(unique(dvmisc::create_qgroups(Ii, groups = 10))))
  a_c <- a[,c(i,5)] 
  names(a_c)[1] <- "coord"
  local_g <- as.matrix(localG(a_c$coord, listw)) 
  local_g <- as.data.frame(cbind(a_c$coord, as.vector(local_g)))
  names(local_g) <- c("cell","gstat")
  #b_S <- cbind(global.inland,local_g)
  
  #local_g$bins <- (cut(local_g$gstat, breaks=c(-1,3,7,11,13)))
  local_g$coord <- coords[i]
  
  local_G <- rbind(local_G,local_g)
}


bon <- format.pval(p.adjustSP(pnorm(2*(abs(local_G$gstat)), lower.tail=FALSE),
                              neighbours, "bonferroni"))

local_G$corrected.p <- bon

