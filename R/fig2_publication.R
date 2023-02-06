library(tidyr)
library(tibb)
library(forcats)
library(tibble)

####biome edges organized II has the actual way of finding the p values of the biome edges in 
#the matrix compared to organized I

# Functions ---------------------------------------------------------------


#1.function prepare dfs
counts.prep.fun <- function(counts.coord.inland2,global.inland){
  require(sf)
  counts.coord.inland2$countss <- rowSums(counts.coord.inland2[,1:4])
  counts.coord.inland.g <- cbind(counts.coord.inland2,global.inland) %>% st_as_sf()
  
  #find which range edges overlap with which biome edges  
  counts.coord.inland.g <- counts.coord.inland.g[!is.nan(counts.coord.inland.g$countss),]
  #counts.coord.inland.g <- counts.coord.inland.g[which(counts.coord.inland.g$countss>0),]
  counts.coord.inland.g$perc <- (counts.coord.inland.g$countss/sum(counts.coord.inland.g$countss)*100)
  
  return(counts.coord.inland.g)
}
Gi.prep.fun <- function(local_G,global.inland){
  local_G.sf <- (local_G)
  local_G.sf$ID <-  rep(seq(nrow(local_G)/4),times =  4)
  local_G.sf = local_G.sf[which(local_G.sf$corrected.p<0.05),]
  
  load("C:/Users/davidle.WISMAIN/Box/lab folder/hotspots/paper/Revision/kernel density/gi hotspots/counts.coord.inland2.RData")
  
  global.inland2 <- global.inland[complete.cases(counts.coord.inland2),]
  
  Gi.cells <- local_G.sf %>% group_by(ID) %>% summarise(countss = sum(gstat))
  st_geometry(Gi.cells) <- st_geometry(global.inland2$geometry[Gi.cells$ID])
  
  return(Gi.cells)
}

#2.function to identify the hexagons at intersection between biomes
steps <- function(b1, b2,coord.intersect.df){
  
  internal.id<- internal %>% mutate(ID = 1:nrow(internal))
  internal2 <- internal.id[which(lengths(internal.id$origins) == 2),]
  internal3 <- internal.id[which(lengths(internal.id$origins) == 3),]
  
  origins <- melt(internal$origins)
  origins2 <- melt(internal2$origins)
  origins2$ID <- rep(internal2$ID, each =2)
  origins3 <- melt(internal3$origins)
  origins3$ID <- rep(internal3$ID, each = 3)
  
  step1 <- origins3$ID[which(origins3$value == b1)]
  step2 <- which(origins3$ID %in% step1)
  step3 <- origins3[step2,]
  step4 <- step3$ID[which(origins3[step2,] == b2)]
  specific3 <- coord.intersect.df$`global cell`[which(coord.intersect.df$`biome combination`%in% step4)]
  
  step1 <- origins2$ID[which(origins2$value == b1)]
  step2 <- which(origins2$ID %in% step1)
  step3 <- origins2[step2,]
  step4 <- step3$ID[which(origins2[step2,] == b2)]
  specific2 <- coord.intersect.df$`global cell`[which(coord.intersect.df$`biome combination` %in% step4)]
  
  specific <- c(specific2,specific3)
  return(specific)
}

#3.function spatially identify intersections between hexagons and biome int
intrsct.fun <- function(RE.df, internal){
  require(sf)
  coord.intersect <- st_intersects(RE.df,internal)
  coord.intersect.df <- as.data.frame(coord.intersect)
  coord.intersect.df$countss <- RE.df$countss[coord.intersect.df$row.id]
  names(coord.intersect.df) <- c('global cell',"biome combination","number RE")
  
  return(coord.intersect.df)
}
#4.find matricies (not normalized)
intrsct.mat.fun <- function(coord.intersect.df,internal,cells.df){
  
  internal.id<- internal %>% mutate(ID = 1:nrow(internal))
  internal2 <- internal.id[which(lengths(internal.id$origins) == 2),]
  internal3 <- internal.id[which(lengths(internal.id$origins) == 3),]
  
  origins <- melt(internal$origins)
  origins2 <- melt(internal2$origins)
  origins2$ID <- rep(internal2$ID, each =2)
  origins3 <- melt(internal3$origins)
  origins3$ID <- rep(internal3$ID, each = 3)
  
  outersct.df <- intrsct.fun(cells.df,wwf_simpl.b)
  coord.intersect.df2 <- coord.intersect.df[!duplicated(coord.intersect.df$`global cell`),]
  sum.mat <- matrix(0, ncol = 14, nrow = 14)
  med.mat <- matrix(0, ncol = 14, nrow = 14)
  perc.mat <- matrix(0, ncol = 14, nrow = 14)
  for (i in 1:14){
    htspINbiome.i <- outersct.df[outersct.df$`biome combination` == i,]
    for (j in 1:14){
      htspINbiome.j <- outersct.df[outersct.df$`biome combination` == j,]
      b1 <- i
      b2 <- j
      specific <- steps(b1,b2,coord.intersect.df)
      sum.inbiome <- sum(htspINbiome.i$`number RE`,htspINbiome.j$`number RE`)
      sum.mat[i,j] <- sum(coord.intersect.df2$`number RE`[which(coord.intersect.df2$`global cell` %in% specific)],na.rm = T)
      med.mat[i,j] <- median(coord.intersect.df2$`number RE`[which(coord.intersect.df2$`global cell` %in% specific)],na.rm = T)
      perc.mat[i,j] <- (sum(coord.intersect.df2$`number RE`[which(coord.intersect.df2$`global cell` %in% specific)],na.rm = T))/sum.inbiome
    }
  }
  
  colnames(sum.mat) <- names_biome$names
  rownames(sum.mat) <- names_biome$names
  
  colnames(med.mat) <- names_biome$names
  rownames(med.mat) <- names_biome$names
  
  colnames(perc.mat) <- names_biome$names
  rownames(perc.mat) <- names_biome$names
  
  #clean the matrix by removing impossible combinations 
  combinations <- data.frame()
  for(i in 1:max(origins$L1)){
    comb.mstep <- origins[origins$L1 == i,]
    if (nrow(comb.mstep) == 0){
      next
    }else{
      ab <- combn(comb.mstep$value,2) %>% as.data.frame()
      combinations <- rbind(combinations,as.data.frame(t(ab)))
    }
  }
  combinations <- unique(combinations)
  
  comb.test <- c(1:14) %>% combn(.,2) %>% as.data.frame() %>% t()
  comb.test <- rbind(combinations,comb.test)
  comb.test <- distinct(comb.test)
  comb.test <- comb.test[(nrow(combinations)+1):nrow(comb.test),]
  
  for (i in 1:nrow(comb.test)){
    a <- comb.test[i,1]
    b <- comb.test[i,2]
    sum.mat[a,b] <- NA
    sum.mat[b,a] <- NA
    med.mat[a,b] <- NA
    med.mat[b,a] <- NA
    perc.mat[a,b] <- NA
    perc.mat[b,a] <- NA
  }
  diag(sum.mat) <- NA
  diag(med.mat) <- NA
  diag(perc.mat) <- NA
  
  sum.mat2 <- (sum.mat/sum(sum.mat,na.rm = T)*100)
  med.mat2 <- (med.mat/sum(med.mat,na.rm = T)*100)
  
  combination.matricies <- list(sum.mat = sum.mat,sum.mat2 = sum.mat2,med.mat = med.mat,med.mat2 = med.mat2, perc.mat = perc.mat)
  
  return(combination.matricies)
}
#5.create a function that allows the removal of nas from the 
dist_no_na <- function(mat) {
  edist <- dist(mat)
  #edist[which(is.na(edist))] <- max(edist, na.rm=TRUE) * 0.0001
  edist[which(is.na(edist))] <- 0
  
  return(edist)
}
#6.numbers and perc of intrsct
calc.int.fun <- function(coord.intersect.df,RE.df){
  a3 <- unique(coord.intersect.df$`global cell`)
  a4.num <- sum(RE.df$countss[a3], na.rm = T)
  
  a4.length <- length(a3)
  a5.length <- round((a4.length/nrow(RE.df))*100,2)
  a5.num <- round((a4.num/sum(RE.df$countss, na.rm = T))*100,2)
  
  df <- data.frame(perc.hex = a5.length, perc.RE = a5.num)
  return(df)
}
#7.bootstrap the biome edges
fct.btstrp.biome.cont.sratified <- function(counts.coord.inland.g,internal){
  
  mat <- intrsct.mat.fun.rand(counts.coord.inland.g,internal)
  
  colnames(mat$sum.mat) <- 1:14
  rownames(mat$sum.mat) <- 1:14
  colnames(mat$perc.mat) <- 1:14
  rownames(mat$perc.mat) <- 1:14
  
  df.sum <- melt(mat$sum.mat)
  
  df.perc <- melt(mat$perc.mat)
  combination.matricies <- list(df.sum = df.sum,df.perc = df.perc)
  
  return(combination.matricies)
}

fct.btstrp.biome.cont.a <- function(counts.coord.inland.g,internal){
  
  rand.counts <- as.data.frame(counts.coord.inland.g$countss[sample(1:nrow(counts.coord.inland.g), nrow(counts.coord.inland.g), replace=F)])
  rand.cell <- as.data.frame(global.inland$geometry[sample(1:nrow(global.inland), nrow(counts.coord.inland.g), replace=F)])
  rand.all <- cbind(rand.counts,rand.cell) %>% st_as_sf() 
  st_crs(rand.all) <- st_crs(global.inland)
  names(rand.all) <- c("countss","geometry")
  
  df <- intrsct.fun(rand.all,internal = internal)
  mat <- intrsct.mat.fun(df,internal, cells.df = counts.coord.inland.g)
  
  mat.sum <- mat$sum.mat
  colnames(mat.sum) <- 1:14
  rownames(mat.sum) <- 1:14
  
  df2 <- melt(mat.sum)
  return(df2)
}

fct.btstrp.biome.cont.b.starified <- function(intrsct.df,df.btsrp,type.test){
  pval.sum <- matrix(ncol = 14, nrow = 14)
  mat.sum <- matrix(ncol = 14, nrow = 14)
  outersct.df <- intrsct.fun(Gi.cells,internal = outernal)
  for(i in 1:14){
    
    htspINbiome.i <- outersct.df[outersct.df$`biome combination` == i,]
    for (j in 1:14){
      b1 <- i
      b2 <- j
      
      htspINbiome.j <- outersct.df[outersct.df$`biome combination` == j,]
      sum.inbiome <- sum(htspINbiome.i$`number RE`,htspINbiome.j$`number RE`)
      
      intrsct.df2 <- intrsct.df[!duplicated(intrsct.df$`global cell`),]
      specific <- steps(b1,b2,intrsct.df)
      
      if (type.test == "sum"){
        RE.num <- sum(intrsct.df2$`number RE`[which(intrsct.df2$`global cell` %in% specific)])
        
      } else if (type.test == "perc"){
        RE.num <- (sum(intrsct.df2$`number RE`[which(intrsct.df2$`global cell` %in% specific)]))/sum.inbiome
      }
      
      p.val.sum <- mean(RE.num <= df.btsrp$sum[which(df.btsrp$ID1==b1&df.btsrp$ID2==b2)])
      
      pval.sum[i,j] <- p.val.sum
      pval.sum[j,i] <- p.val.sum
      
      
      mat.sum[i,j] <- mean(df.btsrp$sum[which(df.btsrp$ID1==b1&df.btsrp$ID2==b2)])
      mat.sum[j,i] <- mean(df.btsrp$sum[which(df.btsrp$ID1==b1&df.btsrp$ID2==b2)])
      
    }
  }
  
  combination.matricies <- list(sum.mat.btsrp = mat.sum, pval.sum = pval.sum)
  
  return(combination.matricies)  
} 
fct.btstrp.biome.cont.b <- function(intrsct.df,df.btsrp){
  pval.sum <- matrix(ncol = 14, nrow = 14)
  mat.sum <- matrix(ncol = 14, nrow = 14)
  for(i in 1:14){
    for (j in 1:14){
      b1 <- i
      b2 <- j
      
      intrsct.df2 <- intrsct.df[!duplicated(intrsct.df$`global cell`),]
      specific <- steps(b1,b2,intrsct.df)
      RE.num <- sum(intrsct.df2$`number RE`[which(intrsct.df2$`global cell` %in% specific)])
      p.val.sum <- mean(RE.num < df.btsrp$sum[which(df.btsrp$ID1==b1&df.btsrp$ID2==b2)])
      
      pval.sum[i,j] <- p.val.sum
      pval.sum[j,i] <- p.val.sum
      
      
      mat.sum[i,j] <- mean(df.btsrp$sum[which(df.btsrp$ID1==b1&df.btsrp$ID2==b2)])
      mat.sum[i,j] <- mean(df.btsrp$sum[which(df.btsrp$ID1==b1&df.btsrp$ID2==b2)])
      
    }
  }
  
  combination.matricies <- list(sum.mat.btsrp = mat.sum, pval.sum = pval.sum)
  
  return(combination.matricies)  
} 

#8 calculate the contribution of each biome
total.cont.fun <- function(comb.mat,alpha){
  comb.df <- melt(comb.mat) 
  comb.df <- comb.df %>% group_by(value) %>% arrange(desc(value))
  comb.df$ID <- c(1:nrow(comb.df))
  comb.df <- comb.df[comb.df$value>0,]
  comb.df$name <- paste(names_biome$names[comb.df$Var1],"vs ",names_biome$names[comb.df$Var2])
  comb.df <- comb.df[complete.cases(comb.df),]
  
  comb.4test <- comb.df[1:max(comb.df$ID,na.rm = T),c(1,3)]
  names(comb.4test)[1] <- "biome"
  num <- data.frame()
  for (i in 1:nrow(comb.4test)){
    num[i,1] <- names_biome$V2[which(names_biome$names==comb.4test$biome[i])]
  }
  
  comb.4test$num <- num$V1
  comb.4test$value <- (comb.4test$value)
  
  anov.4duncan <- aov(log(value) ~ num, data = comb.4test)
  out <- agricolae::duncan.test(anov.4duncan,"num", alpha = alpha)
  cld <- out$groups
  tk <- group_by(comb.4test,biome) %>%
    summarise(mean = mean(value),quant = quantile(value,probs = 0.75)) %>% 
    arrange(desc(mean))
  tk$cld <- cld$groups   
  
  return(list(comb.df = comb.4test, tk = tk))
}
#calculate the bootstrap statified for within each biome 
intrsct.mat.fun.rand <- function(Gi.cells,internal){

  rand.g <- data.frame()
  for (i in 1:14){
    rand.i <- rand[which(rand$`biome combination` == i),]
    biomecells <- biome2globecell.df$col.id[biome2globecell.df$row.id == i]
    rand.cell <- sample(biomecells,size = nrow(rand.i), replace = T)
    rand.i <- rand.i[,c(1,3)]
    names(rand.i) <- c("ID","countss")
    
    rand.i$geometry <- global.inland2$geometry[rand.cell]
    rand.i$ID <- rand.cell
    rand.i <- rand.i %>% st_as_sf()
    rand.g <- rbind(rand.g,rand.i)
  }
  
  intrsct.df.rand <- intrsct.fun(rand.g,internal = internal)      
  intrsct.df.rand2 <- intrsct.df.rand[!duplicated(intrsct.df.rand$`global cell`),]
  
  outersct.df <- intrsct.fun(rand.g,internal = outernal)      
  
  for (i in 1:14){
    htspINbiome.i <- outersct.df[outersct.df$`biome combination` == i,]
    for (j in 1:14) {
      htspINbiome.j <- outersct.df[outersct.df$`biome combination` == j,]
      b1 <- i
      b2 <- j
      specific <- steps(b1,b2,intrsct.df.rand)
      
      sum.inbiome <- sum(htspINbiome.i$`number RE`,htspINbiome.j$`number RE`)
      perc.mat[i,j] <- (sum(intrsct.df.rand2$`number RE`[which(intrsct.df.rand2$`global cell` %in% specific)],na.rm = T))/sum.inbiome
      sum.mat[i,j] <- sum(intrsct.df.rand2$`number RE`[which(intrsct.df.rand2$`global cell` %in% specific)],na.rm = T)
      
    }
  }
  
  colnames(perc.mat) <- names_biome$names
  rownames(perc.mat) <- names_biome$names
  colnames(sum.mat) <- names_biome$names
  rownames(sum.mat) <- names_biome$names
  
  #clean the matrix by removing impossible combinations 
  combinations <- data.frame()
  for(i in 1:max(origins$L1)){
    comb.mstep <- origins[origins$L1 == i,]
    if (nrow(comb.mstep) == 0){
      next
    }else{
      ab <- combn(comb.mstep$value,2) %>% as.data.frame()
      combinations <- rbind(combinations,as.data.frame(t(ab)))
    }
  }
  combinations <- unique(combinations)
  
  comb.test <- c(1:14) %>% combn(.,2) %>% as.data.frame() %>% t()
  comb.test <- rbind(combinations,comb.test)
  comb.test <- distinct(comb.test)
  comb.test <- comb.test[(nrow(combinations)+1):nrow(comb.test),]
  
  for (i in 1:nrow(comb.test)){
    a <- comb.test[i,1]
    b <- comb.test[i,2]
    sum.mat[a,b] <- NA
    sum.mat[b,a] <- NA
    perc.mat[a,b] <- NA
    perc.mat[b,a] <- NA
  }
  diag(sum.mat) <- NA
  diag(perc.mat) <- NA
  
  
  combination.matricies <- list(sum.mat = sum.mat,perc.mat = perc.mat)
  
  return(combination.matricies)
  
}

# script ------------------------------------------------------------------


library(sf)
library(dplyr)
library(ggplot2)
library(reshape2)
library(gplots)
library(RColorBrewer)

sf_use_s2(F)

my_palette <- colorRampPalette(c("white","tomato1", "firebrick2")) (n=20)


load("C:/Users/davidle.WISMAIN//Box/lab folder/hotspots/maps raw data/biomes/wwf_simple_b2.RData")
load("C:/Users/davidle.WISMAIN//Box/lab folder/hotspots/paper/raw datat/4. quadrant freq/global.inland.RData")
load("C:/Users/davidle.WISMAIN/Box/lab folder/hotspots/maps raw data/biomes/names_biome.RData")
load("C:/Users/davidle.WISMAIN/Box/lab folder/hotspots/paper/raw datat/biome.int/interserctions.b2.RData")
load("C:/Users/davidle.WISMAIN/Box/lab folder/hotspots/paper/Revision/kernel density/gi hotspots/local_G.RData")
load("C:/Users/davidle.WISMAIN/Box/lab folder/hotspots/paper/Revision/kernel density/gi hotspots/counts.coord.inland.RData")

counts.coord.inland.NORM <- counts.coord.inland[,2:5]/((counts.in$V2))
#need to increase the buffer of intersection between biomes
intersections.b <- intersections.b %>% st_buffer(0.3)

int <- which(lengths(intersections.b$origins)!=1)
biom <- which(lengths(intersections.b$origins)==1)

internal <- intersections.b[int,]
outernal <- intersections.b[biom,]

#### preparing cells for analysis
Gi.cells <- Gi.prep.fun(local_G,global.inland)
cont.cells <- counts.prep.fun(counts.coord.inland.NORM,global.inland)

#### intersection df and matrix
intrsct.df <- intrsct.fun(Gi.cells,internal = internal)

intrsct.mat <- intrsct.mat.fun(internal = internal, coord.intersect.df = intrsct.df, cells.df = Gi.cells)

################## NOW CHOOSE EITHER BOOTSTRAP "ALL BIOMES" or BOOTSTRAP "STRATIFIED", THEN MOVE ON 
#TO SECTION 3

# sECTION 2A: bootstrap all biomes ----------------------------------------------------

######calculate bootrstrap using global reshuffling (not stratified)

df2 <- data.frame()
for (i in 150:500) {
  df <- data.frame()
  df <- fct.btstrp.biome.cont.a(Gi.cells,internal)
  df2 <- rbind(df2,df)
}

df3.globe.Gi <- df2[complete.cases(df2),]

#load("C:/Users/davidle.WISMAIN/Box/lab folder/hotspots/paper/Revision/bootstrap fig2/df3.globe.RData")

df.use <- df3.globe
names(df.use) <- c("ID1","ID2","sum")

#global bootstrap
btsrp.mat2 <- fct.btstrp.biome.cont.b(intrsct.df,df.use)


# SECTION 2B: bootstrab stratified BIOMES----------------------------------------------------
####calculate bootstrab for each biome intersect (stratified) 


df.perc <- data.frame()
df.sum <- data.frame()
load("C:/Users/davidle.WISMAIN/Box/lab folder/hotspots/paper/Revision/kernel density/gi hotspots/counts.coord.inland2.RData")


global.inland2 <- global.inland[complete.cases(counts.coord.inland2),]
internal.id<- internal %>% mutate(ID = 1:nrow(internal))
internal2 <- internal.id[which(lengths(internal.id$origins) == 2),]
internal3 <- internal.id[which(lengths(internal.id$origins) == 3),]

origins <- melt(internal$origins)
origins2 <- melt(internal2$origins)
origins2$ID <- rep(internal2$ID, each =2)
origins3 <- melt(internal3$origins)
origins3$ID <- rep(internal3$ID, each = 3)

perc.mat <- matrix(0, ncol = 14, nrow = 14)
sum.mat <- matrix(0, ncol = 14, nrow = 14)

rand <- intrsct.fun(Gi.cells,internal = outernal)
rand <- rand[!duplicated(rand$`global cell`),]
biome2globecell <- st_intersects(wwf_simpl.b,global.inland2)
biome2globecell.df <- as.data.frame(biome2globecell)

for (i in 1:500) {
  df <- fct.btstrp.biome.cont.a(Gi.cells,internal)
  df.sum <- rbind(df.sum,df$df.sum)
  df.perc <- rbind(df.perc,df$df.perc)
}

df3.strat.Gi <- list(df.sum = df.sum, df.perc = df.perc)

load("C:/Users/davidle.WISMAIN/Box/lab folder/hotspots/paper/Revision/bootstrap fig2/df3.strat.Gi.RData")

#if using coordinates (rather than hotspots)
df.use <- df3.strat.Gi$df.sum[complete.cases(df3.strat.Gi$df.sum),]
names(df.use) <- c("ID1","ID2","sum")


#stratified biome
btsrp.mat2 <- fct.btstrp.biome.cont.b.starified(intrsct.df,df.use,type.test = "sum")


# -------------------------------------------------------------------------
#Section 3: Create Per-biome contribution matrix

reorder.mat <- c(1,2,3,7,9,14,11,6,12,13,4,5,8,10)

btsrp.mat <- btsrp.mat2$sum.mat.btsrp[reorder.mat,reorder.mat]
comb.mat <- (intrsct.mat$sum.mat)
comb.mat2 <- comb.mat[reorder.mat,reorder.mat]
comb.mat3 <- (comb.mat2)/btsrp.mat

pval.sum <- btsrp.mat2$pval.sum[reorder.mat,reorder.mat]
pval <- matrix(ncol = 14, nrow = 14)
pval[(which(pval.sum<=0.085))] = "*"
diag(pval) <- NA
pval[which(comb.mat2 == 0)] = NA


comb.df = comb.mat3 %>% 
  log() %>%
  as.data.frame() %>%
  rownames_to_column("r.id") %>%
  pivot_longer(-c(r.id), names_to = "c.id", values_to = "counts") %>%
  mutate(c.id= fct_relevel(c.id,colnames(comb.mat2))) %>%
  mutate(r.id= fct_relevel(r.id,colnames(comb.mat2))) 

pval.df <- pval %>% 
  as.data.frame() %>% 
  rownames_to_column("r.id") %>%
  pivot_longer(-c(r.id), names_to = "c.id", values_to = "counts")

comb.df$p.val <- pval.df$counts

ggplot(comb.df,aes(x=r.id, y=c.id, fill=counts)) + 
  geom_tile() +
  geom_text(aes(label = p.val, size = 10)) +
  scale_fill_gradient2(low = "blue",mid = "white",high = "red")
                       
ggsave(filename = "2bi.svg", dpi = 600, height = 5,width = 7)

# SECTION 4 find correlations BETWEEN BIOMES/WITHIN BIOMES -------------------------------------------------------------------------

cont.lst <- total.cont.fun(comb.mat2,0.05)
cont.df <- cont.lst$comb.df
cont.df2 <- cont.df %>% group_by(biome) %>% summarise(summary = mean(value, na.rm = T))

#what is the number of species in each of the biomes
outersct.df <- intrsct.fun(Gi.cells,outernal)
outersct.df2 <- outersct.df %>% group_by(`biome combination`) %>% summarize(sum = mean(`number RE`))
#normalize this number to the number obtained in a randomized bootstrap (just like with the edges)

#outersct.df$NORM <- outersct.df$sum/apply(df.out, 1, mean)
outersct.df2$name <- names_biome$names

reord <- match(cont.df2$biome,outersct.df2$name)
a <- data.frame(edge = cont.df2$summary,biome = outersct.df2$sum[reord], names = names_biome$names[reord])
a[,c(1:2)] <- log(a[,c(1:2)])
m <- lm(edge~biome,data = a)
p1=predict(m,interval="confidence",level=0.95)
a$group = 1
a$group[a$edge<p1[,2] | a$edge>p1[,3]]=0
a$group <- as.character(a$group)
reord2 <- reord[c(6,13,9,14,12,10,3,2,7,1,4,5,8,11)]
a$names <- factor(a$names,levels = names_biome$names[c(1,10:14,2:9)])

p1 <- ggplot(a, aes(x = biome, y = edge, label = names)) + 
  geom_point(aes(fill = names),shape = 21, color = "black", size =5) +
  scale_fill_manual(values = myPalette(14))+
  ggrepel::geom_text_repel(size = 7.5, direction = "both")+
  stat_smooth(method = lm, level = 0.95) + 
  stat_cor(p.accuracy = NULL) + 
  xlab("log(hotspot within biome)") +
  #scale_x_continuous(limits = c(0,10))+
  ylab("log(mean hotspots at biome edge)") + 
  theme(legend.position = "none")


p1<- p1 + theme_test()
ggsave("correlations.log.svg", dpi = 300,width = 17.56, height = 9)

