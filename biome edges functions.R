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
  
  load("C:/Users/davidle.WISMAIN//Box/lab folder/hotspots/paper/raw datat/4. quadrant freq/counts.coord.inland.norm.RData")
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
intrsct.mat.fun <- function(coord.intersect.df,internal){
  
  internal.id<- internal %>% mutate(ID = 1:nrow(internal))
  internal2 <- internal.id[which(lengths(internal.id$origins) == 2),]
  internal3 <- internal.id[which(lengths(internal.id$origins) == 3),]
  
  origins <- melt(internal$origins)
  origins2 <- melt(internal2$origins)
  origins2$ID <- rep(internal2$ID, each =2)
  origins3 <- melt(internal3$origins)
  origins3$ID <- rep(internal3$ID, each = 3)
  
  coord.intersect.df2 <- coord.intersect.df[!duplicated(coord.intersect.df$`global cell`),]
  sum.mat <- matrix(0, ncol = 14, nrow = 14)
  med.mat <- matrix(0, ncol = 14, nrow = 14)
  for (i in 1:14){
    for (j in 1:14){
      b1 <- i
      b2 <- j
      specific <- steps(b1,b2,coord.intersect.df)
      sum.mat[i,j] <- sum(coord.intersect.df2$`number RE`[which(coord.intersect.df2$`global cell` %in% specific)],na.rm = T)
      med.mat[i,j] <- median(coord.intersect.df2$`number RE`[which(coord.intersect.df2$`global cell` %in% specific)],na.rm = T)
    }
  }
  
  colnames(sum.mat) <- names_biome$names
  rownames(sum.mat) <- names_biome$names
  
  colnames(med.mat) <- names_biome$names
  rownames(med.mat) <- names_biome$names
  
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
  }
  diag(sum.mat) <- NA
  diag(med.mat) <- NA
  
  sum.mat2 <- (sum.mat/sum(sum.mat,na.rm = T)*100)
  med.mat2 <- (med.mat/sum(med.mat,na.rm = T)*100)
  
  combination.matricies <- list(sum.mat = sum.mat,sum.mat2 = sum.mat2,med.mat = med.mat,med.mat2 = med.mat2)
  
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
fct.btstrp.biome.cont.a <- function(counts.coord.inland.g,internal){
  
  rand.counts <- as.data.frame(counts.coord.inland.g$countss[sample(1:nrow(counts.coord.inland.g), nrow(counts.coord.inland.g), replace=F)])
  rand.cell <- as.data.frame(global.inland$geometry[sample(1:nrow(global.inland), nrow(counts.coord.inland.g), replace=F)])
  rand.all <- cbind(rand.counts,rand.cell) %>% st_as_sf() 
  st_crs(rand.all) <- st_crs(global.inland)
  names(rand.all) <- c("countss","geometry")
  
  df <- intrsct.fun(rand.all,internal = internal)
  mat <- intrsct.mat.fun(df,internal)
  
  mat.sum <- mat$sum.mat
  colnames(mat.sum) <- 1:14
  rownames(mat.sum) <- 1:14
  
  df2 <- melt(mat.sum)
  return(df2)
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
  
  comb.4test <- comb.df[1:max(comb.df$ID,na.rm = T),c(1,3)]
  names(comb.4test)[1] <- "biome"
  num <- data.frame()
  for (i in 1:nrow(comb.4test)){
    num[i,1] <- names_biome$V2[which(names_biome$names==comb.4test$biome[i])]
  }
  
  comb.4test$num <- num$V1
  comb.4test$value <- log(comb.4test$value)
  anov.4duncan <- aov(value ~ num, data = comb.4test)
  out <- agricolae::duncan.test(anov.4duncan,"num", alpha = alpha)
  cld <- out$groups
  tk <- group_by(comb.4test,biome) %>%
    summarise(mean = mean(value),quant = quantile(value,probs = 0.75)) %>% 
    arrange(desc(mean))
  tk$cld <- cld$groups   
  
  return(list(comb.df = comb.4test, tk = tk))
}

