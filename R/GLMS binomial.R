library(effects) 
library(factoextra) #PCA analyses
library(sjPlot) #plot CI forest plots
library(ggpubr)
library(spdep) #spatial dependency
library(statmod) #tweedie
library(lme4)
library(lmerTest)
#Run models

  


#############################GLM for the different continents

#identify the failiure by removing the number of species (i.e. range edges) in that grid cell - both the coastile and inland

inlandquad <- as.numeric(global.inland$cell)

counts.coord.inland$countss <- rowSums(counts.coord.inland2[,2:5])
counts.coord.inland.unnorm <- counts.coord.inland[inlandquad,]
counts.coord.inland.unnorm$countss <- rowSums(counts.coord.inland.unnorm[,2:5])

counts.coord$countss <- rowSums(counts.coord[,2:5])
counts.failiure <- counts.in
counts.failiure$V2 <- counts.in$V2 - (counts.coord$countss)
counts.failiure$V2[which(counts.failiure$V2<0)] = 0
counts.coord.inland.unnorm$species.t <- counts.failiure$V2[inlandquad]

########################prepare date
meansBioClim2$bio20 <- altitude$m

PV.hex.M2.all$V20 <- altitude.pv[,1]

means.std <- scale(meansBioClim2)
PV.std <- scale(PV.hex.M2.all)

alltog.in2 <- cbind(counts.coord.inland.unnorm,means.std,PV.std)

coords <- st_coordinates(st_centroid(global.inland))
ac <- autocov_dist(as.numeric(cbind(alltog.in2$countss,alltog.in2$species.t)),coords, nbs = 500, zero.policy = T, longlat = T)

m2.df <- data.frame()
aic2 <- data.frame()
r2_LR <- data.frame()


###############test for autocorrelation using Moran's test 
library(MASS)
library(nlme)
  g <- st_coordinates(global.inland$cent) %>% as.data.frame()
  
  plot(g$X, g$Y, col=c("blue","red")[sign(resid(m2))/2+1.5], 
       pch=19,cex=abs(resid(m2))/max(resid(m2))*2, 
       xlab="geographical xcoordinates", 
       ylab="geographical y-coordinates")
  
  
  neighbours <- poly2nb(global.inland2)
  listw2 <- nb2listw(neighbours, zero.policy = T)
  
  moran.test(residuals(m2), listw2, zero.policy = T)
  
###############GLM 
for (i in 1:40)
  {
    alltog2 <- as.data.frame(alltog.in2[,c(6,7,i+7)])
    names(alltog2) <- c("all","failiure","var.d")
    
    #ac <- autocov_dist(all, g, nbs = 10, type = "1", zero.policy = T)
    #alltog2 <- cbind(alltog2,ac)
    #m2 <- glm(cbind(alltog2$all,alltog2$failiure)~ var.d + ac, data = alltog2, family = binomial)
    
    m2 <- glm(cbind(alltog2$all,alltog2$failiure)~ var.d  , data = alltog2, 
              family = binomial)
    
    
    summary(m2)
    
    m2.df[1:4,i]<- (as.data.frame(coef(summary(m2))[2,]))
    aic2[i,1] <-AIC(m2)
    r2_LR[i,1] <- MuMIn::r.squaredLR(m2)
    #r2_mf[i,1] <-((m2$null.deviance-m2$deviance)/m2$null.deviance)*100
  
}
  
  m2.df <- as.data.frame(t(m2.df))
  names(m2.df) <- c("Estimate","std error","Z value","P-value")
  
  r2 <- as.data.frame(r2_LR)
  
  final_table <- as.data.frame(cbind(m2.df,r2,aic2))
  names(final_table)[5:6] <- c("r2","AIC")
  nam = paste0("finaltab",j)
  assign(nam, final_table)
  write.csv(final_table, paste("modelB7,in.csv"))
  # regline <- plot(predictorEffects(m1, residuals=TRUE),axes=list(grid=TRUE, x=list(novelty=seq(1, 2, 0.25))), partial.residuals=list(smooth=TRUE,span=0.75,lty="dashed")) 
  
  r2[,2] <- c(1:nrow(r2))
  r2[,3] <- m2.df$`P-value`
  r2[,4] <- m2.df$Estimate
  
  names.plot <- 1:40
  
  r2$V2 <- factor(r2$V2, labels = names.plot)
  plott1 <- ggplot(r2, aes(x = V2, y = V1, fill = V4,
                           label = ifelse(V3 < 0.05,
                                          ifelse(V3 >0.01,"*","**"),""))) + 
    geom_bar(stat="identity") + 
    theme_test()+
    theme(axis.text.x = element_text(angle = 25, hjust = 1, size = 6),
          legend.text = element_text(angle = -20,size = 6),
          legend.title.align = 0.2,
          legend.title = element_text(size = 14),
          axis.title.x = element_blank()) +
    scale_fill_gradient2(expression(beta),low = "#ff0000",mid = "white", high = "#1f3c88")+
    ylab(bquote(italic(R)^2))+ 
    scale_x_discrete(labels=c(paste0("mean ","BC",names.plot[1:19]),paste0("mean ","elevation"),paste0("PV BC",names.plot[1:19]),paste0("PV ","elevation")))
  
  
  plot(plott1)
