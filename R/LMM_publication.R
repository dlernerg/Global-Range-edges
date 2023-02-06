load("C:/Users/davidle.WISMAIN//Box/lab folder/hotspots/maps raw data/biomes/wwf_simple_b2.RData")
make.continents <- function(){
  library(maptools)
  library(cleangeo) 
  library(rworldmap)
  library(dplyr)
  library(sf)
  
  
  sPDF <- getMap()
  sPDF <- clgeo_Clean(sPDF)  ## Needed to fix up some non-closed polygons 
  cont <-
    sapply(levels(sPDF$REGION),
           FUN = function(i) {
             ## Merge polygons within a continent
             poly <- gUnionCascaded(subset(sPDF, REGION==i))
             ## Give each polygon a unique ID
             poly <- spChFIDs(poly, i)
             ## Make SPDF from SpatialPolygons object
             SpatialPolygonsDataFrame(poly,
                                      data.frame(REGION=i, row.names=i))
           },
           USE.NAMES=TRUE)
  
  ## Bind the 6 continent-level SPDFs into a single SPDF
  cont <- Reduce(spRbind, cont)
  
  cont <- st_as_sf(cont)
  return(cont)
}

continents <- make.continents()
which.cont <- st_intersects(global.inland,continents) %>% data.frame()
which.biom <- st_intersects(global.inland,wwf_simpl.b) %>% data.frame()


which.cont.df.true <- c(1:5851) %>% as.data.frame()
which.cont.df.true$continent <- NA
which.cont.df.true$continent <- which.cont$col.id[match(which.cont.df.true$.,which.cont$row.id)]
which.cont.df.true$biome <- NA
which.cont.df.true$biome <- which.biom$col.id[match(which.cont.df.true$.,which.biom$row.id)]

####not quite the way we should do it. Other clean for repeats or include the repeats appropriately so there is not a bias


# prepare data ------------------------------------------------------------
load("C:/Users/davidle.WISMAIN/Box/lab folder/hotspots/paper/Revision/kernel density/gi hotspots/counts.coord.inland.RData")
load("C:/Users/davidle.WISMAIN/Box/lab folder/hotspots/paper/Revision/kernel density/counts.coord.RData")
load("C:/Users/davidle.WISMAIN/Box/lab folder/hotspots/paper/Revision/kernel density/counts.in.RData")
load("C:/Users/davidle.WISMAIN/Box/lab folder/hotspots/paper/raw datat/5. preparing data for glm/altitude.RData")
load("C:/Users/davidle.WISMAIN/Box/lab folder/hotspots/paper/raw datat/5. preparing data for glm/bio.global.all.RData")
load("C:/Users/davidle.WISMAIN/Box/lab folder/hotspots/paper/raw datat/5. preparing data for glm/elevation.RData")
load("C:/Users/davidle.WISMAIN/Box/lab folder/hotspots/paper/raw datat/5. preparing data for glm/meansBioClim.RData")
load("C:/Users/davidle.WISMAIN/Box/lab folder/hotspots/paper/raw datat/5. preparing data for glm/PV.hex.M2.all.RData")
load("C:/Users/davidle.WISMAIN/Box/lab folder/hotspots/paper/raw datat/5. preparing data for glm/PV.altitude.RData")
load("C:/Users/davidle.WISMAIN/Box/lab folder/hotspots/paper/Revision/envirem/means.envirem.RData")
load("C:/Users/davidle.WISMAIN/Box/lab folder/hotspots/paper/Revision/envirem/PV.envirem.RData")

counts.coord.inland.unnorm <- counts.coord.inland
counts.coord.inland.unnorm$countss <- rowSums(counts.coord.inland.unnorm[,2:5])

counts.coord$countss <- rowSums(counts.coord[,2:5])
counts.failiure <- counts.in
counts.failiure$V2 <- counts.in$V2 - (counts.coord$countss)
counts.failiure$V2[which(counts.failiure$V2<0)] = 0
counts.coord.inland.unnorm$species.t <- counts.failiure$V2

meansBioClim2$bio20 <- altitude$m

PV.hex.M2.all$V20 <- altitude.pv[,1]

means.std <- scale(meansBioClim2)
PV.std <- scale(PV.hex.M2.all)

means.std.e <- scale(means.envirem)
PV.std.e <- scale(PV.envirem)
colnames(means.std.e) <- c(paste0("mean.E", 1:16))
colnames(PV.std.e) <- c(paste0("PV.E", 1:16))

cent <- st_centroid(global.inland) %>% st_coordinates() %>% as.data.frame()

alltog.in2 <- cbind(counts.coord.inland.unnorm[,6:7],means.std,PV.std)

#add ENVIREM variables

#alltog.in2 <- cbind(counts.coord.inland.unnorm[,6:7],means.std.e,PV.std.e)

alltog.in2$cent <- cent$Y
alltog.in2.cont <- alltog.in2
alltog.in2.cont$continent <- which.cont.df.true$continent
alltog.in2.cont$biome <- which.cont.df.true$biome
alltog.in2.cont <- alltog.in2.cont[complete.cases(alltog.in2.cont),]

alltog.in2.cont$continent <- as.factor(alltog.in2.cont$continent)
alltog.in2.cont$biome <- as.factor(alltog.in2.cont$biome)


# LMM ---------------------------------------------------------------------


m2.df <- data.frame()
aic2 <- data.frame()
r2_c <- data.frame()
r2_m <- data.frame()

for (i in 1:(ncol(alltog.in2.cont)-4))
{
  names.i <- names(alltog.in2.cont[i+2])
  alltog2 <- alltog.in2.cont %>% 
    dplyr::select(countss,species.t,continent,biome,all_of(names.i))%>% 
    as.data.frame()
  
  names(alltog2) <- c("all","failiure","continent","biome","var.d")
  
  
  m2 <- lme4::glmer(cbind(alltog2$all,alltog2$failiure)~  var.d + (var.d|continent)+ (var.d|biome), data = alltog2, 
            family = binomial)
  

  a <- summary(m2)
  
  m2.df[1:4,i]<- (as.data.frame(coef(summary(m2))[2,]))
  aic2[i,1] <-AIC(m2)
  r2 <- MuMIn::r.squaredGLMM(m2)
  r2_c[i,1] <- r2[3]
  r2_m[i,1] <- r2[1]
  
}


m2.df.t <- as.data.frame(t(m2.df))
names(m2.df.t) <- c("Estimate","std error","Z value","P-value")

r2.df <- as.data.frame(cbind(r2_m,r2_c))

final_table.all <- as.data.frame(cbind(m2.df.t,r2.df,aic2))
names(final_table.all)[5:7] <- c("r2m","r2c","AIC")


r2.df[,1] <- as.data.frame(r2_m)
r2.df[,2] <- c(1:nrow(r2.df))
r2.df[,3] <- m2.df.t$`P-value`
r2.df[,4] <- m2.df.t$Estimate
names(r2.df) <- c("V1","V1.1","V3","V4")
names.plot <- 1:(ncol(alltog.in2)-7)

r2.df$V1.1 <- factor(r2.df$V1.1, labels = names.plot)
plott1 <- ggplot(r2.df, aes(x = V1.1, y = V1, fill = V4,
                         label = ifelse(V3 < 0.05,
                                        ifelse(V3 >0.01,"*","**"),""))) + 
  geom_bar(stat="identity") + 
  geom_text(vjust = 0)+
  theme_test()+
  theme(axis.text.x = element_text(angle = 25, hjust = 1, size = 6),
        legend.text = element_text(angle = -20,size = 6),
        legend.title.align = 0.2,
        legend.title = element_text(size = 14),
        axis.title.x = element_blank()) +
  scale_fill_gradient2(expression(beta),low = "#ff0000",mid = "white", high = "#1f3c88")+
  ylab(bquote(italic(R)^2))+ 
  scale_x_discrete(labels=c(paste0("mean ","envirem",names.plot[1:16]),paste0("PV envirem",names.plot[1:16])))+
  ggtitle("all")
 

plot(plott1)


ggsave(plott1,filename = "GLM_envirem_all.svg",width = 8.1, height = 3.2)

svg(filename = "GLM_envirem_all.svg",width = 8.1, height = 3.2)
dev.off()

# -------------------------------------------------------------------------

#make global model
var.mod.bioclim <- c(1,5,6,8,10,11,21,22,23,25,26,27,28,29,30,31,32,33,36,40)

var.mean.envirem <- c(1,6,7,8,9,10,11,16)
var.PV.envirem <- c(1,3:16)
var.mod.envirem <- c(var.mean.envirem,(var.PV.envirem+16)) 

var.mod <- c(var.mod.bioclim)
var.mod <- c((var.mod.envirem),33)

var.nam <- names(alltog.in2.cont[var.mod+2])
fmla <- as.formula(paste("cbind(alltog.in2.cont$countss,alltog.in2.cont$species.t) ~ ", 
                         paste(var.nam, collapse= "+"),
                         "+ (1|continent) + (1|biome)"))


global.m <- lme4::glmer(fmla, data = alltog.in2.cont, 
                        family = binomial, na.action = na.fail)

summary(global.m)

# identify correlated variables  ------------------------------------------


is.correlated <- function(i, j, data, conf.level = 0.95, cutoff = 0.7, ...) {
  if(j >= i) return(NA)
  ct <- cor.test(data[, i], data[, j], conf.level = conf.level, ...)
  ct$p.value > (1 - conf.level) || abs(ct$estimate) <= cutoff
}

# Need vectorized function to use with 'outer'
vCorrelated <- Vectorize(is.correlated, c("i", "j"))

# Create logical matrix
nm <- colnames(alltog.in2.cont[,var.mod+2])
smat <- outer(1:length(nm), 1:length(nm), vCorrelated, data = alltog.in2.cont[,var.mod+2])
dimnames(smat) <- list(nm, nm)


# -------------------------------------------------------------------------


library(snow)

cl <- snow::makeCluster(8, type = "SOCK")
clusterEvalQ(cl, library(lme4))
clusterExport(cl,"alltog.in2.cont")

model.dredge <- MuMIn::pdredge(global.m, beta = T, evaluate = T, rank = "AIC" , trace = 2, cluster = cl, subset = smat)


stopCluster(cl)


### 
model.dredge <- model.dredge.bioclim
model.dredge <- model.dredge.envirem


head(model.dredge)
model.f <- subset(model.dredge, delta <= 2)
#summary(model.f)
model.f.avg <- summary(model.avg(model.f))




# plot --------------------------------------------------------------------
m.full <- as.data.frame(model.f.avg$coefmat.subset)
CI <- as.data.frame(confint(model.f.avg, full=F)) # get confidence intervals for full model
m.full <- m.full %>% 
  mutate(CI.min = CI$`2.5 %`) %>% 
  mutate(CI.max = CI$`97.5 %`) %>% 
  mutate(worldclim = row.names(m.full)) %>% 
  filter(`Pr(>|z|)`<0.05) %>%
  mutate(model = "partial")

names(m.full) <- gsub(" ", "", names(m.full)) # remove spaces from column headers

dataset <- "ENVIREM"
forest <- ggplot(data=m.full[-1,], aes(x=worldclim, y=Estimate))+ #again, excluding intercept because estimates so much larger
  geom_hline(yintercept=0, color = "grey",linetype="dashed", lwd=1)+ #add dashed line at zero
  geom_errorbar(aes(ymin=CI.min, ymax=CI.max), colour="black",alpha = 0.7, #adj SE
                width=0, lwd=1.5) +
  #geom_errorbar(aes(ymin=Estimate-AdjustedSE, ymax=Estimate+AdjustedSE), colour="blue", #adj SE
  #             width=0, lwd=2) +
  geom_point(size=4, color = "#FFCCFF")+
  theme_bw(base_size = 14)+ 
  scale_x_discrete(position = "top")+
  coord_flip()+ # flipping x and y axes+
  theme(plot.margin = margin(1, 1, 1, 1, "cm"))+
  ylab("Coefficient") + xlab(paste( dataset,"variables"))

forest
ggsave(paste0("model_sel.svg", dataset, ".png"), dpi = 300)

imprt <- importance(summary(model.avg(subset(model.dredge, delta <= 2)))) %>% 
  as.data.frame()
imprt$variable <- rownames(imprt)

#imprt$variable <- factor(imprt$variable, levels = rownames(imprt))
imprt$n <- 1:nrow(imprt)
ggplot(imprt, aes(x = variable, y = .)) + 
  geom_col(aes(fill = 1/.),show.legend = F) + 
  #ggtitle(paste("Relative Importance", dataset)) + 
  theme_bw()+
  #theme(axis.text.y = element_blank(),
  #     axis.title.y = element_blank()) +
  ylab("Sum of weights") + 
  #scale_x_discrete(limits = rev)+
  coord_flip()

#ggplot2::scale_fill_gradient2(high = "#FF9999")+

ggsave(paste0("Rel.imp.",dataset,".png"), dpi = 300, height = 3.8, width = 2.1)
