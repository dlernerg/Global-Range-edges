#cleaning GBIF data
library(CoordinateCleaner)
library(countrycode)
library(dplyr)

setwd("C:/Users/davidle/Box Sync/lab folder/HotSpots/maps raw data/GBIF")

load("GBIFtree_unfiltered.RData")
load("GBIFtree_f2.RData")

#filter by size 
GBIFtree_f1a = GBIFfinal_australia %>% group_by(species) %>% filter(n()<500000)
GBIFtree_f1a = GBIFtree_f1a %>% group_by(species) %>% filter(n()>300)

count_species_f1a <- GBIFtree_f1a %>% count(species) %>% arrange(desc(n))



#1. Look for Coordinate uncertenties and filter

#plot a cumfreq of the uncertainty coordinates
a = GBIFtree_f1a %>% filter(coordinateUncertaintyInMeters != "NA")
a = a$coordinateUncertaintyInMeters 
a = a/1000

breaks = exp(seq(log(1), log(max(a)), length.out = round(max(a)/100)))
coord_uncer.cut = cut(a, breaks, right=FALSE) 
coord_uncer.freq = table(coord_uncer.cut)
coord_uncer.cumfreq = cumsum(coord_uncer.freq)
coord_uncer.cumfreq = max(coord_uncer.cumfreq) - coord_uncer.cumfreq
coord_uncer.cumfreq = coord_uncer.cumfreq/nrow(GBIFtree_f1a)
coord_uncer.cumfreq0 = c(1, coord_uncer.cumfreq)

plot(breaks, coord_uncer.cumfreq0, xlab="distance in km", ylab = "1-cumfreq", log="xy")
lines(breaks,coord_uncer.cumfreq0)


GBIFtree_f1 = GBIFtree_f1a %>% filter(coordinateUncertaintyInMeters <= 100000| is.na(coordinateUncertaintyInMeters))

count_species_f1 <- GBIFtree_f1 %>% count(species) %>% arrange(desc(n))

#2. Remove fossil and unknown records
table(GBIFtree_f1$basisOfRecord)

#GBIFtree_f2 = filter(GBIFtree_f1, basisOfRecord != "UNKNOWN")
#GBIFtree_f2 = filter(GBIFtree_f2, basisOfRecord != "FOSSIl_SPECIMEN")

GBIFtree_f1small =  GBIFtree_f1[,c(1,10,36)] %>% mutate(filtered = ifelse(basisOfRecord != "FOSSIL_SPECIMEN" & basisOfRecord != "UNKNOWN" ,1,0))
GBIFtree_f2 = GBIFtree_f1[GBIFtree_f1small$filtered == 1,]

count_species_f2 <- GBIFtree_f2 %>% count(species) %>% arrange(species)



#3 Filter using CoordinateCleaner
#convert country code from ISO2c to ISO3c
GBIFtree_f2$countryCode <-  countrycode(GBIFtree_f2$countryCode, origin =  'iso2c', destination = 'iso3c')


#make a smaller tbl_df for the sake of time consuming 
GBIFtree_f2small = GBIFtree_f2[,c(1,10,12,16,22,23,36)]

flags<- clean_coordinates(x = GBIFtree_f2small,
                           lon = "decimalLongitude",
                           lat = "decimalLatitude",
                           countries = "countryCode",
                           species = "species",
                           tests = c("capitals","centroids","equal"))#
flags2<- clean_coordinates(x = GBIFtree_f2small,
                          lon = "decimalLongitude",
                          lat = "decimalLatitude",
                          countries = "countryCode",
                          species = "species",
                          tests = c("gbif","zeros","seas"),
                          )

summary(flags2)

#create tables to know what tree species are the most flagged
{
flag_by_specie.cap = as.data.frame(table(flags$species[flags$.cap == FALSE])) 
flag_by_specie.cap$Var1 = as.character(flag_by_specie.cap$Var1)
#add another row with normalized frequencies
for (i in 1:nrow(flag_by_specie.cap)){flag_by_specie.cap$NormFreq[i] = (flag_by_specie.cap$Freq[i]/(count_species_f2$n[count_species_f2$species==flag_by_specie.cap$Var1[i]]))*100}
flag_by_specie.cap = flag_by_specie.cap%>% arrange(desc(NormFreq))
flag_by_specie.cap$flag = "capital"

#flag centroids                                                   
flag_by_specie.cen = as.data.frame(table(flags$species[flags$.cen == FALSE])) 
flag_by_specie.cen$Var1 = as.character(flag_by_specie.cen$Var1)
for (i in 1:nrow(flag_by_specie.cen)){flag_by_specie.cen$NormFreq[i] = (flag_by_specie.cen$Freq[i]/(count_species_f2$n[count_species_f2$species==flag_by_specie.cen$Var1[i]]))*100}
flag_by_specie.cen = flag_by_specie.cen%>% arrange(desc(NormFreq))
flag_by_specie.cen$flag = "centroid"

#flag equal                                                   
flag_by_specie.equ = as.data.frame(table(flags$species[flags$.equ == FALSE])) 
flag_by_specie.equ$Var1 = as.character(flag_by_specie.equ$Var1)
for (i in 1:nrow(flag_by_specie.equ)){flag_by_specie.equ$NormFreq[i] = (flag_by_specie.equ$Freq[i]/(count_species_f2$n[count_species_f2$species==flag_by_specie.equ$Var1[i]]))*100}
flag_by_specie.equ = flag_by_specie.equ%>% arrange(desc(NormFreq))
flag_by_specie.equ$flag = "equal"

#flag zeros
flag_by_specie.zer = as.data.frame(table(flags2$species[flags2$.zer == FALSE])) 
flag_by_specie.zer$Var1 = as.character(flag_by_specie.zer$Var1)
for (i in 1:nrow(flag_by_specie.zer)){flag_by_specie.zer$NormFreq[i] = (flag_by_specie.zer$Freq[i]/(count_species_f2$n[count_species_f2$species==flag_by_specie.zer$Var1[i]]))*100}
flag_by_specie.zer = flag_by_specie.zer%>% arrange(desc(NormFreq))
flag_by_specie.zer$flag = "zeros"

#flag seas
flag_by_specie.sea = as.data.frame(table(flags2$species[flags2$.sea == FALSE])) 
flag_by_specie.sea$Var1 = as.character(flag_by_specie.sea$Var1)
for (i in 1:nrow(flag_by_specie.sea)){flag_by_specie.sea$NormFreq[i] = (flag_by_specie.sea$Freq[i]/(count_species_f2$n[count_species_f2$species==flag_by_specie.sea$Var1[i]]))*100}
flag_by_specie.sea = flag_by_specie.sea%>% arrange(desc(NormFreq))
flag_by_specie.sea$flag = "sea"

flag_by_specie = bind_rows(flag_by_specie.cap,flag_by_specie.cen,flag_by_specie.zer,flag_by_specie.sea,flag_by_specie.equ)

#flag_by_specie = bind_rows(flag_by_specie.cap,flag_by_specie.cen,flag_by_specie.zer)

}

p = ggplot(flag_by_specie[flag_by_specie$NormFreq>10,], aes(x =  reorder(Var1,NormFreq), y = NormFreq, fill=flag)) + 
  geom_bar(stat = "identity", width=.6, position="dodge") + labs(x = "Specie", y = "Flagged norm Freq") + 
 # scale_y_log10(breaks = trans_breaks("log10", function(y) 10^y),
  #              labels = trans_format("log10", math_format(10^.x)))+
  coord_flip()
p


#filter the flaged
flaggsall <- flags %>% mutate(as.logical(pmin(.summary, flags2$.summary)))

GBIFtree_f3 <- GBIFtree_f2[flaggsall$.summary,]

count_species_f3 <- GBIFtree_f3 %>% count(species) %>% arrange(desc(n))

#4 remove duplicates 

flags<- clean_coordinates(x = GBIFtree_f3,
                          lon = "decimalLongitude",
                          lat = "decimalLatitude",
                          species = "species",
                          tests = c("duplicates"))

#flag duplicates
flag_by_specie.dpl = as.data.frame(table(flags$species[flags$.dpl == FALSE])) 
flag_by_specie.dpl$Var1 = as.character(flag_by_specie.dpl$Var1)
for (i in 1:nrow(flag_by_specie.dpl)){flag_by_specie.dpl$NormFreq[i] = (flag_by_specie.dpl$Freq[i]/(count_species_f3$n[count_species_f3$species==flag_by_specie.dpl$Var1[i]]))*100}
flag_by_specie.dpl = flag_by_specie.dpl%>% arrange(desc(NormFreq))
flag_by_specie.dpl$flag = "duplicates"

p = ggplot(flag_by_specie.dpl, aes(x =  reorder(Var1,NormFreq), y = NormFreq, fill=flag)) + 
  geom_bar(stat = "identity", width=.6, position="dodge") + labs(x = "Specie", y = "Flagged norm Freq") + 
  # scale_y_log10(breaks = trans_breaks("log10", function(y) 10^y),
  #              labels = trans_format("log10", math_format(10^.x)))+
  coord_flip()
p

GBIFtree_f4 <- GBIFtree_f3[flags$.summary,]

count_species_f4 <- GBIFtree_f4 %>% count(species) %>% arrange(desc(n))

#5 Filter by occuerences - only n greater than 500
#count_species <- GBIFtree_f2 %>% count(species) %>% filter(n>500) %>% arrange(desc(n))
#spec_to_pull = as.character(count_species) 
GBIFtree_f5 = GBIFtree_f4 %>% group_by(species) %>% filter(n()>500)

count_species_f5 <- GBIFtree_f5 %>% count(species) %>% arrange(desc(n))

