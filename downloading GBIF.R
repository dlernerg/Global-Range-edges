library(rgbif) # for occ_download
library(CoordinateCleaner)
library(countrycode)
library(dplyr)

a_1 = a%>%
  pull("Taxon name") %>% # use fewer names if you want to just test
  taxize::get_gbifid_(method="backbone") %>% # match names to the GBIF backbone to get taxonkeys
  imap(~ .x %>% mutate(original_sciname = .y)) %>% # add original name back into data.frame
  bind_rows() %T>% # combine all data.frames into one
  readr::write_tsv(path = "all_matches.tsv") %>% # save as side effect for you to inspect if you want
  filter(matchtype == "EXACT" & status == "ACCEPTED") %>% # get only accepted and matched names
  filter(kingdom == "Plantae") %>% # remove anything that might have matched to a non-plant
  pull(usagekey) # get the gbif taxonkeys

taxon_key_10000 = a %>% pull("Taxon name")
taxon_key_10000$taxonkey = c(a_1)

a1_n <- split(a_1, ceiling(seq_along(a_1)/11200))
a_1a <- a1_n[[1]]
a_1b <- a1_n[[2]]
a_1c <- a1_n[[3]]
a_1d <- a1_n[[4]]
a_1e <- a1_n[[5]]


occ_download(
  pred_in("taxonKey", a_1e),
  pred_gte("year", 1980),
  pred("hasGeospatialIssue", FALSE),
  format = "SIMPLE_CSV",
  user=user,pwd=pwd,email=email
)

#check if GBIF has finished preparing the data
occ_download_meta(key = "0317244-200613084148143") #1st half of world
occ_download_meta(key = "0317246-200613084148143") #2nd half of world
occ_download_meta(key = "0317247-200613084148143") #2nd half of world
occ_download_meta(key = "0317263-200613084148143") #2nd half of world
occ_download_meta(key = "0317264-200613084148143") #2nd half of world

#download the file into the working directory
GBIFpre_1a = occ_download_get(key="0317244-200613084148143", overwrite=TRUE) %>% occ_download_import()
GBIFpre_1b = occ_download_get(key="0317246-200613084148143", overwrite=TRUE) %>% occ_download_import()
GBIFpre_1c = occ_download_get(key="0317247-200613084148143", overwrite=TRUE) %>% occ_download_import()
GBIFpre_1d = occ_download_get(key="0317263-200613084148143", overwrite=TRUE) %>% occ_download_import()
GBIFpre_1e = occ_download_get(key="0317264-200613084148143", overwrite=TRUE) %>% occ_download_import()

GBIFkey_all <- c("0317244-200613084148143","0317246-200613084148143","0317247-200613084148143","0317263-200613084148143","0317264-200613084148143")

for (i in 2:5)
  {
  
  newdir <- paste0("Run",i)
  dir.create(newdir)      # should test for error
  cwd <- getwd()          # CURRENT dir
  setwd(newdir) 
  
  #filter by size 
  GBIFkey <- GBIFkey_all[i]
  GBIFpre = occ_download_get(key=GBIFkey, overwrite=TRUE) %>% occ_download_import()
  
  GBIFtree_f1a = GBIFpre %>% group_by(species) %>% filter(n()<500000)
  GBIFtree_f1a = GBIFtree_f1a %>% group_by(species) %>% filter(n()>300)
  
  count_species_f1a <- GBIFtree_f1a %>% count(species) %>% arrange(desc(n))
  
  
  GBIFtree_f1 = GBIFtree_f1a %>% filter(coordinateUncertaintyInMeters <= 100000| is.na(coordinateUncertaintyInMeters))
  count_species_f1 <- GBIFtree_f1 %>% count(species) %>% arrange(desc(n))
  
  
  
  GBIFtree_f1small =  GBIFtree_f1[,c(1,10,36)] %>% mutate(filtered = ifelse(basisOfRecord != "FOSSIL_SPECIMEN" & basisOfRecord != "UNKNOWN" ,1,0))
  GBIFtree_f2 = GBIFtree_f1[GBIFtree_f1small$filtered == 1,]
  remove(GBIFtree_f1small)
  count_species_f2 <- GBIFtree_f2 %>% count(species) %>% arrange(species)
  
  
  #3 Filter using CoordinateCleaner
  #convert country code from ISO2c to ISO3c
  GBIFtree_f2$countryCode <-  countrycode(GBIFtree_f2$countryCode, origin =  'iso2c', destination = 'iso3c')
  
  
  #make a smaller tbl_df for the sake of time consuming 
  GBIFtree_f2b <- GBIFtree_f2[which(!is.na(GBIFtree_f2$decimalLatitude)),]
  GBIFtree_f2small = GBIFtree_f2b[,c(1,10,12,16,22,23,36)]
  
  save(GBIFtree_f2b, file = "GBIFtree_f3.RData")
  
  flags<- clean_coordinates(x = GBIFtree_f2b,
                            lon = "decimalLongitude",
                            lat = "decimalLatitude",
                            countries = "countryCode",
                            species = "species",
                            tests = c("capitals","centroids","equal","gbif","zeros","seas","duplicates"))
  
  
  
  #filter the flaged
  
  GBIFtree_f3 <- GBIFtree_f2b[flags$.summary,]
  count_species_f3 <- GBIFtree_f3 %>% count(species) %>% arrange(desc(n))
  
  save(GBIFtree_f3, file = "GBIFtree_f3.RData")
  remove(GBIFtree_f2)
  remove(GBIFtree_f2small)
  remove(GBIFtree_f2b) 
  
  #4 Filter by occuerences - only n greater than 500
  #count_species <- GBIFtree_f2 %>% count(species) %>% filter(n>500) %>% arrange(desc(n))
  #spec_to_pull = as.character(count_species) 
  GBIFtree_f4 = GBIFtree_f3 %>% group_by(species) %>% filter(n()>300)
  count_species_f4 <- GBIFtree_f4 %>% count(species) %>% arrange(desc(n))
  
  save(GBIFtree_f4, file = "GBIFtree_f4.RData")
  remove(GBIFtree_f3)
  remove(GBIFtree_f4)
  setwd(cwd)
  
}
