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
