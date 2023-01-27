#### Legume Meta Analysis ####
#Authors: Kathryn Bloodworth, Smriti Pehim Limbu, Kim Komatsu

setwd('G:\\.shortcut-targets-by-id\\1w2OXIzBKQqFZ0BCeKP7C9pX36ViGDPBj\\Legume-Meta Analysis') #kim's wd

library(tidyverse)

##### data import #####
plantStatus <- read.csv('Data/legume_strain diversity_meta analysis_plant associations_edited names_presence absence.csv') %>% 
  select(paper_id, genus, species, plant_status)

plantData <- read.csv("Data/legume_strain diversity_meta analysis_plant data.csv") %>% 
  #change all 999 to NAs
  mutate_all(~na_if(., 999)) %>% 
  select(-notes, -plant_status) %>% 
  full_join(plantStatus) %>% 
  mutate(plant_status=ifelse(plant_status=='invasive', 'introduced', plant_status))
plantData[plantData==""]<-NA

#NOTE: There are about 230 observations with mis-matches between the presence/absence datafile and the plant data datafile. Ignoring these for now, since we're just looking for an estimate of whether to go forward with genetic analysis at all.
test <- plantData %>% 
  filter(is.na(plant_status))


##### determine what species have home/away or native/non-native comparisons and genetic data #####
homeAway <- plantData %>% 
  filter(compares_homeaway==1) %>% #76 observations of home and away
  filter(genbank==1) %>% #49 observations of those that have genetic data
  group_by(genus, species, plant_status) %>% 
  summarise(gene_regions=length(unique(genetic_region)), studies=length(unique(paper_id))) %>% 
  ungroup() #12 legume spp with home/away and genetic data, many of which have info for multiple gene regions

nativeInvasive <- plantData %>% 
  filter(compares_natinv==1) %>% #107 observations with native/non-native comparison
  filter(genbank==1) %>%  #88 of those that have genetic data
  group_by(genus, species, plant_status) %>% 
  summarise(gene_regions=length(unique(genetic_region)), studies=length(unique(paper_id))) %>% 
  ungroup() #17 introduced with genetic data that can be compared to at least one native species in the study

allSpp <- plantData %>% 
  select(genus_species, plant_status) %>% 
  filter(!is.na(genus_species)) %>% 
  unique() #790 legume species examined (including native/introduced spp and some duplicates for spp that were sampled both at home and away)

numSpp <- plantData %>% 
  select(genus_species) %>% 
  filter(!is.na(genus_species)) %>% 
  unique() #711 legume species examined (including native/introduced spp and some duplicates for spp that were sampled both at home and away)

introducedSpp <- allSpp %>% 
  filter(plant_status=='introduced') #103 introduced legume spp

anyGenetic <- plantData %>% 
  filter(genbank==1) %>% #2261 observations of those that have genetic data
  select(genus_species, plant_status, genbank, genetic_region) %>% 
  unique() %>% 
  pivot_wider(names_from=genetic_region, values_from=genbank, values_fill=NA) %>% 
  select(-genus_species) %>% 
  group_by(plant_status) %>% 
  summarise_all(funs(sum), na.rm=T) %>% 
  ungroup() #we have a lot of good gene regions that we can use for comparisons across all native vs introduced spp, many fewer for specific comparisons across individual spp in native and introduced ranges


# #Read in strain diversity info and add NA for anywhere that has a blank
# strainDiversity <- read.csv("Data/legume_strain diversity_meta analysis_strain sequences.csv") %>%
#   #change all 999 to NAs
#   mutate_all(~na_if(., 999))
# strainDiversity[strainDiversity==""]<-NA
# strainDiversity <- strainDiversity %>% 
#   select(paper_id, genus, species, plant_status, sample_country, sample_continent, genetic_region, strain, accession_number) %>% 
#   left_join(plantData)
