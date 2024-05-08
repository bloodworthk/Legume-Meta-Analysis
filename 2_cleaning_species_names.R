#### Legume Meta Analysis ####
#Authors: Kathryn Bloodworth, Smriti Pehim Limbu, Kim Komatsu

library(Taxonstand)
#install.packages("WorldFlora")
library(WorldFlora)
library(tidyverse)

setwd('G:\\.shortcut-targets-by-id\\1w2OXIzBKQqFZ0BCeKP7C9pX36ViGDPBj\\Legume-Meta Analysis\\Data')


## In this script the species names are cleaned 
## the full species list is run through the TPL() function from Taxonstand
## to save time so the species list does not always have to be run through TPL()
## intermediate stages are saved to run the rest of the code

#### import data ####

full.list <- read.csv ("legume_old_names.csv") %>% 
  select(-genus_species) %>% 
  unite(col=genus_species, c(genus, species), sep=' ')

genus.species <- read.csv ("legume_old_names.csv") %>% 
  select(genus, species) %>% 
  unite(col=genus_species, c(genus, species), sep=' ')


#### finding recognized names ####

TPL.output <- TPL(genus.species$genus_species)

clean.list <- TPL.output %>% 
  rename(genus_species=Taxon) %>% 
  unite(col=new_name, c(New.Genus, New.Species), sep=' ') %>% 
  select(Plant.Name.Index, genus_species, new_name) %>% 
  #manually cleaning names that were not recognized in World Flora
  mutate(new_name=ifelse(genus_species=='Lespedeza bicolor for. alba', 'Lespedeza bicolor', new_name),
         Plant.Name.Index=ifelse(genus_species %in% c('Acacia mangium x auriculiformis','Acacia tortillis spp. heteracantha','Lespedeza maximowiezzi var. tomentella'), 'TRUE', Plant.Name.Index))

clean.list.annotated <- full.list %>% 
  left_join(clean.list)


#### write data to file ####

write.csv(clean.list.annotated, "legume_clean_names_updated.csv", row.names=F)

