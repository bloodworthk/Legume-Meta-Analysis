#### Legume Meta Analysis ####
#Authors: Kathryn Bloodworth, Smriti Pehim Limbu, Kim Komatsu

library(Taxonstand)
library(WorldFlora)
library(tidyverse)

setwd('G:\\.shortcut-targets-by-id\\1w2OXIzBKQqFZ0BCeKP7C9pX36ViGDPBj\\Legume-Meta Analysis\\Data')


## In this script the species names are cleaned so that they are standardized
## between experiments and to be merged with TRY.
## Two files are read in - the full species list created in the Abundance_merge.R
## file and the list of the corrected names from version 1.0 (note: no code for this
## because a lot of the cleaning for this list was done by hand)
## the full species list is run through the TPL() function from Taxonstand
## to save time so the species list does not always have to be run through TPL()
## intermediate stages are saved to run the rest of the code

##this file creates 1) a full list of species names that are cleaned to merge with TRY, FullList_Nov2021, 2) a list of new species that need categorical traits for newsp2021.csv, and 3) a list of family, and whether it is a tree species for all cleaned names. This needs to then be filled by hand for those that are new and unclear species_families_trees_2021_toclean.csv

#### import data ####

full.list <- read.csv ("legume_old_names.csv") %>% 
  select(-genus_species) %>% 
  unite(col=genus_species, c(genus, species), sep=' ')

genus.species <- full.list %>% 
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

write.csv(clean.list.annotated, "legume_clean_names.csv", row.names=F)

