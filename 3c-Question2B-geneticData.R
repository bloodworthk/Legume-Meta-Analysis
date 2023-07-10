#### Legume Meta Analysis ####
#Authors: Kim Komatsu, Kathryn Bloodworth, Smriti Pehim Limbu 

setwd('G:\\.shortcut-targets-by-id\\1w2OXIzBKQqFZ0BCeKP7C9pX36ViGDPBj\\Legume-Meta Analysis\\Data') #kim's wd

library(lme4)
library(lmerTest)
library(tidyverse)

theme_set(theme_bw())
theme_update(axis.title.x=element_text(size=20, vjust=-0.35, margin=margin(t=15)), axis.text.x=element_text(size=16),
             axis.title.y=element_text(size=20, angle=90, vjust=0.5, margin=margin(r=15)), axis.text.y=element_text(size=16),
             plot.title = element_text(size=24, vjust=2),
             panel.grid.major=element_blank(), panel.grid.minor=element_blank(),
             legend.title=element_blank(), legend.text=element_text(size=20))

##### data import #####
cleanNames <- read.csv('legume_clean_names.csv') %>% 
  separate(old_genus_species,c("genus","species","subspecies","extra"),sep=" ") %>% 
  unite(col=genus_species, c(genus,species,subspecies,extra), sep='_',na.rm=TRUE) %>% 
  rename(clean_name=new_name) %>% 
  select(-notes)

plantStatus <- read.csv('legume_strain diversity_meta analysis_plant associations_edited names_presence absence.csv') %>% 
  left_join(cleanNames) %>% 
  select(paper_id, clean_name, plant_status)

plantData <- read.csv("legume_strain diversity_meta analysis_plant data.csv") %>% 
  #change all 999 to NAs
  mutate_all(~na_if(., 999)) %>% 
  left_join(cleanNames) %>% 
  select(-notes, -genus_species) %>% 
  full_join(plantStatus) %>% 
  mutate(plant_status=ifelse(plant_status=='invasive', 'introduced', plant_status)) %>% 
  filter(cultivation.status!='row crop',
         plant_status!="999",
         !is.na(plant_status),
         species!="sp",
         species!="spp",
         species!="sp.",
         species!="",
         paper_id!=142,
         plant_status!="extinct") %>% 
  filter(!(sample_country %in% c("Japan-China","Japan_China","Kenya-Sudan","Czech_Republic-France-Georgia-Hungary-Italy-Romania-Spain","Senegal-Mauritania-Tunisia-Burundi","Malawi-Zambia-Kenya","Brazil-Venezuela"))) %>% 
  filter(!is.na(strain_richness)) %>% 
  mutate(strain_richness=as.numeric(strain_richness)) %>% #NA introduced for blanks
  mutate(genetic_region=(ifelse(genetic_region %in% c('16S rRNA', 'PCR-RFLP', 'RFLP-16S'), '16S',
                         ifelse(genetic_region=='IGS-RFLP', 'IGS',
                         genetic_region))))

plantData[plantData==""]<-NA


##### determine what species have home/away or native/non-native comparisons and genetic data #####
homeAway <- plantData %>% 
  filter(compares_homeaway==1) %>% #39 observations of home and away
  filter(genbank==1) %>% #20 observations of those that have genetic data
  group_by(clean_name, plant_status) %>% 
  summarise(gene_regions=length(unique(genetic_region)),
            studies=length(unique(paper_id))) %>% 
  ungroup() #3 legume spp with home/away and genetic data, most of which have info for multiple gene regions

nativeInvasive <- plantData %>% 
  filter(compares_natinv==1) %>% #99 observations with native/non-native comparison
  filter(genbank==1) %>%  #84 of those that have genetic data
  group_by(clean_name, plant_status) %>% 
  summarise(gene_regions=length(unique(genetic_region)),
            studies=length(unique(paper_id))) %>% 
  ungroup() #12 introduced with genetic data that can be compared to at least one native species in the study

allSpp <- plantData %>% 
  filter(num_nodules>2) %>%
  select(clean_name, plant_status, compares_homeaway, compares_natinv, genetic_region) %>% 
  filter(!is.na(clean_name),
         compares_homeaway==0) %>% 
  select(clean_name, plant_status, genetic_region) %>% 
  unique() %>%
  mutate(yes=1) %>% 
  pivot_wider(names_from=plant_status, values_from=yes) %>% #figure out which species have the same genetic region from native and introduced ranges
  mutate(both_ranges=native+introduced) %>% 
  filter(both_ranges==2) %>%  #41 species * region combos
  select(-native, -introduced)

#Relating number of nodules collected to strain richness. Below a cutoff of 3, the number of strains is nearly always equal to the number of nodules. At 3 and up, this is less of a problem. Subsetting data to 5+ nodules.
ggplot(data=plantData, aes(x=num_nodules, y=strain_richness)) +
  geom_jitter() +
  geom_abline(slope=1) +
  coord_cartesian(xlim=c(0,10), ylim=c(0,10))

length(unique(allSpp$clean_name)) #19 species that have the same gene region for both native and introduced ranges


##### mixed model for native/invasive using same gene region #####
nativeInvasiveGenetic <- plantData %>% 
  right_join(allSpp) %>% 
  select(paper_id, clean_name, plant_status, annual_perennial, growth_form, habitat_type, sample_country, sample_continent, num_nodules, num_plants, genetic_region, strain_richness, cultivation.status)


#model
nativeInvasiveGeneticModel <- lmer(strain_richness ~ plant_status*annual_perennial +
                                                     (1|num_nodules) +  (1|genetic_region), 
                                   data=nativeInvasiveGenetic)
summary(nativeInvasiveGeneticModel)
anova(nativeInvasiveGeneticModel)

#figure 2b
nativeInvasiveGeneticSummary <- nativeInvasiveGenetic %>% 
  group_by(plant_status, clean_name, genetic_region) %>% 
  summarise(strain_richness=mean(strain_richness)) %>% 
  ungroup()

nativeInvasiveGeneticWide <- nativeInvasiveGeneticSummary %>% 
  pivot_wider(names_from=plant_status, values_from=strain_richness) %>% 
  mutate(strain_diff=native-introduced)

ggplot(data=nativeInvasiveGeneticSummary, 
       aes(x=plant_status, y=strain_richness, group=interaction(clean_name, genetic_region), color=interaction(clean_name, genetic_region))) +
  geom_point() +
  geom_line() +
  ylab('Strain Richness') + xlab('Plant Status') +
  scale_x_discrete(breaks=c('native', 'introduced'),
                   limits=c('native', 'introduced'),
                   labels=c('Native', 'Introduced')) +
  theme(legend.position='none') +
  coord_cartesian(xlim=c(1.35, 1.65))

ggplot(data=nativeInvasiveGeneticSummary, 
       aes(x=plant_status, y=strain_richness, color=clean_name)) +
  geom_point() +
  geom_line() +
  facet_wrap(~genetic_region) +
  ylab('Strain Richness') + xlab('Plant Status') +
  scale_x_discrete(breaks=c('native', 'introduced'),
                   limits=c('native', 'introduced'),
                   labels=c('Native', 'Introduced')) +
  theme(legend.position='none') +
  coord_cartesian(xlim=c(1.35, 1.65))





numSpp <- plantData %>% 
  select(clean_name) %>% 
  filter(!is.na(clean_name)) %>% 
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
