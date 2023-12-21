#### Legume Meta Analysis ####
#Authors: Kim Komatsu, Kathryn Bloodworth, Smriti Pehim Limbu 

setwd('G:\\.shortcut-targets-by-id\\1w2OXIzBKQqFZ0BCeKP7C9pX36ViGDPBj\\Legume-Meta Analysis\\Data') #kim's wd

library(car)
library(lme4)
library(lmerTest)
library(ggsci)
library(ggrepel)
library(tidyverse)

theme_set(theme_bw())
theme_update(axis.title.x=element_text(size=20, vjust=-0.35, margin=margin(t=15)), axis.text.x=element_text(size=16),
             axis.title.y=element_text(size=20, angle=90, vjust=0.5, margin=margin(r=15)), axis.text.y=element_text(size=16),
             plot.title = element_text(size=24, vjust=2),
             panel.grid.major=element_blank(), panel.grid.minor=element_blank(),
             legend.title=element_blank(), legend.text=element_text(size=20))

###bar graph summary statistics function
#barGraphStats(data=, variable="", byFactorNames=c(""))
barGraphStats <- function(data, variable, byFactorNames) {
  count <- length(byFactorNames)
  N <- aggregate(data[[variable]], data[byFactorNames], FUN=length)
  names(N)[1:count] <- byFactorNames
  names(N) <- sub("^x$", "N", names(N))
  mean <- aggregate(data[[variable]], data[byFactorNames], FUN=mean)
  names(mean)[1:count] <- byFactorNames
  names(mean) <- sub("^x$", "mean", names(mean))
  sd <- aggregate(data[[variable]], data[byFactorNames], FUN=sd)
  names(sd)[1:count] <- byFactorNames
  names(sd) <- sub("^x$", "sd", names(sd))
  preSummaryStats <- merge(N, mean, by=byFactorNames)
  finalSummaryStats <- merge(preSummaryStats, sd, by=byFactorNames)
  finalSummaryStats$se <- finalSummaryStats$sd / sqrt(finalSummaryStats$N)
  return(finalSummaryStats)
}  


##### data import #####

# getting clean legume names for standardizing across papers
cleanNames <- read.csv('legume_clean_names.csv') %>% 
  separate(old_genus_species,c("genus","species","subspecies","extra"),sep=" ") %>% 
  unite(col=genus_species, c(genus,species,subspecies,extra), sep='_',na.rm=TRUE) %>% 
  rename(clean_name=new_name) %>% 
  select(-notes)

# getting native/invasive status for each study and legume species
plantStatus <- read.csv('legume_strain diversity_meta analysis_plant associations_edited names_presence absence.csv') %>% 
  left_join(cleanNames) %>%
  rename(paper_id=ï..paper_id) %>% 
  select(paper_id, clean_name, plant_status, sample_country) %>% 
  unique()

# removing species where legume species isn't clear, native/invasive status is unknown, and cultivated species
plantData <- read.csv("legume_strain diversity_meta analysis_plant data.csv") %>% 
  #change all 999 to NAs
  mutate(across(all_of(c('sites_sampled', 'num_nodules', 'strain_richness', 'proportion_novel_strains',
                         'proportion_familiar_strains', 'proportion_overlapping_strains', 'num_spp_compared')), ~na_if(., 999))) %>% 
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
  # filter(!(sample_country %in% c("Japan-China","Japan_China","Kenya-Sudan","Czech_Republic-France-Georgia-Hungary-Italy-Romania-Spain","Senegal-Mauritania-Tunisia-Burundi","Malawi-Zambia-Kenya","Brazil-Venezuela"))) %>% 
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
  ungroup() #8 legume spp with home/away and genetic data, most of which have info for multiple gene regions

nativeInvasive <- plantData %>% 
  filter(compares_natinv==1) %>% #471 observations with native/non-native comparison
  filter(genbank==1) %>%  #325 of those that have genetic data
  group_by(clean_name, plant_status) %>% 
  summarise(gene_regions=length(unique(genetic_region)),
            studies=length(unique(paper_id))) %>% 
  ungroup() #24 introduced with genetic data that can be compared to at least one native species in the study

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
  filter(both_ranges==2) %>%  #35 species * region combos
  select(-native, -introduced)

#Relating number of nodules collected to strain richness. Below a cutoff of 3, the number of strains is nearly always equal to the number of nodules. At 3 and up, this is less of a problem. Subsetting data to 3+ nodules.
ggplot(data=plantData, aes(x=num_nodules, y=strain_richness)) +
  geom_jitter() +
  geom_abline(slope=1) +
  coord_cartesian(xlim=c(0,10), ylim=c(0,10))

length(unique(allSpp$clean_name)) #16 species that have the same gene region for both native and introduced ranges


##### mixed model for native/invasive using same gene region #####
nativeInvasiveGenetic <- plantData %>% 
  right_join(allSpp) %>% 
  select(paper_id, clean_name, plant_status, annual_perennial, growth_form, habitat_type, sample_country, 
         sample_continent, num_nodules, num_plants, genetic_region, strain_richness, cultivation.status) %>%
  mutate(genetic_region2=ifelse(genetic_region %in% c('16S', '16S-23S', '16S-ARDRA', '16S-IGS',
                                                      '16S-RFLP', '16S rDNA', '16S_23S_RFLP', 
                                                      '16S_ARDRA', '16S_BLAST', '16S_PCR_RFLP',
                                                      '16S_rDNA', '16S_RFLP', 'ARDRA', 'PCR',
                                                      'PCR_RFLP', '16S_PCR-RFLP', 'RFLP_16S',
                                                      'RFLP', 'RFLP_PCR'), '16S',
                        ifelse(genetic_region %in% c('23S', '23S_IVS', 'RFLP-23S'), '23S',
                        ifelse(genetic_region %in% c('AFLP', 'AFLP_Pst-A', 'AFLP_Pst-G',
                                                                    'AFLP_Pst-GC'), 'AFLP',
                        ifelse(genetic_region %in% c('BOX', 'BOX-AIR', 'box-PCR', 'BOX-PCR',
                                                     'Box_A1R-PCR', 'BOX_PCR', 'BoxA1R',
                                                     'BOXA1R', 'BOXA1R-PCR', 'BOXAIR'), 'BOX',
                        ifelse(genetic_region %in% c('CLUSTAL_W', 
                                                     'Cluster Analysis'), 'cluster',
                        ifelse(genetic_region %in% c('ERIC', 'ERIC-PCR', 'ERIC_PCR', 
                                                     'RFLP-ERIC'), 'ERIC',
                        ifelse(genetic_region %in% c('IGS', 'IGS_PCR-RFLP', 'RFLP-IGS', 
                                                     'IGSS', 'ITS', 'RFLP_ITS'), 'ITS',
                        ifelse(genetic_region %in% c('nif_KD', 'nifD', 'nifD-K', 'nifh',
                                                     'nifH', 'nifH-nifDK', 'nifHD',
                                                     'RFLP_nifH'), 'nif',
                        ifelse(genetic_region %in% c('nodBC', 'nodC', 'nodC-nodA',
                                                     'nodC-RFLP', 'nodA', 'nodA_PCR_RFLP',
                                                     'nodD', 'nodD1', 'nodD2', 'nodDAB',
                                                     'nodDF', 'nodF', 'nodY/K', 'RFLP_nodb3',
                                                     'RFLP_nodA', 'RFLP_nodb1', 'RFLP_nodb4',
                                                     'RFLP_nodb5', 'RFLP_nodC'), 'nod',
                       ifelse(genetic_region %in% c('rep-PCR', 'REP_PCR', 'REP-PCR', 
                                                    'rep_PCR', 'REP1R-I_REP2-I'), 'REP PCR',
                       ifelse(genetic_region %in% c('PCR-RAPD', 'RAPD'), 'RAPD',
                       ifelse(genetic_region %in% c('recA', 'recA-glnA-dnaK',
                                                    'recA, glnII', 'recA-glnII-atpD'), 'recA',
                       ifelse(genetic_region %in% c('glnA', 'glnB', 'glnII', 'gltA', 'gryB',
                                                    'gyrA', 'gyrB'), 'gln',
                              'other')))))))))))))) %>% 
  filter(genetic_region2=='16S') %>% 
  group_by(plant_status, clean_name) %>% 
  summarise(strain_richness=mean(strain_richness)) %>% 
  ungroup()


#normality
hist(nativeInvasiveGenetic$strain_richness)
qqPlot(nativeInvasiveGenetic$strain_richness)
shapiro.test(nativeInvasiveGenetic$strain_richness)
# W = 0.75401, p-value = 1.058e-05

hist(log(nativeInvasiveGenetic$strain_richness))
qqPlot(log(nativeInvasiveGenetic$strain_richness))
shapiro.test(log(nativeInvasiveGenetic$strain_richness))
# W = 0.9429, p-value = 0.1089


#model
nativeInvasiveGeneticModel <- lmer(log(strain_richness) ~ plant_status + (1|clean_name),
                                   data=nativeInvasiveGenetic)
summary(nativeInvasiveGeneticModel)
anova(nativeInvasiveGeneticModel)

#figure 2b
nativeInvasiveGeneticWide <- nativeInvasiveGenetic %>% 
  pivot_wider(names_from=plant_status, values_from=strain_richness) %>% 
  mutate(strain_diff=native-introduced)

ggplot(data=nativeInvasiveGenetic, 
       aes(x=plant_status, y=strain_richness, group=interaction(clean_name), color=interaction(clean_name), label=clean_name)) +
  geom_line(size=1.5) +
  geom_point(size=5, position=position_jitter(height=0.1, width=0)) +
  geom_text_repel(data = subset(nativeInvasiveGenetic, plant_status == "introduced"),
                  aes(label = clean_name, size=14), 
                  nudge_x = 0.1, direction = "y", hjust = 0,
                  max.overlaps=100
                  # box.padding = 0.5, point.padding = 0.5,
                  # segment.color = 'transparent'
                  ) +
  ylab('Strain Richness') + xlab('Plant Status') +
  scale_x_discrete(breaks=c('native', 'introduced'),
                   limits=c('native', 'introduced'),
                   labels=c('Native\n(N=15)', 'Introduced\n(N=15)'),
                   expand = expansion(mult = 2.5)) +
  scale_color_simpsons() +
  theme(legend.position='none') +
  coord_cartesian(xlim=c(1.35, 1.65))
# ggsave('C:\\Users\\kjkomatsu\\UNCG\\Kathryn Bloodworth - Invasive Legume Meta-Analysis\\Figures\\Fig2b_HomeAway_acrossStudies.png', width=12, height=12, units='in', dpi=300, bg='white')



#### Trends in strain richness by gene region ####
geneRegion <- plantData %>% 
  filter(num_nodules>2) %>%
  filter(!is.na(clean_name),
         compares_homeaway==0) %>% 
  mutate(genetic_region2=ifelse(genetic_region %in% c('16S', '16S-23S', '16S-ARDRA', '16S-IGS',
                                                      '16S-RFLP', '16S rDNA', '16S_23S_RFLP', 
                                                      '16S_ARDRA', '16S_BLAST', '16S_PCR_RFLP',
                                                      '16S_rDNA', '16S_RFLP', 'ARDRA', 'PCR',
                                                      'PCR_RFLP', '16S_PCR-RFLP', 'RFLP_16S',
                                                      'RFLP', 'RFLP_PCR'), '16S',
                          ifelse(genetic_region %in% c('23S', '23S_IVS', 'RFLP-23S'), '23S',
                          ifelse(genetic_region %in% c('AFLP', 'AFLP_Pst-A', 'AFLP_Pst-G',
                                                       'AFLP_Pst-GC'), 'AFLP',
                          ifelse(genetic_region %in% c('BOX', 'BOX-AIR', 'box-PCR', 'BOX-PCR',
                                                       'Box_A1R-PCR', 'BOX_PCR', 'BoxA1R',
                                                       'BOXA1R', 'BOXA1R-PCR', 
                                                       'BOXAIR'), 'BOX',
                          ifelse(genetic_region %in% c('CLUSTAL_W', 
                                                       'Cluster Analysis'), 'cluster',
                          ifelse(genetic_region %in% c('ERIC', 'ERIC-PCR', 'ERIC_PCR', 
                                                       'RFLP-ERIC'), 'ERIC',
                          ifelse(genetic_region %in% c('IGS', 'IGS_PCR-RFLP', 'RFLP-IGS',
                                                       'IGSS', 'ITS', 'RFLP_ITS'), 'ITS',
                          ifelse(genetic_region %in% c('nif_KD', 'nifD', 'nifD-K', 'nifh',
                                                       'nifH', 'nifH-nifDK', 'nifHD',
                                                       'RFLP_nifH'), 'nif',
                          ifelse(genetic_region %in% c('nodBC', 'nodC', 'nodC-nodA', 
                                                       'nodC-RFLP', 'nodA', 'nodA_PCR_RFLP',
                                                       'nodD', 'nodD1', 'nodD2', 'nodDAB',
                                                       'nodDF', 'nodF', 'nodY/K', 'RFLP_nodb3',
                                                       'RFLP_nodA', 'RFLP_nodb1', 'RFLP_nodb4',
                                                       'RFLP_nodb5', 'RFLP_nodC'), 'nod',
                          ifelse(genetic_region %in% c('rep-PCR', 'REP_PCR', 'REP-PCR', 
                                                       'rep_PCR', 'REP1R-I_REP2-I'), 'REP PCR',
                          ifelse(genetic_region %in% c('PCR-RAPD', 'RAPD'), 'RAPD',
                          ifelse(genetic_region %in% c('recA', 'recA-glnA-dnaK', 
                                                       'recA, glnII', 
                                                       'recA-glnII-atpD'), 'recA',
                          ifelse(genetic_region %in% c('glnA', 'glnB', 'glnII', 'gltA', 'gryB',
                                                       'gyrA', 'gyrB'), 'gln',
                                 'other')))))))))))))) %>% 
  select(genetic_region, genetic_region2) %>% 
  unique()
#gene regions do differ from each other, with the fingerprinting methods having more "diversity" than the genotyping methods

geneticData <- plantData %>% 
  right_join(allSpp) %>% 
  select(paper_id, clean_name, plant_status, annual_perennial, growth_form, habitat_type, sample_country, 
         sample_continent, num_nodules, num_plants, genetic_region, strain_richness, cultivation.status) %>% 
  group_by(paper_id, clean_name, plant_status, genetic_region) %>% 
  summarise(strain_richness=mean(strain_richness)) %>% 
  ungroup() %>% 
  filter(genetic_region=='16S') %>% 
  mutate(paper_id=as.character(paper_id))

accessionNumbers <- read.csv('legume_strain diversity_meta analysis_strain sequences.csv') %>% 
  rename(paper_id=ï..paper_id) %>% 
  select(paper_id, author, year, genus_species, genetic_region, strain, accession_number, notes) %>% 
  separate(genus_species,c("genus","species","subspecies","extra"),sep=" ") %>% 
  unite(col=genus_species, c(genus,species,subspecies,extra), sep='_', na.rm=TRUE) %>% 
  left_join(cleanNames) %>% 
  right_join(geneticData) %>% 
  filter(!is.na(accession_number), accession_number!='')

uniqueAccessionNumbers <- accessionNumbers %>% 
  select(accession_number) %>% 
  unique()
# write.csv(uniqueAccessionNumbers, 'accessionNumbers_homeAway.csv', row.names=F)
  




#### Global plant status ####

globalStatus <- read.csv('legume_strain diversity_meta analysis_plant associations_edited names_presence absence.csv') %>% 
  rename(paper_id=ï..paper_id) %>% 
  left_join(cleanNames) %>% 
  mutate(global_plant_status=ifelse((exo_NA+exo_SA+exo_AU+exo_AS+exo_EU+exo_AF)>0,'introduced','native')) %>% 
  select(paper_id, clean_name, global_plant_status) %>% 
  unique()
  
homeAwayAll <- plantData %>% 
  filter(num_nodules>2) %>%
  filter(!is.na(clean_name),
         compares_homeaway==0) %>% 
  left_join(globalStatus) %>% 
  left_join(geneRegion) %>% 
  group_by(global_plant_status, clean_name) %>%  #averaging over papers, gene regions, etc -- plant species are our "replicates" for this question
  summarise(strain_richness=mean(strain_richness)) %>% 
  ungroup() %>% 
  filter(!is.na(strain_richness))
  

#normality
hist(homeAwayAll$strain_richness)
qqPlot(homeAwayAll$strain_richness)
shapiro.test(homeAwayAll$strain_richness)
# W = 0.70895, p-value < 2.2e-16

hist(log(homeAwayAll$strain_richness))
qqPlot(log(homeAwayAll$strain_richness))
shapiro.test(log(homeAwayAll$strain_richness))
# W = 0.9773, p-value = 2.277e-05


#model - all species regardless of whether home and away pairing possible
homeAwayAllModel <- lm(log(strain_richness) ~ global_plant_status,
                         data=homeAwayAll)
summary(homeAwayAllModel)
anova(homeAwayAllModel)

#figure 1
ggplot(data=barGraphStats(data=homeAwayAll, variable="strain_richness", byFactorNames=c("global_plant_status")), aes(x=global_plant_status, y=mean)) +
  geom_bar(stat='identity', fill='white', color='black') +
  geom_errorbar(aes(ymin=mean-se, ymax=mean+se, width=0.2)) +
  ylab('Strain Richness') + xlab('Global Plant Status') +
  scale_x_discrete(breaks=c('native', 'introduced'),
                   limits=c('native', 'introduced'),
                   labels=c('Native\n(N=186)', 'Introduced\n(N=164)')) +
  annotate('text', x=1, y=8, label='a', size=6) +
  annotate('text', x=2, y=10.5, label='b', size=6) +
  theme(legend.position='none')
# ggsave('C:\\Users\\kjkomatsu\\UNCG\\Kathryn Bloodworth - Invasive Legume Meta-Analysis\\Figures\\Fig1_HomeAway_allSpecies.png', width=6, height=6, units='in', dpi=300, bg='white')





numSpp <- plantData %>% 
  select(clean_name) %>% 
  filter(!is.na(clean_name)) %>% 
  unique() #523 legume species examined (including native/introduced spp and some duplicates for spp that were sampled both at home and away)

anyGenetic <- plantData %>% 
  filter(genbank==1) %>% #1717 observations of those that have genetic data
  select(clean_name, plant_status, genbank, genetic_region) %>% 
  unique() %>% 
  pivot_wider(names_from=genetic_region, values_from=genbank, values_fill=NA) %>% 
  select(-clean_name) %>% 
  group_by(plant_status) %>% 
  summarise_all(funs(sum), na.rm=T) %>% 
  ungroup() #we have a lot of good gene regions that we can use for comparisons across all native vs introduced spp, many fewer for specific comparisons across individual spp in native and introduced ranges
