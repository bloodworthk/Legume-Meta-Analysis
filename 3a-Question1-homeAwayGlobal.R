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
  group_by(global_plant_status, clean_name) %>%  #averaging over papers, gene regions, etc -- plant species are our "replicates" for this question
  summarise(strain_richness=mean(strain_richness)) %>% 
  ungroup() %>% 
  filter(!is.na(strain_richness))


##### mixed model for native/invasive globally #####

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