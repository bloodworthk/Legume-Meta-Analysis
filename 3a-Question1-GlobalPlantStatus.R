#### Legume Meta Analysis ####
#Authors: Kim Komatsu, Kathryn Bloodworth, Smriti Pehim Limbu 
#### This script addresses question 1 in manuscript: "Do legume species that are globally native associate with a different number of rhizobial strains than those that are non-native somewhere?"

setwd('G:\\.shortcut-targets-by-id\\1w2OXIzBKQqFZ0BCeKP7C9pX36ViGDPBj\\Legume-Meta Analysis\\Data') #kim's wd

setwd('G:\\.shortcut-targets-by-id\\1w2OXIzBKQqFZ0BCeKP7C9pX36ViGDPBj\\Legume-Meta Analysis\\Data')

#Bloodworth - Mac
setwd("~/Documents/GitHub/Legume-Meta-Analysis/Legume-Meta Analysis_Data_etc/Data")

library(car)
library(lme4)
library(lmerTest)
library(ggsci)
library(ggrepel)
library(tidyverse)


theme_set(theme_bw())
theme_update(axis.title.x=element_text(size=30, vjust=-0.35, margin=margin(t=15)),
             axis.text.x=element_text(size=30), axis.title.y=element_text(size=30, angle=90, vjust=0.5,
                                                                          margin=margin(r=15)), axis.text.y=element_text(size=30), plot.title =
               element_text(size=30, vjust=2), panel.grid.major=element_blank(),
             panel.grid.minor=element_blank(), legend.title=element_blank(),
             legend.text=element_text(size=40))

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
  dplyr::rename(clean_name=new_name) %>% 
  dplyr::select(-notes)

# getting native/invasive status for each study and legume species
plantStatus <- read.csv('legume_strain diversity_meta analysis_plant associations_edited names_presence absence.csv') %>% 
  left_join(cleanNames) %>%
  # rename(paper_id=ï..paper_id) %>%
  dplyr::select(paper_id, clean_name, plant_status, sample_country) %>% 
  unique()

# removing species where legume species isn't clear, native/invasive status is unknown, and cultivated species
plantData <- read.csv("legume_strain diversity_meta analysis_plant data.csv") %>% 
  #change all 999 to NAs
  mutate(across(all_of(c('sites_sampled', 'num_nodules', 'strain_richness', 'proportion_novel_strains',
                         'proportion_familiar_strains', 'proportion_overlapping_strains', 'num_spp_compared')), ~na_if(., 999))) %>% 
  left_join(cleanNames) %>% 
  dplyr::select(-notes, -genus_species) %>% 
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
  # rename(paper_id=ï..paper_id) %>%
  left_join(cleanNames) %>% 
  mutate(global_plant_status=ifelse((exo_NA+exo_SA+exo_AU+exo_AS+exo_EU+exo_AF)>0,'introduced','native')) %>% 
  dplyr::select(paper_id, clean_name, global_plant_status) %>% 
  unique()

homeAwayAll <- plantData %>% 
  filter(num_nodules>2) %>%
  filter(!is.na(clean_name)) %>% 
  left_join(globalStatus) %>% 
  # group_by(global_plant_status, clean_name) %>%  #averaging over papers, gene regions, etc -- plant species are our "replicates" for this question
  # summarise(strain_richness=mean(strain_richness)) %>% 
  # ungroup() %>% 
  filter(!is.na(strain_richness)) %>% 
  filter(!(strain_source %in% c('999', 'rhizosphere_samples', 'NA', 'Look up Vincent (1970)')))

#### Determine number of countries, species, etc. ####

unique(sort(homeAwayAll$sample_country))

NonNative_Species<-homeAwayAll %>% 
  filter(global_plant_status=="introduced")

unique(sort(NonNative_Species$clean_name))
unique(sort(NonNative_Species$sample_country))

Native_Species<-homeAwayAll %>% 
  filter(global_plant_status=="native")

unique(sort(Native_Species$clean_name))
unique(sort(Native_Species$sample_country))

#local plant status
NonNative_Species_Local<-homeAwayAll %>% 
  filter(plant_status=="introduced")

unique(sort(NonNative_Species_Local$clean_name))
unique(sort(NonNative_Species_Local$sample_country))

Native_Species_Local<-homeAwayAll %>% 
  filter(plant_status=="native")

unique(sort(Native_Species_Local$clean_name))
unique(sort(Native_Species_Local$sample_country))

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
homeAwayAllModel <- lmer(strain_richness ~ global_plant_status + (1|clean_name) + (1|genetic_region),
                         data=homeAwayAll)
summary(homeAwayAllModel)
anova(homeAwayAllModel)
res <- resid(homeAwayAllModel)
plot(fitted(homeAwayAllModel), res)
qqnorm(res)
plot(density(res))

# homeAwayAllModel <- lmer(log(strain_richness) ~ global_plant_status + (1|clean_name) + (1|genetic_region),
#                         data=homeAwayAll)
# summary(homeAwayAllModel)
# anova(homeAwayAllModel)

### cannot use poisson(link='log') because data are overdispersed 
#x <- glmer(strain_richness ~ global_plant_status + (1|clean_name) + (1|genetic_region),
                          #data=homeAwayAll, family=poisson(link#='log'))

#check for over dispersion
#library(performance)
#check_overdispersion(x)

library(glmmTMB)
x=glmmTMB(strain_richness ~ global_plant_status + (1|clean_name) + (1|genetic_region),
        data = homeAwayAll, family = nbinom2(link = "log"))

summary(x)
Anova(x)
res <- resid(x)
plot(fitted(x), res)
qqnorm(res)
plot(density(res))


#### Manuscript Figure 2 ####
#Figure 2 boxplot (using in MS)
ggplot(data=homeAwayAll, aes(x=global_plant_status, y=strain_richness, fill=global_plant_status))+
  geom_boxplot(outlier.size=3)+
  ylab('Rhizobial Strain Richness') + xlab('Global Plant Status') +
  scale_x_discrete(breaks=c('native', 'introduced'),
                   limits=c('native', 'introduced'),
                   labels=c('Native\n(n=664)', 'Non-native\n(n=876)')) +
  scale_fill_manual(values=c("#A79371", "#E2E4DE"))+
  theme(legend.position="none") #+
  # expand_limits(y=c(0,65))+
  # scale_y_continuous(breaks=c(0,20,40,60))
#save at 1000x1000


