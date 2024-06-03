#### Legume Meta Analysis ####
#Authors: Kim Komatsu, Kathryn Bloodworth, Smriti Pehim Limbu 
#### This script addresses question 2 in manuscript: Are legumes utilizing the same number of rhizobial strains in their novel range as they do in their native range, either within (2a) or across studies (2b)? 
#### This script also create supplemental figure 1

setwd('G:\\.shortcut-targets-by-id\\1w2OXIzBKQqFZ0BCeKP7C9pX36ViGDPBj\\Legume-Meta Analysis\\Data') #kim's wd
#Bloodworth - Mac
setwd("~/Documents/GitHub/Legume-Meta-Analysis/Legume-Meta Analysis_Data_etc/Data")

library(car)
library(lme4)
library(lmerTest)
library(ggsci)
library(ggrepel)
library(tidyverse)

theme_set(theme_bw())
theme_update(axis.title.x=element_text(size=22, vjust=-0.35, margin=margin(t=15)), axis.text.x=element_text(size=18),
             axis.title.y=element_text(size=22, angle=90, vjust=0.5, margin=margin(r=15)), axis.text.y=element_text(size=18),
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
  rename(paper_id=paper_id) %>% 
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

#### Supplemental Figure S1 ####
#Relating number of nodules collected to strain richness. Below a cutoff of 3, the number of strains is nearly always equal to the number of nodules. At 3 and up, this is less of a problem. Subsetting data to 3+ nodules.
ggplot(data=plantData, aes(x=num_nodules, y=strain_richness)) +
  geom_jitter() +
  geom_abline(slope=1) +
  coord_cartesian(xlim=c(0,10), ylim=c(0,10)) +
  xlab('Number of Nodules Sampled') + ylab('Strain Richness')
# ggsave('C:\\Users\\kjkomatsu\\UNCG\\Kathryn Bloodworth - Invasive Legume Meta-Analysis\\Figures\\SuppFig1_numNodules.png', width=6, height=6, units='in', dpi=300, bg='white')

#### Question 2A data and model ####

paper_plant_status <- plantData %>%
  filter(num_nodules>2) %>%
  filter(!is.na(clean_name),
         compares_homeaway==1) %>% group_by(clean_name) %>% filter(all(c("native", "introduced") %in% plant_status))

str(paper_plant_status)

hist(paper_plant_status$strain_richness)
hist(log(paper_plant_status$strain_richness))


model2a <- lm(log(strain_richness)~plant_status*clean_name, paper_plant_status)
anova(model2a)
plot(model2a)
residuals <- resid(model2a)
qqnorm(residuals)
hist(residuals)

alpha <- 0.05
df_summary <- paper_plant_status %>% group_by(clean_name,plant_status) %>% dplyr::summarize(mean_value = mean(strain_richness))

df_summary$plant_status[df_summary$plant_status == "introduced"] <- "Non-native"
df_summary$plant_status[df_summary$plant_status == "native"] <- "Native"

df_summary1 <- paper_plant_status %>% group_by(clean_name,plant_status, proportion_novel_strains, proportion_familiar_strains) %>% dplyr::summarize(mean_value = mean(strain_richness))


#### Manuscript Figure 3b ####
ggplot(df_summary, aes(x=plant_status, y=mean_value, group = interaction(clean_name), color= interaction(clean_name))) +
  geom_line( size= 1.5) +
  geom_point(aes(shape=clean_name), size = 5) +
  geom_text_repel(data = subset(df_summary, plant_status == "Non-native"),
                  aes(label = clean_name, size=14), 
                  nudge_x = 0.1, direction = "y", hjust = 0,
                  max.overlaps=100,fontface = "italic",
                  size=8
                  # box.padding = 0.5, point.padding = 0.5,
                  # segment.color = 'transparent'  
  ) +
  theme(legend.position = "none")+
  ylab('Rhizobial Strain Richness') + xlab('Local Plant Status') +
  scale_x_discrete(breaks=c('Native', 'Non-native'),
                   limits=c('Native', 'Non-native'),
                   labels=c('Native\n(n=8)', 'Non-native\n(n=13)')
                   # , expand = expansion(mult = 2.5)
  )+
  scale_color_manual(values=c("#414535","#8C4843","#189993","#AB92BF"))+
  scale_shape_manual(values=c(16,17,15,18))+
  coord_cartesian(xlim=c(1.5, 2))+
  expand_limits(y=c(0,15),by=5)
#Save at 1000 x 1000

#### Question 2B data and model####

##### determine what species have home/away or native/non-native comparisons and genetic data #####
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

length(unique(allSpp$clean_name)) #16 species that have the same gene region for both native and introduced ranges


##### mixed model for native/invasive using same gene region #####
nativeInvasiveGenetic <- plantData %>% 
  right_join(allSpp) %>% 
  filter(num_nodules>2) %>% 
  select(paper_id, clean_name, plant_status, annual_perennial, growth_form, habitat_type, sample_country,sample_continent, num_nodules, num_plants, genetic_region, strain_richness, cultivation.status) %>%
  mutate(genetic_region2=ifelse(genetic_region %in% c('16S', '16S-23S', '16S-ARDRA', '16S-IGS','16S-RFLP', '16S rDNA', '16S_23S_RFLP', '16S_ARDRA', '16S_BLAST', '16S_PCR_RFLP', '16S_rDNA', '16S_RFLP', 'ARDRA', 'PCR','PCR_RFLP', '16S_PCR-RFLP', 'RFLP_16S','RFLP', 'RFLP_PCR'), '16S',ifelse(genetic_region %in% c('23S', '23S_IVS', 'RFLP-23S'), '23S',ifelse(genetic_region %in% c('AFLP', 'AFLP_Pst-A', 'AFLP_Pst-G', 'AFLP_Pst-GC'), 'AFLP',ifelse(genetic_region %in% c('BOX', 'BOX-AIR', 'box-PCR', 'BOX-PCR','Box_A1R-PCR', 'BOX_PCR', 'BoxA1R','BOXA1R', 'BOXA1R-PCR', 'BOXAIR'), 'BOX',ifelse(genetic_region %in% c('CLUSTAL_W','Cluster Analysis'), 'cluster',ifelse(genetic_region %in% c('ERIC', 'ERIC-PCR', 'ERIC_PCR','RFLP-ERIC'), 'ERIC',ifelse(genetic_region %in% c('IGS', 'IGS_PCR-RFLP', 'RFLP-IGS','IGSS', 'ITS', 'RFLP_ITS'), 'ITS',ifelse(genetic_region %in% c('nif_KD', 'nifD', 'nifD-K', 'nifh','nifH', 'nifH-nifDK', 'nifHD','RFLP_nifH'), 'nif',ifelse(genetic_region %in% c('nodBC', 'nodC', 'nodC-nodA','nodC-RFLP', 'nodA', 'nodA_PCR_RFLP','nodD', 'nodD1', 'nodD2', 'nodDAB','nodDF', 'nodF', 'nodY/K', 'RFLP_nodb3','RFLP_nodA', 'RFLP_nodb1', 'RFLP_nodb4','RFLP_nodb5', 'RFLP_nodC'), 'nod',ifelse(genetic_region %in% c('rep-PCR', 'REP_PCR', 'REP-PCR', 'rep_PCR', 'REP1R-I_REP2-I'), 'REP PCR',ifelse(genetic_region %in% c('PCR-RAPD', 'RAPD'), 'RAPD',ifelse(genetic_region %in% c('recA', 'recA-glnA-dnaK','recA, glnII', 'recA-glnII-atpD'), 'recA',ifelse(genetic_region %in% c('glnA', 'glnB', 'glnII', 'gltA', 'gryB','gyrA', 'gyrB'), 'gln','other')))))))))))))) %>% 
  # filter(genetic_region2=='16S') %>%
  group_by(plant_status, clean_name) %>% #each species tested across different studies is a replicate (averaged across gene regions for species with multiple)
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

##### Manuscript Figure 4b ####
ggplot(data=nativeInvasiveGenetic, 
       aes(x=plant_status, y=strain_richness, group=interaction(clean_name), color=interaction(clean_name), label=clean_name)) +
  geom_line(size=1.5) +
  geom_point(size=5.5, position=position_jitter(height=0.1, width=0)) +
  geom_text_repel(data = subset(nativeInvasiveGenetic, plant_status == "introduced"),
                  aes(label = clean_name, size=14), 
                  nudge_x = 0.1, direction = "y", hjust = 0,
                  max.overlaps=100,fontface = "italic",
                  size=8
                  # box.padding = 0.5, point.padding = 0.5,
                  # segment.color = 'transparent'
                  ) +
  ylab('Rhizobial Strain Richness') + xlab('Local Plant Status') +
  scale_x_discrete(breaks=c('native', 'introduced'),
                   limits=c('native', 'introduced'),
                   labels=c('Native\n(n=16)', 'Non-native\n(n=16)')
                   # , expand = expansion(mult = 2.5)
                   ) +
  scale_color_simpsons() +
  theme(legend.position='none') +
  coord_cartesian(xlim=c(1.5, 2))+
  expand_limits(y=30)
#save at 1200 x 1000
