#### Legume Meta Analysis ####
#Authors: Kim Komatsu, Kathryn Bloodworth, Smriti Pehim Limbu 

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
  rename(clean_name=new_name) %>% 
  select(-notes)

# getting native/invasive status for each study and legume species
plantStatus <- read.csv('legume_strain diversity_meta analysis_plant associations_edited names_presence absence.csv') %>% 
  left_join(cleanNames) %>%
  #rename(paper_id=ï..paper_id) %>% 
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
  #rename(paper_id=ï..paper_id) %>% 
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

#Figure 1 boxplot (using in MS)
ggplot(data=homeAwayAll,aes(x=global_plant_status,y=strain_richness,fill=global_plant_status))+
  geom_boxplot()+
  ylab('Rhizobial Strain Richness') + xlab('Global Plant Status') +
  scale_x_discrete(breaks=c('native', 'introduced'),
                   limits=c('native', 'introduced'),
                   labels=c('Native\n(n=186)', 'Non-native\n(n=164)')) +
  scale_fill_manual(values=c("#A79371","#E2E4DE"))+
  theme(legend.position="none")+
  expand_limits(y=c(0,65))+
  scale_y_continuous(breaks=c(0,20,40,60))


# Question 2A

paper_plant_status <- plantData %>% 
  filter(num_nodules>2) %>%
  filter(!is.na(clean_name),
         compares_homeaway==0)

hist(paper_plant_status$strain_richness)
hist(log(paper_plant_status$strain_richness))


model2a.1 <- lm(log(strain_richness)~plant_status*clean_name, paper_plant_status)
anova(model2a.1)
plot(model2a.1)
residuals <- resid(model2a.1)
qqnorm(residuals)
hist(residuals)

#Question 2A figure

df_summary_root <- emm_summary.below %>%
  group_by(Carbon, Mutualist) %>%
  dplyr::summarize(mean_value = mean(Root_biomass), sd=sd(Root_biomass), se = sqrt(var(Root_biomass)/length(Root_biomass)), lower = mean(Root_biomass) - qt(1- alpha/2, (n() - 1))*sd(Root_biomass)/sqrt(n()),upper = mean(Root_biomass) + qt(1- alpha/2, (n() - 1))*sd(Root_biomass)/sqrt(n()))
##root figure
x.expression <- expression(CO[2] ~ Treatment~(ppm))

q <- ggplot(emm_summary.below, aes(x=Carbon, y=emmean, group = Mutualist, color=Mutualist))+
  geom_line(linetype="dashed", size= 1.5) +
  geom_point(aes(shape=Mutualist), size = 5)+
  geom_ribbon(aes(ymin=lower.CL,ymax=upper.CL, fill=Mutualist), line_type = "dashed", alpha=0.1)  + theme_bw()+ scale_x_discrete(limits = c("400", "600", "1000"), expand = c(0.05, 0.05)) +
  labs(x = x.expression, y = "Root biomass (g)") + ggplot2::theme(axis.text.x = element_text(color = "black", size = 16, angle = 0, hjust = 1, vjust = 0, face = "plain"), axis.text.y = element_text(color = "black", size = 16, angle = 0, hjust = 1, vjust = 0, face = "plain"),  axis.title.x = element_text(color = "black", size = 16, angle = 0, hjust = .5, vjust = 0, face = "plain"), axis.title.y = element_text(color = "black", size = 16, angle = 90, hjust = .5, vjust = 1.5, face = "plain"))+
  ggplot2::theme(legend.text=element_text(size=14)) +
  ggplot2::theme (legend.title = element_text(size=16)) +
  ggplot2::theme(panel.border = element_rect(colour = "black", fill=NA, size=2)) +
  ggplot2::theme(axis.ticks.length=unit(0.15,"cm"), axis.ticks = element_line(size = 1)) +
  ggplot2::theme(panel.grid = element_blank(), panel.background = element_rect(fill="white")) +
  ggplot2::theme(legend.position = "top") + ggplot2::scale_color_manual(values = c("#FF7F50","#999999","#009E73", "#1f77b4"))+
  scale_shape_manual(values = c(15, 16, 17, 18))

q

#### Creating CSV to add to Supplemental Table 1 ####
## creating dataframe using homeAwayAll to make csv with paper information ##


CSV_Papers_Fig2 <- plantData %>% 
  filter(num_nodules>2) %>%
  filter(!is.na(clean_name),
         compares_homeaway==0) %>% 
  left_join(globalStatus) %>% 
  filter(!is.na(strain_richness)) %>% 
  select(paper_id,sample_country,clean_name) %>% 
  left_join(Paper_Information) %>% 
  unique()

#######note -- need to resave with species names and do the same in 4-Figure4 and update SupplementalTable Document in 3_dataAnalysis 

write.csv(CSV_Papers_Fig2,"Fig2_Papers.csv")
