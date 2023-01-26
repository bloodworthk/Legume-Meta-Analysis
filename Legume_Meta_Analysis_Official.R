#### Legume Meta Analysis ####
#Authors: Kathryn Bloodworth, Smriti Pehim Limbu, Kim Komatsu


#### Install and Load Libraries #### 

#install.packages("ggplot2")
library(ggplot2)
#install.packages("lme4")
library(lme4)
#install.packages("lmerTest")
library(lmerTest)
#install.packages("nlme")
library(nlme)
#install.packages("tidyverse")
library(tidyverse)

#### Set Working Directory ####

#Bloodworth - Mac
setwd("/Volumes/GoogleDrive/My Drive/Projects/Legume-Meta Analysis")


#### Read in Data ####

#Read in paper information and add NA for anywhere that has a blank
Manuscript_Info<-read.csv("Data/legume_strain diversity_meta analysis_paper information.csv") %>% 
  #change all 999 to NAs
  mutate_all(~na_if(., 999))
Manuscript_Info[Manuscript_Info==""]<-NA

#Read in plant association information and add NA for anywhere that has a blank
Plant_Associations<-read.csv("Data/legume_strain diversity_meta analysis_plant associations_edited names_presence absence.csv") %>% 
  #change all 999 to NAs
  mutate_all(~na_if(., 999))
Plant_Associations[Plant_Associations==""]<-NA

#Read in plant data info and add NA for anywhere that has a blank
Plant_Data<-read.csv("Data/legume_strain diversity_meta analysis_plant data.csv") %>% 
  #change all 999 to NAs
  mutate_all(~na_if(., 999))
Plant_Data[Plant_Data==""]<-NA

#Read in strain diversity info and add NA for anywhere that has a blank
Strain_Diversity<-read.csv("Data/legume_strain diversity_meta analysis_strain sequences.csv") %>%
  #change all 999 to NAs
  mutate_all(~na_if(., 999))
Strain_Diversity[Strain_Diversity==""]<-NA

#Read in strain sequence information
accession_numbers<-read.csv("Data/legume_strain diversity_meta analysis_strain sequences.csv",stringsAsFactors = FALSE,na.strings="") %>% 
  #change all 999 to NAs
  mutate_all(~na_if(., 999))

#### Set ggplot2 theme ####
theme_set(theme_bw())
#Update ggplot2 theme - make box around the x-axis title size 30, vertically justify x-axis title to 0.35, Place a margin of 15 around the x-axis title.  Make the x-axis title size 30. For y-axis title, make the box size 30, put the writing at a 90 degree angle, and vertically justify the title to 0.5.  Add a margin of 15 and make the y-axis text size 25. Make the plot title size 30 and vertically justify it to 2.  Do not add any grid lines.  Do not add a legend title, and make the legend size 20
theme_update(axis.title.x=element_text(size=40, vjust=-0.35, margin=margin(t=15)),
             axis.text.x=element_text(size=40), axis.title.y=element_text(size=40, angle=90, vjust=0.5,
                                                                          margin=margin(r=15)), axis.text.y=element_text(size=40), plot.title =
               element_text(size=40, vjust=2), panel.grid.major=element_blank(), 
             panel.grid.minor=element_blank(),
             legend.text=element_text(size=30))

#### Basic Info about Study ####

## Determine number of papers assessed (327)
Manuscript_Info%>%
  #Removing papers that had comments saying it was a review or irrelevent or didnt exist
  filter(paper_id!=313 & paper_id!=312 & paper_id!=317) %>% 
  group_by(paper_id, author, year)%>%
  unique()%>%
  summarise(num_papers=length(paper_id))%>%
  ungroup() %>% 
  summarise(papers=length(paper_id))%>%
  summarise(num_papers=sum(papers))

## Determine number plant species (696)
Plant_Associations%>%
  group_by(genus_species)%>%
  unique()%>%
  summarise(num_species=length(genus_species))%>%
  ungroup()%>%
  summarise(num_species=length(genus_species))

## Determine number of countries assessed (85)
Plant_Associations%>%
  group_by(sample_country)%>%
  unique()%>%
  summarise(num_countries=length(sample_country))%>%
  ungroup()%>% 
  summarise(num_countries=length(sample_country))

#### Basic Data Analyses ####

#### Introduced vs. Native ####

# Create data frame with native vs. introduced species
Introduced_native_Number<-Plant_Associations%>%
  mutate(Plant_status=ifelse(plant_status=="introduced","non-native",ifelse(plant_status=="invasive","non-native",ifelse(plant_status=="native","native", ifelse(plant_status=="extinct","extinct",ifelse(plant_status=="introduced-native","non-native",ifelse(plant_status==" native","native",plant_status)))))))%>%
  group_by(Plant_status)%>%
  unique()%>%
  summarise(plant_status_counts=length(Plant_status))%>%
  ungroup()

# Make a graph for non-native vs native
Introduced_Native_Number_Graph<-Introduced_native_Number%>%
  na.omit()%>%
  filter(Plant_status!="extinct") 

#Make graph showing number of species and plant species status
ggplot(Introduced_Native_Number_Graph,aes(x=Plant_status,y=plant_status_counts,fill=Plant_status))+
  #Make a bar graph where the height of the bars is equal to the data (stat=identity) and you preserve the vertical position while adjusting the horizontal(position_dodge), and outline the bars with the color black.
  geom_bar(stat="identity", position=position_dodge(),color="black")+
  #make error bars using the Standard error from the mean and place them at 0.9 with a width of 0.2
  #Label the x-axis "Watershed"
  xlab("Plant Species Status")+
  #Label the y-axis "Species Richness"
  ylab("Number of Species")+
  expand_limits(y=1000)+
  #Fill the bar graphs with grey and white according and label them "Inside Fence" and "Outside Fence"
  scale_fill_manual(values=c("darkslategrey","cadetblue4"), labels=c("Native", "Non-Native"))  

#### Rhizobial Associates by plant status ####

#Make a new data frame called Wide_rhizobial_associations 
Rhizobial_Associations<-Plant_Associations%>%
  mutate(Plant_status=ifelse(plant_status=="introduced","non_native",ifelse(plant_status=="invasive","non_native",ifelse(plant_status=="native","native",ifelse(plant_status=="extinct","extinct",ifelse(plant_status=="introduced-native","non_native",ifelse(plant_status==" native","native",plant_status)))))))

#Sum across all columns with presence absence data for each row, giving us the number of assiciations for each species from each paper
Rhizobial_Associations$diversity_sum <- as.numeric(apply(Rhizobial_Associations[,27:312], 1, sum))

#Make a new condensed dataframe with only necessary information
Rhizobial_Associations_Condensed<-Rhizobial_Associations %>% 
  select(paper_id, genus_species, sample_country,diversity_sum,Plant_status, num_nodules) %>% 
  filter(Plant_status!="extinct") %>% 
  na.omit()

#Make a data frame calculating the mean, standard error, and number of observations for native vs. non-native to make graph
Average_rhizobial_symbionts_genus_sp<-Rhizobial_Associations_Condensed%>%
  group_by(Plant_status,genus_species)%>%
  summarize(symbionts_Std=sd(diversity_sum),symbionts_Mean=mean(diversity_sum),species_status_n=length(diversity_sum))%>%
  #Make a new column called "Richness_St_Error" and divide "Richness_Std" by the square root of "Richness_n"
  mutate(symbionts_St_Error=symbionts_Std/sqrt(species_status_n)) %>% 
  ungroup()

#Make a data frame calculating the mean, standard error, and number of observations for native vs. non-native to make graph
Average_rhizobial_symbionts_Graph<-Average_rhizobial_symbionts_genus_sp%>%
  group_by(Plant_status)%>%
  summarize(symbionts_Std=sd(symbionts_Mean),symbionts_Mean2=mean(symbionts_Mean),status_n=length(symbionts_Mean))%>%
  #Make a new column called "Richness_St_Error" and divide "Richness_Std" by the square root of "Richness_n"
  mutate(symbionts_St_Error=symbionts_Std/sqrt(status_n)) %>% 
  ungroup()

#Make graph showing number of introduced vs. native
ggplot(Average_rhizobial_symbionts_Graph,aes(x=Plant_status,y=symbionts_Mean2,fill=Plant_status))+
  #Make a bar graph where the height of the bars is equal to the data (stat=identity) and you preserve the vertical position while adjusting the horizontal(position_dodge), and outline the bars with the color black.
  geom_bar(stat="identity", position=position_dodge(),color="black")+
  geom_errorbar(aes(ymin=symbionts_Mean2-symbionts_St_Error,ymax=symbionts_Mean2+symbionts_St_Error),position=position_dodge(0.9),width=0.2)+
  #make error bars using the Standard error from the mean and place them at 0.9 with a width of 0.2
  #Label the x-axis "Watershed"
  xlab("Plant Species Status")+
  #Label the y-axis "Species Richness"
  ylab("Average Number of Symbionts")+
  #Fill the bar graphs with grey and white according and label them "Inside Fence" and "Outside Fence"
  scale_fill_manual(values=c("darkslategrey","cadetblue4"), labels=c("Native","Non-native"))  

#histogram for Rhizobial Associations
hist(Rhizobial_Associations_Condensed$diversity_sum)

#### K question -- which one? ####

#Run anova comparing # of symbionts by plant status accounting for how many studies per species 
summary (Mixed_Model_Rhiz_Status<- lmer(symbionts_Mean ~ Plant_status + (1 | species_status_n), data = Average_rhizobial_symbionts_genus_sp))
#anova accounting for type II error - for unbalanced data or interaction terms that are significant 
anova(Mixed_Model_Rhiz_Status,type = 2)

#random effect -- i think this is correct
summary(Mixed_Model_Rhiz_Status_NoMean_Random<- lmer(diversity_sum ~ Plant_status + (1 | genus_species) + (1 | num_nodules), data = Rhizobial_Associations_Condensed)) 
#anova accounting for type II error - for unbalanced data or interaction terms that are significant 
anova(Mixed_Model_Rhiz_Status_NoMean_Random,type = 2) 
