#### Legume Meta Analysis ####
#Authors: Kathryn Bloodworth, Smriti Pehim Limbu, Kim Komatsu


#### Install and Load Libraries #### 

#install.packages("ggplot2")
library(ggplot2)
#install.packages("lmerTest")
library(lmerTest)
#install.packages("tidyverse")
library(tidyverse)

#### Set Working Directory ####

#Bloodworth - Mac
setwd("~/Documents/GitHub/Legume-Meta-Analysis/Legume-Meta Analysis_Data_etc")

#### Read in Data ####

#Read in paper information and add NA for anywhere that has a blank
Manuscript_Info<-read.csv("Data/legume_strain diversity_meta analysis_paper information.csv") %>% 
  #change all 999 to NAs
  mutate_all(~na_if(., 999))
Manuscript_Info[Manuscript_Info==""]<-NA

#Read in clean species file and remove extra rows and place NAs anywhere that has no information
Clean_Species<-read.csv("Data/legume_clean_names.csv") %>% 
  filter(old_genus_species!="") %>% 
  separate(old_genus_species,c("genus","species","subspecies","extra"),sep=" ") %>% 
  unite(col=genus_species, c(genus,species,subspecies,extra), sep='_',na.rm=TRUE) %>% 
  mutate(genus_species=ifelse(genus_species=="Acacia_mangium_x_auriculiformis","Acacia_mangium x auriculiformis",ifelse(genus_species=="Acacia_nilotica_adansonii","Acacia_nilotica adansonii",ifelse(genus_species=="Acacia_nilotica_adansonii","Acacia_nilotica adansonii",ifelse(genus_species== "Acacia_nilotica_tomentosa", "Acacia_nilotica tomentosa",ifelse(genus_species=="Acacia_tortillis_spp._heteracantha","Acacia_tortillis spp. heteracantha",ifelse(genus_species=="Camptotheca_acuminata_Decne","Camptotheca_acuminata Decne",ifelse(genus_species=="Hedysarum_confertum_Desf.","Hedysarum_confertum Desf.",ifelse(genus_species=="Hedysarum_coronarium_L._Leguminosae","Hedysarum_coronarium L. Leguminosae",ifelse(genus_species=="Lespedeza_bicolor_for._alba","Lespedeza_bicolor for. alba",ifelse(genus_species=="Lespedeza_maximowiezzi_var._tomentella","Lespedeza_maximowiezzi var. tomentella",ifelse(genus_species=="Oxytropis_ochrocephala_Bunge","Oxytropis_ochrocephala Bunge",ifelse(genus_species=="Robinia_pseudoacacia_L.","Robinia_pseudoacacia L.",ifelse(genus_species=="phaseolus_vulgaris","Phaseolus_vulgaris",ifelse(genus_species=="Stylosanthes_ viscosa","Stylosanthes_viscosa",ifelse(genus_species=="Teline__canariensis","Teline_canariensis",ifelse(genus_species=="Teline__stenopetala","Teline_stenopetala",ifelse(genus_species=="Desmodium__multiflorum","Desmodium_multiflorum",genus_species)))))))))))))))))) %>% 
  rename(Clean_Species_notes=notes)
Clean_Species[Clean_Species==""]<-NA

#Read in plant association information and add NA for anywhere that has a blank
Plant_Associations<-read.csv("Data/legume_strain diversity_meta analysis_plant associations_edited names_presence absence.csv") %>% 
  #change all 999 to NAs
  mutate_all(~na_if(., 999)) %>% 
  left_join(Clean_Species) %>% 
  select(genus_species,new_name,1:316)
Plant_Associations[Plant_Associations==""]<-NA

#Read in plant data info and add NA for anywhere that has a blank
Plant_Data<-read.csv("Data/legume_strain diversity_meta analysis_plant data.csv") %>% 
  #change all 999 to NAs
  mutate_all(~na_if(., 999)) %>% 
  left_join(Clean_Species) %>% 
  select(genus_species,new_name,1:29)
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
  #Removing papers that had comments saying it was a review or irrelevant or didnt exist
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
  ylab("Number of Observations")+
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

#### Rhizobial Associates by Growth Form  ####

Plant_Data_simplified<-Plant_Data %>% 
  select(paper_id,author,year,genus_species,annual_perennial,growth_form)

#Make a new data frame called rhizobial_associates_growth_form
Rhizobial_Associations_growth_form<-Plant_Associations%>%
  mutate(Plant_status=ifelse(plant_status=="introduced","non_native",ifelse(plant_status=="invasive","non_native",ifelse(plant_status=="native","native",ifelse(plant_status=="extinct","extinct",ifelse(plant_status=="introduced-native","non_native",ifelse(plant_status==" native","native",plant_status))))))) %>% 
  left_join(Plant_Data_simplified)

#Sum across all columns with presence absence data for each row, giving us the number of assiciations for each species from each paper
Rhizobial_Associations_growth_form$diversity_sum <- as.numeric(apply(Rhizobial_Associations_growth_form[,27:312], 1, sum))

#Make a new condensed dataframe with only necessary information
Rhizobial_Associations_GF_Condensed<-Rhizobial_Associations_growth_form %>% 
  select(paper_id, genus_species, sample_country,diversity_sum,annual_perennial,growth_form,num_nodules,Plant_status) %>% 
  na.omit() %>% 
  filter(Plant_status %in% c("native","non_native"))

#Make a data frame calculating the mean, standard error, and number of observations for native vs. non-native to make graph
Average_growth_form_Graph<-Rhizobial_Associations_GF_Condensed%>%
  group_by(growth_form)%>%
  summarize(symbionts_Std=sd(diversity_sum),symbionts_Mean=mean(diversity_sum),symbionts_n=length(diversity_sum))%>%
  #Make a new column called "Richness_St_Error" and divide "Richness_Std" by the square root of "Richness_n"
  mutate(symbionts_St_Error=symbionts_Std/sqrt(symbionts_n))

#Make graph showing number of introduced vs. native
ggplot(Average_growth_form_Graph,aes(x=growth_form,y=symbionts_Mean,fill=growth_form))+
  #Make a bar graph where the height of the bars is equal to the data (stat=identity) and you preserve the vertical position while adjusting the horizontal(position_dodge), and outline the bars with the color black.
  geom_bar(stat="identity", position=position_dodge(),color="black")+
  geom_errorbar(aes(ymin=symbionts_Mean-symbionts_St_Error,ymax=symbionts_Mean+symbionts_St_Error),position=position_dodge(0.9),width=0.2)+
  #make error bars using the Standard error from the mean and place them at 0.9 with a width of 0.2
  #Label the x-axis "Watershed"
  xlab("Growth Form")+
  #Label the y-axis "Species Richness"
  ylab("Average Number of Symbionts")#+
  #scale_fill_manual(values=c("darkslategrey","cadetblue4","cadetblue3"))  


####Anova for growth form ####
summary(Mixed_Model_Rhiz_Status_Growth_Form<- lmer(diversity_sum ~ growth_form + (1|num_nodules), data = Rhizobial_Associations_GF_Condensed)) #
anova(Mixed_Model_Rhiz_Status_Growth_Form,type=2) 


#### Rhizobial Associates by Annual/Perennial  ####

#Make a data frame calculating the mean, standard error, and number of observations for native vs. non-native to make graph
Average_annual_perennial_Graph<-Rhizobial_Associations_GF_Condensed%>%
  group_by(annual_perennial,Plant_status)%>%
  summarize(symbionts_Std=sd(diversity_sum),symbionts_Mean=mean(diversity_sum),symbionts_n=length(diversity_sum))%>%
  #Make a new column called "Richness_St_Error" and divide "Richness_Std" by the square root of "Richness_n"
  mutate(symbionts_St_Error=symbionts_Std/sqrt(symbionts_n))

#Make graph showing number of introduced vs. native and perennial / annual
ggplot(Average_annual_perennial_Graph,aes(x=annual_perennial,y=symbionts_Mean,fill=Plant_status))+
  #Make a bar graph where the height of the bars is equal to the data (stat=identity) and you preserve the vertical position while adjusting the horizontal(position_dodge), and outline the bars with the color black.
  geom_bar(stat="identity", position=position_dodge(),color="black")+
  geom_errorbar(aes(ymin=symbionts_Mean-symbionts_St_Error,ymax=symbionts_Mean+symbionts_St_Error),position=position_dodge(0.9),width=0.2)+
  #make error bars using the Standard error from the mean and place them at 0.9 with a width of 0.2
  #Label the x-axis "Watershed"
  xlab("Growth Timing")+
  #Label the y-axis "Species Richness"
  ylab("Average Number of Symbionts")
#scale_fill_manual(values=c("darkslategrey","cadetblue4"))  

#Make a data frame calculating the mean, standard error, and number of observations for native vs. non-native to make graph
Average_status_herb_Graph<-Rhizobial_Associations_GF_Condensed %>%
  mutate(Growth_Status=paste(growth_form,annual_perennial,sep = "::")) %>% 
  group_by(Growth_Status,Plant_status) %>%
  summarize(symbionts_Std=sd(diversity_sum),symbionts_Mean=mean(diversity_sum),symbionts_n=length(diversity_sum)) %>%
  #Make a new column called "Richness_St_Error" and divide "Richness_Std" by the square root of "Richness_n"
  mutate(symbionts_St_Error=symbionts_Std/sqrt(symbionts_n))

#Make graph showing number of introduced vs. native and herb/woody and annual/perennial
ggplot(subset(Average_status_herb_Graph,Growth_Status!="woody::annual"),aes(x=Growth_Status,y=symbionts_Mean,fill=Plant_status))+
  #Make a bar graph where the height of the bars is equal to the data (stat=identity) and you preserve the vertical position while adjusting the horizontal(position_dodge), and outline the bars with the color black.
  geom_bar(stat="identity", position=position_dodge(),color="black")+
  geom_errorbar(aes(ymin=symbionts_Mean-symbionts_St_Error,ymax=symbionts_Mean+symbionts_St_Error),position=position_dodge(0.9),width=0.2)+
  #make error bars using the Standard error from the mean and place them at 0.9 with a width of 0.2
  #Label the x-axis "Watershed"
  xlab("Growth Timing")+
  #Label the y-axis "Species Richness"
  ylab("Average Number of Symbionts")#+
#scale_fill_manual(values=c("darkslategrey","cadetblue4"))  

####Anova for annual/perennial ####
summary(Mixed_Model_Rhiz_Status_A_P<- lmer(diversity_sum ~ annual_perennial + (1|num_nodules), data = Rhizobial_Associations_GF_Condensed)) #
anova(Mixed_Model_Rhiz_Status_A_P,type=2) 

#t.test(Rhizobial_Associations_GF_Condensed$diversity_sum~Rhizobial_Associations_GF_Condensed$annual_perennial) 

####Anova for annual/perennial*native/nonnative ####
summary(Mixed_Model_Rhiz_Status_A_P_Plant_status<- lmer(diversity_sum ~ annual_perennial*Plant_status + (1 | genus_species) + (1|num_nodules), data = Rhizobial_Associations_GF_Condensed))
#type 3 is used because we have an interaction
anova(Mixed_Model_Rhiz_Status_A_P_Plant_status,type = 3) 

####Anova for annual/perennial*native/nonnative*growth form ####
summary(Mixed_Model_Rhiz_Status_A_P_Status_GF<- lmer(diversity_sum ~ annual_perennial*Plant_status + growth_form*Plant_status + (1 | genus_species) + (1|num_nodules), data = Rhizobial_Associations_GF_Condensed)) #
anova(Mixed_Model_Rhiz_Status_A_P_Status_GF,type = 3) 

#### Number Nodules (sample number) by strain richness ####

Nod_by_strain<-Plant_Data %>%
  filter(strain_richness!="NEED SUPPLEMENTAL PAPER")%>%
  filter(!is.na(num_nodules),!is.na(strain_richness))%>% 
  mutate(strain_richness=as.numeric(strain_richness)) %>% 
  filter(num_nodules>10) %>% 
  group_by(genus_species) %>% 
  summarize(num_nod_Mean=mean(num_nodules),num_nod_n=length(num_nodules),strain_rich_Mean=mean(strain_richness),strain_rich_n=length(strain_richness))%>%
  ungroup()

hist(Nod_by_strain$num_nod_Mean) 

#graph looking at strain richness by nodule number
ggplot(Nod_by_strain,aes(y=strain_rich_Mean,x=num_nod_Mean))+
  geom_point()+
  geom_smooth()+
  #keep entire graph but zoom in on certain area
  coord_cartesian(xlim=c(0,100))

#merge plant status with number of nodules
Rhiz_symbionts_Nod_Num<-Nod_by_strain %>% 
  left_join(Average_rhizobial_symbionts_genus_sp) %>% 
  select(-symbionts_Std,-symbionts_St_Error) %>% 
  filter(!is.na(num_nod_Mean),!is.na(Plant_status))

#look at difference between strain richness number and symbionts number
ggplot(Rhiz_symbionts_Nod_Num,aes(y=strain_rich_Mean,x=symbionts_Mean))+
  geom_point()+
  geom_abline()

#make dataframe for graph
Strain_rich_graph<-Rhiz_symbionts_Nod_Num %>% 
  group_by(Plant_status) %>% 
  summarize(strain_rich_Mean2=mean(strain_rich_Mean),strain_rich_n2=length(strain_rich_Mean),strain_rich_Std=sd(strain_rich_Mean))%>%
  mutate(strain_rich_St_Error=strain_rich_Std/sqrt(strain_rich_n2)) %>% 
  ungroup()

#Make graph showing average strain richness by introduced vs. native
ggplot(Strain_rich_graph,aes(x=Plant_status,y=strain_rich_Mean2,fill=Plant_status))+
  #Make a bar graph where the height of the bars is equal to the data (stat=identity) and you preserve the vertical position while adjusting the horizontal(position_dodge), and outline the bars with the color black.
  geom_bar(stat="identity", position=position_dodge(),color="black")+
  geom_errorbar(aes(ymin=strain_rich_Mean2-strain_rich_St_Error,ymax=strain_rich_Mean2+strain_rich_St_Error),position=position_dodge(0.9),width=0.2)+
  #make error bars using the Standard error from the mean and place them at 0.9 with a width of 0.2
  #Label the x-axis "Watershed"
  xlab("Plant Species Status")+
  #Label the y-axis "Species Richness"
  ylab("Average Strain Richness")#+
  #scale_fill_manual(values=c("darkslategrey","cadetblue4"), labels=c("Native","Non-native"))  

hist(Rhizobial_Associations_Condensed$diversity_sum)

#Run anova accounting for how many species and how many nodules they sample
#### run this again without mean and accounting for nodule number #### did this above
summary (Mixed_Model_Rhiz_Status_Strain_rich<- lmer(strain_rich_Mean ~ Plant_status + (1 | species_status_n) + (1 | num_nod_Mean), data = Rhiz_symbionts_Nod_Num))
anova(Mixed_Model_Rhiz_Status_Strain_rich,type = 2)


#### Create Supplemental Table 1 with paper information ####
Fig2_Papers<-read.csv("Fig2_Papers.csv")
Fig4_Papers<-read.csv("Fig4_Papers.csv")

SupplementalTable1<-Fig2_Papers %>%
  rbind(Fig4_Papers) %>% 
  select(paper_id,author,year,title,clean_name,sample_country) %>% 
  unique()

#save supplemental table 1#
write.csv(SupplementalTable1,"SupplementalTable1.csv")
