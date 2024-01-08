#Original title: "3a-Qtitle: "3a-Question1-R-legutitle: "3a-Question1-R-legu
#date: "2023-03-17"
#Authors: Kathryn Bloodworth, Smriti Pehim Limbu, Kimberly Komatsu

library(ggplot2)
library(lme4)
library(lmerTest)
library(tidyverse)
library(dplyr)
library(readr)
library(olsrr)

#set working directory
setwd("~/Documents/GitHub/Legume-Meta-Analysis/Legume-Meta Analysis_Data_etc/Data")


#Set ggplot2 theme to black and white
theme_set(theme_bw())
#Update ggplot2 theme - make box around the x-axis title size 30, vertically justify x-axis title to 0.35, Place a margin of 15 around the x-axis title.  Make the x-axis title size 30. For y-axis title, make the box size 30, put the writing at a 90 degree angle, and vertically justify the title to 0.5.  Add a margin of 15 and make the y-axis text size 25. Make the plot title size 30 and vertically justify it to 2.  Do not add any grid lines.  Do not add a legend title, and make the legend size 20
theme_update(axis.title.x=element_text(size=30, vjust=-0.35, margin=margin(t=15)),
             axis.text.x=element_text(size=30), axis.title.y=element_text(size=30, angle=90, vjust=0.5,
                                                                          margin=margin(r=15)), axis.text.y=element_text(size=30), plot.title =
               element_text(size=30, vjust=2), panel.grid.major=element_blank(),
             panel.grid.minor=element_blank(), legend.title=element_blank(),
             legend.text=element_text(size=40))

####Read in Dataframes and Clean Data ####
### Starting a new analysis and this chunk is for loading files, Plant_Associations is the cleaned file
Clean_Species<-read_csv("~/Documents/GitHub/Legume-Meta-Analysis/Legume-Meta Analysis_Data_etc/Data/legume_clean_names.csv") %>% 
  filter(old_genus_species!="") %>% 
  separate(old_genus_species,c("genus","species","subspecies","extra"),sep=" ") %>% 
  unite(col=genus_species, c(genus,species,subspecies,extra), sep='_',na.rm=TRUE)

Plant_Associations<-read.csv("~/Documents/GitHub/Legume-Meta-Analysis/Legume-Meta Analysis_Data_etc/Data/legume_strain diversity_meta analysis_plant associations_edited names_presence absence.csv") %>% pivot_longer(names_to = "rhizobia_sp", values_to = "presence_absence", cols=Bradyrhizobium_nitroreducens:Paraburkholderia_caribensis)

##Now filtering out row crop and 999, 0 and 1 presence and absence and joining clean_species
Plant_Associations_Clean <- Plant_Associations %>%
  left_join(Clean_Species) %>%
  filter(cultivation.status!="row crop", genus_species!="999", presence_absence=="1", plant_status!="999", species!="sp",species!="spp",species!="sp.",paper_id!=142, plant_status!="extinct") %>% 
  #filter(!(sample_country %in% c("Japan-China","Japan_China","Kenya-Sudan","Czech_Republic-France-Georgia-Hungary-Italy-Romania-Spain","Senegal-Mauritania-Tunisia-Burundi","Malawi-Zambia-Kenya","Brazil-Venezuela"))) %>% 
  select(-introduced,-invasive,-cultivated) %>% 
  rename(paper_plant_status=plant_status) %>% 
  mutate(global_plant_status=ifelse((exo_NA+exo_SA+exo_AU+exo_AS+exo_EU+exo_AF)>0,1,0)) %>% 
  filter(num_nodules>4) %>% 
  mutate(global_plant_status=ifelse(global_plant_status==1,"nonnative","native")) %>% 
  mutate(paper_plant_status=ifelse(paper_plant_status=="invasive","nonnative",ifelse(paper_plant_status=="introduced","nonnative", paper_plant_status)))

colnames(Plant_Associations_Clean)


Plant_Data<-read.csv("~/Documents/GitHub/Legume-Meta-Analysis/Legume-Meta Analysis_Data_etc/Data/legume_strain diversity_meta analysis_plant data.csv")

##Now filtering out row crop and 999, 0 and 1 presence and absence and joining clean_species
Plant_Data_Clean <- Plant_Data %>%
  rename(notes_plantdata=notes) %>% 
  left_join(Clean_Species) %>%
  filter(cultivation.status!="row crop", genus_species!="999", species!="sp",species!="spp",species!="sp.",paper_id!=142) %>% 
  filter(num_nodules>4)

colnames(Plant_Data_Clean)


#### Figure 4 Dataframe Creation ####
#r Create dataframe with just papers where native/non natives were studied within a given location

#Create a new dataframe that has only plant associations with species 
Plant_Status<-Plant_Associations_Clean %>% 
  select(paper_id,genus_species,new_name,cultivation.status,sample_country, global_plant_status,paper_plant_status) %>% 
  unique()

Native_NonNative_Fig1 <- Plant_Data_Clean %>%
  #merge plant associations so that we can determine plant status, etc.
  merge(Plant_Status,by = c("paper_id", "genus_species", "new_name", "cultivation.status","sample_country"),all=TRUE) %>% 
  select(paper_id,new_name,paper_plant_status,annual_perennial,growth_form,sample_country,sample_continent,num_nodules,genetic_region,strain_richness,compares_natinv) %>% 
  #only keep papers that looked at native and non native species within the same location
  filter(compares_natinv==1) %>% 
  na.omit(plant_paper_status) %>% 
  mutate(Graph_x=ifelse(paper_plant_status=="native","Native (n=47)","Non Native (n=19)")) %>% 
  #remove papers that no longer have a native or non-native partner after filtering steps
  filter(!(paper_id %in% c(289,325))) %>% 
  #average across plant ID, species, and status
  group_by(paper_id,new_name,sample_country,paper_plant_status,annual_perennial,growth_form,sample_continent,Graph_x) %>% 
  summarise(avg_strain_richness=mean(as.numeric(strain_richness))) %>% 
  ungroup()

#Calculate N number
Native_NonNative_N<-Native_NonNative_Fig1 %>% 
  group_by(Graph_x) %>% 
  summarise(status=length(new_name)) #47 native, 19 non native

#### Figure 4 Creation ####
ggplot(Native_NonNative_Fig1,aes(x=Graph_x,y=avg_strain_richness,fill=Graph_x))+
  geom_boxplot()+
  scale_x_discrete(labels = function(x) str_wrap(x, width = 10))+
  xlab("Local Plant Status")+
  ylab("Rhizobial Strain Richness")+
  scale_fill_manual(values=c("#E2E4DE","#A79371"))+
  theme(legend.position="none")+
  expand_limits(y=c(0,25))


#### Figure 4 Stats ####
#LMER native/non natives within a given location

#Data normality
Native_NonNative_Norm <- lm(data = Native_NonNative_Fig1, avg_strain_richness ~ paper_plant_status)
ols_plot_resid_hist(Native_NonNative_Norm) 
ols_test_normality(Native_NonNative_Norm) #left skewed

Native_NonNative_Fig1<-Native_NonNative_Fig1 %>% 
  mutate(strain_richness_TF=log(avg_strain_richness))

Native_NonNative_Norm_TF <- lm(data = Native_NonNative_Fig1, strain_richness_TF ~ paper_plant_status)
ols_plot_resid_hist(Native_NonNative_Norm_TF) 
ols_test_normality(Native_NonNative_Norm_TF) #normal

#lmer
Native_NonNative_lmer_2 <- lmerTest::lmer(data = Native_NonNative_Fig1, strain_richness_TF ~ paper_plant_status + (1|paper_id))
anova(Native_NonNative_lmer_2) #NS

#### Pie Chart Data Frame####

#Create a new dataframe that has only plant associations with species 
PieChart<-Native_NonNative_Fig1 %>% 
  left_join(Plant_Associations_Clean) %>% 
  select(paper_id,new_name,sample_country,sample_continent,paper_plant_status,Graph_x,avg_strain_richness,rhizobia_sp,presence_absence) %>% 
  #remove all information except for local plant status and rhizobia species 
  select(Graph_x,rhizobia_sp,presence_absence) %>% 
  #remove duplicates so we're left with a list of just the rhizobia that are found with any native legume anywhere, and the rhizobia that are found wiht any non-native legume anywhere
  unique() %>% 
  #fix error where "R..cellulosilyticum" should be Rhizobium_cellulosilyticum",etc.
  mutate(rhizobia_sp=ifelse(rhizobia_sp=="R..cellulosilyticum","Rhizobium_cellulosilyticum",ifelse(rhizobia_sp=="E_kummerowiae","Ensifer_kummerowiae",rhizobia_sp))) %>% 
  #create wide format dataframe
  spread(key=rhizobia_sp,value=presence_absence,fill=0) 

write.csv(PieChart,"Figure3_PieChart.csv")

#in excel: made a new row where the rhizobial strain recieves a 1 if it's novel in the non-native range, a 2 if its is found in both the native and non native range, and a 0 if it is only found in the native range 
PieChart_Updated <- read_csv("Figure3_PieChart_updated.csv")

PieChart_Updated[,2:58] <- sapply(PieChart_Updated[,2:58],as.numeric)

PieChart_NovelNotNovel <- PieChart_Updated %>% 
  #select only the row called Novel_NonNovel
  filter(Graph_x=="Novel_NonNovel") %>% 
  #Make dataframe long instead of wide
  gather(key="Rhizobial_Strain",value="Novel_Status",2:58) %>% 
#filter out any 0s because this means they were not found in non-native ranges
  dplyr::filter(Novel_Status>0) %>% 
  mutate(Graph_x=ifelse(Novel_Status==2,"NotNovel","Novel"))


#count the number of strains that are novel
PieChart_Novel<- PieChart_NovelNotNovel %>% 
  filter(Graph_x=="Novel") %>% 
  mutate(CountNovelNonNovel=length(Novel_Status)) %>% 
  select(-Rhizobial_Strain) %>% 
  unique()

#count the number of strains that are not novel
PieChart_NotNovel<- PieChart_NovelNotNovel %>% 
  filter(Graph_x=="NotNovel") %>% 
  mutate(CountNovelNonNovel=length(Novel_Status))%>% 
  select(-Rhizobial_Strain) %>% 
  unique()


PieChart_Graph<-PieChart_Novel %>% 
rbind(PieChart_NotNovel) %>% 
select(-Novel_Status) %>% 
mutate(Sum=29) %>% 
mutate(Percent=CountNovelNonNovel/Sum*100)###

#Write CSV so I can create pie chart in excel
write.csv(PieChart_Graph,"PieChart_Graph.csv")


#not working
ggplot(PieChart_Graph, aes(x=Graph_x, y=Percent, fill=Graph_x)) +
  geom_bar(stat="identity", width=1) +
  coord_polar("y", start=0)


























































































































