#### Legume Meta Analysis ####
#Authors: Kim Komatsu, Kathryn Bloodworth, Smriti Pehim Limbu 
#### This script creates the basis for Supp table 1 and relies on data from other code


##### data import #####

# getting clean legume names for standardizing across papers
cleanNames <- read.csv('legume_clean_names.csv') %>% 
  separate(old_genus_species,c("genus","species","subspecies","extra"),sep=" ") %>% 
  unite(col=genus_species, c(genus,species,subspecies,extra), sep='_',na.rm=TRUE) %>% 
  rename(clean_name=new_name) %>% 
  select(-notes)

#read in paper information dataframe 
Paper_Information<-read.csv("legume_strain diversity_meta analysis_paper information.csv")

Plant_Associations<-read.csv("~/Documents/GitHub/Legume-Meta-Analysis/Legume-Meta Analysis_Data_etc/Data/legume_strain diversity_meta analysis_plant associations_edited names_presence absence.csv") %>% pivot_longer(names_to = "rhizobia_sp", values_to = "presence_absence", cols=Bradyrhizobium_nitroreducens:Paraburkholderia_caribensis)

##Now filtering out row crop and 999, 0 and 1 presence and absence and joining clean_species
Plant_Associations_Clean <- Plant_Associations %>%
  left_join(cleanNames) %>%
  filter(cultivation.status!="row crop", genus_species!="999", presence_absence=="1", plant_status!="999", species!="sp",species!="spp",species!="sp.",paper_id!=142, plant_status!="extinct") %>% 
  #filter(!(sample_country %in% c("Japan-China","Japan_China","Kenya-Sudan","Czech_Republic-France-Georgia-Hungary-Italy-Romania-Spain","Senegal-Mauritania-Tunisia-Burundi","Malawi-Zambia-Kenya","Brazil-Venezuela"))) %>% 
  select(-introduced,-invasive,-cultivated,-author,-year) %>% 
  rename(paper_plant_status=plant_status) %>% 
  mutate(global_plant_status=ifelse((exo_NA+exo_SA+exo_AU+exo_AS+exo_EU+exo_AF)>0,1,0)) %>% 
  filter(num_nodules>2) %>% 
  mutate(global_plant_status=ifelse(global_plant_status==1,"nonnative","native")) %>% 
  mutate(paper_plant_status=ifelse(paper_plant_status=="invasive","nonnative",ifelse(paper_plant_status=="introduced","nonnative", paper_plant_status))) %>% 
  select(paper_id,clean_name,sample_country,paper_plant_status,global_plant_status)


Table_File<-Plant_Associations_Clean %>% 
  left_join(Paper_Information) %>% 
  select(paper_id,author,year,journal,genus_species,sample_country,paper_plant_status,global_plant_status) %>% 
  unique()

write.csv(Table_File,"~/Documents/GitHub/Legume-Meta-Analysis/Legume-Meta Analysis_Data_etc/Data/SupportingTable1_Help.csv")



Question1<-homeAwayAll %>% #Question 1 (Figure 2)
  select(paper_id,clean_name,sample_country) %>% 
  mutate(Figure="Figure 2")
  
Question2a<-paper_plant_status %>% #Question 2 (Figure 3)
  select(paper_id,clean_name,sample_country) %>% 
  mutate(Figure="Figure 3")

Question2b<-nativeInvasiveGenetic %>% #Question 2 (Figure 4)
  ungroup() %>% 
  select(paper_id,clean_name,sample_country)%>% 
  mutate(Figure="Figure 4")

Question3<-Native_NonNative_Fig1 %>%  #Question 3 (Figure 5)
  ungroup() %>% 
  select(paper_id,new_name,sample_country) %>% 
  rename(clean_name=new_name)%>% 
  mutate(Figure="Figure 5")

Table_File <- Question1 %>% 
  rbind(Question2a) %>% 
  rbind(Question2b) %>% 
  rbind(Question3) %>% 
  unique() %>% 
  left_join(Plant_Associations_Clean) %>% 
  left_join(Paper_Information) %>%
  select(paper_id,author,year,journal,clean_name,sample_country,paper_plant_status,global_plant_status,Figure) %>%
  mutate(FigurePresence=1) %>% 
  unique() %>% 
  spread(key=Figure,value=FigurePresence,fill=NA)

  
write.csv(Table_File,"~/Documents/GitHub/Legume-Meta-Analysis/Legume-Meta Analysis_Data_etc/Data/SupportingTable1_Help2.csv")

