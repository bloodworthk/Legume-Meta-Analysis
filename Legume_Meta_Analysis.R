#### Legume Meta Analysis ####
#Authors: Kathryn Bloodworth and Kim Komatsu


#### Install and Load Libraries #### 

#package allows for extracting sequences using accession numbers from genbank
#install.packages("ape",type="source") 
library(ape)
#install.packages("seqinr")
library(seqinr)
options(repos="https://CRAN.R-project.org")
#install.packages("rentrez")
packages = c("ape", "rentrez")
package.check <- lapply(
  packages,
  FUN = function(x) {
    if (!require(x, character.only = TRUE)) {
      install.packages(x, dependencies = TRUE)
      library(x, character.only = TRUE)
    }
  }
)
#install.packages("ggplot2")
library(ggplot2)
#install.packages("raster")
library(raster)
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
setwd("/Users/kathrynbloodworth/Dropbox (Smithsonian)/Projects/Invasive Legume Meta-Analysis/data")

#### Read in Data ####
accession_numbers<-read.csv("legume_strain diversity_meta analysis_strain sequences.csv",stringsAsFactors = FALSE,na.strings="")

Manuscript_Info<-read.csv("legume_strain diversity_meta analysis_paper information.csv")
Manuscript_Info[Manuscript_Info==""]<-NA

Plant_Associations<-read.csv("legume_strain diversity_meta analysis_plant associations_edited names_presence absence.csv")
Plant_Associations[Plant_Associations==""]<-NA

Plant_Data<-read.csv("legume_strain diversity_meta analysis_plant data.csv")
Plant_Data[Plant_Data==""]<-NA

Strain_Diversity<-read.csv("legume_strain diversity_meta analysis_strain sequences.csv")
Strain_Diversity[Strain_Diversity==""]<-NA

#Set ggplot2 theme to black and white
theme_set(theme_bw())
#Update ggplot2 theme - make box around the x-axis title size 30, vertically justify x-axis title to 0.35, Place a margin of 15 around the x-axis title.  Make the x-axis title size 30. For y-axis title, make the box size 30, put the writing at a 90 degree angle, and vertically justify the title to 0.5.  Add a margin of 15 and make the y-axis text size 25. Make the plot title size 30 and vertically justify it to 2.  Do not add any grid lines.  Do not add a legend title, and make the legend size 20
theme_update(axis.title.x=element_text(size=40, vjust=-0.35, margin=margin(t=15)),
             axis.text.x=element_text(size=40), axis.title.y=element_text(size=40, angle=90, vjust=0.5,
                                                                          margin=margin(r=15)), axis.text.y=element_text(size=40), plot.title =
               element_text(size=40, vjust=2), panel.grid.major=element_blank(), 
             panel.grid.minor=element_blank(),
             legend.text=element_text(size=30))


#### Get GenBank data ####
#https://www.vikram-baliga.com/blog/2016/1/23/using-r-to-download-genbank-sequences

#Remove any NAs from accession_numbers
accession_numbers_nNA <- accession_numbers %>% 
  select(paper_id,author,title,year,genus,species,genus_species,plant_status,sample_country,sample_continent,genetic_region,strain,accession_number) %>% 
  drop_na(accession_number)

#convert column of accession numbers into list
accession_number_list<-(accession_numbers_nNA$accession_number)
str(accession_number_list)

#split list into multiple
accession_number_list1<-accession_number_list[1:100]
accession_number_list2<-accession_number_list[101:200]
accession_number_list3<-accession_number_list[201:300]
accession_number_list4<-accession_number_list[301:400]
accession_number_list5<-accession_number_list[401:500]
accession_number_list6<-accession_number_list[501:600]
accession_number_list7<-accession_number_list[601:700]
accession_number_list8<-accession_number_list[701:800]
accession_number_list9<-accession_number_list[801:900]
accession_number_list10<-accession_number_list[901:1000]
accession_number_list11<-accession_number_list[1001:1100]
accession_number_list12<-accession_number_list[1101:1200]
accession_number_list13<-accession_number_list[1201:1300]
accession_number_list14<-accession_number_list[1301:1400]
accession_number_list15<-accession_number_list[1401:1500]
accession_number_list16<-accession_number_list[1501:1600]
accession_number_list17<-accession_number_list[1601:1700]
accession_number_list18<-accession_number_list[1701:1800]
accession_number_list19<-accession_number_list[1801:1900]
accession_number_list20<-accession_number_list[1901:2000]
accession_number_list21<-accession_number_list[2001:2100]
accession_number_list22<-accession_number_list[2101:2200]
accession_number_list23<-accession_number_list[2201:2300]
accession_number_list24<-accession_number_list[2301:2400]
accession_number_list25<-accession_number_list[2401:2500]
accession_number_list26<-accession_number_list[2501:2600]
accession_number_list27<-accession_number_list[2601:2700]
accession_number_list28<-accession_number_list[2701:2800]
accession_number_list29<-accession_number_list[2801:2900]
accession_number_list30<-accession_number_list[2901:3000]
accession_number_list31<-accession_number_list[3001:3100]
accession_number_list32<-accession_number_list[3101:3200]
accession_number_list33<-accession_number_list[3201:3300]
accession_number_list34<-accession_number_list[3301:3400]
accession_number_list35<-accession_number_list[3401:3500]
accession_number_list36<-accession_number_list[3501:3600]
accession_number_list37<-accession_number_list[3601:3700]
accession_number_list38<-accession_number_list[3701:3800]
accession_number_list39<-accession_number_list[3801:3900]
accession_number_list40<-accession_number_list[3901:4000]
accession_number_list41<-accession_number_list[4001:4100]
accession_number_list42<-accession_number_list[4101:4200]
accession_number_list43<-accession_number_list[4201:4300]
accession_number_list44<-accession_number_list[4301:4400]
accession_number_list45<-accession_number_list[4401:4500]
accession_number_list46<-accession_number_list[4501:4600]
accession_number_list47<-accession_number_list[4601:4700]
accession_number_list48<-accession_number_list[4701:4800]
accession_number_list49<-accession_number_list[4801:4900]
accession_number_list50<-accession_number_list[4901:5000]
accession_number_list51<-accession_number_list[5001:5100]
accession_number_list52<-accession_number_list[5101:5200]
accession_number_list53<-accession_number_list[5201:5300]
accession_number_list54<-accession_number_list[5301:5400]
accession_number_list55<-accession_number_list[5401:5500]
accession_number_list56<-accession_number_list[5501:5600]
accession_number_list57<-accession_number_list[5601:5700]
accession_number_list58<-accession_number_list[5701:5800]
accession_number_list59<-accession_number_list[5801:5900]
accession_number_list60<-accession_number_list[5901:6000]
accession_number_list61<-accession_number_list[6001:6100]
accession_number_list62<-accession_number_list[6101:6200]
accession_number_list63<-accession_number_list[6201:6300]
accession_number_list64<-accession_number_list[6301:6400]
accession_number_list65<-accession_number_list[6401:6500]
accession_number_list66<-accession_number_list[6501:6600]
accession_number_list67<-accession_number_list[6601:6700]
accession_number_list68<-accession_number_list[6701:6800]
accession_number_list69<-accession_number_list[6801:6900]
accession_number_list70<-accession_number_list[6901:7000]
accession_number_list71<-accession_number_list[7001:7100]
accession_number_list72<-accession_number_list[7101:7200]
accession_number_list73<-accession_number_list[7201:7300]
accession_number_list74<-accession_number_list[7301:7400]
accession_number_list75<-accession_number_list[7401:7500]
accession_number_list76<-accession_number_list[7501:7558]

#connect to GenBank database and download the sequences

read.GenBank(d,species.name=TRUE)



sequences1 <- entrez_fetch(id = accession_number_list1,
                         db = "nuccore", 
                         rettype = "fasta")
sequences2 <- entrez_fetch(id = accession_number_list1,
                           db = "nuccore", 
                           rettype = "fasta")
sequences3 <- entrez_fetch(id = accession_number_list1,
                           db = "nuccore", 
                           rettype = "fasta")
sequences4 <- entrez_fetch(id = accession_number_list1,
                           db = "nuccore", 
                           rettype = "fasta")
sequences5 <- entrez_fetch(id = accession_number_list1,
                           db = "nuccore", 
                           rettype = "fasta")
sequences6 <- entrez_fetch(id = accession_number_list1,
                           db = "nuccore", 
                           rettype = "fasta")
sequences7 <- entrez_fetch(id = accession_number_list1,
                           db = "nuccore", 
                           rettype = "fasta")
sequences8 <- entrez_fetch(id = accession_number_list1,
                           db = "nuccore", 
                           rettype = "fasta")
sequences9 <- entrez_fetch(id = accession_number_list1,
                           db = "nuccore", 
                           rettype = "fasta")
sequences10 <- entrez_fetch(id = accession_number_list1,
                           db = "nuccore", 
                           rettype = "fasta")

sequences11 <- entrez_fetch(id = accession_number_list1,
                           db = "nuccore", 
                           rettype = "fasta")
sequences12 <- entrez_fetch(id = accession_number_list1,
                           db = "nuccore", 
                           rettype = "fasta")
sequences13 <- entrez_fetch(id = accession_number_list1,
                           db = "nuccore", 
                           rettype = "fasta")
sequences14 <- entrez_fetch(id = accession_number_list1,
                           db = "nuccore", 
                           rettype = "fasta")
sequences15 <- entrez_fetch(id = accession_number_list1,
                           db = "nuccore", 
                           rettype = "fasta")
sequences16 <- entrez_fetch(id = accession_number_list1,
                           db = "nuccore", 
                           rettype = "fasta")
sequences17 <- entrez_fetch(id = accession_number_list1,
                           db = "nuccore", 
                           rettype = "fasta")
sequences18 <- entrez_fetch(id = accession_number_list1,
                           db = "nuccore", 
                           rettype = "fasta")
sequences19 <- entrez_fetch(id = accession_number_list1,
                           db = "nuccore", 
                           rettype = "fasta")
sequences20 <- entrez_fetch(id = accession_number_list1,
                            db = "nuccore", 
                            rettype = "fasta")

sequences21 <- entrez_fetch(id = accession_number_list1,
                           db = "nuccore", 
                           rettype = "fasta")
sequences22 <- entrez_fetch(id = accession_number_list1,
                           db = "nuccore", 
                           rettype = "fasta")
sequences23 <- entrez_fetch(id = accession_number_list1,
                           db = "nuccore", 
                           rettype = "fasta")
sequences24 <- entrez_fetch(id = accession_number_list1,
                           db = "nuccore", 
                           rettype = "fasta")
sequences25 <- entrez_fetch(id = accession_number_list1,
                           db = "nuccore", 
                           rettype = "fasta")
sequences26 <- entrez_fetch(id = accession_number_list1,
                           db = "nuccore", 
                           rettype = "fasta")
sequences27 <- entrez_fetch(id = accession_number_list1,
                           db = "nuccore", 
                           rettype = "fasta")
sequences28 <- entrez_fetch(id = accession_number_list1,
                           db = "nuccore", 
                           rettype = "fasta")
sequences29 <- entrez_fetch(id = accession_number_list1,
                           db = "nuccore", 
                           rettype = "fasta")
sequences30 <- entrez_fetch(id = accession_number_list1,
                            db = "nuccore", 
                            rettype = "fasta")

sequences31 <- entrez_fetch(id = accession_number_list1,
                           db = "nuccore", 
                           rettype = "fasta")
sequences32 <- entrez_fetch(id = accession_number_list1,
                           db = "nuccore", 
                           rettype = "fasta")
sequences33 <- entrez_fetch(id = accession_number_list1,
                           db = "nuccore", 
                           rettype = "fasta")
sequences34 <- entrez_fetch(id = accession_number_list1,
                           db = "nuccore", 
                           rettype = "fasta")
sequences35 <- entrez_fetch(id = accession_number_list1,
                           db = "nuccore", 
                           rettype = "fasta")
sequences36 <- entrez_fetch(id = accession_number_list1,
                           db = "nuccore", 
                           rettype = "fasta")
sequences37 <- entrez_fetch(id = accession_number_list1,
                           db = "nuccore", 
                           rettype = "fasta")
sequences38 <- entrez_fetch(id = accession_number_list1,
                           db = "nuccore", 
                           rettype = "fasta")
sequences39 <- entrez_fetch(id = accession_number_list1,
                           db = "nuccore", 
                           rettype = "fasta")
sequences40 <- entrez_fetch(id = accession_number_list1,
                            db = "nuccore", 
                            rettype = "fasta")

sequences41 <- entrez_fetch(id = accession_number_list1,
                           db = "nuccore", 
                           rettype = "fasta")
sequences42 <- entrez_fetch(id = accession_number_list1,
                           db = "nuccore", 
                           rettype = "fasta")
sequences43 <- entrez_fetch(id = accession_number_list1,
                           db = "nuccore", 
                           rettype = "fasta")
sequences44 <- entrez_fetch(id = accession_number_list1,
                           db = "nuccore", 
                           rettype = "fasta")
sequences45 <- entrez_fetch(id = accession_number_list1,
                           db = "nuccore", 
                           rettype = "fasta")
sequences46 <- entrez_fetch(id = accession_number_list1,
                           db = "nuccore", 
                           rettype = "fasta")
sequences47 <- entrez_fetch(id = accession_number_list1,
                           db = "nuccore", 
                           rettype = "fasta")
sequences48 <- entrez_fetch(id = accession_number_list1,
                           db = "nuccore", 
                           rettype = "fasta")
sequences49 <- entrez_fetch(id = accession_number_list1,
                           db = "nuccore", 
                           rettype = "fasta")
sequences50 <- entrez_fetch(id = accession_number_list1,
                            db = "nuccore", 
                            rettype = "fasta")

sequences51 <- entrez_fetch(id = accession_number_list1,
                           db = "nuccore", 
                           rettype = "fasta")
sequences52 <- entrez_fetch(id = accession_number_list1,
                           db = "nuccore", 
                           rettype = "fasta")
sequences53 <- entrez_fetch(id = accession_number_list1,
                           db = "nuccore", 
                           rettype = "fasta")
sequences54 <- entrez_fetch(id = accession_number_list1,
                           db = "nuccore", 
                           rettype = "fasta")
sequences55 <- entrez_fetch(id = accession_number_list1,
                           db = "nuccore", 
                           rettype = "fasta")
sequences56 <- entrez_fetch(id = accession_number_list1,
                           db = "nuccore", 
                           rettype = "fasta")
sequences57 <- entrez_fetch(id = accession_number_list1,
                           db = "nuccore", 
                           rettype = "fasta")
sequences58 <- entrez_fetch(id = accession_number_list1,
                           db = "nuccore", 
                           rettype = "fasta")
sequences59 <- entrez_fetch(id = accession_number_list1,
                           db = "nuccore", 
                           rettype = "fasta")
sequences60 <- entrez_fetch(id = accession_number_list1,
                            db = "nuccore", 
                            rettype = "fasta")

sequences61 <- entrez_fetch(id = accession_number_list1,
                           db = "nuccore", 
                           rettype = "fasta")
sequences62 <- entrez_fetch(id = accession_number_list1,
                           db = "nuccore", 
                           rettype = "fasta")
sequences63 <- entrez_fetch(id = accession_number_list1,
                           db = "nuccore", 
                           rettype = "fasta")
sequences64 <- entrez_fetch(id = accession_number_list1,
                           db = "nuccore", 
                           rettype = "fasta")
sequences65 <- entrez_fetch(id = accession_number_list1,
                           db = "nuccore", 
                           rettype = "fasta")
sequences66 <- entrez_fetch(id = accession_number_list1,
                           db = "nuccore", 
                           rettype = "fasta")
sequences67 <- entrez_fetch(id = accession_number_list1,
                           db = "nuccore", 
                           rettype = "fasta")
sequences68 <- entrez_fetch(id = accession_number_list1,
                           db = "nuccore", 
                           rettype = "fasta")
sequences69 <- entrez_fetch(id = accession_number_list1,
                           db = "nuccore", 
                           rettype = "fasta")
sequences70 <- entrez_fetch(id = accession_number_list1,
                            db = "nuccore", 
                            rettype = "fasta")


sequences71 <- entrez_fetch(id = accession_number_list1,
                           db = "nuccore", 
                           rettype = "fasta")
sequences72 <- entrez_fetch(id = accession_number_list1,
                           db = "nuccore", 
                           rettype = "fasta")
sequences73 <- entrez_fetch(id = accession_number_list1,
                           db = "nuccore", 
                           rettype = "fasta")
sequences74 <- entrez_fetch(id = accession_number_list1,
                           db = "nuccore", 
                           rettype = "fasta")
sequences75 <- entrez_fetch(id = accession_number_list1,
                           db = "nuccore", 
                           rettype = "fasta")
sequences76 <- entrez_fetch(id = accession_number_list1,
                           db = "nuccore", 
                           rettype = "fasta")
write(sequences76,file="sequences76.fasta")


#### Counting number of papers,species, etc ####

## Determine number of papers assessed
Paper_Number<-Manuscript_Info%>%
  group_by(paper_id, author, year)%>%
  unique()%>%
  summarise(num_papers=length(paper_id))%>%
  ungroup() %>% 
  summarise(papers=length(paper_id))%>%
  summarise(num_papers=sum(papers))

## Determine number plant species
Species_Number<-Plant_Associations%>%
group_by(genus_species)%>%
  unique()%>%
  summarise(num_species=length(genus_species))%>%
  ungroup()%>%
  filter(genus_species!=999) %>% 
  summarise(num_species=length(genus_species))

## Determine number of countries assessed
Country_Number<-Plant_Associations%>%
  group_by(sample_country)%>%
  unique()%>%
  summarise(num_countries=length(sample_country))%>%
  ungroup()%>%
  filter(sample_country!=999) %>% 
  summarise(num_countries=length(sample_country))

#### Introduced vs. Native ####

# Create data frame with native vs. introduced species
Introduced_native_Number<-Plant_Associations%>%
  mutate(Plant_status=ifelse(plant_status=="introduced","non-native",ifelse(plant_status=="invasive","non-native",ifelse(plant_status=="native","native", ifelse(plant_status=="extinct","extinct",ifelse(plant_status=="introduced-native","introduced-native", ifelse(plant_status==999,999,ifelse(plant_status==" native","native",plant_status))))))))%>%
  group_by(Plant_status)%>%
  unique()%>%
  summarise(plant_status_counts=length(Plant_status))%>%
  ungroup()

# Make a graph for non-native vs native
Introduced_Native_Number_Graph<-Introduced_native_Number%>%
  filter(Plant_status!=999)%>%
  filter(Plant_status!=" native")%>%
  filter(Plant_status!="introduced-native")%>%
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
  mutate(Plant_status=ifelse(plant_status=="introduced","non_native",ifelse(plant_status=="invasive","non_native",ifelse(plant_status=="native","native",ifelse(plant_status=="extinct","extinct",ifelse(plant_status=="introduced-native","introduced-native",ifelse(plant_status==999,999,ifelse(plant_status==" native","native",plant_status)))))))) 

#Change all 999 to NAs for next addition step
Rhizobial_Associations[ , 27:290 ][Rhizobial_Associations[,27:290]==999]<-NA

#Sum across all columns with presence absence data for each row, giving us the number of assiciations for each species from each paper
Rhizobial_Associations$diversity_sum <- as.numeric(apply(Rhizobial_Associations[,27:290], 1, sum))

#Make a new condensed dataframe with only necessary information
Rhizobial_Associations_Condensed<-Rhizobial_Associations %>% 
  select(paper_id, genus_species, sample_country,diversity_sum,Plant_status) %>% 
  filter(Plant_status!=999)%>%
  filter(Plant_status!="introduced-native")%>%
  filter(Plant_status!="extinct") %>% 
  na.omit()

#Make a data frame calculating the mean, standard error, and number of observations for native vs. non-native to make graph
Average_rhizobial_symbionts_Graph<-Rhizobial_Associations_Condensed%>%
  filter(diversity_sum!=999)%>%
  group_by(Plant_status)%>%
  summarize(symbionts_Std=sd(diversity_sum),symbionts_Mean=mean(diversity_sum),symbionts_n=length(diversity_sum))%>%
  #Make a new column called "Richness_St_Error" and divide "Richness_Std" by the square root of "Richness_n"
  mutate(symbionts_St_Error=symbionts_Std/sqrt(symbionts_n))

#Make graph showing number of introduced vs. native
ggplot(Average_rhizobial_symbionts_Graph,aes(x=Plant_status,y=symbionts_Mean,fill=Plant_status))+
  #Make a bar graph where the height of the bars is equal to the data (stat=identity) and you preserve the vertical position while adjusting the horizontal(position_dodge), and outline the bars with the color black.
  geom_bar(stat="identity", position=position_dodge(),color="black")+
  geom_errorbar(aes(ymin=symbionts_Mean-symbionts_St_Error,ymax=symbionts_Mean+symbionts_St_Error),position=position_dodge(0.9),width=0.2)+
  #make error bars using the Standard error from the mean and place them at 0.9 with a width of 0.2
  #Label the x-axis "Watershed"
  xlab("Plant Species Status")+
  #Label the y-axis "Species Richness"
  ylab("Average Number of Symbionts")+
  #Fill the bar graphs with grey and white according and label them "Inside Fence" and "Outside Fence"
  scale_fill_manual(values=c("darkslategrey","cadetblue4"), labels=c("Native","Non-native"))  

t.test(Rhizobial_Associations_Condensed$diversity_sum~Rhizobial_Associations_Condensed$Plant_status)  

#### Rhizobial Associates by Growth Form  ####

Plant_Data_simplified<-Plant_Data %>% 
  select(paper_id,author,year,genus_species,annual_perennial,growth_form)
  
#Make a new data frame called rhizobial_associates_growth_form
Rhizobial_Associations_growth_form<-Plant_Associations%>%
  mutate(Plant_status=ifelse(plant_status=="introduced","non_native",ifelse(plant_status=="invasive","non_native",ifelse(plant_status=="native","native",ifelse(plant_status=="extinct","extinct",ifelse(plant_status=="introduced-native","introduced-native",ifelse(plant_status==999,999,ifelse(plant_status==" native","native",plant_status)))))))) %>% 
  left_join(Plant_Data_simplified)

#Change all 999 to NAs for next addition step
Rhizobial_Associations_growth_form[ , 27:290 ][Rhizobial_Associations_growth_form[,27:290]==999]<-NA

#Sum across all columns with presence absence data for each row, giving us the number of assiciations for each species from each paper
Rhizobial_Associations_growth_form$diversity_sum <- as.numeric(apply(Rhizobial_Associations_growth_form[,27:290], 1, sum))
  

#Make a new condensed dataframe with only necessary information
Rhizobial_Associations_GF_Condensed<-Rhizobial_Associations_growth_form %>% 
  select(paper_id, genus_species, sample_country,diversity_sum,annual_perennial,growth_form) %>% 
  na.omit()


#Make a data frame calculating the mean, standard error, and number of observations for native vs. non-native to make graph
Average_growth_form_Graph<-Rhizobial_Associations_GF_Condensed%>%
  filter(diversity_sum!=999)%>%
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
  ylab("Average Number of Symbionts")+
  scale_fill_manual(values=c("darkslategrey","cadetblue4"))  

t.test(Rhizobial_Associations_GF_Condensed$diversity_sum~Rhizobial_Associations_GF_Condensed$growth_form)

#### Rhizobial Associates by Annual/Perennial  ####

#Make a data frame calculating the mean, standard error, and number of observations for native vs. non-native to make graph
Average_annual_perennial_Graph<-Rhizobial_Associations_GF_Condensed%>%
  filter(diversity_sum!=999)%>%
  group_by(annual_perennial)%>%
  summarize(symbionts_Std=sd(diversity_sum),symbionts_Mean=mean(diversity_sum),symbionts_n=length(diversity_sum))%>%
  #Make a new column called "Richness_St_Error" and divide "Richness_Std" by the square root of "Richness_n"
  mutate(symbionts_St_Error=symbionts_Std/sqrt(symbionts_n))

#Make graph showing number of introduced vs. native
ggplot(Average_annual_perennial_Graph,aes(x=annual_perennial,y=symbionts_Mean,fill=annual_perennial))+
  #Make a bar graph where the height of the bars is equal to the data (stat=identity) and you preserve the vertical position while adjusting the horizontal(position_dodge), and outline the bars with the color black.
  geom_bar(stat="identity", position=position_dodge(),color="black")+
  geom_errorbar(aes(ymin=symbionts_Mean-symbionts_St_Error,ymax=symbionts_Mean+symbionts_St_Error),position=position_dodge(0.9),width=0.2)+
  #make error bars using the Standard error from the mean and place them at 0.9 with a width of 0.2
  #Label the x-axis "Watershed"
  xlab("Growth Timing")+
  #Label the y-axis "Species Richness"
  ylab("Average Number of Symbionts")+
  scale_fill_manual(values=c("darkslategrey","cadetblue4"))  

t.test(Rhizobial_Associations_GF_Condensed$diversity_sum~Rhizobial_Associations_GF_Condensed$annual_perennial) 
