#### Legume Meta Analysis ####
#Authors: Kathryn Bloodworth, Smriti Pehim Limbu, Kim Komatsu


#### Install and Load Libraries #### 

#install.packages("ggplot2")
library(ggplot2)
#install.packages("tidyverse")
library(tidyverse)

#### Set Working Directory ####

#Bloodworth - Mac
setwd("/Volumes/GoogleDrive/My Drive/Projects/Legume-Meta Analysis")


#### Read in Data ####

#Read in paper information and add NA for anywhere that has a blank
Manuscript_Info<-read.csv("Data/legume_strain diversity_meta analysis_paper information.csv")
Manuscript_Info[Manuscript_Info==""]<-NA

#Read in plant association information and add NA for anywhere that has a blank
Plant_Associations<-read.csv("Data/legume_strain diversity_meta analysis_plant associations_edited names_presence absence.csv")
Plant_Associations[Plant_Associations==""]<-NA

#Read in plant data info and add NA for anywhere that has a blank
Plant_Data<-read.csv("Data/legume_strain diversity_meta analysis_plant data.csv")
Plant_Data[Plant_Data==""]<-NA

#Read in strain diversity info and add NA for anywhere that has a blank
Strain_Diversity<-read.csv("Data/legume_strain diversity_meta analysis_strain sequences.csv")
Strain_Diversity[Strain_Diversity==""]<-NA

#Read in strain sequence information
accession_numbers<-read.csv("Data/legume_strain diversity_meta analysis_strain sequences.csv",stringsAsFactors = FALSE,na.strings="")
