#### Legume Meta Analysis ####
#Authors: Kathryn Bloodworth, Smriti Pehim Limbu, Kim Komatsu


#### Install and Load Libraries #### 

#install.packages("ggplot2")
library(ggplot2)
#install.packages("tidyverse")
library(tidyverse)

#### Set Working Directory ####

#Bloodworth - Mac
setwd("~/Library/CloudStorage/Box-Box/Projects/Invasive Legume Meta-Analysis/data")

#Bloodworth - PC
setwd("/Users/kjbloodw/Box-Box/Projects/Invasive Legume Meta-Analysis/data")

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