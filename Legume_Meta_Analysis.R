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
#install.packages("car")
library(car)


#### Set Working Directory ####
#Bloodworth - Mac
setwd("~/Library/CloudStorage/Box-Box/Projects/Invasive Legume Meta-Analysis/data")

#pc
setwd("/Users/kjbloodw/Dropbox (Smithsonian)/Projects//Invasive Legume Meta-Analysis/data/")

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
accession_numbers_nNA_multiple <- accession_numbers %>% 
  select(paper_id,author,title,year,genus,species,genus_species,plant_status,sample_country,sample_continent,genetic_region,strain,accession_number) %>% 
  drop_na(accession_number)

#Remove any accession numbers for papers where only native species were studied. Keep all studies where introduced species were studied and keep all studies where native AND introduced species were studied either within a paper or across papers
accession_numbers_nNA<-accession_numbers_nNA_multiple %>% 
  filter(!paper_id %in% c(3,7,10,13,16,17,18,23,25,28,29,30,32,33,34,37,40,41,42,43,47,54,60,61,64,66,67,68,70,71,72,73,74,75,81,82,84,88,90,92,94,95,96,98,99,100,101,102,103,104,105,106,107,108,109,110,112,116,118,120,121,123,124,127,130,131,132,135,136,140,141,142,143,145,146,149,151,152,153,155,157,161,163,165,166,167,168,169,170,171,172,174,175,178,181,183,189,190,192,193,194,195,196,198,199,200,201,202,204,205,206,207,208,209,212,213,214,215,220,221,222,223,224,225,227,228,229,231,232,233,234,239,241,242,243,244,245,246,247,248,250,251,252,255,256,258,259,260,261,262,263,264,268,269,271,272,273,275,276,278,279,281,282,283,284,285,286,288,290,292,293,294,295,296,300,301))

#remove any repeated accession numbers 
accession_numbers_Filtered <- distinct(accession_numbers_nNA,accession_number, .keep_all = TRUE) #removed 180 accession numbers

#make subsetted data frames for each genetic region and convert into just a list of accession numbers

accession_numbers_16S<-accession_numbers_Filtered %>% 
  filter(genetic_region %in% c("16S","16S-RFLP"))
#make it a list
accession_numbers_16S_list<-(accession_numbers_16S$accession_number)
str(accession_numbers_16S_list)

accession_numbers_nodC<-accession_numbers_Filtered %>% 
  filter(genetic_region %in% c("nodC","NodC")) 
#make it a list
accession_numbers_nodC_list<-(accession_numbers_nodC$accession_number)
  
accession_numbers_nifH<-accession_numbers_Filtered %>% 
  filter(genetic_region=="nifH")
#make it a list
accession_numbers_nifH_list<-(accession_numbers_nifH$accession_number)

accession_numbers_nifHD<-accession_numbers_Filtered %>% 
  filter(genetic_region=="nifHD") 
#make it a list
accession_numbers_nifHD_list<-(accession_numbers_nifHD$accession_number)

accession_numbers_nodA<-accession_numbers_Filtered %>% 
  filter(genetic_region=="nodA") 
#make it a list
accession_numbers_nodA_list<-(accession_numbers_nodA$accession_number)


accession_numbers_nodZ<-accession_numbers_Filtered %>% 
  filter(genetic_region=="nodZ") 
#make it a list
accession_numbers_nodZ_list<-(accession_numbers_nodZ$accession_number)

accession_numbers_glnII<-accession_numbers_Filtered %>% 
  filter(genetic_region=="glnII")
#make it a list
accession_numbers_glnII_list<-(accession_numbers_glnII$accession_number)

accession_numbers_atpD<-accession_numbers_Filtered %>% 
  filter(genetic_region %in% c("atpd","atpD")) 
#make it a list
accession_numbers_atpD_list<-(accession_numbers_atpD$accession_number)

accession_numbers_recA<-accession_numbers_Filtered %>% 
  filter(genetic_region %in% c("recA","RecA","rec A"))
#make it a list
accession_numbers_recA_list<-(accession_numbers_recA$accession_number)


accession_numbers_dnaK<-accession_numbers_Filtered %>% 
  filter(genetic_region %in% c("dnaK","dnaKÊ")) 
#make it a list
accession_numbers_dnaK_list<-(accession_numbers_dnaK$accession_number)

accession_numbers_noeI<-accession_numbers_Filtered %>% 
  filter(genetic_region=="noeI")
#make it a list
accession_numbers_noeI_list<-(accession_numbers_noeI$accession_number)

accession_numbers_nodD<-accession_numbers_Filtered %>% 
  filter(genetic_region=="nodD") 
#make it a list
accession_numbers_nodD_list<-(accession_numbers_nodD$accession_number)

accession_numbers_ITS<-accession_numbers_Filtered %>% 
  filter(genetic_region %in% c("ITS","IGS")) 
#make it a list
accession_numbers_ITS_list<-(accession_numbers_ITS$accession_number)

accession_numbers_23S<-accession_numbers_Filtered %>% 
  filter(genetic_region=="23S") 
#make it a list
accession_numbers_23S_list<-(accession_numbers_23S$accession_number)

accession_numbers_nifD<-accession_numbers_Filtered %>% 
  filter(genetic_region=="nifD") 
#make it a list
accession_numbers_nifD_list<-(accession_numbers_nifD$accession_number)

accession_numbers_thrC<-accession_numbers_Filtered %>% 
  filter(genetic_region=="thrC") 
#make it a list
accession_numbers_thrC_list<-(accession_numbers_thrC$accession_number)

accession_numbers_glnA<-accession_numbers_Filtered %>% 
  filter(genetic_region=="glnA")
#make it a list
accession_numbers_glnA_list<-(accession_numbers_glnA$accession_number)

accession_numbers_glnB<-accession_numbers_Filtered %>% 
  filter(genetic_region=="glnB") 
#make it a list
accession_numbers_glnB_list<-(accession_numbers_glnB$accession_number)

accession_numbers_gyrB<-accession_numbers_Filtered %>% 
  filter(genetic_region=="gyrB")
#make it a list
accession_numbers_gyrB_list<-(accession_numbers_gyrB$accession_number)

accession_numbers_rrs<-accession_numbers_Filtered %>% 
  filter(genetic_region=="rrs") 
#make it a list
accession_numbers_rrs_list<-(accession_numbers_rrs$accession_number)

accession_numbers_cpn60<-accession_numbers_Filtered %>% 
  filter(genetic_region=="cpn60") 
#make it a list
accession_numbers_cpn60_list<-(accession_numbers_cpn60$accession_number)

accession_numbers_gltA<-accession_numbers_Filtered %>% 
  filter(genetic_region=="gltA") 
#make it a list
accession_numbers_gltA_list<-(accession_numbers_gltA$accession_number)

accession_numbers_rpoB<-accession_numbers_Filtered %>% 
  filter(genetic_region=="rpoB") 
#make it a list
accession_numbers_rpoB_list<-(accession_numbers_rpoB$accession_number)




#split 16S list into multiple
accession_numbers_16S_1<-accession_numbers_16S_list[1:100]
accession_numbers_16S_2<-accession_numbers_16S_list[101:200]
accession_numbers_16S_3<-accession_numbers_16S_list[201:300]
accession_numbers_16S_4<-accession_numbers_16S_list[301:400]
accession_numbers_16S_5<-accession_numbers_16S_list[401:500]
accession_numbers_16S_6<-accession_numbers_16S_list[501:600]
accession_numbers_16S_7<-accession_numbers_16S_list[601:700]
accession_numbers_16S_8<-accession_numbers_16S_list[701:800]
accession_numbers_16S_9<-accession_numbers_16S_list[801:900]
accession_numbers_16S_10<-accession_numbers_16S_list[901:1000]
accession_numbers_16S_11<-accession_numbers_16S_list[1001:1100]
accession_numbers_16S_12<-accession_numbers_16S_list[1101:1200]
accession_numbers_16S_13<-accession_numbers_16S_list[1201:1265]

#split atpD list into multiple
accession_numbers_atpD_1<-accession_numbers_atpD_list[1:100]
accession_numbers_atpD_2<-accession_numbers_atpD_list[101:200]
accession_numbers_atpD_3<-accession_numbers_atpD_list[201:220]

#split glnII list into multiple
accession_numbers_glnII_1<-accession_numbers_glnII_list[1:100]
accession_numbers_glnII_2<-accession_numbers_glnII_list[101:106]

#split ITS list into multiple
accession_numbers_ITS_1<-accession_numbers_ITS_list[1:100]
accession_numbers_ITS_2<-accession_numbers_ITS_list[101:200]
accession_numbers_ITS_3<-accession_numbers_ITS_list[201:236]


#split nodA list into multiple
accession_numbers_nodA_1<-accession_numbers_nodA_list[1:100]
accession_numbers_nodA_2<-accession_numbers_nodA_list[101:200]
accession_numbers_nodA_3<-accession_numbers_nodA_list[201:247]

#split nodC list into multiple
accession_numbers_nodC_1<-accession_numbers_nodC_list[1:100]
accession_numbers_nodC_2<-accession_numbers_nodC_list[101:200]
accession_numbers_nodC_3<-accession_numbers_nodC_list[201:300]
accession_numbers_nodC_4<-accession_numbers_nodC_list[301:389]

#split recA list into multiple
accession_numbers_recA_1<-accession_numbers_recA_list[1:100]
accession_numbers_recA_2<-accession_numbers_recA_list[101:200]
accession_numbers_recA_3<-accession_numbers_recA_list[201:300]
accession_numbers_recA_4<-accession_numbers_recA_list[301:400]
accession_numbers_recA_5<-accession_numbers_recA_list[401:500]
accession_numbers_recA_6<-accession_numbers_recA_list[501:557]

#split nifH list into multiple
accession_numbers_nifH_1<-accession_numbers_nifH_list[1:100]
accession_numbers_nifH_2<-accession_numbers_nifH_list[101:200]
accession_numbers_nifH_3<-accession_numbers_nifH_list[201:300]
accession_numbers_nifH_4<-accession_numbers_nifH_list[301:324]

#16S 1:100
#connect to GenBank database and download the sequence
read.GenBank(accession_numbers_16S_1,species.name=TRUE)
#put sequence into new dataframe using fasta format
sequences_16S_1 <- entrez_fetch(id = accession_numbers_16S_1,
                         db = "nuccore", 
                         rettype = "fasta")
#write fasta file to working directory
write(sequences_16S_1,file="Updated_sequences_16S_1-100.fasta")

#16S 101:200
#connect to GenBank database and download the sequence
read.GenBank(accession_numbers_16S_2,species.name=TRUE)
#put sequence into new dataframe using fasta format
sequences_16S_2 <- entrez_fetch(id = accession_numbers_16S_2,
                                    db = "nuccore", 
                                    rettype = "fasta")
#write fasta file to working directory
write(sequences_16S_2,file="Updated_sequences_16S_101-200.fasta")

#16S 201:300
#connect to GenBank database and download the sequence
read.GenBank(accession_numbers_16S_3,species.name=TRUE)
#put sequence into new dataframe using fasta format
sequences_16S_3 <- entrez_fetch(id = accession_numbers_16S_3,
                                db = "nuccore", 
                                rettype = "fasta")
#write fasta file to working directory
write(sequences_16S_3,file="Updated_sequences_16S_201-300.fasta")

#16S 301:400
#connect to GenBank database and download the sequence
read.GenBank(accession_numbers_16S_4,species.name=TRUE)
#put sequence into new dataframe using fasta format
sequences_16S_4 <- entrez_fetch(id = accession_numbers_16S_4,
                                db = "nuccore", 
                                rettype = "fasta")
#write fasta file to working directory
write(sequences_16S_4,file="Updated_sequences_16S_301-400.fasta")

#16S 401:500
#connect to GenBank database and download the sequence
read.GenBank(accession_numbers_16S_5,species.name=TRUE)
#put sequence into new dataframe using fasta format
sequences_16S_5 <- entrez_fetch(id = accession_numbers_16S_5,
                                db = "nuccore", 
                                rettype = "fasta")
#write fasta file to working directory
write(sequences_16S_5,file="Updated_sequences_16S_401-500.fasta")

#16S 501:600
#connect to GenBank database and download the sequence
read.GenBank(accession_numbers_16S_6,species.name=TRUE)
#put sequence into new dataframe using fasta format
sequences_16S_6 <- entrez_fetch(id = accession_numbers_16S_6,
                                db = "nuccore", 
                                rettype = "fasta")
#write fasta file to working directory
write(sequences_16S_6,file="Updated_sequences_16S_501-600.fasta")

#16S 601:700
#connect to GenBank database and download the sequence
read.GenBank(accession_numbers_16S_7,species.name=TRUE)
#put sequence into new dataframe using fasta format
sequences_16S_7 <- entrez_fetch(id = accession_numbers_16S_7,
                                db = "nuccore", 
                                rettype = "fasta")
#write fasta file to working directory
write(sequences_16S_7,file="Updated_sequences_16S_601-700.fasta")

#16S 701:800
#connect to GenBank database and download the sequence
read.GenBank(accession_numbers_16S_8,species.name=TRUE)
#put sequence into new dataframe using fasta format
sequences_16S_8 <- entrez_fetch(id = accession_numbers_16S_8,
                                db = "nuccore", 
                                rettype = "fasta")
#write fasta file to working directory
write(sequences_16S_8,file="Updated_sequences_16S_701-800.fasta")

#16S 801:900
#connect to GenBank database and download the sequence
read.GenBank(accession_numbers_16S_9,species.name=TRUE)
#put sequence into new dataframe using fasta format
sequences_16S_9 <- entrez_fetch(id = accession_numbers_16S_9,
                                db = "nuccore", 
                                rettype = "fasta")
#write fasta file to working directory
write(sequences_16S_9,file="Updated_sequences_16S_801-900.fasta")

#16S 901:1000
#connect to GenBank database and download the sequence
read.GenBank(accession_numbers_16S_10,species.name=TRUE)
#put sequence into new dataframe using fasta format
sequences_16S_10 <- entrez_fetch(id = accession_numbers_16S_10,
                                db = "nuccore", 
                                rettype = "fasta")
#write fasta file to working directory
write(sequences_16S_10,file="Updated_sequences_16S_901-1000.fasta")

#16S 1001:1100
#connect to GenBank database and download the sequence
read.GenBank(accession_numbers_16S_11,species.name=TRUE)
#put sequence into new dataframe using fasta format
sequences_16S_11 <- entrez_fetch(id = accession_numbers_16S_11,
                                db = "nuccore", 
                                rettype = "fasta")
#write fasta file to working directory
write(sequences_16S_11,file="Updated_sequences_16S_1001-1100.fasta")

#16S 1101:1200
#connect to GenBank database and download the sequence
read.GenBank(accession_numbers_16S_12,species.name=TRUE)
#put sequence into new dataframe using fasta format
sequences_16S_12 <- entrez_fetch(id = accession_numbers_16S_12,
                                db = "nuccore", 
                                rettype = "fasta")
#write fasta file to working directory
write(sequences_16S_12,file="Updated_sequences_16S_1101-1200.fasta")

#16S 1201:1300
#connect to GenBank database and download the sequence
read.GenBank(accession_numbers_16S_13,species.name=TRUE)
#put sequence into new dataframe using fasta format
sequences_16S_13 <- entrez_fetch(id = accession_numbers_16S_13,
                                db = "nuccore", 
                                rettype = "fasta")
#write fasta file to working directory
write(sequences_16S_13,file="Updated_sequences_16S_1201-1300.fasta")


#23S - 1:30
#connect to GenBank database and download the sequence
read.GenBank(accession_numbers_23S_list,species.name=TRUE)
#put sequence into new dataframe using fasta format
sequences_23S_1 <- entrez_fetch(id = accession_numbers_23S_list,
                                 db = "nuccore", 
                                 rettype = "fasta")
#write fasta file to working directory
write(sequences_23S_1,file="Updated_sequences_23S_1-30.fasta")


#atpD 1:100
#connect to GenBank database and download the sequence
read.GenBank(accession_numbers_atpD_1,species.name=TRUE)
#put sequence into new dataframe using fasta format
sequences_atpD_1 <- entrez_fetch(id = accession_numbers_atpD_1,
                                db = "nuccore", 
                                rettype = "fasta")
#write fasta file to working directory
write(sequences_atpD_1,file="Updated_sequences_atpD_1-100.fasta")

#atpD 101:200
#connect to GenBank database and download the sequence
read.GenBank(accession_numbers_atpD_2,species.name=TRUE)
#put sequence into new dataframe using fasta format
sequences_atpD_2 <- entrez_fetch(id = accession_numbers_atpD_2,
                                db = "nuccore", 
                                rettype = "fasta")
#write fasta file to working directory
write(sequences_atpD_2,file="Updated_sequences_atpD_101-200.fasta")

#atpD 201:300
#connect to GenBank database and download the sequence
read.GenBank(accession_numbers_atpD_3,species.name=TRUE)
#put sequence into new dataframe using fasta format
sequences_atpD_3 <- entrez_fetch(id = accession_numbers_atpD_3,
                                db = "nuccore", 
                                rettype = "fasta")
#write fasta file to working directory
write(sequences_atpD_3,file="Updated_sequences_atpD_201-300.fasta")

#cpn60 - 1:24
#connect to GenBank database and download the sequence
read.GenBank(accession_numbers_cpn60_list,species.name=TRUE)
#put sequence into new dataframe using fasta format
sequences_cpn60_1 <- entrez_fetch(id = accession_numbers_cpn60_list,
                                db = "nuccore", 
                                rettype = "fasta")
#write fasta file to working directory
write(sequences_cpn60_1,file="Updated_sequences_cpn60_1-24.fasta")

#dnaK - 1:86
#connect to GenBank database and download the sequence
read.GenBank(accession_numbers_dnaK_list,species.name=TRUE)
#put sequence into new dataframe using fasta format
sequences_dnaK_1 <- entrez_fetch(id = accession_numbers_dnaK_list,
                                  db = "nuccore", 
                                  rettype = "fasta")
#write fasta file to working directory
write(sequences_dnaK_1,file="Updated_sequences_dnaK_1-86.fasta")

#glnA - 1:36
#connect to GenBank database and download the sequence
read.GenBank(accession_numbers_glnA_list,species.name=TRUE)
#put sequence into new dataframe using fasta format
sequences_glnA_1 <- entrez_fetch(id = accession_numbers_glnA_list,
                                 db = "nuccore", 
                                 rettype = "fasta")
#write fasta file to working directory
write(sequences_glnA_1,file="Updated_sequences_glnA_1-36.fasta")

#glnB - 1:37
#connect to GenBank database and download the sequence
read.GenBank(accession_numbers_glnB_list,species.name=TRUE)
#put sequence into new dataframe using fasta format
sequences_glnB_1 <- entrez_fetch(id = accession_numbers_glnB_list,
                                 db = "nuccore", 
                                 rettype = "fasta")
#write fasta file to working directory
write(sequences_glnB_1,file="Updated_sequences_glnB_1-37.fasta")

#glnII - 1:100
#connect to GenBank database and download the sequence
read.GenBank(accession_numbers_glnII_1,species.name=TRUE)
#put sequence into new dataframe using fasta format
sequences_glnII_1 <- entrez_fetch(id = accession_numbers_glnII_1,
                                 db = "nuccore", 
                                 rettype = "fasta")
#write fasta file to working directory
write(sequences_glnII_1,file="Updated_sequences_glnII_1-100.fasta")

#glnII - 101:106
#connect to GenBank database and download the sequence
read.GenBank(accession_numbers_glnII_2,species.name=TRUE)
#put sequence into new dataframe using fasta format
sequences_glnII_2 <- entrez_fetch(id = accession_numbers_glnII_2,
                                  db = "nuccore", 
                                  rettype = "fasta")
#write fasta file to working directory
write(sequences_glnII_2,file="Updated_sequences_glnII_100-106.fasta")

#gltA - 1:36
#connect to GenBank database and download the sequence
read.GenBank(accession_numbers_gltA_list,species.name=TRUE)
#put sequence into new dataframe using fasta format
sequences_gltA_1 <- entrez_fetch(id = accession_numbers_gltA_list,
                                 db = "nuccore", 
                                 rettype = "fasta")
#write fasta file to working directory
write(sequences_gltA_1,file="Updated_sequences_gltA_1-36.fasta")

#gyrB - 1:67
#connect to GenBank database and download the sequence
read.GenBank(accession_numbers_gyrB_list,species.name=TRUE)
#put sequence into new dataframe using fasta format
sequences_gyrB_1 <- entrez_fetch(id = accession_numbers_gyrB_list,
                                 db = "nuccore", 
                                 rettype = "fasta")
#write fasta file to working directory
write(sequences_gyrB_1,file="Updated_sequences_gyrB_1-67.fasta")

#ITS 1:100
#connect to GenBank database and download the sequence
read.GenBank(accession_numbers_ITS_1,species.name=TRUE)
#put sequence into new dataframe using fasta format
sequences_ITS_1 <- entrez_fetch(id = accession_numbers_ITS_1,
                                db = "nuccore", 
                                rettype = "fasta")
#write fasta file to working directory
write(sequences_ITS_1,file="Updated_sequences_ITS_1-100.fasta")

#ITS 101:200
#connect to GenBank database and download the sequence
read.GenBank(accession_numbers_ITS_2,species.name=TRUE)
#put sequence into new dataframe using fasta format
sequences_ITS_2 <- entrez_fetch(id = accession_numbers_ITS_2,
                                db = "nuccore", 
                                rettype = "fasta")
#write fasta file to working directory
write(sequences_ITS_2,file="Updated_sequences_ITS_101-200.fasta")

#ITS 201:244
#connect to GenBank database and download the sequence
read.GenBank(accession_numbers_ITS_3,species.name=TRUE)
#put sequence into new dataframe using fasta format
sequences_ITS_3 <- entrez_fetch(id = accession_numbers_ITS_3,
                                db = "nuccore", 
                                rettype = "fasta")
#write fasta file to working directory
write(sequences_ITS_3,file="Updated_sequences_ITS_201-244.fasta")

#nifD - 1:41
#connect to GenBank database and download the sequence
read.GenBank(accession_numbers_nifD_list,species.name=TRUE)
#put sequence into new dataframe using fasta format
sequences_nifD_1 <- entrez_fetch(id = accession_numbers_nifD_list,
                                 db = "nuccore", 
                                 rettype = "fasta")
#write fasta file to working directory
write(sequences_nifD_1,file="Updated_sequences_nifD_1-40.fasta")

#nifH 1:100
#connect to GenBank database and download the sequence
read.GenBank(accession_numbers_nifH_1,species.name=TRUE)
#put sequence into new dataframe using fasta format
sequences_nifH_1 <- entrez_fetch(id = accession_numbers_nifH_1,
                                db = "nuccore", 
                                rettype = "fasta")
#write fasta file to working directory
write(sequences_nifH_1,file="Updated_sequences_nifH_1-100.fasta")

#nifH 101:200
#connect to GenBank database and download the sequence
read.GenBank(accession_numbers_nifH_2,species.name=TRUE)
#put sequence into new dataframe using fasta format
sequences_nifH_2 <- entrez_fetch(id = accession_numbers_nifH_2,
                                db = "nuccore", 
                                rettype = "fasta")
#write fasta file to working directory
write(sequences_nifH_2,file="Updated_sequences_nifH_101-200.fasta")

#nifH 201:300
#connect to GenBank database and download the sequence
read.GenBank(accession_numbers_nifH_3,species.name=TRUE)
#put sequence into new dataframe using fasta format
sequences_nifH_3 <- entrez_fetch(id = accession_numbers_nifH_3,
                                db = "nuccore", 
                                rettype = "fasta")
#write fasta file to working directory
write(sequences_nifH_3,file="Updated_sequences_nifH_201-300.fasta")

#nifH 301:339
#connect to GenBank database and download the sequence
read.GenBank(accession_numbers_nifH_4,species.name=TRUE)
#put sequence into new dataframe using fasta format
sequences_nifH_4 <- entrez_fetch(id = accession_numbers_nifH_4,
                                db = "nuccore", 
                                rettype = "fasta")
#write fasta file to working directory
write(sequences_nifH_4,file="Updated_sequences_nifH_301-339.fasta")

#nifHD - 1:1
#connect to GenBank database and download the sequence
read.GenBank(accession_numbers_nifHD_list,species.name=TRUE)
#put sequence into new dataframe using fasta format
sequences_nifHD_1 <- entrez_fetch(id = accession_numbers_nifHD_list,
                                 db = "nuccore", 
                                 rettype = "fasta")
#write fasta file to working directory
write(sequences_nifHD_1,file="Updated_sequences_nifHD_1.fasta")

#nifD - 1:41
#connect to GenBank database and download the sequence
read.GenBank(accession_numbers_nifD_list,species.name=TRUE)
#put sequence into new dataframe using fasta format
sequences_nifD_1 <- entrez_fetch(id = accession_numbers_nifD_list,
                                 db = "nuccore", 
                                 rettype = "fasta")
#write fasta file to working directory
write(sequences_nifD_1,file="Updated_sequences_nifD_1-40.fasta")

#nodA 1:100
#connect to GenBank database and download the sequence
read.GenBank(accession_numbers_nodA_1,species.name=TRUE)
#put sequence into new dataframe using fasta format
sequences_nodA_1 <- entrez_fetch(id = accession_numbers_nodA_1,
                                 db = "nuccore", 
                                 rettype = "fasta")
#write fasta file to working directory
write(sequences_nodA_1,file="Updated_sequences_nodA_1-100.fasta")

#nodA 101:200
#connect to GenBank database and download the sequence
read.GenBank(accession_numbers_nodA_2,species.name=TRUE)
#put sequence into new dataframe using fasta format
sequences_nodA_2 <- entrez_fetch(id = accession_numbers_nodA_2,
                                 db = "nuccore", 
                                 rettype = "fasta")
#write fasta file to working directory
write(sequences_nodA_2,file="Updated_sequences_nodA_101-200.fasta")

#nodA 201:254
#connect to GenBank database and download the sequence
read.GenBank(accession_numbers_nodA_3,species.name=TRUE)
#put sequence into new dataframe using fasta format
sequences_nodA_3 <- entrez_fetch(id = accession_numbers_nodA_3,
                                 db = "nuccore", 
                                 rettype = "fasta")
#write fasta file to working directory
write(sequences_nodA_3,file="Updated_sequences_nodA_201-254.fasta")

#nodC 1:100
#connect to GenBank database and download the sequence
read.GenBank(accession_numbers_nodC_1,species.name=TRUE)
#put sequence into new dataframe using fasta format
sequences_nodC_1 <- entrez_fetch(id = accession_numbers_nodC_1,
                                 db = "nuccore", 
                                 rettype = "fasta")
#write fasta file to working directory
write(sequences_nodC_1,file="Updated_sequences_nodC_1-100.fasta")

#nodC 101:200
#connect to GenBank database and download the sequence
read.GenBank(accession_numbers_nodC_2,species.name=TRUE)
#put sequence into new dataframe using fasta format
sequences_nodC_2 <- entrez_fetch(id = accession_numbers_nodC_2,
                                 db = "nuccore", 
                                 rettype = "fasta")
#write fasta file to working directory
write(sequences_nodC_2,file="Updated_sequences_nodC_101-200.fasta")

#nodC 201:300
#connect to GenBank database and download the sequence
read.GenBank(accession_numbers_nodC_3,species.name=TRUE)
#put sequence into new dataframe using fasta format
sequences_nodC_3 <- entrez_fetch(id = accession_numbers_nodC_3,
                                 db = "nuccore", 
                                 rettype = "fasta")
#write fasta file to working directory
write(sequences_nodC_3,file="Updated_sequences_nodC_201-300.fasta")

#nodC 301:400
#connect to GenBank database and download the sequence
read.GenBank(accession_numbers_nodC_4,species.name=TRUE)
#put sequence into new dataframe using fasta format
sequences_nodC_4 <- entrez_fetch(id = accession_numbers_nodC_4,
                                db = "nuccore", 
                                rettype = "fasta")
#write fasta file to working directory
write(sequences_nodC_4,file="Updated_sequences_nodC_301-400.fasta")

#nodC 401:500
#connect to GenBank database and download the sequence
read.GenBank(accession_numbers_nodC_5,species.name=TRUE)
#put sequence into new dataframe using fasta format
sequences_nodC_5 <- entrez_fetch(id = accession_numbers_nodC_5,
                                db = "nuccore", 
                                rettype = "fasta")
#write fasta file to working directory
write(sequences_nodC_5,file="Updated_sequences_nodC_401-500.fasta")

#nodD - 1:14
#connect to GenBank database and download the sequence
read.GenBank(accession_numbers_nodD_list,species.name=TRUE)
#put sequence into new dataframe using fasta format
sequences_nodD_1 <- entrez_fetch(id = accession_numbers_nodD_list,
                                 db = "nuccore", 
                                 rettype = "fasta")
#write fasta file to working directory
write(sequences_nodD_1,file="Updated_sequences_nodD_1-14.fasta")

#nodZ - 1:17
#connect to GenBank database and download the sequence
read.GenBank(accession_numbers_nodZ_list,species.name=TRUE)
#put sequence into new dataframe using fasta format
sequences_nodZ_1 <- entrez_fetch(id = accession_numbers_nodZ_list,
                                 db = "nuccore", 
                                 rettype = "fasta")
#write fasta file to working directory
write(sequences_nodZ_1,file="Updated_sequences_nodZ_1-17.fasta")

#noeI - 1:18
#connect to GenBank database and download the sequence
read.GenBank(accession_numbers_noeI_list,species.name=TRUE)
#put sequence into new dataframe using fasta format
sequences_noeI_1 <- entrez_fetch(id = accession_numbers_noeI_list,
                                 db = "nuccore", 
                                 rettype = "fasta")
#write fasta file to working directory
write(sequences_noeI_1,file="Updated_sequences_noeI_1-18.fasta")

#recA 1:100
#connect to GenBank database and download the sequence
read.GenBank(accession_numbers_recA_1,species.name=TRUE)
#put sequence into new dataframe using fasta format
sequences_recA_1 <- entrez_fetch(id = accession_numbers_recA_1,
                                 db = "nuccore", 
                                 rettype = "fasta")
#write fasta file to working directory
write(sequences_recA_1,file="Updated_sequences_recA_1-100.fasta")

#recA 101:200
#connect to GenBank database and download the sequence
read.GenBank(accession_numbers_recA_2,species.name=TRUE)
#put sequence into new dataframe using fasta format
sequences_recA_2 <- entrez_fetch(id = accession_numbers_recA_2,
                                 db = "nuccore", 
                                 rettype = "fasta")
#write fasta file to working directory
write(sequences_recA_2,file="Updated_sequences_recA_101-200.fasta")

#recA 201:300
#connect to GenBank database and download the sequence
read.GenBank(accession_numbers_recA_3,species.name=TRUE)
#put sequence into new dataframe using fasta format
sequences_recA_3 <- entrez_fetch(id = accession_numbers_recA_3,
                                 db = "nuccore", 
                                 rettype = "fasta")
#write fasta file to working directory
write(sequences_recA_3,file="Updated_sequences_recA_201-300.fasta")

#recA 301:400
#connect to GenBank database and download the sequence
read.GenBank(accession_numbers_recA_4,species.name=TRUE)
#put sequence into new dataframe using fasta format
sequences_recA_4 <- entrez_fetch(id = accession_numbers_recA_4,
                                 db = "nuccore", 
                                 rettype = "fasta")
#write fasta file to working directory
write(sequences_recA_4,file="Updated_sequences_recA_301-400.fasta")

#nodC 401:500
#connect to GenBank database and download the sequence
read.GenBank(accession_numbers_recA_5,species.name=TRUE)
#put sequence into new dataframe using fasta format
sequences_recA_5 <- entrez_fetch(id = accession_numbers_recA_5,
                                 db = "nuccore", 
                                 rettype = "fasta")
#write fasta file to working directory
write(sequences_recA_5,file="Updated_sequences_recA_401-500.fasta")

#recA 501:573
#connect to GenBank database and download the sequence
read.GenBank(accession_numbers_recA_6,species.name=TRUE)
#put sequence into new dataframe using fasta format
sequences_recA_6 <- entrez_fetch(id = accession_numbers_recA_6,
                                db = "nuccore", 
                                rettype = "fasta")
#write fasta file to working directory
write(sequences_recA_6,file="Updated_sequences_recA_501-573.fasta")

#rpoB - 1:39
#connect to GenBank database and download the sequence
read.GenBank(accession_numbers_rpoB_list,species.name=TRUE)
#put sequence into new dataframe using fasta format
sequences_rpoB_1 <- entrez_fetch(id = accession_numbers_rpoB_list,
                                 db = "nuccore", 
                                 rettype = "fasta")
#write fasta file to working directory
write(sequences_rpoB_1,file="Updated_sequences_rpoB_1-39.fasta")

#rrs - 1:9
#connect to GenBank database and download the sequence
read.GenBank(accession_numbers_rrs_list,species.name=TRUE)
#put sequence into new dataframe using fasta format
sequences_rrs_1 <- entrez_fetch(id = accession_numbers_rrs_list,
                                 db = "nuccore", 
                                 rettype = "fasta")
#write fasta file to working directory
write(sequences_rrs_1,file="Updated_sequences_rrs_1-9.fasta")

#thrC - 1:36
#connect to GenBank database and download the sequence
read.GenBank(accession_numbers_thrC_list,species.name=TRUE)
#put sequence into new dataframe using fasta format
sequences_thrC_1 <- entrez_fetch(id = accession_numbers_thrC_list,
                                db = "nuccore", 
                                rettype = "fasta")
#write fasta file to working directory
write(sequences_thrC_1,file="Updated_sequences_thrC_1-36.fasta")

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
  select(paper_id, genus_species, sample_country,diversity_sum,Plant_status, num_nodules) %>% 
  filter(Plant_status!=999)%>%
  filter(Plant_status!="introduced-native")%>%
  filter(Plant_status!="extinct") %>% 
  na.omit()

#Make a data frame calculating the mean, standard error, and number of observations for native vs. non-native to make graph
Average_rhizobial_symbionts_genus_sp<-Rhizobial_Associations_Condensed%>%
  filter(diversity_sum!=999)%>%
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

hist(Rhizobial_Associations_Condensed$diversity_sum)

#Run anova accounting for how many species using averaged number 
summary (Mixed_Model_Rhiz_Status<- lmer(symbionts_Mean ~ Plant_status + (1 | species_status_n), data = Average_rhizobial_symbionts_genus_sp))
anova(Mixed_Model_Rhiz_Status,type = 2)

###Run ANOVA using non-averaged numbers + species and nod number as a fixed effect 
#### Question: is this how you do multiple fixed effects? #### 
#summary(Mixed_Model_Rhiz_Status_NoMean<- lmer(diversity_sum ~ genus_species*num_nodules + (1 | Plant_status), data = Rhizobial_Associations_Condensed)) #seems incorrect to do species as fixed effect? putting that as random effect below
#Anova(Mixed_Model_Rhiz_Status_NoMean,type = 2) 

#random effect -- i think this is correct?! get a warning 
summary(Mixed_Model_Rhiz_Status_NoMean_Random<- lmer(diversity_sum ~ Plant_status + (1 | genus_species) + (1 | num_nodules), data = Rhizobial_Associations_Condensed)) 

#anova accounting for type II error - for unbalanced data or interaction terms that are significant 
anova(Mixed_Model_Rhiz_Status_NoMean_Random,type = 2) 



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
  select(paper_id, genus_species, sample_country,diversity_sum,annual_perennial,growth_form,num_nodules,Plant_status) %>% 
  na.omit() %>% 
  filter(Plant_status %in% c("native","non_native"))


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


####Anova for growth form ####
summary(Mixed_Model_Rhiz_Status_Growth_Form<- lmer(diversity_sum ~ growth_form + (1|num_nodules), data = Rhizobial_Associations_GF_Condensed)) #
anova(Mixed_Model_Rhiz_Status_Growth_Form,type=2) 

#t.test(Rhizobial_Associations_GF_Condensed$diversity_sum~Rhizobial_Associations_GF_Condensed$growth_form)

#### Rhizobial Associates by Annual/Perennial  ####

#Make a data frame calculating the mean, standard error, and number of observations for native vs. non-native to make graph
Average_annual_perennial_Graph<-Rhizobial_Associations_GF_Condensed%>%
  filter(diversity_sum!=999)%>%
  group_by(annual_perennial,Plant_status)%>%
  summarize(symbionts_Std=sd(diversity_sum),symbionts_Mean=mean(diversity_sum),symbionts_n=length(diversity_sum))%>%
  #Make a new column called "Richness_St_Error" and divide "Richness_Std" by the square root of "Richness_n"
  mutate(symbionts_St_Error=symbionts_Std/sqrt(symbionts_n))

#Make graph showing number of introduced vs. native
#### switch with graph below #####
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


####Anova for annual/perennial ####
summary(Mixed_Model_Rhiz_Status_A_P<- lmer(diversity_sum ~ annual_perennial + (1|num_nodules), data = Rhizobial_Associations_GF_Condensed)) #
anova(Mixed_Model_Rhiz_Status_A_P,type=2) 

#t.test(Rhizobial_Associations_GF_Condensed$diversity_sum~Rhizobial_Associations_GF_Condensed$annual_perennial) 

####Anova for annual/perennial*native/nonnative ####
summary(Mixed_Model_Rhiz_Status_A_P_Plant_status<- lmer(diversity_sum ~ annual_perennial*Plant_status + (1 | genus_species) + (1|num_nodules), data = Rhizobial_Associations_GF_Condensed))
#type 3 is used because we have an interaction
anova(Mixed_Model_Rhiz_Status_A_P_Plant_status,type = 3) 

#figure for above mixed model
#Make graph showing number of introduced vs. native and annual perennial
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



####Anova for annual/perennial*native/nonnative*growth form ####
summary(Mixed_Model_Rhiz_Status_A_P_Status_GF<- lmer(diversity_sum ~ annual_perennial*Plant_status + growth_form*Plant_status + (1 | genus_species) + (1|num_nodules), data = Rhizobial_Associations_GF_Condensed)) #
anova(Mixed_Model_Rhiz_Status_A_P_Status_GF,type = 3) 

#Make a data frame calculating the mean, standard error, and number of observations for native vs. non-native to make graph
Average_status_herb_Graph<-Rhizobial_Associations_GF_Condensed %>%
  filter(diversity_sum!=999) %>%
  mutate(Growth_Status=paste(growth_form,annual_perennial,sep = "::")) %>% 
  group_by(Growth_Status,Plant_status) %>%
  summarize(symbionts_Std=sd(diversity_sum),symbionts_Mean=mean(diversity_sum),symbionts_n=length(diversity_sum)) %>%
  #Make a new column called "Richness_St_Error" and divide "Richness_Std" by the square root of "Richness_n"
  mutate(symbionts_St_Error=symbionts_Std/sqrt(symbionts_n))

#Make graph showing number of introduced vs. native and herb/woody
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




#### Number Nodules (sample number) by strain richness ####

Nod_by_strain<-Plant_Data %>%
  filter(num_nodules!=999)%>%
  filter(strain_richness!=999,strain_richness!="NEED SUPPLEMENTAL PAPER")%>%
  filter(genus_species!=999)%>%
  filter(!is.na(num_nodules),!is.na(strain_richness))%>% 
  mutate(strain_richness=as.numeric(strain_richness)) %>% 
  filter(num_nodules>10) %>% 
  group_by(genus_species) %>% 
  summarize(num_nod_Mean=mean(num_nodules),num_nod_n=length(num_nodules),strain_rich_Mean=mean(strain_richness),strain_rich_n=length(strain_richness))%>%
  ungroup()

hist(subset_nod_by_strain$num_nod_Mean) 

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

#Make graph showing number of introduced vs. native
ggplot(Strain_rich_graph,aes(x=Plant_status,y=strain_rich_Mean2,fill=Plant_status))+
  #Make a bar graph where the height of the bars is equal to the data (stat=identity) and you preserve the vertical position while adjusting the horizontal(position_dodge), and outline the bars with the color black.
  geom_bar(stat="identity", position=position_dodge(),color="black")+
  geom_errorbar(aes(ymin=strain_rich_Mean2-strain_rich_St_Error,ymax=strain_rich_Mean2+strain_rich_St_Error),position=position_dodge(0.9),width=0.2)+
  #make error bars using the Standard error from the mean and place them at 0.9 with a width of 0.2
  #Label the x-axis "Watershed"
  xlab("Plant Species Status")+
  #Label the y-axis "Species Richness"
  ylab("Average Strain Richness")+
  #Fill the bar graphs with grey and white according and label them "Inside Fence" and "Outside Fence"
  scale_fill_manual(values=c("darkslategrey","cadetblue4"), labels=c("Native","Non-native"))  

hist(Rhizobial_Associations_Condensed$diversity_sum)

#Run anova accounting for how many species and how many nodules they sample
#### run this again without mean and accounting for nodule number #### did this above
summary (Mixed_Model_Rhiz_Status_Strain_rich<- lmer(strain_rich_Mean ~ Plant_status + (1 | species_status_n) + (1 | num_nod_Mean), data = Rhiz_symbionts_Nod_Num))
anova(Mixed_Model_Rhiz_Status_Strain_rich,type = 2)

#### run this again without mean and accounting for nodule number #### still need to merge symbionts to data frame where mean isn't taken to change this
summary (Mixed_Model_Rhiz_Status_symbionts<- lmer(symbionts_Mean ~ Plant_status + (1 | species_status_n) + (1 | num_nod_Mean), data = (Rhiz_symbionts_Nod_Num)))
anova(Mixed_Model_Rhiz_Status_symbionts,type = 2)

