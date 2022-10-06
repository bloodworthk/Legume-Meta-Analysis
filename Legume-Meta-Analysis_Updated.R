#### Legume Meta Analysis ####
#Created by Kathryn Bloodworth, Smriti Pehim Limbu, and Kimberly Komatsu

#### Set working directory #### 
#should be one drive folder

#Kathryn Bloodworth - PC
setwd("C:/Users/kjbloodw/OneDrive - UNCG/Invasive Legume Meta-Analysis")

#### Install and load in packages ####
#If you have a PC just installing and loading metagear should allow you to run everything needed if you allow it to install needed packages
#If you have a mac, this will be different and I have a series of other packages that will need to be installed and run. Let me know if you need this.

#install.packages("metagear")
library(metagear)


#### Load in Data ####

#read in data frame with the Web of Science search results from Oct. 2022
WoS_Abstracts_2022<-read.csv("Methods/Web_of_Science_10_05_Search_Results.csv")

#### Setting up Screening Tasks - DO NOT RUN ANY OF THIS AGAIN ####

#This step adds 4 new columns (Study ID, Reviewer, 2 columns for include/or dont include)
References_Screening_Ready<-effort_initialize(WoS_Abstracts_2022)

#Create dataframe with names of reviewers
theTeam<-c("Kathryn","Kim","Smriti")

#Randomly delegate screening efforts between the team of reviewers evenly and then save these as seperate files to the directory given
References_unscreened<-effort_distribute(References_Screening_Ready,reviewers=theTeam, save_split = TRUE,directory = "C:/Users/kjbloodw/OneDrive - UNCG/Invasive Legume Meta-Analysis/Methods/2022_Reference_Screening")

#### Screen Abstracts - START HERE ####
#KIM AND SMRITI START HERE: Run only the line that applies to you (has your csv file in it)

abstract_screener(file=file.choose("C:/Users/kjbloodw/OneDrive - UNCG/Invasive Legume Meta-Analysis/Methods/2022_Reference_Screening/effort_Kathryn.csv"),
                  aReviewer = "Kathryn",
                  reviewerColumnName = "REVIEWERS",
                  unscreenedColumnName = "INCLUDE",
                  unscreenedValue = "not vetted",
                  abstractColumnName = "Abstract",
                  titleColumnName = "Article.Title",
                  browserSearch = "https://www.google.com/search?q=",
                  fontSize = 13,
                  windowWidth = 70,
                  windowHeight = 16,
                  theButtons = c("YES","maybe","NO"),
                  keyBindingToButtons = c("y","m","n"),
                  buttonSize = 10,
                  highlightColor = "powderblue",
                  highlightKeywords = c("legume","native","invasive","non-native","non native","nonnative","rhizobia")
)

