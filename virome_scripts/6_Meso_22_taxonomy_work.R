#Mesocosm 2022 Taxonomy anaysis 
#For Genera Taxonomic level

#The goal is to compare how genera differ across populations and time points and to assess diversity metrics

#Load libraries
library(ggplot2)
library(ggVennDiagram)
library(RVenn)

#set working directory
setwd("~/Library/CloudStorage/GoogleDrive-sbirch1@charlotte.edu/My Drive/Postdoc_work/Mesocosom_exps/Mesocosom_Summer_2022/Transcriptomics/updated_Blast-Taxonomy")

#venndiagram resources: 
#https://r-charts.com/part-whole/ggvenndiagram/
#https://gaospecial.github.io/ggVennDiagram/
#https://cran.r-project.org/web/packages/RVenn/vignettes/vignette.html
#https://peat-clark.github.io/BIO381/veganTutorial.html


#import full dataset 
genus_dat <- read.csv("~/Library/CloudStorage/GoogleDrive-sbirch1@charlotte.edu/My Drive/Postdoc_work/Mesocosom_exps/Mesocosom_Summer_2022/Transcriptomics/updated_Blast-Taxonomy/genus_dat_v2.csv", stringsAsFactors=TRUE)


#### Location comparisons T0 vs T14 + Control (Field T14) ####

#NH T0 vs T14
NS_field <- list(T0 = genus_dat$NS.T0, T14=genus_dat$NS.T14, ctrl_Field = genus_dat$FIELD.T14)
str(NS_field)
ggVennDiagram(NS_field,category.names=c("NS T0", "NS T14", "Field")) + scale_fill_gradient(low="#3D405B",high = "#E07A5F") + ggtitle("Comparison of Genera present in NS T0 vs NS T14 vs Field")

#ME T0 vs T14
ME_field <- list(T0 = genus_dat$ME.T0, T14=genus_dat$ME.T14, ctrl_Field = genus_dat$FIELD.T14)
str(ME_field)
ggVennDiagram(ME_field,category.names=c("ME T0", "ME T14", "Field")) + scale_fill_gradient(low="#3D405B",high = "#E07A5F") + ggtitle("Comparison of Genera present in ME T0 vs ME T14 vs Field")

#NH T0 vs T14
NH_field <- list(T0 = genus_dat$NH.T0, T14=genus_dat$NH.T14, ctrl_Field = genus_dat$FIELD.T14)
str(NH_field)
ggVennDiagram(NH_field,category.names=c("NH T0", "NH T14", "Field")) + scale_fill_gradient(low="#3D405B",high = "#E07A5F") + ggtitle("Comparison of Genera present in NH T0 vs NH T14 vs Field")

#MA T0 vs T14
MA_field <- list(T0 = genus_dat$MA.T0, T14=genus_dat$MA.T14, ctrl_Field = genus_dat$FIELD.T14)
str(MA_field)
ggVennDiagram(MA_field,category.names=c("MA T0", "MA T14", "Field")) + scale_fill_gradient(low="#3D405B",high = "#E07A5F") + ggtitle("Comparison of Genera present in MA T0 vs MA T14 vs Field")

#SC T0 vs T14
SC_field <- list(T0 = genus_dat$SC.T0, T14=genus_dat$SC.T14, ctrl_Field = genus_dat$FIELD.T14)
str(SC_field)
ggVennDiagram(SC_field,category.names=c("SC T0", "SC T14", "Field")) + scale_fill_gradient(low="#3D405B",high = "#E07A5F") + ggtitle("Comparison of Genera present in SC T0 vs SC T14 vs Field")

#FL T0 vs T14
FL_field <- list(T0 = genus_dat$FL.T0, T14=genus_dat$FL.T14, ctrl_Field = genus_dat$FIELD.T14)
str(FL_field)
ggVennDiagram(FL_field,category.names=c("FL T0", "FL T14", "Field")) + scale_fill_gradient(low= "#3D405B",high = "#E07A5F") + ggtitle("Comparison of Genera present in FL T0 vs FL T14 vs Field")




#### Get members of intersections ####
library(RVenn)

genus_dat_v2 <- read.csv("~/Library/CloudStorage/GoogleDrive-sbirch1@charlotte.edu/My Drive/Postdoc_work/Mesocosom_exps/Mesocosom_Summer_2022/Transcriptomics/updated_Blast-Taxonomy/genus_dat_v2.csv")
str(genus_dat_v2)

#make lists of each of the three-way comparisons 
NS_dat = list(NS_T0 = genus_dat_v2$NS.T0,NS_T14 = genus_dat_v2$NS.T14, Field = genus_dat_v2$FIELD.T14) 
ME_dat = list(ME_T0 = genus_dat_v2$ME.T0,ME_T14 = genus_dat_v2$ME.T14, Field = genus_dat_v2$FIELD.T14)
NH_dat = list(NH_T0 = genus_dat_v2$NH.T0,NH_T14 = genus_dat_v2$NH.T14, Field = genus_dat_v2$FIELD.T14)
MA_dat = list(MA_T0 = genus_dat_v2$MA.T0,MA_T14 = genus_dat_v2$MA.T14, Field = genus_dat_v2$FIELD.T14)
SC_dat = list(SC_T0 = genus_dat_v2$SC.T0,SC_T14 = genus_dat_v2$SC.T14, Field = genus_dat_v2$FIELD.T14)
FL_dat = list(FL_T0 = genus_dat_v2$FL.T0,FL_T14 = genus_dat_v2$FL.T14, Field = genus_dat_v2$FIELD.T14)
All_T14_dat = list(NS_T14 = genus_dat_v2$NS.T14, ME_T14 = genus_dat_v2$ME.T14, NH_T14 = genus_dat_v2$NH.T14, MA_T14 = genus_dat_v2$MA.T14, SC_T14 = genus_dat_v2$SC.T14, FL_T14 = genus_dat_v2$FL.T14, Field = genus_dat_v2$FIELD.T14)
All_T0_dat = list(NS_T0 = genus_dat_v2$NS.T0, ME_T0 = genus_dat_v2$ME.T0, NH_T0 = genus_dat_v2$NH.T0, MA_T0 = genus_dat_v2$MA.T0, SC_T0 = genus_dat_v2$SC.T0, FL_T0 = genus_dat_v2$FL.T0, Field = genus_dat_v2$FIELD.T14)
str(NS_dat)

#to get rid of empty spaces: 
library(vctrs)
ME_dat_empty <-list_drop_empty(ME_dat)


#make venn object
All_14_venn<-Venn(All_T14_dat) #For all T14 locations
All_T0_venn<-Venn(All_T0_dat) #For all T0 locations + Field
Total_venn<-Venn(genus_dat_v2) #For all locations
NS_venn<-Venn(NS_dat) 
ME_venn<-Venn(ME_dat) 
NH_venn<-Venn(NH_dat) 
MA_venn<-Venn(MA_dat) 
SC_venn<-Venn(SC_dat) 
FL_venn<-Venn(FL_dat) 


# make a venn diagram with this program 
#ggvenn(test_venn) 
ggvenn(NS_venn) 
ggvenn(ME_venn)
ggvenn(NH_venn)
ggvenn(MA_venn)
ggvenn(SC_venn)
ggvenn(FL_venn)


#get overlap of all sets
tot_overlap<-overlap(Total_venn) #total overlap for all locations (13 genera: "g__Alphabaculovirus","g__Betabaculovirus","g__unclassified Caudoviricetes genus" ,"g__Avipoxvirus" ,"g__Chlorovirus","g__Prasinovirus"  "g__unclassified Phycodnaviridae genus"  "g__Rheavirus"  "g__unclassified Marseilleviridae genus" "g__Betaentomopoxvirus"  "g__Marseillevirus"   "g__Oceanusvirus"  "g__Tupanvirus")
tot_overlap_T0s<-overlap(All_T0_venn) #tot overlap for all T0 locations (13 genera)
tot_overlap_T14s<-overlap(All_14_venn) #tot overlap for all T14 locations (97 genera)
tot_overlap_NS<-overlap(NS_venn)
tot_overlap_ME<-overlap(ME_venn)
tot_overlap_NH<-overlap(NH_venn)
tot_overlap_MA<-overlap(MA_venn)
tot_overlap_SC<-overlap(SC_venn)
tot_overlap_FL<-overlap(FL_venn)


#get total pairwise overlap of each sample
NS_T0vT14_tot <-overlap(NS_venn, c("NS_T0","NS_T14"))
NS_T0vField_tot<-overlap(NS_venn,c("NS_T0","Field"))
NS_T14vField_tot<-overlap(NS_venn,c("NS_T14","Field"))

ME_T0vT14_tot <-overlap(ME_venn, c("ME_T0","ME_T14"))
ME_T0vField_tot<-overlap(ME_venn,c("ME_T0","Field"))
ME_T14vField_tot<-overlap(ME_venn,c("ME_T14","Field"))

NH_T0vT14_tot <-overlap(NH_venn, c("NH_T0","NH_T14"))
NH_T0vField_tot<-overlap(NH_venn,c("NH_T0","Field"))
NH_T14vField_tot<-overlap(NH_venn,c("NH_T14","Field"))

MA_T0vT14_tot <-overlap(MA_venn, c("MA_T0","MA_T14"))
MA_T0vField_tot<-overlap(MA_venn,c("MA_T0","Field"))
MA_T14vField_tot<-overlap(MA_venn,c("MA_T14","Field"))

SC_T0vT14_tot <-overlap(SC_venn, c("SC_T0","SC_T14"))
SC_T0vField_tot<-overlap(SC_venn,c("SC_T0","Field"))
SC_T14vField_tot<-overlap(SC_venn,c("SC_T14","Field"))

FL_T0vT14_tot <-overlap(FL_venn, c("FL_T0","FL_T14"))
FL_T0vField_tot<-overlap(FL_venn,c("FL_T0","Field"))
FL_T14vField_tot<-overlap(FL_venn,c("FL_T14","Field"))


#make new list of total overlap and each pairwise overlap
NS_list = list(tot = tot_overlap_NS, T0v14 = NS_T0vT14_tot, T0vF = NS_T0vField_tot, T14vF = NS_T14vField_tot)
NS_list_venn <- Venn(NS_list)

ME_list = list(tot = tot_overlap_ME, T0v14 = ME_T0vT14_tot, T0vF = ME_T0vField_tot, T14vF = ME_T14vField_tot)
ME_list_venn <- Venn(ME_list)

NH_list = list(tot = tot_overlap_NH, T0v14 = NH_T0vT14_tot, T0vF = NH_T0vField_tot, T14vF = NH_T14vField_tot)
NH_list_venn <- Venn(NH_list)

MA_list = list(tot = tot_overlap_MA, T0v14 = MA_T0vT14_tot, T0vF = MA_T0vField_tot, T14vF = MA_T14vField_tot)
MA_list_venn <- Venn(MA_list)

SC_list = list(tot = tot_overlap_SC, T0v14 = SC_T0vT14_tot, T0vF = SC_T0vField_tot, T14vF = SC_T14vField_tot)
SC_list_venn <- Venn(SC_list)

FL_list = list(tot = tot_overlap_FL, T0v14 = FL_T0vT14_tot, T0vF = FL_T0vField_tot, T14vF = FL_T14vField_tot)
FL_list_venn <- Venn(FL_list)


#get the differences
NS_T0vT14_actual<- discern(NS_list_venn, c("T0v14","tot"))
NS_T0vField_actual<- discern(NS_list_venn, c("T0vF","tot"))
NS_T14vField_actual<- discern(NS_list_venn, c("T14vF","tot"))

ME_T0vT14_actual<- discern(ME_list_venn, c("T0v14","tot"))
ME_T0vField_actual<- discern(ME_list_venn, c("T0vF","tot"))
ME_T14vField_actual<- discern(ME_list_venn, c("T14vF","tot"))

NH_T0vT14_actual<- discern(NH_list_venn, c("T0v14","tot"))
NH_T0vField_actual<- discern(NH_list_venn, c("T0vF","tot"))
NH_T14vField_actual<- discern(NH_list_venn, c("T14vF","tot"))

MA_T0vT14_actual<- discern(MA_list_venn, c("T0v14","tot"))
MA_T0vField_actual<- discern(MA_list_venn, c("T0vF","tot"))
MA_T14vField_actual<- discern(MA_list_venn, c("T14vF","tot"))

SC_T0vT14_actual<- discern(SC_list_venn, c("T0v14","tot"))
SC_T0vField_actual<- discern(SC_list_venn, c("T0vF","tot"))
SC_T14vField_actual<- discern(SC_list_venn, c("T14vF","tot"))

FL_T0vT14_actual<- discern(FL_list_venn, c("T0v14","tot"))
FL_T0vField_actual<- discern(FL_list_venn, c("T0vF","tot"))
FL_T14vField_actual<- discern(FL_list_venn, c("T14vF","tot"))


#make presence/absence heatmap 
setmap(Total_venn, title= "Presence/Absence of Genera across all groups")
setmap(All_T0_venn, title= "Presence/Absence of Genera across all T0 groups")
setmap(All_14_venn, title= "Presence/Absence of Genera across all T14 groups")
setmap(NS_venn, title = "Presence/Absence of Genera in Field x NS_T0 x NS_T14")
setmap(ME_venn, title = "Presence/Absence of Genera in Field x ME_T0 x ME_T14")
setmap(NH_venn, title = "Presence/Absence of Genera in Field x NH_T0 x NH_T14")
setmap(MA_venn, title = "Presence/Absence of Genera in Field x MA_T0 x MA_T14")
setmap(SC_venn, title = "Presence/Absence of Genera in Field x SC_T0 x SC_T14")
setmap(FL_venn, title = "Presence/Absence of Genera in Field x FL_T0 x FL_T14")


#north  #54 shared across all (163 total)
north_14s_list = list(NS_Field_T14 = NS_T14vField_actual, ME_Field_T14 = ME_T14vField_actual, NH_Field_T14 = NH_T14vField_actual)
north_venn<- Venn(north_14s_list)
north_overlap<-overlap(north_venn)
ggvenn(north_venn)
ggVennDiagram(north_14s_list,category.names=c("NS T14-Field", "ME T14-Field", "NH T14-Field")) + scale_fill_gradient(low= "#3D405B",high = "#E07A5F") + ggtitle("Comparison of Genera present in Field vs North 14s")
setmap(north_venn, title= "Presence/Absence of Genera Field vs North T14s")


#south #71 shared across all (119 total)
south_14s_list = list(MA_Field_T14 = MA_T14vField_actual, SC_Field_T14 = SC_T14vField_actual, FL_Field_T14 = FL_T14vField_actual)
south_venn<- Venn(south_14s_list)
south_overlap<-overlap(south_venn)
ggvenn(south_venn)
ggVennDiagram(south_14s_list,category.names=c("MA T14-Field", "SC T14-Field", "FL T14-Field")) + scale_fill_gradient(low= "#3D405B",high = "#E07A5F") + ggtitle("Comparison of Genera present in Field vs South 14s")
setmap(south_venn, title= "Presence/Absence of Genera Field vs South T14s")

#make a list of field and overlaps with T14
Field_T14_list = list(NS_Field_T14 = NS_T14vField_actual, ME_Field_T14 = ME_T14vField_actual, NH_Field_T14 = NH_T14vField_actual, MA_Field_T14 = MA_T14vField_actual, SC_Field_T14 = SC_T14vField_actual, FL_Field_T14 = FL_T14vField_actual)
Field_T14_venn <- Venn(Field_T14_list)
Field_T14_overlap<-overlap(Field_T14_venn) #49 shared across all T14s
setmap(Field_T14_venn, title= "Presence/Absence of Genera Field vs T14s")



#### Write out files for members of interactions ####

# 1) Write comparison output to file - not working figure out later 
write_list_NS = list(Total_overlap = tot_overlap_NS, T0vT14_overlap = NS_T0vT14_actual, T0vField_overlap = NS_T0vField_actual, T14vField_overlap = NS_T14vField_actual)
write_list_ME = list(Total_overlap = tot_overlap_ME, T0vT14_overlap = ME_T0vT14_actual, T0vField_overlap = ME_T0vField_actual, T14vField_overlap = ME_T14vField_actual)
write_list_NH = list(Total_overlap = tot_overlap_NH, T0vT14_overlap = NH_T0vT14_actual, T0vField_overlap = NH_T0vField_actual, T14vField_overlap = NH_T14vField_actual)
write_list_MA = list(Total_overlap = tot_overlap_MA, T0vT14_overlap = MA_T0vT14_actual, T0vField_overlap = MA_T0vField_actual, T14vField_overlap = MA_T14vField_actual)
write_list_SC = list(Total_overlap = tot_overlap_SC, T0vT14_overlap = SC_T0vT14_actual, T0vField_overlap = SC_T0vField_actual, T14vField_overlap = SC_T14vField_actual)
write_list_FL = list(Total_overlap = tot_overlap_FL, T0vT14_overlap = FL_T0vT14_actual, T0vField_overlap = FL_T0vField_actual, T14vField_overlap = FL_T14vField_actual)

# 2) Add NAs to make into an even list/dataframe - make all list the same as the longest by adding NA 
#NS
length(write_list_NS$Total_overlap) #28
length(write_list_NS$T0vT14_overlap) #16
length(write_list_NS$T0vField_overlap) #0
length(write_list_NS$T14vField_overlap) #91

# list to adjust <- longest list
length(write_list_NS$Total_overlap) <- length(write_list_NS$T14vField_overlap)
length(write_list_NS$T0vT14_overlap) <- length(write_list_NS$T14vField_overlap)
length(write_list_NS$T0vField_overlap) <- length(write_list_NS$T14vField_overlap)

#ME
length(write_list_ME$Total_overlap) #37
length(write_list_ME$T0vT14_overlap) #24
length(write_list_ME$T0vField_overlap) #2
length(write_list_ME$T14vField_overlap) #77

# list to adjust <- longest list
length(write_list_ME$Total_overlap) <- length(write_list_ME$T14vField_overlap)
length(write_list_ME$T0vT14_overlap) <- length(write_list_ME$T14vField_overlap)
length(write_list_ME$T0vField_overlap) <- length(write_list_ME$T14vField_overlap)

#NH
length(write_list_NH$Total_overlap) #46
length(write_list_NH$T0vT14_overlap) #52
length(write_list_NH$T0vField_overlap) #15
length(write_list_NH$T14vField_overlap) #73

# list to adjust <- longest list
length(write_list_NH$Total_overlap) <- length(write_list_NH$T14vField_overlap)
length(write_list_NH$T0vT14_overlap) <- length(write_list_NH$T14vField_overlap)
length(write_list_NH$T0vField_overlap) <- length(write_list_NH$T14vField_overlap)

#MA
length(write_list_MA$Total_overlap) #28
length(write_list_MA$T0vT14_overlap) #26
length(write_list_MA$T0vField_overlap) #11
length(write_list_MA$T14vField_overlap) #94

# list to adjust <- longest list
length(write_list_MA$Total_overlap) <- length(write_list_MA$T14vField_overlap)
length(write_list_MA$T0vT14_overlap) <- length(write_list_MA$T14vField_overlap)
length(write_list_MA$T0vField_overlap) <- length(write_list_MA$T14vField_overlap)

#SC
length(write_list_SC$Total_overlap) #31
length(write_list_SC$T0vT14_overlap) #25
length(write_list_SC$T0vField_overlap) #11
length(write_list_SC$T14vField_overlap) #93

# list to adjust <- longest list
length(write_list_SC$Total_overlap) <- length(write_list_SC$T14vField_overlap)
length(write_list_SC$T0vT14_overlap) <- length(write_list_SC$T14vField_overlap)
length(write_list_SC$T0vField_overlap) <- length(write_list_SC$T14vField_overlap)

#FL
length(write_list_FL$Total_overlap) #31
length(write_list_FL$T0vT14_overlap) #37
length(write_list_FL$T0vField_overlap) #4
length(write_list_FL$T14vField_overlap) #89

# list to adjust <- longest list
length(write_list_FL$Total_overlap) <- length(write_list_FL$T14vField_overlap)
length(write_list_FL$T0vT14_overlap) <- length(write_list_FL$T14vField_overlap)
length(write_list_FL$T0vField_overlap) <- length(write_list_FL$T14vField_overlap)


# 3) Write out file
write.csv(write_list_NS, file= "NS_total_overlap_members.csv")
write.csv(write_list_ME, file= "ME_total_overlap_members.csv")
write.csv(write_list_NH, file= "NH_total_overlap_members.csv")
write.csv(write_list_MA, file= "MA_total_overlap_members.csv")
write.csv(write_list_SC, file= "SC_total_overlap_members.csv")
write.csv(write_list_FL, file= "FL_total_overlap_members.csv")




#### Visualize Genomic Composition of viral taxa ####

library(ggplot2)
library(dplyr)
library(ggsignif)
setwd("~/Library/CloudStorage/GoogleDrive-sbirch1@charlotte.edu/My Drive/Postdoc_work/Mesocosom_exps/Mesocosom_Summer_2022/Transcriptomics/updated_Blast-Taxonomy")

#Read in data
gen_comp_dat <- read.csv("~/Library/CloudStorage/GoogleDrive-sbirch1@uncc.edu/My Drive/Postdoc_work/Mesocosom_exps/Mesocosom_Summer_2022/Transcriptomics/Blast-taxonomy/genomic_content_refseq_dat.csv", stringsAsFactors=TRUE)

#colors
["#780000","#c1121f","#fdf0d5","#003049","#669bbc"]

["#81b29a","#f4f1de","#e07a5f","#3d405b","#f2cc8f"]



#All samples - With Legend (5 colors needed)
gen_comp_dat %>%
  arrange(count) %>%
  mutate(Location=factor(Location, levels=c("FL-T14","FL-T0","SC-T14","SC-T0","MA-T14","MA-T0","FIELD-T14","FIELD-T0","NH-T14","NH-T0","ME-T14","ME-T0","NS-T14","NS-T0"))) %>%
  ggplot( aes(x=count, y=Location, fill = gen_comp)) + geom_bar(stat= "identity") +theme_minimal() + scale_fill_manual(values=c("#81b29a","#95190C","#e07a5f","#3d405b","#f2cc8f" ))


#T0 only
gen_comp_T0_dat <- read.csv("~/Library/CloudStorage/GoogleDrive-sbirch1@charlotte.edu/My Drive/Postdoc_work/Mesocosom_exps/Mesocosom_Summer_2022/Transcriptomics/updated_Blast-Taxonomy/gen_comp_T0_dat.csv")

#With Legend (15 colors needed)
gen_comp_T0_dat %>%
  arrange(Count) %>%
  mutate(Location=factor(Location, levels=c("FL-T14","FL-T0","SC-T14","SC-T0","MA-T14","MA-T0","FIELD-T14","FIELD-T0","NH-T14","NH-T0","ME-T14","ME-T0","NS-T14","NS-T0"))) %>%
  ggplot( aes(x=Count, y=Location, fill = Gen_comp)) + geom_bar(stat= "identity") +theme_minimal() + scale_fill_manual(values=c("#81b29a","#95190C","#e07a5f","#f2cc8f","#3d405b","#B497D6" )) +
  xlim(NA,1200)

#With Legend (15 colors needed) - filled rel abund
gen_comp_T0_dat %>%
  arrange(Count) %>%
  mutate(Location=factor(Location, levels=c("FL-T14","FL-T0","SC-T14","SC-T0","MA-T14","MA-T0","FIELD-T14","FIELD-T0","NH-T14","NH-T0","ME-T14","ME-T0","NS-T14","NS-T0"))) %>%
  ggplot( aes(x=Count, y=Location, fill = Gen_comp)) + geom_bar(stat= "identity",position = "fill") +theme_minimal() + scale_fill_manual(values=c("#81b29a","#95190C","#e07a5f","#f2cc8f","#3d405b","#B497D6" )) 




#T14 only
gen_comp_T14_dat <- read.csv("~/Library/CloudStorage/GoogleDrive-sbirch1@charlotte.edu/My Drive/Postdoc_work/Mesocosom_exps/Mesocosom_Summer_2022/Transcriptomics/updated_Blast-Taxonomy/gen_comp_T14_dat.csv")

#With Legend (15 colors needed)
gen_comp_T14_dat %>%
  arrange(Count) %>%
  mutate(Location=factor(Location, levels=c("FL-T14","FL-T0","SC-T14","SC-T0","MA-T14","MA-T0","FIELD-T14","FIELD-T0","NH-T14","NH-T0","ME-T14","ME-T0","NS-T14","NS-T0"))) %>%
  ggplot( aes(x=Count, y=Location, fill = Gen_comp)) + geom_bar(stat= "identity") +theme_minimal() + scale_fill_manual(values=c("#81b29a","#95190C","#e07a5f","#f2cc8f","#4A8FE7","#3d405b","#B497D6" )) +
  xlim(NA,1200)


#With Legend (15 colors needed) - filled rel abund
gen_comp_T14_dat %>%
  arrange(Count) %>%
  mutate(Location=factor(Location, levels=c("FL-T14","FL-T0","SC-T14","SC-T0","MA-T14","MA-T0","FIELD-T14","FIELD-T0","NH-T14","NH-T0","ME-T14","ME-T0","NS-T14","NS-T0"))) %>%
  ggplot( aes(x=Count, y=Location, fill = Gen_comp)) + geom_bar(stat= "identity",position = "fill") +theme_minimal() + scale_fill_manual(values=c("#81b29a","#95190C","#e07a5f","#f2cc8f","#4A8FE7","#3d405b","#B497D6" )) 


#####################




##### Diversity Stats: FIRST MAKE A MATRIX OF YOUR DATA - ALL GENERA IN FIRST COL THEN EACH SAMPLE COUNT IN REST OF COLUMNS #### 

## Load in genes data
genus_dat <- read.csv("~/Library/CloudStorage/GoogleDrive-sbirch1@charlotte.edu/My Drive/Postdoc_work/Mesocosom_exps/Mesocosom_Summer_2022/Transcriptomics/updated_Blast-Taxonomy/genus_dat_v2.csv", stringsAsFactors=TRUE)

library(dplyr)

#using dplyr get counts/number of observations for each genus in particular location/TP

#Make DF for each location/timepoint
ns_T0_counts<-genus_dat%>% count(NS.T0)
ns_T14_counts<-genus_dat%>% count(NS.T14)

me_T0_counts<-genus_dat%>% count(ME.T0)
me_T14_counts<-genus_dat%>% count(ME.T14)

nh_T0_counts<-genus_dat%>% count(NH.T0)
nh_T14_counts<-genus_dat%>% count(NH.T14)

field_T0_counts<-genus_dat%>% count(FIELD.T0)
field_T14_counts<-genus_dat%>% count(FIELD.T14)

ma_T0_counts<-genus_dat%>% count(MA.T0)
ma_T14_counts<-genus_dat%>% count(MA.T14)

sc_T0_counts<-genus_dat%>% count(SC.T0)
sc_T14_counts<-genus_dat%>% count(SC.T14)

fl_T0_counts<-genus_dat%>% count(FL.T0)
fl_T14_counts<-genus_dat%>% count(FL.T14)


#remove first row in DF (has total - don't want there)
ns_T0_counts<-ns_T0_counts[-1,]
ns_T14_counts<-ns_T14_counts[-1,]

me_T0_counts<-me_T0_counts[-1,]
me_T14_counts<-me_T14_counts[-1,]

nh_T0_counts<-nh_T0_counts[-1,]
nh_T14_counts<-nh_T14_counts[-1,]

field_T0_counts<-field_T0_counts[-1,]
field_T14_counts<-field_T14_counts[-1,]

ma_T0_counts<-ma_T0_counts[-1,]
ma_T14_counts<-ma_T14_counts[-1,]

sc_T0_counts<-sc_T0_counts[-1,]
sc_T14_counts<-sc_T14_counts[-1,]

fl_T0_counts<-fl_T0_counts[-1,]
fl_T14_counts<-fl_T14_counts[-1,]


#Write out csv file of the dataframe for each locaiton (optional) 
write.csv(ns_T0_counts, file= "NS_T0_genera_counts.csv")
write.csv(ns_T14_counts, file= "NS_T14_genera_counts.csv")

write.csv(me_T0_counts, file= "ME_T0_genera_counts.csv")
write.csv(me_T14_counts, file= "ME_T14_genera_counts.csv")

write.csv(nh_T0_counts, file= "NH_T0_genera_counts.csv")
write.csv(nh_T14_counts, file= "NH_T14_genera_counts.csv")

write.csv(field_T0_counts, file= "FIELD_T0_genera_counts.csv")
write.csv(field_T14_counts, file= "FIELD_T14_genera_counts.csv")

write.csv(ma_T0_counts, file= "MA_T0_genera_counts.csv")
write.csv(ma_T14_counts, file= "MA_T14_genera_counts.csv")

write.csv(sc_T0_counts, file= "SC_T0_genera_counts.csv")
write.csv(sc_T14_counts, file= "SC_T14_genera_counts.csv")

write.csv(fl_T0_counts, file= "FL_T0_genera_counts.csv")
write.csv(fl_T14_counts, file= "FL_T14_genera_counts.csv")

#Make one large matrix with ID first column then each column is the sample with its count

#rename count column
names(ns_T0_counts)[names(ns_T0_counts) == "n"]<-"NS.T0_count"
names(ns_T14_counts)[names(ns_T14_counts) == "n"]<-"NS.T14_count"

names(me_T0_counts)[names(me_T0_counts) == "n"]<-"ME.T0_count"
names(me_T14_counts)[names(me_T14_counts) == "n"]<-"ME.T14_count"

names(nh_T0_counts)[names(nh_T0_counts) == "n"]<-"NH.T0_count"
names(nh_T14_counts)[names(nh_T14_counts) == "n"]<-"NH.T14_count"

names(ma_T0_counts)[names(ma_T0_counts) == "n"]<-"MA.T0_count"
names(ma_T14_counts)[names(ma_T14_counts) == "n"]<-"MA.T14_count"

names(sc_T0_counts)[names(sc_T0_counts) == "n"]<-"SC.T0_count"
names(sc_T14_counts)[names(sc_T14_counts) == "n"]<-"SC.T14_count"

names(fl_T0_counts)[names(fl_T0_counts) == "n"]<-"FL.T0_count"
names(fl_T14_counts)[names(fl_T14_counts) == "n"]<-"FL.T14_count"

names(field_T0_counts)[names(field_T0_counts) == "n"]<-"FIELD.T0_count"
names(field_T14_counts)[names(field_T14_counts) == "n"]<-"FIELD.T14_count"

#rename ID (genera) column 
names(ns_T0_counts)[names(ns_T0_counts) == "NS.T0"]<-"ID"
names(ns_T14_counts)[names(ns_T14_counts) == "NS.T14"]<-"ID"

names(me_T0_counts)[names(me_T0_counts) == "ME.T0"]<-"ID"
names(me_T14_counts)[names(me_T14_counts) == "ME.T14"]<-"ID"

names(nh_T0_counts)[names(nh_T0_counts) == "NH.T0"]<-"ID"
names(nh_T14_counts)[names(nh_T14_counts) == "NH.T14"]<-"ID"

names(ma_T0_counts)[names(ma_T0_counts) == "MA.T0"]<-"ID"
names(ma_T14_counts)[names(ma_T14_counts) == "MA.T14"]<-"ID"

names(sc_T0_counts)[names(sc_T0_counts) == "SC.T0"]<-"ID"
names(sc_T14_counts)[names(sc_T14_counts) == "SC.T14"]<-"ID"

names(fl_T0_counts)[names(fl_T0_counts) == "FL.T0"]<-"ID"
names(fl_T14_counts)[names(fl_T14_counts) == "FL.T14"]<-"ID"

names(field_T0_counts)[names(field_T0_counts) == "FIELD.T0"]<-"ID"
names(field_T14_counts)[names(field_T14_counts) == "FIELD.T14"]<-"ID"


#Run a function to merge all dfs at the same time - make list first of df then merge together:
df_list<-list(ns_T0_counts,ns_T14_counts,me_T0_counts,me_T14_counts,nh_T0_counts,nh_T14_counts,ma_T0_counts,ma_T14_counts,sc_T0_counts,sc_T14_counts,fl_T0_counts,fl_T14_counts,field_T0_counts,field_T14_counts)
genera_counts_merge<-Reduce(function(x,y) merge(x,y, by = "ID", all=T), df_list)

#replace NAs with 0
genera_counts_merge <- genera_counts_merge %>% replace(is.na(.), 0)



##### diversity stats ##### 

#Maybe later try Vegan tutorial to get diversity
#https://peat-clark.github.io/BIO381/veganTutorial.html

#load library
library(vegan)

#data from tutorial:
#data(package="vegan")
#data(dune)
#str(dune)


#look at structure of data - transpose
str(genera_counts_merge)

#try changing row names to genera 
rownames<-genera_counts_merge$ID

#change the rownames to your genera and remove first column 
row.names(genera_counts_merge)<-rownames
genera_counts_merge<-genera_counts_merge[,-1]
str(genera_counts_merge)

#try transposing now
genera_counts_trans<-t(genera_counts_merge)
str(genera_counts_trans)

#change from list to dataframe ## Both methods work - after you transpose it turns from df to a list - need to convert back
genera_counts_trans_df<-data.frame(genera_counts_trans)
str(genera_counts_trans_df)

#calc simpson's 1-D index of diversity for each location #closer to 1 = greater diversity
simpson<-diversity(genera_counts_trans_df,index="simpson")


#calc shannons #Typically ranges from 1.5 - 3.4, higher = more diverse 
shannon<-diversity(genera_counts_trans_df)

# lets compare the two
par(mfrow = c(1, 2))  # use par to generate panels with 1 row of 2 graphs
hist(simpson)
hist(shannon)


#Now calc pairwise dissimilaryt (distance) measures between sites based on their species composition 
#Vegdist computes dissimilarity indices. We are using gower and bray-curtis which are good in detecting underlying ecological gradients
#Both indices are used to quantify the compositional dissimilarity between two different sites. They are bounded between 0 and 1, where 0 = same composition, 1 = maximally dissimilar.

par(mfrow = c(1, 2))
bray = vegdist(genera_counts_trans_df, "bray") 
gower = vegdist(genera_counts_trans_df, "gower")
hist(bray, xlim = range(0.0,1.0))
hist(gower, xlim = range(0.0,1.0))


#Rarefaction - assess expected species richness
spAbund <- rowSums(genera_counts_trans_df)  #gives the number of individuals found in each plot
spAbund # view observations per plot

raremin <- min(rowSums(genera_counts_trans_df))  #rarefaction uses the smallest number of observations per sample to extrapolate the expected number if all other samples only had that number of observations
raremin # view smallest # of obs (site 17)

sRare <- rarefy(genera_counts_trans_df, raremin) # now use function rarefy
sRare #gives an "expected"rarefied" number of species (not obs) if only 15 individuals were present

par(mfrow=c(1,1))
rarecurve(genera_counts_trans_df, col = "blue") # produces rarefaction curves # squares are site numbers positioned at observed space. To "rarefy" a larger site, follow the rarefaction curve until the curve corresponds with the lesser site obs. This gives you rarefied species richness



#### Non-metric multidemensional scaling (NMDS,MDS,NMS) ####

#NMDS does not use the absolute abundances of species in communities, but rather their rank orders and as a result is a flexible 
#technique that accepts a variety of types of data. (It’s also where the “non-metric” part of the name comes from.)

#To run the NMDS, we will use the function metaMDS. The function requires only a community-by-species matrix 
#(which we will create randomly).

#set.seed(2) # random no. generator / way to specify seeds, 2=no. of integers?
#community_matrix=matrix(
#  sample(1:100,300,replace=T),nrow=10, # counts up to 100, 300 cells
#  dimnames=list(paste("community",1:10,sep=""),paste("sp",1:30,sep="")))
#head(community_matrix)


#example_NMDS=metaMDS(community_matrix, # Our community-by-species matrix
#                     k=2) # The number of reduced dimensions. Increase if high stress is problem. 

NMDS=metaMDS(genera_counts_trans_df, # Our community-by-species matrix
             k=2)

#"The stress, or the disagreement between 2-D configuration and predicted values from the regression"

#A good rule of thumb: stress > 0.05 provides an excellent representation in reduced dimensions, > 0.1 is great, >0.2 is good/ok, and stress > 0.3 provides a poor representation

plot(NMDS)

ordiplot(NMDS,type="n") #Ordination plot function especially for congested plots
orditorp(NMDS,display="species",col="red",air=0.01) #The function adds text or points to ordination plots
orditorp(NMDS,display="sites",cex=1.25,air=0.01)

treat=c("T0","T14","T0","T14","T0","T14","T0","T14","T0","T14","T0","T14","T0","T14")
ordiplot(NMDS,type="n")
ordihull(NMDS,groups=treat,draw="polygon",col=c("#BB9D9E","#D5B381"),label=F)
orditorp(NMDS,display="species",col="#6596A4",air=0.01)
orditorp(NMDS,display="sites",col=c("#540B0E","#E09F3E","#540B0E","#E09F3E","#540B0E","#E09F3E","#540B0E","#E09F3E","#540B0E","#E09F3E","#540B0E","#E09F3E","#540B0E","#E09F3E"),
         air=0.01,cex=1.25)



#Following a different vegan tutorial: https://grunwaldlab.github.io/metacoder_documentation/workshop--07--diversity_stats.html

library(vegan)
genera_counts_trans_df$alpha <-diversity(genera_counts_trans_df, index="invsimpson")
genera_counts_trans_df$Location <-c("NS","NS"	,"ME",	"ME",	"NH",	"NH",	"FIELD"	,"FIELD",	"MA",	"MA"	,"SC",	"SC",	"FL",	"FL")
genera_counts_trans_df$Location <-as.factor(genera_counts_trans_df$Location)
genera_counts_trans_df$Timepoint <- c("T0","T14","T0","T14","T0","T14","T0","T14","T0","T14","T0","T14","T0","T14")
genera_counts_trans_df$Timepoint <-as.factor(genera_counts_trans_df$Timepoint)
genera_counts_trans_df$loc.time <- c("NS.T0","NS.T14","ME.T0","ME.T14","NH.T0","NH.T14","FIELD.T0","FIELD.T14","MA.T0","MA.T14","SC.T0","SC.T14","FL.T0","FL.T14")
genera_counts_trans_df$loc.time <-as.factor(genera_counts_trans_df$loc.time)
genera_counts_trans_df$region <- c("North","North","North","North","North","North","Field","Field","South","South","South","South","South","South")
genera_counts_trans_df$region <-as.factor(genera_counts_trans_df$region)
genera_counts_trans_df$region.time <- c("North.T0","North.T14","North.T0","North.T14","North.T0","North.T14","Field","Field","South.T0","South.T14","South.T0","South.T14","South.T0","South.T14")
genera_counts_trans_df$region.time <-as.factor(genera_counts_trans_df$region.time)

T0 vs T14 colors: '#83C5BE','#006D77'

pal <- c("#8ecae6","#2093B0","#023047","#ffb703","#fb8500","#d36582","#ac92a6","#a80874","#f991cc","#e5446d","#BA1200","#636940")
pal <- c("#83C5BE","#006D77")

library(ggplot2)
library(dplyr)
library(ggsignif)
library(agricolae)

region_sub<-subset(genera_counts_trans_df, genera_counts_trans_df$region != "Field")
region_time_sub<-subset(genera_counts_trans_df, genera_counts_trans_df$region.time != "Field")

south_sub<-subset(genera_counts_trans_df, genera_counts_trans_df$region.time == "South.T14")
north_sub<-subset(genera_counts_trans_df, genera_counts_trans_df$region.time == "North.T14")

anova_result <- aov(alpha ~ region.time, region_time_sub)
summary(anova_result)
#region.time  3 1792.4   597.5   5.969 0.0194 *
#Residuals    8  800.7   100.1 

tukey_result <- HSD.test(anova_result, "region.time", group = TRUE)
tukey_result

## Timepoint
time_inv_simp<-ggplot(genera_counts_trans_df, aes(x = Timepoint, y = alpha, fill = Timepoint)) + 
  geom_boxplot() + scale_fill_manual(values = c("#83C5BE","#006D77")) +
  theme(legend.position = "none",
        plot.background = element_rect("white"),
        panel.background = element_rect("gray95"),
        axis.line = element_line("gray25"),
        axis.text = element_text(size = 12, color = "gray25"),
        axis.title = element_text(color = "gray25"),
        legend.text = element_text(size = 12)) + 
  labs(x = "Timepoint",
       y = "Alpha",
       title = "Alpha Diversity - Inverse Simpson")

metadata3 %>%
  filter(!is.na(shannon_entropy)) %>%
  arrange(shannon_entropy) %>%
  mutate(Location=factor(Location, levels=c("NS","ME","NH","Field","MA","SC","FL"))) %>%
  ggplot(aes(x=Location, y=shannon_entropy, color=Timepoint)) +
  geom_boxplot() + geom_jitter(shape=16, size= 2, position = position_jitter(0.2))  +
  xlab("Location") +
  ylab("Shannon Diversity") +
  theme_bw() + # try other themes like theme_bw() or theme_classic()
  scale_color_manual(values=c("#669BBC", "#003049"))

anova_result <- aov(alpha ~ Timepoint, genera_counts_trans_df)
summary(anova_result)

## Location
ggplot(genera_counts_trans_df, aes(x = Location, y = alpha, fill = Location)) + 
  geom_boxplot() + scale_fill_manual(values = c("#8ecae6","#2093B0","#023047","#ffb703","#fb8500","#d36582","#ac92a6","#a80874","#f991cc","#e5446d","#BA1200","#636940")) +
  theme(legend.position = "none",
        plot.background = element_rect("white"),
        panel.background = element_rect("gray95"),
        axis.line = element_line("gray25"),
        axis.text = element_text(size = 12, color = "gray25"),
        axis.title = element_text(color = "gray25"),
        legend.text = element_text(size = 12)) + 
  labs(x = "Location",
       y = "Alpha",
       title = "Alpha Diversity - Inverse Simpson")

anova_result <- aov(alpha ~ Location, genera_counts_trans_df)
summary(anova_result)
#Location     6  431.7    72.0   0.165  0.978
#Residuals    7 3055.5   436.5 

genera_counts_trans_df %>%
  arrange(alpha) %>%
  mutate(Location=factor(Location, levels=c("NS","ME","NH","FIELD","MA","SC","FL"))) %>%
  ggplot( aes(x=Location,y=alpha,fill=Location)) +
  geom_boxplot() + scale_fill_manual(values = c("#8ecae6","#2093B0","#023047","#ffb703","#fb8500","#d36582","#ac92a6","#a80874","#f991cc","#e5446d","#BA1200","#636940")) +
  theme(legend.position = "none",
        plot.background = element_rect("white"),
        panel.background = element_rect("gray95"),
        axis.line = element_line("gray25"),
        axis.text = element_text(size = 12, color = "gray25"),
        axis.title = element_text(color = "gray25"),
        legend.text = element_text(size = 12)) + 
  labs(x = "Location",
       y = "Alpha",
       title = "Alpha Diversity - Inverse Simpson")

## Location
ggplot(genera_counts_trans_df, aes(x = loc.time, y = alpha, color = loc.time)) + 
  geom_point(size=3) + scale_color_manual(values = c("#8ecae6","#2093B0","#023047","#ffb703","#fb8500","#d36582","#ac92a6","#a80874","#f991cc","#e5446d","#BA1200","#636940","#023047","#ffb703")) +
  theme(legend.position = "none",
        plot.background = element_rect("white"),
        panel.background = element_rect("gray95"),
        axis.line = element_line("gray25"),
        axis.text = element_text(size = 12, color = "gray25"),
        axis.title = element_text(color = "gray25"),
        axis.text.x=element_text(angle=45,hjust=1,vjust=1),
        legend.text = element_text(size = 12)) + 
  labs(x = "Location * Timepoint",
       y = "Alpha",
       title = "Alpha Diversity - Inverse Simpson")

genera_counts_trans_df %>%
  arrange(alpha) %>%
  mutate(loc.time=factor(loc.time, levels=c("NS.T0","NS.T14","ME.T0","ME.T14","NH.T0","NH.T14","FIELD.T0","FIELD.T14","MA.T0","MA.T14","SC.T0","SC.T14","FL.T0","FL.T14"))) %>%
  ggplot( aes(x=loc.time,y=alpha,color=loc.time)) +
  geom_point(size=5) + scale_color_manual(values = c("#8ecae6","#023047","#8ecae6","#023047","#8ecae6","#023047","#d36582","#a80874","#ffb703","#fb8500","#ffb703","#fb8500","#ffb703","#fb8500")) +
  theme(legend.position = "none",
        plot.background = element_rect("white"),
        panel.background = element_rect("gray95"),
        axis.line = element_line("gray25"),
        axis.text = element_text(size = 12, color = "gray25"),
        axis.title = element_text(color = "gray25"),
        axis.text.x=element_text(angle=45,hjust=1,vjust=1),
        legend.text = element_text(size = 12)) + 
  labs(x = "Location",
       y = "Alpha",
       title = "Alpha Diversity - Inverse Simpson")



## Region
ggplot(genera_counts_trans_df, aes(x = region, y = alpha, fill = region)) + 
  geom_boxplot() + scale_fill_manual(values = c("#8ecae6","#2093B0","#023047","#ffb703","#fb8500","#d36582","#ac92a6","#a80874","#f991cc","#e5446d","#BA1200","#636940")) +
  theme(legend.position = "none",
        plot.background = element_rect("white"),
        panel.background = element_rect("gray95"),
        axis.line = element_line("gray25"),
        axis.text = element_text(size = 12, color = "gray25"),
        axis.title = element_text(color = "gray25"),
        legend.text = element_text(size = 12)) + 
  labs(x = "Region",
       y = "Alpha",
       title = "Alpha Diversity - Inverse Simpson")

anova_result <- aov(alpha ~ region, genera_counts_trans_df)
summary(anova_result)
#region       2    150   75.07   0.247  0.785
#Residuals   11   3337  303.37

genera_counts_trans_df %>%
  arrange(alpha) %>%
  mutate(Location=factor(Location, levels=c("NS","ME","NH","FIELD","MA","SC","FL"))) %>%
  ggplot( aes(x=Location,y=alpha,fill=Location)) +
  geom_boxplot() + scale_fill_manual(values = c("#8ecae6","#2093B0","#023047","#ffb703","#fb8500","#d36582","#ac92a6","#a80874","#f991cc","#e5446d","#BA1200","#636940")) +
  theme(legend.position = "none",
        plot.background = element_rect("white"),
        panel.background = element_rect("gray95"),
        axis.line = element_line("gray25"),
        axis.text = element_text(size = 12, color = "gray25"),
        axis.title = element_text(color = "gray25"),
        legend.text = element_text(size = 12)) + 
  labs(x = "Location",
       y = "Alpha",
       title = "Alpha Diversity - Inverse Simpson")
#"North" = '#6B0504', "Field" = "#D64045","South" = "#F19C79"
## Region Time
ggplot(genera_counts_trans_df, aes(x = region.time, y = alpha, fill = region.time)) + 
  geom_boxplot() + scale_fill_manual(values = c('#6B0504','#6B0504',"#D64045","#F19C79","#F19C79")) +
  theme(legend.position = "none",
        plot.background = element_rect("white"),
        panel.background = element_rect("gray95"),
        axis.line = element_line("gray25"),
        axis.text = element_text(size = 12, color = "gray25"),
        axis.title = element_text(color = "gray25"),
        legend.text = element_text(size = 12)) + 
  labs(x = "Region X Timepoint",
       y = "Alpha",
       title = "Alpha Diversity - Inverse Simpson")

anova_result <- aov(alpha ~ region.time, genera_counts_trans_df)
summary(anova_result)
#region.time  4   1871   467.7   2.604  0.107
#Residuals    9   1616   179.6 

n_s_inv_simp_a<-genera_counts_trans_df %>%
  arrange(alpha) %>%
  mutate(region.time=factor(region.time, levels=c("North.T0","North.T14","Field","South.T0","South.T14"))) %>%
  ggplot( aes(x=region.time,y=alpha,fill=region.time)) +
  geom_boxplot() + scale_fill_manual(values = c('#6B0504','#6B0504',"#D64045","#F19C79","#F19C79")) +
  theme(legend.position = "none",
        plot.background = element_rect("white"),
        panel.background = element_rect("gray95"),
        axis.line = element_line("gray25"),
        axis.text = element_text(size = 12, color = "gray25"),
        axis.title = element_text(color = "gray25"),
        legend.text = element_text(size = 12)) + 
  labs(x = "Location",
       y = "Alpha",
       title = "Alpha Diversity - Inverse Simpson")

n_s_inv_simp_b<-genera_counts_trans_df %>%
  arrange(alpha) %>%
  mutate(region.time=factor(region.time, levels=c("North.T0","North.T14","Field","South.T0","South.T14"))) %>%
  ggplot( aes(x=region.time,y=alpha,fill=region.time)) +
  geom_boxplot() + scale_fill_manual(values = c("#83C5BE","#006D77","#D64045","#83C5BE","#006D77")) +
  theme(legend.position = "none",
        plot.background = element_rect("white"),
        panel.background = element_rect("gray95"),
        axis.line = element_line("gray25"),
        axis.text = element_text(size = 12, color = "gray25"),
        axis.title = element_text(color = "gray25"),
        legend.text = element_text(size = 12)) + 
  labs(x = "Location",
       y = "Alpha",
       title = "Alpha Diversity - Inverse Simpson")

#convert to a matrix
genera_counts_trans_matrix<-data.matrix(genera_counts_trans_df)



## Beta
beta_dist <- vegdist(genera_counts_trans_matrix,index = "bray")
mds <- metaMDS(beta_dist)
mds_data <- as.data.frame(mds$points)
mds_data$SampleID <- rownames(mds_data)
genera_counts_trans_df_3<-genera_counts_trans_df
genera_counts_trans_df_3$SampleID <-rownames(genera_counts_trans_df_3)
mds_data <- dplyr::left_join(mds_data, genera_counts_trans_df_3)

genera_counts_trans_df_3$Timepoint <- c("T0","T14","T0","T14","T0","T14","T0","T14","T0","T14","T0","T14","T0","T14")
mds_data$Location <-c("NS","NS"	,"ME",	"ME",	"NH",	"NH",	"FIELD"	,"FIELD",	"MA",	"MA"	,"SC",	"SC",	"FL",	"FL")
mds_data$Timepoint <-c("T0","T14","T0","T14","T0","T14","T0","T14","T0","T14","T0","T14","T0","T14")
mds_data$region <-c("North","North","North","North","North","North","Field","Field","South","South","South","South","South","South")
mds_data$region.time <-c("North.T0","North.T14","North.T0","North.T14","North.T0","North.T14","Field","Field","South.T0","South.T14","South.T0","South.T14","South.T0","South.T14")


library(ggplot2)
ggplot(mds_data, aes(x = MDS1, y = MDS2, color = Timepoint)) +
  geom_point(size = 3, aes(shape=Location)) +
  scale_shape_manual(values=c(19,1,17,15,4,8,3))

ggplot(mds_data, aes(x = MDS1, y = MDS2, color = Timepoint)) +
  geom_point(size = 4, aes(shape=Location)) +
  scale_shape_manual(values=c(19,1,17,15,4,8,3))  +stat_ellipse(lwd=1.2) + scale_color_manual(values= c("#8ecae6","#fb8500","#2093B0","#023047","#ffb703","#fb8500","#d36582","#ac92a6","#a80874","#f991cc"))

ggplot(mds_data, aes(x = MDS1, y = MDS2, color = region)) +
  geom_point(size = 4, aes(shape=region)) +
  scale_shape_manual(values=c(19,1,17,15,4,8,3))  +stat_ellipse(lwd=1.2) + scale_color_manual(values= c("#8ecae6","#fb8500","#d36582","#2093B0","#023047","#ffb703","#fb8500","#d36582","#ac92a6","#a80874","#f991cc"))

ggplot(mds_data, aes(x = MDS1, y = MDS2, color = Timepoint)) +
  geom_point(size = 4, aes(shape=region.time)) +
  scale_shape_manual(values=c(19,1,17,15,4,8,3))  +stat_ellipse(lwd=1.2) + scale_color_manual(values= c("#8ecae6","#fb8500","#d36582","#2093B0","#023047","#ffb703","#fb8500","#d36582","#ac92a6","#a80874","#f991cc"))

ggplot(mds_data, aes(x = MDS1, y = MDS2, color = Timepoint)) +
  geom_point(size = 5, aes(shape=region)) +
  scale_shape_manual(values=c(19,1,17,15,4,8,3))  +stat_ellipse(lwd=1.2) + scale_color_manual(values= c("#8ecae6","#fb8500","#d36582","#2093B0","#023047","#ffb703","#fb8500","#d36582","#ac92a6","#a80874","#f991cc"))



#### GOOD INFO #### --> https://www.davidzeleny.net/anadat-r/doku.php/en:rarefaction_r
#https://www.randigriffin.com/2017/05/23/mosquito-community-ecology-in-vegan.html
#https://grunwaldlab.github.io/metacoder_documentation/workshop--07--diversity_stats.html
#https://peat-clark.github.io/BIO381/veganTutorial.html



###### Follwoing a different vegan tutorial: #####
#https://www.rpubs.com/an-bui/vegan-cheat-sheet

# libraries
library(tidyverse)
library(vegan)
library(ggvegan)
#install.packages("remotes")
#remotes::install_github("gavinsimpson/ggvegan")

#data: 
#use genera_counts_trans_df made above - df of your genera counts
#write.csv(genera_counts_trans_df, file= "genera_counts_trans_df.csv")
genera_df <- read_csv(here::here( "genera_counts_trans_df.csv")) 

genera_counts_trans_df_subset<-subset(genera_counts_trans_df,select=-c(Location, Timepoint,region,region.time,loc.time,alpha))
genera_df <- genera_counts_trans_df_subset

# environmental variables
#env <- read_csv(here::here("env_meta_data.csv"))
env <- read.csv("env_meta_data.csv")

# set up a "metadata" frame - will be useful for plotting later!
site_type <- env %>% 
  select(Site, Timepoint,Sample_Name,Region,Location)

head(genera_df,3)
str(genera_df)

#How many species are in each sample
species_num<-specnumber(genera_df)
str(species_num)

# analysis of variance takes the same form as the usual models you'd see in R
# response ~ dependent, data = environmental grouping
sppr_aov <- aov(species_num ~ Timepoint, data = site_type)
summary(sppr_aov)

#            Df Sum Sq Mean Sq F value  Pr(>F)   
#Timepoint    1  66930   66930   16.57 0.00155 **
#Residuals   12  48461    4038                   
#---
#Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1


sppr_aov_2 <- aov(species_num ~ Region*Timepoint, data = site_type)
summary(sppr_aov_2)

species_num<-specnumber(genera_df)
str(species_num)


#This little bit works
species_num<-specnumber(genera_df)
str(species_num)
pal <- c("lightsalmon1", "palegreen4")

boxplot(species_num ~ Timepoint, data = site_type, col=pal,ylab="Richness")


## Timepoint
time_rich<-ggplot(genera_counts_trans_df, aes(x = Timepoint, y = species_num, fill = Timepoint)) + 
  geom_boxplot() + scale_fill_manual(values = c("#83C5BE","#006D77")) +
  theme(legend.position = "none",
        plot.background = element_rect("white"),
        panel.background = element_rect("gray95"),
        axis.line = element_line("gray25"),
        axis.text = element_text(size = 12, color = "gray25"),
        axis.title = element_text(color = "gray25"),
        legend.text = element_text(size = 12)) + 
  labs(x = "Timepoint",
       y = "Richness",
       title = "Richness")

# Region x time
n_s_rich_a<-genera_counts_trans_df %>%
  arrange(species_num) %>%
  mutate(region.time=factor(region.time, levels=c("North.T0","North.T14","Field","South.T0","South.T14"))) %>%
  ggplot( aes(x=region.time,y=species_num,fill=region.time)) +
  geom_boxplot() + scale_fill_manual(values = c('#6B0504','#6B0504',"#D64045","#F19C79","#F19C79")) +
  theme(legend.position = "none",
        plot.background = element_rect("white"),
        panel.background = element_rect("gray95"),
        axis.line = element_line("gray25"),
        axis.text = element_text(size = 12, color = "gray25"),
        axis.title = element_text(color = "gray25"),
        legend.text = element_text(size = 12)) + 
  labs(x = "Location",
       y = "Richness",
       title = "Alpha Diversity - Richness")

n_s_rich_b<-genera_counts_trans_df %>%
  arrange(species_num) %>%
  mutate(region.time=factor(region.time, levels=c("North.T0","North.T14","Field","South.T0","South.T14"))) %>%
  ggplot( aes(x=region.time,y=species_num,fill=region.time)) +
  geom_boxplot() + scale_fill_manual(values = c("#83C5BE","#006D77","#D64045","#83C5BE","#006D77")) +
  theme(legend.position = "none",
        plot.background = element_rect("white"),
        panel.background = element_rect("gray95"),
        axis.line = element_line("gray25"),
        axis.text = element_text(size = 12, color = "gray25"),
        axis.title = element_text(color = "gray25"),
        legend.text = element_text(size = 12)) + 
  labs(x = "Location",
       y = "Richness",
       title = "Alpha Diversity - Richness")

# location
genera_counts_trans_df %>%
  arrange(species_num) %>%
  mutate(Location=factor(Location, levels=c("NS","ME","NH","FIELD","MA","SC","FL"))) %>%
  ggplot( aes(x=Location,y=species_num,fill=Location)) +
  geom_boxplot() + scale_fill_manual(values = c("#8ecae6","#2093B0","#023047","#ffb703","#fb8500","#d36582","#ac92a6","#a80874","#f991cc","#e5446d","#BA1200","#636940")) +
  theme(legend.position = "none",
        plot.background = element_rect("white"),
        panel.background = element_rect("gray95"),
        axis.line = element_line("gray25"),
        axis.text = element_text(size = 12, color = "gray25"),
        axis.title = element_text(color = "gray25"),
        legend.text = element_text(size = 12)) + 
  labs(x = "Location",
       y = "Richness",
       title = "Alpha Diversity - Richness")



###

boxplot(sRare ~ Timepoint, data = site_type, col=pal,ylab="Rarefied richness")
mod.rarefied<-lm(sRare~Timepoint,data=site_type)
anova(mod.rarefied)
#0.0004529 ***


plot_sppr <- ggplot(sppr_df, aes(x = Timepoint, y = value, fill = Timepoint)) +
  geom_boxplot() +
  scale_fill_manual(values = pal) +
  scale_x_discrete(labels = c("T0", "T14")) +
  theme(legend.position = "none",
        plot.background = element_rect("white"),
        panel.background = element_rect("white"),
        panel.grid = element_line("grey90"),
        axis.line = element_line("gray25"),
        axis.text = element_text(size = 12, color = "gray25"),
        axis.title = element_text(color = "gray25"),
        legend.text = element_text(size = 12)) + 
  labs(x = "Ecological landtype",
       y = "Number of species per site",
       title = "Species richness")
plot_sppr




## Shannobn - how diverse are my communities
shannondiv<-diversity(genera_df)

#anova on shannon diversity
sppdiv_aov <- aov(shannondiv ~ Timepoint, data = site_type)
summary(sppdiv_aov)

#            Df Sum Sq Mean Sq F value   Pr(>F)    
#Timepoint    1  2.936  2.9363   22.45 0.000481 ***
#Residuals   12  1.569  0.1308                     
#---
#Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1


sppdiv_aov <- aov(shannondiv ~ Location, data = site_type)
summary(sppdiv_aov)

sppdiv_aov <- aov(shannondiv ~ Region, data = site_type)
summary(sppdiv_aov)


## Timepoint
time_shan<-ggplot(genera_counts_trans_df, aes(x = Timepoint, y = shannondiv, fill = Timepoint)) + 
  geom_boxplot() + scale_fill_manual(values = c("#83C5BE","#006D77")) +
  theme(legend.position = "none",
        plot.background = element_rect("white"),
        panel.background = element_rect("gray95"),
        axis.line = element_line("gray25"),
        axis.text = element_text(size = 12, color = "gray25"),
        axis.title = element_text(color = "gray25"),
        legend.text = element_text(size = 12)) + 
  labs(x = "Timepoint",
       y = "Shannon",
       title = "Alpha Diversity - Shannon")

# Region x time
n_s_shan_a<-genera_counts_trans_df %>%
  arrange(shannondiv) %>%
  mutate(region.time=factor(region.time, levels=c("North.T0","North.T14","Field","South.T0","South.T14"))) %>%
  ggplot( aes(x=region.time,y=shannondiv,fill=region.time)) +
  geom_boxplot() + scale_fill_manual(values = c('#6B0504','#6B0504',"#D64045","#F19C79","#F19C79")) +
  theme(legend.position = "none",
        plot.background = element_rect("white"),
        panel.background = element_rect("gray95"),
        axis.line = element_line("gray25"),
        axis.text = element_text(size = 12, color = "gray25"),
        axis.title = element_text(color = "gray25"),
        legend.text = element_text(size = 12)) + 
  labs(x = "Region x Timepoint",
       y = "Shannon",
       title = "Alpha Diversity - Shannon")

n_s_shan_b<-genera_counts_trans_df %>%
  arrange(shannondiv) %>%
  mutate(region.time=factor(region.time, levels=c("North.T0","North.T14","Field","South.T0","South.T14"))) %>%
  ggplot( aes(x=region.time,y=shannondiv,fill=region.time)) +
  geom_boxplot() + scale_fill_manual(values = c("#83C5BE","#006D77","#D64045","#83C5BE","#006D77")) +
  theme(legend.position = "none",
        plot.background = element_rect("white"),
        panel.background = element_rect("gray95"),
        axis.line = element_line("gray25"),
        axis.text = element_text(size = 12, color = "gray25"),
        axis.title = element_text(color = "gray25"),
        legend.text = element_text(size = 12),axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1)) + 
  labs(x = "Region x Timepoint",
       y = "Shannon",
       title = "Alpha Diversity - Shannon")

# location
genera_counts_trans_df %>%
  arrange(shannondiv) %>%
  mutate(Location=factor(Location, levels=c("NS","ME","NH","FIELD","MA","SC","FL"))) %>%
  ggplot( aes(x=Location,y=shannondiv,fill=Location)) +
  geom_boxplot() + scale_fill_manual(values = c("#8ecae6","#2093B0","#023047","#ffb703","#fb8500","#d36582","#ac92a6","#a80874","#f991cc","#e5446d","#BA1200","#636940")) +
  theme(legend.position = "none",
        plot.background = element_rect("white"),
        panel.background = element_rect("gray95"),
        axis.line = element_line("gray25"),
        axis.text = element_text(size = 12, color = "gray25"),
        axis.title = element_text(color = "gray25"),
        legend.text = element_text(size = 12)) + 
  labs(x = "Location",
       y = "Shannon",
       title = "Alpha Diversity - Shannon")





shandiv_df <- shannondiv %>% 
  # put all those calculations into a data frame
  enframe() %>% 
  # rename columns for ease of joining
  rename(Site = name,
         shan_div = value)

div_plot_df <- shandiv_df %>% 
  # join with site_type
  full_join(site_type, ., by = "Site") %>% 
  # group by landtype
  group_by(Timepoint) %>% 
  # calculate mean and standard error of diversity
  summarize(mean = round(mean(shan_div), 2),
            err = sd(shan_div)/sqrt(length(shan_div))) %>% 
  dplyr::mutate(label = "mean") %>% 
  unite("mean_label", label, mean, sep = " = ", remove = FALSE)

clean_background <- theme(plot.background = element_rect("white"),
                          panel.background = element_rect("gray95"),
                          panel.grid = element_line("white"),
                          axis.line = element_line("gray25"),
                          axis.text = element_text(size = 12, color = "gray25"),
                          axis.title = element_text(color = "gray25"),
                          legend.text = element_text(size = 12),
                          legend.key = element_rect("white"))

plot_shandiv <- ggplot(div_plot_df, aes(x = Timepoint, y = mean, fill = Timepoint)) +
  geom_col(color = "black") +
  scale_fill_manual(values = pal) +
  geom_errorbar(aes(ymin = mean - err, ymax = mean + err), width = 0.5) +
  geom_text(aes(x = Timepoint, y = mean + err + 0.07, label = mean_label)) +
  scale_x_discrete(labels = c("T0","T14")) +
  clean_background + 
  theme(legend.position = "none") +
  labs(x = "Timepoint",
       y = "Mean Shannon diversity",
       title = "Shannon diversity")
plot_shandiv


# How different are my communities in species composition (PerMANOVA)

genera_perm <- adonis2(genera_counts_trans_df_subset ~ Timepoint, data = env)
genera_perm

#Permutation test for adonis under reduced model
#Permutation: free
#Number of permutations: 999

#adonis2(formula = genera_counts_trans_df ~ Timepoint, data = env)
#         Df SumOfSqs      R2      F Pr(>F)   
#Model     1  0.52536 0.24278 3.8475  0.003 **
#Residual 12  1.63856 0.75722                 
#Total    13  2.16392 1.00000                 
#---
#Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

genera_perm <- adonis2(genera_counts_trans_df_subset ~ Location, data = env)
genera_perm


#PCA
gen_PCA<-rda(genera_counts_trans_df)
gen_PCA  

PCAscores <- scores(gen_PCA, display = "sites") %>% 
  as.data.frame() %>% 
  rownames_to_column("Site") %>% 
  full_join(site_type, by = "Site")

PCAvect <- scores(gen_PCA, display = "species") %>% 
  as.data.frame()

plot_PCA <- ggplot() +
  geom_point(data = PCAscores, aes(x = PC1, y = PC2, color = Timepoint)) +
  scale_color_manual(values = pal) +
  geom_vline(xintercept = c(0), color = "grey70", linetype = 2) +
  geom_hline(yintercept = c(0), color = "grey70", linetype = 2) +
  geom_segment(data = PCAvect, aes(x = 0, y = 0, xend = PC1, yend = PC2), arrow = arrow(length = unit(0.2, "cm"))) +
  geom_text(data = PCAvect, aes(x = PC1, y = PC2, label = rownames(PCAvect))) +
  clean_background +
  labs(x = "PC1 (39.35%)",
       y = "PC2 (17.52%)",
       title = "Principal Components Analysis") 
plot_PCA  


plot_PCA <- ggplot() +
  geom_point(data = PCAscores, aes(x = PC1, y = PC2, color = Site)) +
  scale_color_manual(values = c("#DE3C00","#FFA07A","#08548F","aquamarine3","#FFC125","#8B4500","#8B668B","#FFB6C1","#556B2F","#A2CD5A","#CDC1C5","#8B8386","#FF82AB","#8B475D")) +
  geom_vline(xintercept = c(0), color = "grey70", linetype = 2) +
  geom_hline(yintercept = c(0), color = "grey70", linetype = 2) +
  geom_segment(data = PCAvect, aes(x = 0, y = 0, xend = PC1, yend = PC2), arrow = arrow(length = unit(0.2, "cm"))) +
  geom_text(data = PCAvect, aes(x = PC1, y = PC2, label = rownames(PCAvect))) +
  clean_background +
  labs(x = "PC1 (39.35%)",
       y = "PC2 (17.52%)",
       title = "Principal Components Analysis") 
plot_PCA 


PCA_biplot <- autoplot(gen_PCA)
PCA_biplot



#do nmds
pal<-c("#D64045",'#6B0504',"#F19C79")
pal<-c("#83C5BE","#006D77")
pal <- c("#2093B0","#fb8500","#023047","#d36582","#8ecae6","#ffb703","#fb8500","#d36582","#ac92a6","#a80874","#f991cc","#e5446d","#BA1200","#636940")

bird_NMDS <- metaMDS(genera_counts_trans_df_subset)
stressplot(bird_NMDS)
plot(bird_NMDS)

plot_df <- scores(bird_NMDS, display = "sites") %>% 
  as.data.frame() %>% 
  rownames_to_column("Site") %>% 
  full_join(site_type, by = "Site")

plot_nmds_time <- ggplot(plot_df, aes(x = NMDS1, y = NMDS2, color = Timepoint, shape = Timepoint)) +
  geom_point(size = 3, alpha = 0.8) +
  scale_color_manual(values = pal) +
  stat_ellipse(linetype = 2, linewidth = 1) +
  clean_background +
  labs(title = "NMDS")
plot_nmds_time

plot_nmds_reg <- ggplot(plot_df, aes(x = NMDS1, y = NMDS2, color = Region, shape = Timepoint)) +
  geom_point(size = 3, alpha = 0.8) +
  scale_color_manual(values = pal) +
  stat_ellipse(linetype = 2, linewidth = 1) +
  clean_background +
  labs(title = "NMDS")
plot_nmds_reg

plot_nmds_reg <- ggplot(plot_df, aes(x = NMDS1, y = NMDS2, color = Region, shape = Timepoint)) +
  geom_point(size = 3, alpha = 0.8) +
  scale_color_manual(values = pal) +
  stat_ellipse(linetype = 2, linewidth = 1) +
  clean_background +
  labs(title = "NMDS")
plot_nmds_reg

clean_background <- theme(plot.background = element_rect("white"),
                          panel.background = element_rect("gray95"),
                          panel.grid = element_line("white"),
                          axis.line = element_line("gray25"),
                          axis.text = element_text(size = 12, color = "gray25"),
                          axis.title = element_text(color = "gray25"),
                          legend.position = "none")

legend.position = "none"
plot_nmds <- ggplot(plot_df, aes(x = NMDS1, y = NMDS2, color = Timepoint, shape = Region)) +
  geom_point(size = 3, alpha = 0.8) +
  scale_color_manual(values = pal) +
  stat_ellipse(linetype = 2, linewidth = 1) +
  clean_background +
  labs(title = "NMDS")
plot_nmds



#Make cowplots of diversity data for supp

#three panel figure
library(cowplot)

#in a square
plot_grid(p,v,ph)
#in single row
plot_grid(p,v,ph, nrow=1,ncol=3, rel_widths = 3)
plot_grid(time_inv_simp,n_s_inv_simp_b, time_rich, n_s_rich_b, time_shan, n_s_shan_b,plot_nmds_reg, nrow=4,ncol=2, rel_heights = 3)
plot_grid(time_inv_simp,n_s_inv_simp_a, time_rich, n_s_rich_a, time_shan, n_s_shan_a,plot_nmds_reg, nrow=4,ncol=2, rel_heights = 3)


