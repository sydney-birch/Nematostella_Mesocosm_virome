## Mesocosm 2022 - Metacerberus work

# We statistically examined the Metacerberus output (counts across Term databases), visualized the terms, and ran statistics on specific terms


setwd("/Users/Sydney/Library/CloudStorage/GoogleDrive-sbirch1@charlotte.edu/My Drive/Postdoc_work/Mesocosom_exps/Mesocosom_Summer_2022/Transcriptomics/updated_Blast-Taxonomy")

#Import data
stats_db <- read.csv("~/Library/CloudStorage/GoogleDrive-sbirch1@charlotte.edu/My Drive/Postdoc_work/Mesocosom_exps/Mesocosom_Summer_2022/Transcriptomics/updated_Blast-Taxonomy/prot_count_comparison_v2.csv", stringsAsFactors=TRUE)
str(stats_db)

#Ordered Bar plot
library(dplyr)
library(ggsignif)
stats_db %>%
  arrange(Prot_Count) %>%
  mutate(location_tp=factor(location_tp, levels=c("NS_T0","NS_T14","ME_T0","ME_T14","NH_T0","NH_T14","FIELD_T0","FIELD_T14","MA_T0","MA_T14","SC_T0","SC_T14","FL_T0","FL_T14"))) %>%
  ggplot( aes(x=location_tp,y=Prot_Count,fill=progam)) +
  geom_bar(stat="identity", color="black", position=position_dodge())+theme_minimal() +
  scale_fill_manual(values=c('#264653','#CB8589','#2A9D8F')) + ggtitle("Protein Counts Comparison") 


#VOG subset data
stats_vog_dat <- read.csv("~/Library/CloudStorage/GoogleDrive-sbirch1@uncc.edu/My Drive/Postdoc_work/Mesocosom_exps/Mesocosom_Summer_2022/Transcriptomics/metacerberus_1st_run/stats_vog.csv", stringsAsFactors=TRUE)
str(stats_vog_dat)

#VOG ID Counts
library(dplyr)
library(ggsignif)
stats_vog_dat %>%
  arrange(VOG_count) %>%
  mutate(Location=factor(Location, levels=c("NS","ME","NH","FIELD","MA","SC","FL"))) %>%
  ggplot( aes(x=Location,y=VOG_count,fill=time_point)) +
  geom_bar(stat="identity", color="black", position=position_dodge())+theme_minimal() +
  scale_fill_manual(values=c('#BBFFFF','#668B8B')) + ylim(NA,2000) + ggtitle("VOG ID Counts") 


#PHROG subset data
stats_phrog_dat <- read.csv("~/Library/CloudStorage/GoogleDrive-sbirch1@uncc.edu/My Drive/Postdoc_work/Mesocosom_exps/Mesocosom_Summer_2022/Transcriptomics/metacerberus_1st_run/stats_phrog.csv", stringsAsFactors=TRUE)
str(stats_phrog_dat)

stats_phrog_dat %>%
  arrange(PHROG_count) %>%
  mutate(Location=factor(Location, levels=c("NS","ME","NH","FIELD","MA","SC","FL"))) %>%
  ggplot( aes(x=Location,y=PHROG_count,fill=time_point)) +
  geom_bar(stat="identity", color="black", position=position_dodge())+theme_minimal() +
  scale_fill_manual(values=c('#B4EEB4','#698B69')) + ylim(NA,1200)+ ggtitle("PHROG ID Counts")


#Protein subset data
stats_protien_dat <- read.csv("~/Library/CloudStorage/GoogleDrive-sbirch1@uncc.edu/My Drive/Postdoc_work/Mesocosom_exps/Mesocosom_Summer_2022/Transcriptomics/metacerberus_1st_run/stats_protien.csv", stringsAsFactors=TRUE)
str(stats_protien_dat)

stats_protien_dat %>%
  arrange(Protien_count) %>%
  mutate(Location=factor(Location, levels=c("NS","ME","NH","FIELD","MA","SC","FL"))) %>%
  ggplot(aes(x=Location,y=Protien_count,fill=time_point)) +
  geom_bar(stat="identity", color="black", position=position_dodge())+theme_minimal() +
  scale_fill_manual(values=c('#EE6A50','#8B3E2F')) + ggtitle("Protein Counts")


#Kegg subset data
stats_Kegg_dat <- read.csv("~/Library/CloudStorage/GoogleDrive-sbirch1@charlotte.edu/My Drive/Postdoc_work/Mesocosom_exps/Mesocosom_Summer_2022/Transcriptomics/metacerberus_1st_run/stats_kegg.csv", stringsAsFactors=TRUE)
str(stats_Kegg_dat)

stats_Kegg_dat %>%
  arrange(KEGG_count) %>%
  mutate(Location=factor(Location, levels=c("NS","ME","NH","FIELD","MA","SC","FL"))) %>%
  ggplot(aes(x=Location,y=KEGG_count,fill=time_point)) +
  geom_bar(stat="identity", color="black", position=position_dodge())+theme_minimal() +
  scale_fill_manual(values=c('#83C5BE','#006D77')) + ggtitle("KEGG Counts")




#### Make stacked bar plot of VOG Level 1 terms ####

setwd("/Users/Sydney/Library/CloudStorage/GoogleDrive-sbirch1@uncc.edu/My Drive/Postdoc_work/Mesocosom_exps/Mesocosom_Summer_2022/Transcriptomics/metacerberus_1st_run")

#Read in data
vog_barplot_dat <- read.csv("~/Library/CloudStorage/GoogleDrive-sbirch1@uncc.edu/My Drive/Postdoc_work/Mesocosom_exps/Mesocosom_Summer_2022/Transcriptomics/metacerberus_1st_run/vog_level_1_barplot_dat.csv", stringsAsFactors=TRUE)
str(vog_barplot_dat)

library(ggplot2)
library(dplyr)
library(ggsignif)


#With Legend (15 colors needed)
vog_barplot_dat %>%
  arrange(Counts) %>%
  mutate(Location=factor(Location, levels=c("FL_T14","FL_T0","SC_T14","SC_T0","MA_T14","MA_T0","FIELD_T14","FIELD_T0","NH_T14","NH_T0","ME_T14","ME_T0","NS_T14","NS_T0"))) %>%
  ggplot( aes(x=Counts, y=Location, fill = Function)) + geom_bar(stat= "identity", color = "black") +theme_minimal() + scale_fill_manual(values=c("#127475","#8cb369","#f4e285","#f4a259","#5b8e7d","#bc4b51","#2d2d2a","#0c4767","#820933","#34623f","#546a76","#562c2c","#f2542d","#f5dfbb","#0e9594" ))

#no legend
vog_barplot_dat %>%
  arrange(Counts) %>%
  mutate(Location=factor(Location, levels=c("FL_T14","FL_T0","SC_T14","SC_T0","MA_T14","MA_T0","FIELD_T14","FIELD_T0","NH_T14","NH_T0","ME_T14","ME_T0","NS_T14","NS_T0"))) %>%
  ggplot( aes(x=Counts, y=Location, fill = Function)) + geom_bar(stat= "identity", show.legend = FALSE) +theme_minimal() + scale_fill_manual(values=c("#127475","#8cb369","#f4e285","#f4a259","#5b8e7d","#bc4b51","#2d2d2a","#0c4767","#820933","#34623f","#546a76","#562c2c","#f2542d","#f5dfbb","#0e9594" ))

#Stacked barplot VOG with Function unknown - T0 and T14 Seperated

#T0 only
vog_dat_T0 <- read.csv("~/Library/CloudStorage/GoogleDrive-sbirch1@uncc.edu/My Drive/Postdoc_work/Mesocosom_exps/Mesocosom_Summer_2022/Transcriptomics/metacerberus_1st_run/vog_level_1_barplot_dat_T0.csv", stringsAsFactors=TRUE)

#With Legend (15 colors needed)
vog_dat_T0 %>%
  arrange(Counts) %>%
  mutate(Location=factor(Location, levels=c("FL_T14","FL_T0","SC_T14","SC_T0","MA_T14","MA_T0","FIELD_T14","FIELD_T0","NH_T14","NH_T0","ME_T14","ME_T0","NS_T14","NS_T0"))) %>%
  ggplot( aes(x=Counts, y=Location, fill = Function)) + geom_bar(stat= "identity", color = "black") +theme_minimal() + scale_fill_manual(values=c("#127475","#8cb369","#f4e285","#f4a259","#5b8e7d","#bc4b51","#2d2d2a","#0c4767","#820933","#34623f","#546a76","#562c2c","#f2542d","#f5dfbb","#0e9594" ))

#no legend
vog_dat_T0 %>%
  arrange(Counts) %>%
  mutate(Location=factor(Location, levels=c("FL_T14","FL_T0","SC_T14","SC_T0","MA_T14","MA_T0","FIELD_T14","FIELD_T0","NH_T14","NH_T0","ME_T14","ME_T0","NS_T14","NS_T0"))) %>%
  ggplot( aes(x=Counts, y=Location, fill = Function)) + geom_bar(stat= "identity", show.legend = FALSE) +theme_minimal() + scale_fill_manual(values=c("#127475","#8cb369","#f4e285","#f4a259","#5b8e7d","#bc4b51","#2d2d2a","#0c4767","#820933","#34623f","#546a76","#562c2c","#f2542d","#f5dfbb","#0e9594" )) + 
  xlim(NA,2000)



#T14 only
vog_dat_T14 <- read.csv("~/Library/CloudStorage/GoogleDrive-sbirch1@uncc.edu/My Drive/Postdoc_work/Mesocosom_exps/Mesocosom_Summer_2022/Transcriptomics/metacerberus_1st_run/vog_level_1_barplot_dat_T14.csv", stringsAsFactors=TRUE)

#With Legend (15 colors needed)
vog_dat_T14 %>%
  arrange(Counts) %>%
  mutate(Location=factor(Location, levels=c("FL_T14","FL_T0","SC_T14","SC_T0","MA_T14","MA_T0","FIELD_T14","FIELD_T0","NH_T14","NH_T0","ME_T14","ME_T0","NS_T14","NS_T0"))) %>%
  ggplot( aes(x=Counts, y=Location, fill = Function)) + geom_bar(stat= "identity", color = "black") +theme_minimal() + scale_fill_manual(values=c("#127475","#8cb369","#f4e285","#f4a259","#5b8e7d","#bc4b51","#2d2d2a","#0c4767","#820933","#34623f","#546a76","#562c2c","#f2542d","#f5dfbb","#0e9594" ))

#no legend
vog_dat_T14 %>%
  arrange(Counts) %>%
  mutate(Location=factor(Location, levels=c("FL_T14","FL_T0","SC_T14","SC_T0","MA_T14","MA_T0","FIELD_T14","FIELD_T0","NH_T14","NH_T0","ME_T14","ME_T0","NS_T14","NS_T0"))) %>%
  ggplot( aes(x=Counts, y=Location, fill = Function)) + geom_bar(stat= "identity", show.legend = FALSE) +theme_minimal() + scale_fill_manual(values=c("#127475","#8cb369","#f4e285","#f4a259","#5b8e7d","#bc4b51","#2d2d2a","#0c4767","#820933","#34623f","#546a76","#562c2c","#f2542d","#f5dfbb","#0e9594" )) +
  xlim(NA,2000)





#Stacked barplot VOG without Function unknown
vog_2_dat <- read.csv("~/Library/CloudStorage/GoogleDrive-sbirch1@uncc.edu/My Drive/Postdoc_work/Mesocosom_exps/Mesocosom_Summer_2022/Transcriptomics/metacerberus_1st_run/vog_level_1_barplot_dat_2.csv", stringsAsFactors=TRUE)
str(vog_2_dat)

#With Legend
vog_2_dat %>%
  arrange(Counts) %>%
  mutate(Location=factor(Location, levels=c("FL_T14","FL_T0","SC_T14","SC_T0","MA_T14","MA_T0","FIELD_T14","FIELD_T0","NH_T14","NH_T0","ME_T14","ME_T0","NS_T14","NS_T0"))) %>%
  ggplot( aes(x=Counts, y=Location, fill = Function)) + geom_bar(stat= "identity", color = "black") +theme_minimal() + scale_fill_manual(values=c("#8cb369","#f4e285","#f4a259","#5b8e7d","#bc4b51","#2d2d2a","#0c4767","#820933","#34623f","#546a76","#562c2c","#f2542d","#f5dfbb","#0e9594" ))

#no legend
vog_2_dat %>%
  arrange(Counts) %>%
  mutate(Location=factor(Location, levels=c("FL_T14","FL_T0","SC_T14","SC_T0","MA_T14","MA_T0","FIELD_T14","FIELD_T0","NH_T14","NH_T0","ME_T14","ME_T0","NS_T14","NS_T0"))) %>%
  ggplot( aes(x=Counts, y=Location, fill = Function)) + geom_bar(stat= "identity", show.legend = FALSE) +theme_minimal() + scale_fill_manual(values=c("#8cb369","#f4e285","#f4a259","#5b8e7d","#bc4b51","#2d2d2a","#0c4767","#820933","#34623f","#546a76","#562c2c","#f2542d","#f5dfbb","#0e9594" ))



# vog - without function unknown - T0 only
vog_2_t0_dat <- read.csv("~/Library/CloudStorage/GoogleDrive-sbirch1@uncc.edu/My Drive/Postdoc_work/Mesocosom_exps/Mesocosom_Summer_2022/Transcriptomics/metacerberus_1st_run/vog_level_1_barplot_dat_2_T0.csv", stringsAsFactors=TRUE)

#with legend
vog_2_t0_dat %>%
  arrange(Counts) %>%
  mutate(Location=factor(Location, levels=c("FL_T14","FL_T0","SC_T14","SC_T0","MA_T14","MA_T0","FIELD_T14","FIELD_T0","NH_T14","NH_T0","ME_T14","ME_T0","NS_T14","NS_T0"))) %>%
  ggplot( aes(x=Counts, y=Location, fill = Function)) + geom_bar(stat= "identity", color = "black") +theme_minimal() + scale_fill_manual(values=c("#8cb369","#f4e285","#f4a259","#5b8e7d","#bc4b51","#2d2d2a","#0c4767","#820933","#34623f","#546a76","#562c2c","#f2542d","#f5dfbb","#0e9594" ))

#no legend
vog_2_t0_dat %>%
  arrange(Counts) %>%
  mutate(Location=factor(Location, levels=c("FL_T14","FL_T0","SC_T14","SC_T0","MA_T14","MA_T0","FIELD_T14","FIELD_T0","NH_T14","NH_T0","ME_T14","ME_T0","NS_T14","NS_T0"))) %>%
  ggplot( aes(x=Counts, y=Location, fill = Function)) + geom_bar(stat= "identity", show.legend = FALSE) +theme_minimal() + scale_fill_manual(values=c("#8cb369","#f4e285","#f4a259","#5b8e7d","#bc4b51","#2d2d2a","#0c4767","#820933","#34623f","#546a76","#562c2c","#f2542d","#f5dfbb","#0e9594" )) +
  xlim(NA,600)


# vog - without function unknown - T14 only
vog_2_t14_dat <- read.csv("~/Library/CloudStorage/GoogleDrive-sbirch1@uncc.edu/My Drive/Postdoc_work/Mesocosom_exps/Mesocosom_Summer_2022/Transcriptomics/metacerberus_1st_run/vog_level_1_barplot_dat_2_T14.csv", stringsAsFactors=TRUE)

#with legend #14 colors
vog_2_t14_dat %>%
  arrange(Counts) %>%
  mutate(Location=factor(Location, levels=c("FL_T14","FL_T0","SC_T14","SC_T0","MA_T14","MA_T0","FIELD_T14","FIELD_T0","NH_T14","NH_T0","ME_T14","ME_T0","NS_T14","NS_T0"))) %>%
  ggplot( aes(x=Counts, y=Location, fill = Function)) + geom_bar(stat= "identity", color = "black") +theme_minimal() + scale_fill_manual(values=c("#8cb369","#f4e285","#f4a259","#5b8e7d","#bc4b51","#2d2d2a","#0c4767","#820933","#34623f","#546a76","#562c2c","#f2542d","#f5dfbb","#0e9594" ))

#no legend #14 colors - love this - took out the black boarder between bars 
vog_2_t14_dat %>%
  arrange(Counts) %>%
  mutate(Location=factor(Location, levels=c("FL_T14","FL_T0","SC_T14","SC_T0","MA_T14","MA_T0","FIELD_T14","FIELD_T0","NH_T14","NH_T0","ME_T14","ME_T0","NS_T14","NS_T0"))) %>%
  ggplot( aes(x=Counts, y=Location, fill = Function)) + geom_bar(stat= "identity", show.legend = FALSE) +theme_minimal() + scale_fill_manual(values=c("#8cb369","#f4e285","#f4a259","#5b8e7d","#bc4b51","#2d2d2a","#0c4767","#820933","#34623f","#546a76","#562c2c","#f2542d","#f5dfbb","#0e9594" )) +
  xlim(NA,600)




#### Make stacked bar plot of PHROG Level 1 terms ####

setwd("/Users/Sydney/Library/CloudStorage/GoogleDrive-sbirch1@uncc.edu/My Drive/Postdoc_work/Mesocosom_exps/Mesocosom_Summer_2022/Transcriptomics/metacerberus_1st_run")

#Read in data (10 colors needed)
phrog_dat <- read.csv("~/Library/CloudStorage/GoogleDrive-sbirch1@uncc.edu/My Drive/Postdoc_work/Mesocosom_exps/Mesocosom_Summer_2022/Transcriptomics/metacerberus_1st_run/PHROG_level_1_dat.csv", stringsAsFactors=TRUE)
str(phrog_dat)

library(ggplot2)
library(dplyr)
library(ggsignif)

#with legend
phrog_dat %>%
  arrange(Counts) %>%
  mutate(Location=factor(Location, levels=c("FL_T14","FL_T0","SC_T14","SC_T0","MA_T14","MA_T0","FIELD_T14","FIELD_T0","NH_T14","NH_T0","ME_T14","ME_T0","NS_T14","NS_T0"))) %>%
  ggplot( aes(x=Counts, y=Location, fill = Function)) + geom_bar(stat= "identity", color = "black") +theme_minimal()

#no legend
phrog_dat %>%
  arrange(Counts) %>%
  mutate(Location=factor(Location, levels=c("FL_T14","FL_T0","SC_T14","SC_T0","MA_T14","MA_T0","FIELD_T14","FIELD_T0","NH_T14","NH_T0","ME_T14","ME_T0","NS_T14","NS_T0"))) %>%
  ggplot( aes(x=Counts, y=Location, fill = Function)) + geom_bar(stat= "identity", show.legend = FALSE) +theme_minimal() + scale_fill_manual(values=c("#c33c54","#34a0a4","#522b47","#6f6866","#85cbd0","#184e77","#1a759f","#a9de8e","#c6878f","#52b600"))



#Phrog dat mirrored - T0
PHROG_dat_T0 <- read.csv("~/Library/CloudStorage/GoogleDrive-sbirch1@uncc.edu/My Drive/Postdoc_work/Mesocosom_exps/Mesocosom_Summer_2022/Transcriptomics/metacerberus_1st_run/PHROG_level_1_dat_T0.csv", stringsAsFactors=TRUE)

#with legend
PHROG_dat_T0 %>%
  arrange(Counts) %>%
  mutate(Location=factor(Location, levels=c("FL_T14","FL_T0","SC_T14","SC_T0","MA_T14","MA_T0","FIELD_T14","FIELD_T0","NH_T14","NH_T0","ME_T14","ME_T0","NS_T14","NS_T0"))) %>%
  ggplot( aes(x=Counts, y=Location, fill = Function)) + geom_bar(stat= "identity", color = "black") +theme_minimal() + scale_fill_manual(values=c("#c33c54","#34a0a4","#522b47","#6f6866","#85cbd0","#184e77","#1a759f","#a9de8e","#c6878f","#52b600"))

#no legend
PHROG_dat_T0 %>%
  arrange(Counts) %>%
  mutate(Location=factor(Location, levels=c("FL_T14","FL_T0","SC_T14","SC_T0","MA_T14","MA_T0","FIELD_T14","FIELD_T0","NH_T14","NH_T0","ME_T14","ME_T0","NS_T14","NS_T0"))) %>%
  ggplot( aes(x=Counts, y=Location, fill = Function)) + geom_bar(stat= "identity", show.legend = FALSE) +theme_minimal() + scale_fill_manual(values=c("#c33c54","#34a0a4","#522b47","#6f6866","#85cbd0","#184e77","#1a759f","#a9de8e","#c6878f","#52b600")) +
  xlim(NA,1200)


#Phrog dat mirrored - T14
PHROG_dat_T14 <- read.csv("~/Library/CloudStorage/GoogleDrive-sbirch1@uncc.edu/My Drive/Postdoc_work/Mesocosom_exps/Mesocosom_Summer_2022/Transcriptomics/metacerberus_1st_run/PHROG_level_1_dat_T14.csv", stringsAsFactors=TRUE)

#with legend
PHROG_dat_T14 %>%
  arrange(Counts) %>%
  mutate(Location=factor(Location, levels=c("FL_T14","FL_T0","SC_T14","SC_T0","MA_T14","MA_T0","FIELD_T14","FIELD_T0","NH_T14","NH_T0","ME_T14","ME_T0","NS_T14","NS_T0"))) %>%
  ggplot( aes(x=Counts, y=Location, fill = Function)) + geom_bar(stat= "identity", color = "black") +theme_minimal() + scale_fill_manual(values=c("#c33c54","#34a0a4","#522b47","#6f6866","#85cbd0","#184e77","#1a759f","#a9de8e","#c6878f","#52b600"))

#no legend
PHROG_dat_T14 %>%
  arrange(Counts) %>%
  mutate(Location=factor(Location, levels=c("FL_T14","FL_T0","SC_T14","SC_T0","MA_T14","MA_T0","FIELD_T14","FIELD_T0","NH_T14","NH_T0","ME_T14","ME_T0","NS_T14","NS_T0"))) %>%
  ggplot( aes(x=Counts, y=Location, fill = Function)) + geom_bar(stat= "identity", show.legend = FALSE) +theme_minimal() + scale_fill_manual(values=c("#c33c54","#34a0a4","#522b47","#6f6866","#85cbd0","#184e77","#1a759f","#a9de8e","#c6878f","#52b600")) +
  xlim(NA,1200)




#### Make stacked bar plot of COG Level 1 terms ####

setwd("/Users/Sydney/Library/CloudStorage/GoogleDrive-sbirch1@charlotte.edu/My Drive/Postdoc_work/Mesocosom_exps/Mesocosom_Summer_2022/Transcriptomics/metacerberus_1st_run")

#Read in data (25 colors needed)
cog_dat <- read.csv("~/Library/CloudStorage/GoogleDrive-sbirch1@uncc.edu/My Drive/Postdoc_work/Mesocosom_exps/Mesocosom_Summer_2022/Transcriptomics/metacerberus_1st_run/COG_level_1_dat.csv", stringsAsFactors=TRUE)
str(cog_dat)

library(ggplot2)
library(dplyr)
library(ggsignif)

#Location x axis
#ggplot(cog_dat, aes(x=Location, y=Counts, fill = Function)) + geom_bar(stat= "identity", color = "black") 

#Counts x axis
#ggplot(cog_dat, aes(x=Counts, y=Location, fill = Function)) + geom_bar(stat= "identity", color = "black")

#with legend
cog_dat %>%
  arrange(Counts) %>%
  mutate(Location=factor(Location, levels=c("FL_T14","FL_T0","SC_T14","SC_T0","MA_T14","MA_T0","FIELD_T14","FIELD_T0","NH_T14","NH_T0","ME_T14","ME_T0","NS_T14","NS_T0"))) %>%
  ggplot( aes(x=Counts, y=Location, fill = Function)) + geom_bar(stat= "identity", color = "black") +theme_minimal()

#no legend
cog_dat %>%
  arrange(Counts) %>%
  mutate(Location=factor(Location, levels=c("FL_T14","FL_T0","SC_T14","SC_T0","MA_T14","MA_T0","FIELD_T14","FIELD_T0","NH_T14","NH_T0","ME_T14","ME_T0","NS_T14","NS_T0"))) %>%
  ggplot( aes(x=Counts, y=Location, fill = Function)) + geom_bar(stat= "identity", show.legend = FALSE) +theme_minimal() + scale_fill_manual(values=c("#0496ff","#8f2d56","#ffbc42","#d81159","#ff499e","#611c35","#ce2d4f","#9c528b","#006ba6","#ec4e20","#1dd3b0","#3c1642","#b2ff9e","#e2e4f6","#fc440f","#086375","#affc41","#cc76a1","#d6a2ad","#ffb7c3","#f0ead2","#dde5b6","#adc178","#a98467","#6c584c"))



#COG dat mirrored - T0 - 25 colors 
COG_dat_T0 <- read.csv("~/Library/CloudStorage/GoogleDrive-sbirch1@uncc.edu/My Drive/Postdoc_work/Mesocosom_exps/Mesocosom_Summer_2022/Transcriptomics/metacerberus_1st_run/COG_level_1_dat_T0.csv", stringsAsFactors=TRUE)

#with legend
COG_dat_T0 %>%
  arrange(Counts) %>%
  mutate(Location=factor(Location, levels=c("FL_T14","FL_T0","SC_T14","SC_T0","MA_T14","MA_T0","FIELD_T14","FIELD_T0","NH_T14","NH_T0","ME_T14","ME_T0","NS_T14","NS_T0"))) %>%
  ggplot( aes(x=Counts, y=Location, fill = Function)) + geom_bar(stat= "identity", color = "black") +theme_minimal() + scale_fill_manual(values=c("#0496ff","#8f2d56","#ffbc42","#d81159","#ff499e","#611c35","#ce2d4f","#9c528b","#006ba6","#ec4e20","#1dd3b0","#3c1642","#b2ff9e","#e2e4f6","#fc440f","#086375","#affc41","#cc76a1","#d6a2ad","#ffb7c3","#f0ead2","#dde5b6","#adc178","#a98467","#6c584c")) 

#no legend
COG_dat_T0 %>%
  arrange(Counts) %>%
  mutate(Location=factor(Location, levels=c("FL_T14","FL_T0","SC_T14","SC_T0","MA_T14","MA_T0","FIELD_T14","FIELD_T0","NH_T14","NH_T0","ME_T14","ME_T0","NS_T14","NS_T0"))) %>%
  ggplot( aes(x=Counts, y=Location, fill = Function)) + geom_bar(stat= "identity", show.legend = FALSE) +theme_minimal() + scale_fill_manual(values=c("#0496ff","#8f2d56","#ffbc42","#d81159","#ff499e","#611c35","#ce2d4f","#9c528b","#006ba6","#ec4e20","#1dd3b0","#3c1642","#b2ff9e","#e2e4f6","#fc440f","#086375","#affc41","#cc76a1","#d6a2ad","#ffb7c3","#f0ead2","#dde5b6","#adc178","#a98467","#6c584c")) +
  xlim(NA,3000)





#COG dat mirrored - T14
COG_dat_T14 <- read.csv("~/Library/CloudStorage/GoogleDrive-sbirch1@uncc.edu/My Drive/Postdoc_work/Mesocosom_exps/Mesocosom_Summer_2022/Transcriptomics/metacerberus_1st_run/COG_level_1_dat_T14.csv", stringsAsFactors=TRUE)

#with legend
COG_dat_T14 %>%
  arrange(Counts) %>%
  mutate(Location=factor(Location, levels=c("FL_T14","FL_T0","SC_T14","SC_T0","MA_T14","MA_T0","FIELD_T14","FIELD_T0","NH_T14","NH_T0","ME_T14","ME_T0","NS_T14","NS_T0"))) %>%
  ggplot( aes(x=Counts, y=Location, fill = Function)) + geom_bar(stat= "identity", color = "black") +theme_minimal() + scale_fill_manual(values=c("#0496ff","#8f2d56","#ffbc42","#d81159","#ff499e","#611c35","#ce2d4f","#9c528b","#006ba6","#ec4e20","#1dd3b0","#3c1642","#b2ff9e","#e2e4f6","#fc440f","#086375","#affc41","#cc76a1","#d6a2ad","#ffb7c3","#f0ead2","#dde5b6","#adc178","#a98467","#6c584c"))

#no legend
COG_dat_T14 %>%
  arrange(Counts) %>%
  mutate(Location=factor(Location, levels=c("FL_T14","FL_T0","SC_T14","SC_T0","MA_T14","MA_T0","FIELD_T14","FIELD_T0","NH_T14","NH_T0","ME_T14","ME_T0","NS_T14","NS_T0"))) %>%
  ggplot( aes(x=Counts, y=Location, fill = Function)) + geom_bar(stat= "identity", show.legend = FALSE) +theme_minimal() + scale_fill_manual(values=c("#0496ff","#8f2d56","#ffbc42","#d81159","#ff499e","#611c35","#ce2d4f","#9c528b","#006ba6","#ec4e20","#1dd3b0","#3c1642","#b2ff9e","#e2e4f6","#fc440f","#086375","#affc41","#cc76a1","#d6a2ad","#ffb7c3","#f0ead2","#dde5b6","#adc178","#a98467","#6c584c")) +
  xlim(NA,3000)




#### Make stacked bar plot of KEGG Level 2 terms of interest ####

setwd("/Users/Sydney/Library/CloudStorage/GoogleDrive-sbirch1@charlotte.edu/My Drive/Postdoc_work/Mesocosom_exps/Mesocosom_Summer_2022/Transcriptomics/metacerberus_1st_run")

#Read in data (25 colors needed)
kegg_dat <- read.csv("~/Library/CloudStorage/GoogleDrive-sbirch1@charlotte.edu/My Drive/Postdoc_work/Mesocosom_exps/Mesocosom_Summer_2022/Transcriptomics/metacerberus_1st_run/KEGG-barplot-All.csv", stringsAsFactors=TRUE)
str(kegg_dat)

library(ggplot2)
library(dplyr)
library(ggsignif)

#Location x axis
#ggplot(cog_dat, aes(x=Location, y=Counts, fill = Function)) + geom_bar(stat= "identity", color = "black") 

#Counts x axis
#ggplot(cog_dat, aes(x=Counts, y=Location, fill = Function)) + geom_bar(stat= "identity", color = "black")

colors: "#264653","#2a9d8f","#e9c46a","#f4a261","#e76f51","#996888","#eac4d5","#4c2b36","#aba9bf","#db6c79","#BA5624","#721817","#FF70A6","#7C90DB"


#with legend
kegg_dat %>%
  arrange(Counts) %>%
  mutate(Location=factor(Location, levels=c("FL_T14","FL_T0","SC_T14","SC_T0","MA_T14","MA_T0","FIELD_T14","FIELD_T0","NH_T14","NH_T0","ME_T14","ME_T0","NS_T14","NS_T0"))) %>%
  ggplot( aes(x=Counts, y=Location, fill = Function)) + geom_bar(stat= "identity", color = "black") +theme_minimal() + scale_fill_manual(values=c("#264653","#2a9d8f","#e9c46a","#f4a261","#e76f51","#996888","#eac4d5","#4c2b36","#aba9bf","#db6c79","#BA5624","#721817","#FF70A6","#7C90DB"))


#no legend
kegg_dat %>%
  arrange(Counts) %>%
  mutate(Location=factor(Location, levels=c("FL_T14","FL_T0","SC_T14","SC_T0","MA_T14","MA_T0","FIELD_T14","FIELD_T0","NH_T14","NH_T0","ME_T14","ME_T0","NS_T14","NS_T0"))) %>%
  ggplot( aes(x=Counts, y=Location, fill = Function)) + geom_bar(stat= "identity", show.legend = FALSE) +theme_minimal() + scale_fill_manual(values=c("#264653","#2a9d8f","#e9c46a","#f4a261","#e76f51","#996888","#eac4d5","#4c2b36","#aba9bf","#db6c79","#BA5624","#721817","#FF70A6","#7C90DB"))



#kegg dat mirrored - T0 - 14 colors 
kegg_dat_T0 <- read.csv("~/Library/CloudStorage/GoogleDrive-sbirch1@charlotte.edu/My Drive/Postdoc_work/Mesocosom_exps/Mesocosom_Summer_2022/Transcriptomics/metacerberus_1st_run/KEGG-barplot-T0.csv", stringsAsFactors=TRUE)

#with legend
kegg_dat_T0 %>%
  arrange(Counts) %>%
  mutate(Location=factor(Location, levels=c("FL_T14","FL_T0","SC_T14","SC_T0","MA_T14","MA_T0","FIELD_T14","FIELD_T0","NH_T14","NH_T0","ME_T14","ME_T0","NS_T14","NS_T0"))) %>%
  ggplot( aes(x=Counts, y=Location, fill = Function)) + geom_bar(stat= "identity") +theme_minimal() + scale_fill_manual(values=c("#264653","#2a9d8f","#e9c46a","#f4a261","#e76f51","#996888","#eac4d5","#4c2b36","#aba9bf","#db6c79","#BA5624","#721817","#FF70A6","#7C90DB")) 

#no legend
kegg_dat_T0 %>%
  arrange(Counts) %>%
  mutate(Location=factor(Location, levels=c("FL_T14","FL_T0","SC_T14","SC_T0","MA_T14","MA_T0","FIELD_T14","FIELD_T0","NH_T14","NH_T0","ME_T14","ME_T0","NS_T14","NS_T0"))) %>%
  ggplot( aes(x=Counts, y=Location, fill = Function)) + geom_bar(stat= "identity", show.legend = FALSE) +theme_minimal() + scale_fill_manual(values=c("#264653","#2a9d8f","#e9c46a","#f4a261","#e76f51","#996888","#eac4d5","#4c2b36","#aba9bf","#db6c79","#BA5624","#721817","#FF70A6","#7C90DB")) +
  xlim(NA,6500)





#kegg dat mirrored - T14
kegg_dat_T14 <- read.csv("~/Library/CloudStorage/GoogleDrive-sbirch1@charlotte.edu/My Drive/Postdoc_work/Mesocosom_exps/Mesocosom_Summer_2022/Transcriptomics/metacerberus_1st_run/KEGG-barplot-T14.csv", stringsAsFactors=TRUE)
#with legend
kegg_dat_T14 %>%
  arrange(Counts) %>%
  mutate(Location=factor(Location, levels=c("FL_T14","FL_T0","SC_T14","SC_T0","MA_T14","MA_T0","FIELD_T14","FIELD_T0","NH_T14","NH_T0","ME_T14","ME_T0","NS_T14","NS_T0"))) %>%
  ggplot( aes(x=Counts, y=Location, fill = Function)) + geom_bar(stat= "identity", color = "black") +theme_minimal() + scale_fill_manual(values=c("#264653","#2a9d8f","#e9c46a","#f4a261","#e76f51","#996888","#eac4d5","#4c2b36","#aba9bf","#db6c79","#BA5624","#721817","#FF70A6","#7C90DB"))

#no legend
kegg_dat_T14 %>%
  arrange(Counts) %>%
  mutate(Location=factor(Location, levels=c("FL_T14","FL_T0","SC_T14","SC_T0","MA_T14","MA_T0","FIELD_T14","FIELD_T0","NH_T14","NH_T0","ME_T14","ME_T0","NS_T14","NS_T0"))) %>%
  ggplot( aes(x=Counts, y=Location, fill = Function)) + geom_bar(stat= "identity", show.legend = FALSE) +theme_minimal() + scale_fill_manual(values=c("#264653","#2a9d8f","#e9c46a","#f4a261","#e76f51","#996888","#eac4d5","#4c2b36","#aba9bf","#db6c79","#BA5624","#721817","#FF70A6","#7C90DB")) +
  xlim(NA,6500)




#### Make stacked bar plot of FOAM level 1 ####

setwd("/Users/Sydney/Library/CloudStorage/GoogleDrive-sbirch1@charlotte.edu/My Drive/Postdoc_work/Mesocosom_exps/Mesocosom_Summer_2022/Transcriptomics/metacerberus_1st_run")

#Read in data (25 colors needed)
foam_dat <- read.csv("~/Library/CloudStorage/GoogleDrive-sbirch1@charlotte.edu/My Drive/Postdoc_work/Mesocosom_exps/Mesocosom_Summer_2022/Transcriptomics/metacerberus_1st_run/FOAM-barplot_all.csv", stringsAsFactors=TRUE)
str(foam_dat)

library(ggplot2)
library(dplyr)
library(ggsignif)


#with legend
foam_dat %>%
  arrange(Counts) %>%
  mutate(Location=factor(Location, levels=c("FL_T14","FL_T0","SC_T14","SC_T0","MA_T14","MA_T0","FIELD_T14","FIELD_T0","NH_T14","NH_T0","ME_T14","ME_T0","NS_T14","NS_T0"))) %>%
  ggplot( aes(x=Counts, y=Location, fill = Function)) + geom_bar(stat= "identity", color = "black") +theme_minimal() + scale_fill_manual(values=c("#ffbc42","#d81159","#611c35","#ce2d4f","#9c528b","#006ba6","#ec4e20","#1dd3b0","#3c1642","#b2ff9e","#e2e4f6","#fc440f","#086375","#affc41","#cc76a1","#d6a2ad","#ffb7c3","#f0ead2","#adc178","#a98467","#6c584c"))


#no legend
foam_dat %>%
  arrange(Counts) %>%
  mutate(Location=factor(Location, levels=c("FL_T14","FL_T0","SC_T14","SC_T0","MA_T14","MA_T0","FIELD_T14","FIELD_T0","NH_T14","NH_T0","ME_T14","ME_T0","NS_T14","NS_T0"))) %>%
  ggplot( aes(x=Counts, y=Location, fill = Function)) + geom_bar(stat= "identity", show.legend = FALSE) +theme_minimal() + scale_fill_manual(values=c("#ffbc42","#d81159","#611c35","#ce2d4f","#9c528b","#006ba6","#ec4e20","#1dd3b0","#3c1642","#b2ff9e","#e2e4f6","#fc440f","#086375","#affc41","#cc76a1","#d6a2ad","#ffb7c3","#f0ead2","#adc178","#a98467","#6c584c"))



#kegg dat mirrored - T0 - 14 colors 
foam_dat_T0 <- read.csv("~/Library/CloudStorage/GoogleDrive-sbirch1@charlotte.edu/My Drive/Postdoc_work/Mesocosom_exps/Mesocosom_Summer_2022/Transcriptomics/metacerberus_1st_run/FOAM-barplot-T0.csv", stringsAsFactors=TRUE)

#with legend
foam_dat_T0 %>%
  arrange(Counts) %>%
  mutate(Location=factor(Location, levels=c("FL_T14","FL_T0","SC_T14","SC_T0","MA_T14","MA_T0","FIELD_T14","FIELD_T0","NH_T14","NH_T0","ME_T14","ME_T0","NS_T14","NS_T0"))) %>%
  ggplot( aes(x=Counts, y=Location, fill = Function)) + geom_bar(stat= "identity") +theme_minimal() + scale_fill_manual(values=c("#ffbc42","#d81159","#611c35","#ce2d4f","#9c528b","#006ba6","#ec4e20","#1dd3b0","#3c1642","#b2ff9e","#e2e4f6","#fc440f","#086375","#affc41","#cc76a1","#d6a2ad","#ffb7c3","#f0ead2","#adc178","#a98467","#6c584c")) 

#no legend
foam_dat_T0 %>%
  arrange(Counts) %>%
  mutate(Location=factor(Location, levels=c("FL_T14","FL_T0","SC_T14","SC_T0","MA_T14","MA_T0","FIELD_T14","FIELD_T0","NH_T14","NH_T0","ME_T14","ME_T0","NS_T14","NS_T0"))) %>%
  ggplot( aes(x=Counts, y=Location, fill = Function)) + geom_bar(stat= "identity", show.legend = FALSE) +theme_minimal() + scale_fill_manual(values=c("#ffbc42","#d81159","#611c35","#ce2d4f","#9c528b","#006ba6","#ec4e20","#1dd3b0","#3c1642","#b2ff9e","#e2e4f6","#fc440f","#086375","#affc41","#cc76a1","#d6a2ad","#ffb7c3","#f0ead2","#adc178","#a98467","#6c584c")) +
  xlim(NA,2500)



#kegg dat mirrored - T14
foam_dat_T14 <- read.csv("~/Library/CloudStorage/GoogleDrive-sbirch1@charlotte.edu/My Drive/Postdoc_work/Mesocosom_exps/Mesocosom_Summer_2022/Transcriptomics/metacerberus_1st_run/FOAM-barplot-T14.csv", stringsAsFactors=TRUE)

#with legend
foam_dat_T14 %>%
  arrange(Counts) %>%
  mutate(Location=factor(Location, levels=c("FL_T14","FL_T0","SC_T14","SC_T0","MA_T14","MA_T0","FIELD_T14","FIELD_T0","NH_T14","NH_T0","ME_T14","ME_T0","NS_T14","NS_T0"))) %>%
  ggplot( aes(x=Counts, y=Location, fill = Function)) + geom_bar(stat= "identity", color = "black") +theme_minimal() + scale_fill_manual(values=c("#ffbc42","#d81159","#611c35","#ce2d4f","#9c528b","#006ba6","#ec4e20","#1dd3b0","#3c1642","#b2ff9e","#e2e4f6","#fc440f","#086375","#affc41","#cc76a1","#d6a2ad","#ffb7c3","#f0ead2","#adc178","#a98467","#6c584c"))

#no legend
foam_dat_T14 %>%
  arrange(Counts) %>%
  mutate(Location=factor(Location, levels=c("FL_T14","FL_T0","SC_T14","SC_T0","MA_T14","MA_T0","FIELD_T14","FIELD_T0","NH_T14","NH_T0","ME_T14","ME_T0","NS_T14","NS_T0"))) %>%
  ggplot( aes(x=Counts, y=Location, fill = Function)) + geom_bar(stat= "identity", show.legend = FALSE) +theme_minimal() + scale_fill_manual(values=c("#ffbc42","#d81159","#611c35","#ce2d4f","#9c528b","#006ba6","#ec4e20","#1dd3b0","#3c1642","#b2ff9e","#e2e4f6","#fc440f","#086375","#affc41","#cc76a1","#d6a2ad","#ffb7c3","#f0ead2","#adc178","#a98467","#6c584c")) +
  xlim(NA,2500)




## ANOVAs ##

setwd("/Users/Sydney/Library/CloudStorage/GoogleDrive-sbirch1@charlotte.edu/My Drive/Postdoc_work/Mesocosom_exps/Mesocosom_Summer_2022/Transcriptomics/metacerberus_1st_run")

Kegg_anova_dat <- read.csv("~/Library/CloudStorage/GoogleDrive-sbirch1@charlotte.edu/My Drive/Postdoc_work/Mesocosom_exps/Mesocosom_Summer_2022/Transcriptomics/metacerberus_1st_run/Kegg_anova_dat.csv", stringsAsFactors=TRUE)
str(Kegg_anova_dat)


stats_mod1<-lm(Environmental.adaptation ~ Time_point*comparison, Kegg_anova_dat) 
anova(stats_mod1)

stats_mod2<-lm(Immune.system ~ Time_point*comparison, Kegg_anova_dat) 
anova(stats_mod2)

stats_mod3<-lm(Viral.protein.families ~ Time_point*comparison, Kegg_anova_dat) 
anova(stats_mod3)

stats_mod4<-lm(Viral.protein.families ~ Time_point*comparison, Kegg_anova_dat) 
anova(stats_mod4)

stats_mod5<-lm(Signaling.molecules.and.interaction ~ Time_point*comparison, Kegg_anova_dat) 
anova(stats_mod5)

stats_mod6<-lm(Protein.families..genetic.information.processing   ~ Time_point*comparison, Kegg_anova_dat) 
anova(stats_mod6)

stats_mod7<-lm(Protein.families..metabolism ~ Time_point*comparison, Kegg_anova_dat) 
anova(stats_mod7)

stats_mod8<-lm(Protein.families..signaling.and.cellular.processes   ~ Time_point*comparison, Kegg_anova_dat) 
anova(stats_mod8)





#North vs South 
boxplot(Immune.system ~ comparison, Kegg_anova_dat, main = "Effect of North vs South on env adapt",
        xlab="comparison", ylab= "Counts")

#Time point
boxplot(Immune.system ~ Time_point, Kegg_anova_dat, main = "The Effect of time pint on env adapt",
        xlab="time point", ylab= "Counts")

#Location
boxplot(Immune.system ~ Location, Kegg_anova_dat, main = "The Effect of Location on env adapt",
        xlab="Location", ylab= "Counts")
