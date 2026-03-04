## Metacerberus MAPPED NEMATOSTELLA DATA 

#This is all looking at the stats output from metacerberus 


setwd("/Users/Sydney/Library/CloudStorage/GoogleDrive-sbirch1@charlotte.edu/My Drive/Postdoc_work/Mesocosom_exps/Mesocosom_Summer_2022/Transcriptomics/Metacerberus_2nd_run/MAPPED_All_comparisons")

#Import data
stats_db <- read.csv("~/Library/CloudStorage/GoogleDrive-sbirch1@charlotte.edu/My Drive/Postdoc_work/Mesocosom_exps/Mesocosom_Summer_2022/Transcriptomics/Metacerberus_2nd_run/MAPPED_All_comparisons/stats_db_mapped.csv", stringsAsFactors=TRUE)
str(stats_db)


#one way anova --> comparison, KEGG
stats_mod<-lm(KEGG_count ~ comparison, stats_db) 
anova(stats_mod)
#Analysis of Variance Table
#Response: KEGG_count
#Df   Sum Sq  Mean Sq F value    Pr(>F)    
#comparison  2 52168123 26084061  46.341 2.207e-07 ***
#Residuals  16  9006030   562877


#one way anova --> comparison, FOAM
stats_mod<-lm(FOAM_count ~ comparison, stats_db) 
anova(stats_mod)

#Analysis of Variance Table
#Response: FOAM_count
#           Df Sum Sq Mean Sq F value    Pr(>F)    
#comparison  2 569159  284579  60.314 3.537e-08 ***
#Residuals  16  75493    4718                      
---
  #Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
  
  
  #one way anova --> comparison, Protien
  stats_mod<-lm(Protien_count ~ comparison, stats_db) 
anova(stats_mod)

#Analysis of Variance Table
#Response: Protien_count
#.          Df   Sum Sq  Mean Sq F value    Pr(>F)    
#comparison  2 34704323 17352161  33.057 2.078e-06 ***
#Residuals  16  8398769   524923 



## Making Barplots with Locations Ordered

#KEGG subset data
stats_kegg <- read.csv("~/Library/CloudStorage/GoogleDrive-sbirch1@charlotte.edu/My Drive/Postdoc_work/Mesocosom_exps/Mesocosom_Summer_2022/Transcriptomics/Metacerberus_2nd_run/MAPPED_All_comparisons/stats_kegg-maped.csv", stringsAsFactors=TRUE)
str(stats_kegg)

#VOG ID Counts
library(dplyr)
library(ggsignif)
stats_kegg %>%
  arrange(KEGG_count) %>%
  mutate(Location=factor(Location, levels=c("NS","ME","NH","FIELD","MA","SC","FL"))) %>%
  ggplot( aes(x=Location,y=KEGG_count,fill=comparison)) +
  geom_bar(stat="identity", color="black", position=position_dodge())+theme_minimal() +
  scale_fill_manual(values=c("#A8DADC","#457B9D", "#1D3557")) + ggtitle("KEGG ID Counts") 


#FOAM subset data
stats_foam <- read.csv("~/Library/CloudStorage/GoogleDrive-sbirch1@charlotte.edu/My Drive/Postdoc_work/Mesocosom_exps/Mesocosom_Summer_2022/Transcriptomics/Metacerberus_2nd_run/MAPPED_All_comparisons/stats_foam-mapped.csv", stringsAsFactors=TRUE)
str(stats_foam)

stats_foam %>%
  arrange(FOAM_count) %>%
  mutate(Location=factor(Location, levels=c("NS","ME","NH","FIELD","MA","SC","FL"))) %>%
  ggplot( aes(x=Location,y=FOAM_count,fill=comparison)) +
  geom_bar(stat="identity", color="black", position=position_dodge())+theme_minimal() +
  scale_fill_manual(values=c("#E5989B","#B5838D","#6D6875" )) + ggtitle("FOAM ID Counts")


#Protein subset data
stats_prot <- read.csv("~/Library/CloudStorage/GoogleDrive-sbirch1@charlotte.edu/My Drive/Postdoc_work/Mesocosom_exps/Mesocosom_Summer_2022/Transcriptomics/Metacerberus_2nd_run/MAPPED_All_comparisons/stats_protien-mapped.csv", stringsAsFactors=TRUE)
str(stats_prot)

stats_prot %>%
  arrange(Protien_count) %>%
  mutate(Location=factor(Location, levels=c("NS","ME","NH","FIELD","MA","SC","FL"))) %>%
  ggplot(aes(x=Location,y=Protien_count,fill=comparison)) +
  geom_bar(stat="identity", color="black", position=position_dodge())+theme_minimal() +
  scale_fill_manual(values=c("#A3B18A","#588157","#3A5A40")) + ggtitle("Protein Counts")




#### Make stacked bar plot of Foam Level 1 terms ####

setwd("/Users/Sydney/Library/CloudStorage/GoogleDrive-sbirch1@charlotte.edu/My Drive/Postdoc_work/Mesocosom_exps/Mesocosom_Summer_2022/Transcriptomics/Metacerberus_2nd_run/MAPPED_All_comparisons")

## FOAM T0 vs CTRL ## 

#Read in data
foam_T0_vs_CTRL <- read.csv("~/Library/CloudStorage/GoogleDrive-sbirch1@charlotte.edu/My Drive/Postdoc_work/Mesocosom_exps/Mesocosom_Summer_2022/Transcriptomics/Metacerberus_2nd_run/MAPPED_All_comparisons/Mapped_FOAM-T0_vs_CTRL_dat.csv", stringsAsFactors=TRUE)
str(foam_T0_vs_CTRL)

library(ggplot2)
library(dplyr)
library(ggsignif)


#With Legend (21 colors needed)
foam_T0_vs_CTRL %>%
  arrange(Counts) %>%
  mutate(Location=factor(Location, levels=c("FL_T0_vs_CTRL","SC_T0_vs_CTRL","MA_T0_vs_CTRL","NH_T0_vs_CTRL","ME_T0_vs_CTRL","NS_T0_vs_CTRL"))) %>%
  ggplot( aes(x=Counts, y=Location, fill = Function)) + geom_bar(stat= "identity", color = "black") +theme_minimal() + scale_fill_manual(values=c("#ffbc42","#d81159","#611c35","#ce2d4f","#9c528b","#006ba6","#ec4e20","#1dd3b0","#3c1642","#b2ff9e","#e2e4f6","#fc440f","#086375","#affc41","#cc76a1","#d6a2ad","#ffb7c3","#f0ead2","#adc178","#a98467","#6c584c" ))

#no legend
foam_T0_vs_CTRL %>%
  arrange(Counts) %>%
  mutate(Location=factor(Location, levels=c("FL_T0_vs_CTRL","SC_T0_vs_CTRL","MA_T0_vs_CTRL","NH_T0_vs_CTRL","ME_T0_vs_CTRL","NS_T0_vs_CTRL"))) %>%
  ggplot( aes(x=Counts, y=Location, fill = Function)) + geom_bar(stat= "identity", show.legend = FALSE) +theme_minimal() + scale_fill_manual(values=c("#ffbc42","#d81159","#611c35","#ce2d4f","#9c528b","#006ba6","#ec4e20","#1dd3b0","#3c1642","#b2ff9e","#e2e4f6","#fc440f","#086375","#affc41","#cc76a1","#d6a2ad","#ffb7c3","#f0ead2","#adc178","#a98467","#6c584c" )) +
  xlim(NA,500)



## FOAM T14 vs CTRL ## 

#Read in data
foam_T14_vs_CTRL <- read.csv("~/Library/CloudStorage/GoogleDrive-sbirch1@charlotte.edu/My Drive/Postdoc_work/Mesocosom_exps/Mesocosom_Summer_2022/Transcriptomics/Metacerberus_2nd_run/MAPPED_All_comparisons/Mapped_FOAM-T14_vs_CTRL_dat.csv", stringsAsFactors=TRUE)
str(foam_T14_vs_CTRL)

library(ggplot2)
library(dplyr)
library(ggsignif)


#With Legend (21 colors needed)
foam_T14_vs_CTRL %>%
  arrange(Counts) %>%
  mutate(Location=factor(Location, levels=c("FL_T14_vs_CTRL","SC_T14_vs_CTRL","MA_T14_vs_CTRL","NH_T14_vs_CTRL","ME_T14_vs_CTRL","NS_T14_vs_CTRL"))) %>%
  ggplot( aes(x=Counts, y=Location, fill = Function)) + geom_bar(stat= "identity", color = "black") +theme_minimal() + scale_fill_manual(values=c("#ffbc42","#d81159","#611c35","#ce2d4f","#9c528b","#006ba6","#ec4e20","#1dd3b0","#3c1642","#b2ff9e","#e2e4f6","#fc440f","#086375","#affc41","#cc76a1","#d6a2ad","#ffb7c3","#f0ead2","#adc178","#a98467","#6c584c" ))

#no legend
foam_T14_vs_CTRL %>%
  arrange(Counts) %>%
  mutate(Location=factor(Location, levels=c("FL_T14_vs_CTRL","SC_T14_vs_CTRL","MA_T14_vs_CTRL","NH_T14_vs_CTRL","ME_T14_vs_CTRL","NS_T14_vs_CTRL"))) %>%
  ggplot( aes(x=Counts, y=Location, fill = Function)) + geom_bar(stat= "identity", show.legend = FALSE) +theme_minimal() + scale_fill_manual(values=c("#ffbc42","#d81159","#611c35","#ce2d4f","#9c528b","#006ba6","#ec4e20","#1dd3b0","#3c1642","#b2ff9e","#e2e4f6","#fc440f","#086375","#affc41","#cc76a1","#d6a2ad","#ffb7c3","#f0ead2","#adc178","#a98467","#6c584c" )) +
  xlim(NA,500)



## FOAM T0 vs T14s in each location ## 

#Read in data
foam_T0_vs_T14s <- read.csv("~/Library/CloudStorage/GoogleDrive-sbirch1@charlotte.edu/My Drive/Postdoc_work/Mesocosom_exps/Mesocosom_Summer_2022/Transcriptomics/Metacerberus_2nd_run/MAPPED_All_comparisons/Mapped_FOAM_T0_vs_T14s.csv", stringsAsFactors=TRUE)
str(foam_T0_vs_T14s)

library(ggplot2)
library(dplyr)
library(ggsignif)


#With Legend (21 colors needed)
foam_T0_vs_T14s %>%
  arrange(Counts) %>%
  mutate(Location=factor(Location, levels=c("FL_T0_vs_T14","SC_T0_vs_T14","MA_T0_vs_T14","FIELD_T0_vs_T14","NH_T0_vs_T14","ME_T0_vs_T14","NS_T0_vs_T14"))) %>%
  ggplot( aes(x=Counts, y=Location, fill = Function)) + geom_bar(stat= "identity", color = "black") +theme_minimal() + scale_fill_manual(values=c("#ffbc42","#d81159","#611c35","#ce2d4f","#9c528b","#006ba6","#ec4e20","#1dd3b0","#3c1642","#b2ff9e","#e2e4f6","#fc440f","#086375","#affc41","#cc76a1","#d6a2ad","#ffb7c3","#f0ead2","#adc178","#a98467","#6c584c" ))

#no legend
foam_T0_vs_T14s %>%
  arrange(Counts) %>%
  mutate(Location=factor(Location, levels=c("FL_T0_vs_T14","SC_T0_vs_T14","MA_T0_vs_T14","FIELD_T0_vs_T14","NH_T0_vs_T14","ME_T0_vs_T14","NS_T0_vs_T14"))) %>%
  ggplot( aes(x=Counts, y=Location, fill = Function)) + geom_bar(stat= "identity", show.legend = FALSE) +theme_minimal() + scale_fill_manual(values=c("#ffbc42","#d81159","#611c35","#ce2d4f","#9c528b","#006ba6","#ec4e20","#1dd3b0","#3c1642","#b2ff9e","#e2e4f6","#fc440f","#086375","#affc41","#cc76a1","#d6a2ad","#ffb7c3","#f0ead2","#adc178","#a98467","#6c584c" )) 





#### Make stacked bar plot of KEGG Level 2 terms ####

setwd("/Users/Sydney/Library/CloudStorage/GoogleDrive-sbirch1@charlotte.edu/My Drive/Postdoc_work/Mesocosom_exps/Mesocosom_Summer_2022/Transcriptomics/Metacerberus_2nd_run/MAPPED_All_comparisons")

## KEGG T0 vs T14 Viral and Gen ## 

#Read in data
kegg_t0vt14 <- read.csv("~/Library/CloudStorage/GoogleDrive-sbirch1@charlotte.edu/My Drive/Postdoc_work/Mesocosom_exps/Mesocosom_Summer_2022/Transcriptomics/Metacerberus_2nd_run/MAPPED_All_comparisons/Mapped_Kegg-virall_and_gen.csv")
str(kegg_t0vt14)

library(ggplot2)
library(dplyr)
library(ggsignif)


#With Legend (21 colors needed)
kegg_t0vt14 %>%
  arrange(Counts) %>%
  mutate(Comparison=factor(Comparison, levels=c("FL_T0_V_T14","SC_T0_V_T14","MA_T0_V_T14","NH_T0_V_T14","ME_T0_V_T14","NS_T0_V_T14"))) %>%
  ggplot( aes(x=Counts, y=Comparison, fill = Function)) + geom_bar(stat= "identity", color = "black") +theme_minimal() + scale_fill_manual(values=c("#ffbc42","#d81159","#611c35","#ce2d4f","#9c528b","#006ba6","#ec4e20","#1dd3b0","#3c1642","#b2ff9e","#e2e4f6","#fc440f" ))

#no legend
kegg_t0vt14 %>%
  arrange(Counts) %>%
  mutate(Comparison=factor(Comparison, levels=c("FL_T0_V_T14","SC_T0_V_T14","MA_T0_V_T14","FIELD_T0_V_T14", "NH_T0_V_T14","ME_T0_V_T14","NS_T0_V_T14"))) %>%
  ggplot( aes(x=Counts, y=Comparison, fill = Function)) + geom_bar(stat= "identity", show.legend = FALSE) +theme_minimal() + scale_fill_manual(values=c("#ffbc42","#d81159","#611c35","#ce2d4f","#9c528b","#006ba6","#ec4e20","#1dd3b0","#3c1642","#b2ff9e","#e2e4f6","#fc440f" )) 



## KEGG T0 vs T14 Cellular ## 

#Read in data
kegg_t0vt14_cell <- read.csv("~/Library/CloudStorage/GoogleDrive-sbirch1@charlotte.edu/My Drive/Postdoc_work/Mesocosom_exps/Mesocosom_Summer_2022/Transcriptomics/Metacerberus_2nd_run/MAPPED_All_comparisons/mapped_Kegg-cellular.csv")
str(kegg_t0vt14_cell)

library(ggplot2)
library(dplyr)
library(ggsignif)


#With Legend (21 colors needed)
kegg_t0vt14_cell %>%
  arrange(Counts) %>%
  mutate(Comparison=factor(Comparison, levels=c("FL_T0_V_T14","SC_T0_V_T14","MA_T0_V_T14","NH_T0_V_T14","ME_T0_V_T14","NS_T0_V_T14"))) %>%
  ggplot( aes(x=Counts, y=Comparison, fill = Function)) + geom_bar(stat= "identity", color = "black") +theme_minimal() + scale_fill_manual(values=c("#ffbc42","#d81159","#611c35","#ffb7c3","#9c528b","#006ba6","#ec4e20","#1dd3b0","#3c1642","#b2ff9e","#e2e4f6","#fc440f" ))

#no legend
kegg_t0vt14_cell %>%
  arrange(Counts) %>%
  mutate(Comparison=factor(Comparison, levels=c("FL_T0_V_T14","SC_T0_V_T14","MA_T0_V_T14","FIELD_T0_V_T14", "NH_T0_V_T14","ME_T0_V_T14","NS_T0_V_T14"))) %>%
  ggplot( aes(x=Counts, y=Comparison, fill = Function)) + geom_bar(stat= "identity", show.legend = FALSE) +theme_minimal() + scale_fill_manual(values=c("#ffbc42","#d81159","#611c35","#ffb7c3","#9c528b","#006ba6","#ec4e20","#1dd3b0","#3c1642","#b2ff9e","#e2e4f6","#fc440f" )) 




## KEGG T0 vs CTRL Viral and gen ## 

#Read in data
kegg_t0vCtrl_v <- read.csv("~/Library/CloudStorage/GoogleDrive-sbirch1@charlotte.edu/My Drive/Postdoc_work/Mesocosom_exps/Mesocosom_Summer_2022/Transcriptomics/Metacerberus_2nd_run/MAPPED_All_comparisons/Mapped_kegg_T0vs_Ctrl_viral_gen.csv")
str(kegg_t0vCtrl_v)

library(ggplot2)
library(dplyr)
library(ggsignif)


#With Legend (13 colors needed)
kegg_t0vCtrl_v %>%
  arrange(Counts) %>%
  mutate(Comparison=factor(Comparison, levels=c("FL_T0_V_CTRL","SC_T0_V_CTRL","MA_T0_V_CTRL","NH_T0_V_CTRL","ME_T0_V_CTRL","NS_T0_V_CTRL"))) %>%
  ggplot( aes(x=Counts, y=Comparison, fill = Function)) + geom_bar(stat= "identity", color = "black") +theme_minimal() + scale_fill_manual(values=c("#ffbc42","#d81159","#611c35","#ffb7c3","#9c528b","#006ba6","#ec4e20","#1dd3b0","#3c1642","#b2ff9e","#e2e4f6","#fc440f","#BDC696" ))

#no legend
kegg_t0vCtrl_v %>%
  arrange(Counts) %>%
  mutate(Comparison=factor(Comparison, levels=c("FL_T0_V_CTRL","SC_T0_V_CTRL","MA_T0_V_CTRL","NH_T0_V_CTRL","ME_T0_V_CTRL","NS_T0_V_CTRL"))) %>%
  ggplot( aes(x=Counts, y=Comparison, fill = Function)) + geom_bar(stat= "identity", show.legend = FALSE) +theme_minimal() + scale_fill_manual(values=c("#ffbc42","#d81159","#611c35","#ffb7c3","#9c528b","#006ba6","#ec4e20","#1dd3b0","#3c1642","#b2ff9e","#e2e4f6","#fc440f","#BDC696" )) +
  xlim(NA,3000)



## KEGG T0 vs CTRL Cellular ## 

#Read in data
kegg_t0vCtrl_c <- read.csv("~/Library/CloudStorage/GoogleDrive-sbirch1@charlotte.edu/My Drive/Postdoc_work/Mesocosom_exps/Mesocosom_Summer_2022/Transcriptomics/Metacerberus_2nd_run/MAPPED_All_comparisons/Mapped_Kegg_T0_vs_Ctrl_celluar.csv")
str(kegg_t0vCtrl_c)

library(ggplot2)
library(dplyr)
library(ggsignif)


#With Legend (13 colors needed)
kegg_t0vCtrl_c %>%
  arrange(Counts) %>%
  mutate(Comparison=factor(Comparison, levels=c("FL_T0_V_CTRL","SC_T0_V_CTRL","MA_T0_V_CTRL","NH_T0_V_CTRL","ME_T0_V_CTRL","NS_T0_V_CTRL"))) %>%
  ggplot( aes(x=Counts, y=Comparison, fill = Function)) + geom_bar(stat= "identity", color = "black") +theme_minimal() + scale_fill_manual(values=c("#ffbc42","#d81159","#611c35","#ffb7c3","#9c528b","#006ba6","#ec4e20","#1dd3b0","#3c1642","#b2ff9e","#e2e4f6","#fc440f","#BDC696" ))

#no legend
kegg_t0vCtrl_c %>%
  arrange(Counts) %>%
  mutate(Comparison=factor(Comparison, levels=c("FL_T0_V_CTRL","SC_T0_V_CTRL","MA_T0_V_CTRL","NH_T0_V_CTRL","ME_T0_V_CTRL","NS_T0_V_CTRL"))) %>%
  ggplot( aes(x=Counts, y=Comparison, fill = Function)) + geom_bar(stat= "identity", show.legend = FALSE) +theme_minimal() + scale_fill_manual(values=c("#ffbc42","#d81159","#611c35","#ffb7c3","#9c528b","#006ba6","#ec4e20","#1dd3b0","#3c1642","#b2ff9e","#e2e4f6","#fc440f","#BDC696" )) +
  xlim(NA,1000)






## KEGG T14 vs CTRL Viral and gen ## 

#Read in data
kegg_t14vCtrl_v <- read.csv("~/Library/CloudStorage/GoogleDrive-sbirch1@charlotte.edu/My Drive/Postdoc_work/Mesocosom_exps/Mesocosom_Summer_2022/Transcriptomics/Metacerberus_2nd_run/MAPPED_All_comparisons/Mapped_kegg_T14_ctrl_viral_gen.csv")
str(kegg_t14vCtrl_v)

library(ggplot2)
library(dplyr)
library(ggsignif)


#With Legend (13 colors needed)
kegg_t14vCtrl_v %>%
  arrange(Counts) %>%
  mutate(Comparison=factor(Comparison, levels=c("FL_T14_V_CTRL","SC_T14_V_CTRL","MA_T14_V_CTRL","NH_T14_V_CTRL","ME_T14_V_CTRL","NS_T14_V_CTRL"))) %>%
  ggplot( aes(x=Counts, y=Comparison, fill = Function)) + geom_bar(stat= "identity", color = "black") +theme_minimal() + scale_fill_manual(values=c("#ffbc42","#d81159","#611c35","#ffb7c3","#9c528b","#006ba6","#ec4e20","#1dd3b0","#3c1642","#b2ff9e","#e2e4f6","#fc440f","#BDC696" ))

#no legend
kegg_t14vCtrl_v %>%
  arrange(Counts) %>%
  mutate(Comparison=factor(Comparison, levels=c("FL_T14_V_CTRL","SC_T14_V_CTRL","MA_T14_V_CTRL","NH_T14_V_CTRL","ME_T14_V_CTRL","NS_T14_V_CTRL"))) %>%
  ggplot( aes(x=Counts, y=Comparison, fill = Function)) + geom_bar(stat= "identity", show.legend = FALSE) +theme_minimal() + scale_fill_manual(values=c("#ffbc42","#d81159","#611c35","#ffb7c3","#9c528b","#006ba6","#ec4e20","#1dd3b0","#3c1642","#b2ff9e","#e2e4f6","#fc440f","#BDC696" )) +
  xlim(NA,3000)



## KEGG T14 vs CTRL Cellular ## 

#Read in data
kegg_t14vCtrl_c <- read.csv("~/Library/CloudStorage/GoogleDrive-sbirch1@charlotte.edu/My Drive/Postdoc_work/Mesocosom_exps/Mesocosom_Summer_2022/Transcriptomics/Metacerberus_2nd_run/MAPPED_All_comparisons/Mapped_Kegg_T14_vs_Ctrl_cellular.csv")
str(kegg_t0vCtrl_c)

library(ggplot2)
library(dplyr)
library(ggsignif)


#With Legend (13 colors needed)
kegg_t14vCtrl_c %>%
  arrange(Counts) %>%
  mutate(Comparison=factor(Comparison, levels=c("FL_T14_V_CTRL","SC_T14_V_CTRL","MA_T14_V_CTRL","NH_T14_V_CTRL","ME_T14_V_CTRL","NS_T14_V_CTRL"))) %>%
  ggplot( aes(x=Counts, y=Comparison, fill = Function)) + geom_bar(stat= "identity", color = "black") +theme_minimal() + scale_fill_manual(values=c("#ffbc42","#d81159","#611c35","#ffb7c3","#9c528b","#006ba6","#ec4e20","#1dd3b0","#3c1642","#b2ff9e","#e2e4f6","#fc440f","#BDC696" ))

#no legend
kegg_t14vCtrl_c %>%
  arrange(Counts) %>%
  mutate(Comparison=factor(Comparison, levels=c("FL_T14_V_CTRL","SC_T14_V_CTRL","MA_T14_V_CTRL","NH_T14_V_CTRL","ME_T14_V_CTRL","NS_T14_V_CTRL"))) %>%
  ggplot( aes(x=Counts, y=Comparison, fill = Function)) + geom_bar(stat= "identity", show.legend = FALSE) +theme_minimal() + scale_fill_manual(values=c("#ffbc42","#d81159","#611c35","#ffb7c3","#9c528b","#006ba6","#ec4e20","#1dd3b0","#3c1642","#b2ff9e","#e2e4f6","#fc440f","#BDC696" )) +
  xlim(NA,1000)



