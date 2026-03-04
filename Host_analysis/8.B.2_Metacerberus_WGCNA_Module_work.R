## Metacerberus WGCNA Module stacked bar plots and stats 

#This is all looking at the stats output from metacerberus 

setwd("~/Library/CloudStorage/GoogleDrive-sbirch1@charlotte.edu/My Drive/Postdoc_work/Mesocosom_exps/Mesocosom_Summer_2022/Transcriptomics/EdgeR_MAPPED_4-29-24/WGCNA_Seqs")


#Import data
metecerb_module_stats_2 <- read.csv("~/Library/CloudStorage/GoogleDrive-sbirch1@charlotte.edu/My Drive/Postdoc_work/Mesocosom_exps/Mesocosom_Summer_2022/Transcriptomics/EdgeR_MAPPED_4-29-24/WGCNA_Seqs/metecerb_module_stats_2.csv", stringsAsFactors=TRUE)
str(metecerb_module_stats_2)

#Make a barplot

library(ggplot2)
#Protien stats Counts
ggplot(data=metecerb_module_stats_2, aes(x=Module,y=Count,fill=Type)) +
  geom_bar(stat="identity", color="black", position=position_dodge())+theme_minimal() +
  scale_fill_manual(values=c('#335C67','#E09F3E','#9E2A2B', '#540B0E')) +  ggtitle("Module Prot Counts")


#### Make stacked bar plots KEGG level 2 ####

library(ggplot2)
library(dplyr)
library(ggsignif)


#Stacked barplot Kegg Level-2 combined
kegg_L2 <- read.csv("~/Library/CloudStorage/GoogleDrive-sbirch1@charlotte.edu/My Drive/Postdoc_work/Mesocosom_exps/Mesocosom_Summer_2022/Transcriptomics/EdgeR_MAPPED_4-29-24/WGCNA_Seqs/combined_kegg_level2-V2-subset.csv", stringsAsFactors=TRUE)
str(kegg_L2)

#With Legend (12 colors needed)
kegg_L2 %>%
  arrange(Count) %>%
  mutate(Module=factor(Module, levels=c("ME_33","ME_21","ME_60","ME_17","ME_96"))) %>%
  ggplot( aes(x=Count, y=Module, fill = Type)) + geom_bar(stat= "identity", color = "black") +theme_minimal() + scale_fill_manual(values=c("#127475","#8cb369","#f4e285","#f4a259","#5b8e7d","#bc4b51","#2d2d2a","#0c4767","#820933","#34623f","#546a76","#562c2c"))

#no legend
kegg_L2 %>%
  arrange(Count) %>%
  mutate(Module=factor(Module, levels=c("ME_33","ME_21","ME_60","ME_17","ME_96"))) %>%
  ggplot( aes(x=Count, y=Module, fill = Type)) + geom_bar(stat= "identity", show.legend = FALSE) +theme_minimal() + scale_fill_manual(values=c("#127475","#8cb369","#f4e285","#f4a259","#5b8e7d","#bc4b51","#2d2d2a","#0c4767","#820933","#34623f","#546a76","#562c2c"))



#Stacked barplot FOAM Level-1 combined

foam_L1 <- read.csv("~/Library/CloudStorage/GoogleDrive-sbirch1@charlotte.edu/My Drive/Postdoc_work/Mesocosom_exps/Mesocosom_Summer_2022/Transcriptomics/EdgeR_MAPPED_4-29-24/WGCNA_Seqs/FOAM_Level_1_combined_V2.csv", stringsAsFactors=TRUE)
str(foam_L1)

#With Legend (16 colors needed)
foam_L1 %>%
  arrange(Count) %>%
  mutate(Module=factor(Module, levels=c("ME_33","ME_21","ME_60","ME_17","ME_96"))) %>%
  ggplot( aes(x=Count, y=Module, fill = Type)) + geom_bar(stat= "identity", color = "black") +theme_minimal() + scale_fill_manual(values=c("#127475","#8cb369","#f4e285","#f4a259","#5b8e7d","#bc4b51","#2d2d2a","#0c4767","#820933","#34623f","#546a76","#562c2c","#A09ABC","#669BBC","#FDF0D5",'#C1121F'))

#no legend
foam_L1 %>%
  arrange(Count) %>%
  mutate(Module=factor(Module, levels=c("ME_33","ME_21","ME_60","ME_17","ME_96"))) %>%
  ggplot( aes(x=Count, y=Module, fill = Type)) + geom_bar(stat= "identity", show.legend = FALSE) +theme_minimal() + scale_fill_manual(values=c("#127475","#8cb369","#f4e285","#f4a259","#5b8e7d","#bc4b51","#2d2d2a","#0c4767","#820933","#34623f","#546a76","#562c2c","#A09ABC","#669BBC","#FDF0D5",'#C1121F'))



### Individual barchart - rel abundance
"#001219", "#005f73", "#0a9396", "#94d2bd", "#e9d8a6", "#ee9b00", "#ca6702", "#bb3e03", "#ae2012", "#9b2226","#A09ABC"





## ME 33

#import data
ME_33_rel_abund_dat <- read.csv("~/Library/CloudStorage/GoogleDrive-sbirch1@charlotte.edu/My Drive/Postdoc_work/Mesocosom_exps/Mesocosom_Summer_2022/Transcriptomics/EdgeR_MAPPED_4-29-24/WGCNA_Seqs/ME_33_rel_abund_dat.csv", stringsAsFactors=TRUE)
str(ME_33_rel_abund_dat)

#plot counts
ggplot(ME_33_rel_abund_dat, aes(x=ME_33_Count, y=Name, fill = Name)) + geom_bar(stat= "identity") +theme_minimal() + scale_fill_manual(values=c("#001219", "#005f73", "#0a9396", "#94d2bd", "#e9d8a6", "#ee9b00", "#ca6702", "#bb3e03", "#ae2012", "#9b2226","#A09ABC"))
ggplot(ME_33_rel_abund_dat, aes(x=ME_33_Count, y=Module, fill = Name)) + geom_bar(stat= "identity") +theme_minimal() + scale_fill_manual(values=c("#001219", "#005f73", "#0a9396", "#94d2bd", "#e9d8a6", "#ee9b00", "#ca6702", "#bb3e03", "#ae2012", "#9b2226","#A09ABC"))

ME_33_rel_abund_dat %>%
  arrange(ME_33_Count) %>%
  mutate(Name=factor(Name, levels=c("Protein families: genetic information processing","Signal transduction","Protein families: signaling and cellular processes","Neurodegenerative disease","Immune system","Infectious disease: viral","Cancer: overview","Infectious disease: bacterial","Cancer: specific types","Endocrine system","Other Functions(35)"))) %>%
  ggplot( aes(x=ME_33_Count, y=Name, fill = Name)) + geom_bar(stat= "identity") +theme_minimal() + scale_fill_manual(values=c("#001219", "#005f73", "#0a9396", "#94d2bd", "#e9d8a6", "#ee9b00", "#ca6702","#bb3e03", "#ae2012", "#9b2226","#A09ABC"))
#Rplot-ME_33-top_10_counts-indiv

ME_33_rel_abund_dat %>%
  arrange(ME_33_Count) %>%
  mutate(Name=factor(Name, levels=c("Other Functions(35)","Endocrine system","Cancer: specific types","Infectious disease: bacterial","Cancer: overview","Infectious disease: viral","Immune system","Neurodegenerative disease","Protein families: signaling and cellular processes","Signal transduction","Protein families: genetic information processing"))) %>%
  ggplot( aes(x=ME_33_Count, y=Name, fill = Name)) + geom_bar(stat= "identity") +theme_minimal() + scale_fill_manual(values=c("#A09ABC","#001219", "#005f73", "#0a9396", "#94d2bd", "#e9d8a6", "#ee9b00", "#ca6702","#bb3e03", "#ae2012", "#9b2226"))
#Rplot-ME_33-top_10_counts-indiv-high-low



#plot rel abundance
ggplot(ME_33_rel_abund_dat, aes(x=rel_abund, y=Name, fill = Name)) + geom_bar(stat= "identity") +theme_minimal() + scale_fill_manual(values=c("#001219", "#005f73", "#0a9396", "#94d2bd", "#e9d8a6", "#ee9b00", "#ca6702", "#bb3e03", "#ae2012", "#9b2226","#A09ABC"))
ggplot(ME_33_rel_abund_dat, aes(x=rel_abund, y=Module, fill = Name)) + geom_bar(stat= "identity") +theme_minimal() + scale_fill_manual(values=c("#001219", "#005f73", "#0a9396", "#94d2bd", "#e9d8a6", "#ee9b00", "#ca6702", "#A09ABC","#bb3e03", "#ae2012", "#9b2226"))

ME_33_rel_abund_dat %>%
  arrange(rel_abund) %>%
  mutate(Name=factor(Name, levels=c("Protein families: genetic information processing","Signal transduction","Protein families: signaling and cellular processes","Neurodegenerative disease","Immune system","Infectious disease: viral","Cancer: overview","Infectious disease: bacterial","Cancer: specific types","Endocrine system","Other Functions(35)"))) %>%
  ggplot( aes(x=rel_abund, y=Name, fill = Name)) + geom_bar(stat= "identity") +theme_minimal() + scale_fill_manual(values=c("#001219", "#005f73", "#0a9396", "#94d2bd", "#e9d8a6", "#ee9b00", "#ca6702","#bb3e03", "#ae2012", "#9b2226","#A09ABC"))
#Rplot-ME_33-top_10_rel_abund-indiv


#plot stacked bart plot rel abundance
ME_33_rel_abund_dat %>%
  arrange(rel_abund) %>%
  mutate(Name=factor(Name, levels=c("Protein families: genetic information processing","Signal transduction","Protein families: signaling and cellular processes","Neurodegenerative disease","Immune system","Infectious disease: viral","Cancer: overview","Infectious disease: bacterial","Cancer: specific types","Endocrine system","Other Functions(35)"))) %>%
  ggplot( aes(x=rel_abund, y=Module, fill = Name)) + geom_bar(stat= "identity") +theme_minimal() + scale_fill_manual(values=c("#001219", "#005f73", "#0a9396", "#94d2bd", "#e9d8a6", "#ee9b00", "#ca6702","#bb3e03", "#ae2012", "#9b2226","#A09ABC"))
#Rplot-ME_33-top_10_rel_abund_stacked

#no legend
ME_33_rel_abund_dat %>%
  arrange(rel_abund) %>%
  mutate(Name=factor(Name, levels=c("Protein families: genetic information processing","Signal transduction","Protein families: signaling and cellular processes","Neurodegenerative disease","Immune system","Infectious disease: viral","Cancer: overview","Infectious disease: bacterial","Cancer: specific types","Endocrine system","Other Functions(35)"))) %>%
  ggplot( aes(x=rel_abund, y=Module, fill = Name)) + geom_bar(stat= "identity", show.legend = FALSE) +theme_minimal() + scale_fill_manual(values=c("#001219", "#005f73", "#0a9396", "#94d2bd", "#e9d8a6", "#ee9b00", "#ca6702","#bb3e03", "#ae2012", "#9b2226","#A09ABC"))
#Rplot-ME_33-top_10_rel_abund_stacked_NO_LEGEND





### FOAM LEVEL 1 16 total 

### Individual barchart - rel abundance
"#001219", "#005f73", "#0a9396", "#94d2bd", "#e9d8a6", "#ee9b00", "#ca6702", "#bb3e03", "#ae2012", "#9b2226","#A09ABC"

'#ecc8af', '#e7ad99', '#ce796b', '#c18c5d', '#495867','#006989', '#4281a4', '#274c77', '#697a21',"#A7E2E3" ,'#475841','#504746',"#ae2012","#4F345A","#C33149","#0B3948"



## ME 33

#import data
ME_33_foam_dat <- read.delim("~/Library/CloudStorage/GoogleDrive-sbirch1@charlotte.edu/My Drive/Postdoc_work/Mesocosom_exps/Mesocosom_Summer_2022/Transcriptomics/EdgeR_MAPPED_4-29-24/WGCNA_Seqs/ME_33_KOFam_eukaryote_FOAM_level-1.tsv", stringsAsFactors=TRUE)
str(ME_33_foam_dat)

#plot counts

ME_33_foam_dat %>%
  arrange(ME_33_Count) %>%
  mutate(Name=factor(Name, levels=c("01_Fermentation","02_Homoacetogenesis","05_Fatty acid oxidation","06_Amino acid utilization biosynthesis metabolism","07_Nucleic acid metabolism","09_Carbohydrate Active enzyme - CAZy","15_Methylotrophy","16_Embden Meyerhof-Parnos (EMP)","17_Gluconeogenesis","19_Saccharide and derivated synthesis","20_Hydrolysis of polymers","21_Cellular response to stress"))) %>%
  ggplot( aes(x=ME_33_Count, y=Name, fill = Name)) + geom_bar(stat= "identity") +theme_minimal() + scale_fill_manual(values=c("#0B3948","#C33149","#ae2012",'#504746','#475841','#697a21','#006989','#495867','#c18c5d','#ce796b','#e7ad99','#ecc8af'))
#Rplot-FOAM_ME_33-indiv

ME_33_foam_dat %>%
  arrange(ME_33_Count) %>%
  mutate(Name=factor(Name, levels=c("21_Cellular response to stress","20_Hydrolysis of polymers","19_Saccharide and derivated synthesis","17_Gluconeogenesis","16_Embden Meyerhof-Parnos (EMP)","15_Methylotrophy","09_Carbohydrate Active enzyme - CAZy","07_Nucleic acid metabolism","06_Amino acid utilization biosynthesis metabolism","05_Fatty acid oxidation","02_Homoacetogenesis","01_Fermentation"))) %>%
  ggplot( aes(x=ME_33_Count, y=Name, fill = Name)) + geom_bar(stat= "identity") +theme_minimal() + scale_fill_manual(values=c('#ecc8af','#e7ad99','#ce796b','#c18c5d','#495867','#006989','#697a21','#475841','#504746',"#ae2012","#C33149","#0B3948"))
#Rplot-FOAM_ME_33-indiv-high-low



#plot rel abundance

ME_33_foam_dat %>%
  arrange(rel_abund) %>%
  mutate(Name=factor(Name, levels=c("01_Fermentation","02_Homoacetogenesis","05_Fatty acid oxidation","06_Amino acid utilization biosynthesis metabolism","07_Nucleic acid metabolism","09_Carbohydrate Active enzyme - CAZy","15_Methylotrophy","16_Embden Meyerhof-Parnos (EMP)","17_Gluconeogenesis","19_Saccharide and derivated synthesis","20_Hydrolysis of polymers","21_Cellular response to stress"))) %>%
  ggplot( aes(x=rel_abund, y=Name, fill = Name)) + geom_bar(stat= "identity") +theme_minimal() + scale_fill_manual(values=c("#0B3948","#C33149","#ae2012",'#504746','#475841','#697a21','#006989','#495867','#c18c5d','#ce796b','#e7ad99','#ecc8af'))
#Rplot-FOAM_ME_33-top_10_rel_abund-indiv


#plot stacked bart plot rel abundance
ME_33_foam_dat %>%
  arrange(rel_abund) %>%
  mutate(Name=factor(Name, levels=c("21_Cellular response to stress","20_Hydrolysis of polymers","19_Saccharide and derivated synthesis","17_Gluconeogenesis","16_Embden Meyerhof-Parnos (EMP)","15_Methylotrophy","09_Carbohydrate Active enzyme - CAZy","07_Nucleic acid metabolism","06_Amino acid utilization biosynthesis metabolism","05_Fatty acid oxidation","02_Homoacetogenesis","01_Fermentation"))) %>%
  ggplot( aes(x=rel_abund, y=Name, fill = Name)) + geom_bar(stat= "identity") +theme_minimal() + scale_fill_manual(values=c('#ecc8af','#e7ad99','#ce796b','#c18c5d','#495867','#006989','#697a21','#475841','#504746',"#ae2012","#C33149","#0B3948"))
#Rplot-FOAM_ME_33-top_10_rel_abund_stacked



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



#### North vs South #### 



#### Make stacked bar plots KEGG level 2 ####

library(ggplot2)
library(dplyr)
library(ggsignif)

setwd("~/Library/CloudStorage/GoogleDrive-sbirch1@charlotte.edu/My Drive/Postdoc_work/Mesocosom_exps/Mesocosom_Summer_2022/Transcriptomics/EdgeR_MAPPED_4-29-24/WGCNA/WGCNA_North-vs-South_outputs")

#Stacked barplot Kegg Level-2 combined
kegg_combined <- read.csv("~/Library/CloudStorage/GoogleDrive-sbirch1@charlotte.edu/My Drive/Postdoc_work/Mesocosom_exps/Mesocosom_Summer_2022/Transcriptomics/EdgeR_MAPPED_4-29-24/WGCNA/WGCNA_North-vs-South_outputs/north_south_kegg_subset.csv", stringsAsFactors=TRUE)
str(kegg_combined)

#With Legend (12 colors needed)
kegg_combined %>%
  arrange(Count) %>%
  mutate(Module=factor(Module, levels=c("ME_35","ME_15","ME_5","ME_55","ME_32"))) %>%
  ggplot( aes(x=Count, y=Module, fill = Type)) + geom_bar(stat= "identity", color = "black") +theme_minimal() + scale_fill_manual(values=c("#127475","#8cb369","#f4e285","#f4a259","#5b8e7d","#bc4b51","#2d2d2a","#0c4767","#820933","#34623f","#546a76","#562c2c",'#ce796b'))

#no legend
kegg_combined %>%
  arrange(Count) %>%
  mutate(Module=factor(Module, levels=c("ME_35","ME_15","ME_5","ME_55","ME_32"))) %>%
  ggplot( aes(x=Count, y=Module, fill = Type)) + geom_bar(stat= "identity", show.legend = FALSE,position="fill") +theme_minimal() + scale_fill_manual(values=c("#127475","#8cb369","#f4e285","#f4a259","#5b8e7d","#bc4b51","#2d2d2a","#0c4767","#820933","#34623f","#546a76","#562c2c",'#ce796b'))

kegg_combined %>%
  arrange(Count) %>%
  mutate(Module=factor(Module, levels=c("ME_35","ME_15","ME_5","ME_55","ME_32"))) %>%
  ggplot( aes(x=Count, y=Module, fill = Type)) + geom_bar(stat= "identity", show.legend = FALSE) +theme_minimal() + scale_fill_manual(values=c("#127475","#8cb369","#f4e285","#f4a259","#5b8e7d","#bc4b51","#2d2d2a","#0c4767","#820933","#34623f","#546a76","#562c2c",'#ce796b'))




## Plot individual modules  

## ME 35

#import data
ME_35 <- read.csv("~/Library/CloudStorage/GoogleDrive-sbirch1@charlotte.edu/My Drive/Postdoc_work/Mesocosom_exps/Mesocosom_Summer_2022/Transcriptomics/EdgeR_MAPPED_4-29-24/WGCNA/WGCNA_North-vs-South_outputs/ME_35.csv", stringsAsFactors=TRUE)
str(ME_35)

#plot counts
ggplot(ME_35, aes(x=ME_35_Count, y=Name, fill = Name)) + geom_bar(stat= "identity") +theme_minimal() + scale_fill_manual(values=c("#001219", "#005f73", "#0a9396", "#94d2bd", "#e9d8a6", "#ee9b00", "#ca6702", "#bb3e03", "#ae2012", "#9b2226","#A09ABC"))
ggplot(ME_35, aes(x=ME_35_Count, y=Module, fill = Name)) + geom_bar(stat= "identity") +theme_minimal() + scale_fill_manual(values=c("#001219", "#005f73", "#0a9396", "#94d2bd", "#e9d8a6", "#ee9b00", "#ca6702", "#bb3e03", "#ae2012", "#9b2226","#A09ABC"))

ME_35 %>%
  arrange(ME_35_Count) %>%
  mutate(Name=factor(Name, levels=c("Protein families: genetic information processing","Protein families: signaling and cellular processes","Protein families: metabolism","Signal transduction","Cancer: overview","Signaling molecules and interaction","Infectious disease: bacterial","Infectious disease: viral","Transport and catabolism","Cardiovascular disease","Other Functions (22)"))) %>%
  ggplot( aes(x=ME_35_Count, y=Name, fill = Name)) + geom_bar(stat= "identity") +theme_minimal() + scale_fill_manual(values=c("#001219", "#005f73", "#0a9396", "#94d2bd", "#e9d8a6", "#ee9b00", "#ca6702","#bb3e03", "#ae2012", "#9b2226","#A09ABC"))
#Rplot-ME_21-top_10_counts-indiv

ME_35 %>%
  arrange(ME_35_Count) %>%
  mutate(Name=factor(Name, levels=c("Other Functions (22)","Cardiovascular disease","Transport and catabolism","Infectious disease: viral","Infectious disease: bacterial","Signaling molecules and interaction","Cancer: overview","Signal transduction","Protein families: metabolism","Protein families: signaling and cellular processes","Protein families: genetic information processing"))) %>%
  ggplot( aes(x=ME_35_Count, y=Name, fill = Name)) + geom_bar(stat= "identity") +theme_minimal() + scale_fill_manual(values=c("#A09ABC","#001219", "#005f73", "#0a9396", "#94d2bd", "#e9d8a6", "#ee9b00", "#ca6702","#bb3e03", "#ae2012", "#9b2226"))
#Rplot-ME_35-top_10_counts-indiv-high-low

ME_35 %>%
  arrange(ME_35_Count) %>%
  mutate(Name=factor(Name, levels=c("Other Functions (22)","Cardiovascular disease","Transport and catabolism","Infectious disease: viral","Infectious disease: bacterial","Signaling molecules and interaction","Cancer: overview","Signal transduction","Protein families: metabolism","Protein families: signaling and cellular processes","Protein families: genetic information processing"))) %>%
  ggplot( aes(x=ME_35_Count, y=Name, fill = Name)) + geom_bar(stat= "identity",show.legend = FALSE) +theme_minimal() + scale_fill_manual(values=c("#A09ABC","#001219", "#005f73", "#0a9396", "#94d2bd", "#e9d8a6", "#ee9b00", "#ca6702","#bb3e03", "#ae2012", "#9b2226"))
#Rplot-ME_35-top_10_counts-indiv-high-low-no_key


