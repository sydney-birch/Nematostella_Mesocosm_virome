### 16s QIIME2 Visualizations 

## Ran the analysis on the HPC - making figures here


#Following: https://forum.qiime2.org/t/tutorial-integrating-qiime2-and-r-for-data-visualization-and-analysis-using-qiime2r/4121

#Set working directory: 
setwd("~/Library/CloudStorage/GoogleDrive-sbirch1@charlotte.edu/My Drive/Postdoc_work/Mesocosom_exps/Mesocosom_Summer_2022/16s seq/QIIME Output/R")

#Load libraries 
library(ape)
library(Biostrings)
library(biomformat)
library(phyloseq)
library(Hmisc)
library(yaml)
library(tidyr)
library(dplyr)
library(stats)
library(utils)
library(ggplot2)
library(qiime2R)
library(tidyverse)

#if (!require("BiocManager", quietly = TRUE))
#  install.packages("BiocManager")
#BiocManager::install("biomformat")
#BiocManager::install("phyloseq")

if (!requireNamespace("devtools", quietly = TRUE)){install.packages("devtools")}
devtools::install_github("jbisanz/qiime2R")



#read in table of sequence variants (SVs)
SVs<-read_qza("table.qza")

#read in the metadata file
metadata<-read_q2metadata("metadata-file-2022_v2.tsv")

#read in taxonomy
taxonomy<-read_qza("taxonomy_v2.qza")
head(taxonomy$data)

#create a phyloseq object
physeq<-qza_to_phyloseq(
  features="table.qza",
  tree="rooted-tree.qza","taxonomy_v2.qza",
  metadata = "metadata-file-2022_v2.tsv"
)
physeq


# Alpha diversity over time

#Shannon
metadata<-read_q2metadata("metadata-file-2022_v2.tsv")
shannon<-read_qza("shannon_vector.qza")

shannon<-shannon$data %>% rownames_to_column("sample-id") # this moves the sample names to a new column that matches the metadata and allows them to be merged
#rename sample id column in shannon file
names(shannon)[1]<-"SampleID"
#see how many samples shared between shannon data and metadata
gplots::venn(list(metadata=metadata$SampleID, shannon=shannon$SampleID))


#add shannon data to metadata file
metadata2<-
  metadata2 %>% 
  left_join(shannon,by="SampleID")
head(metadata2)

#remove blank
metadata3<-metadata2[!(row.names(metadata2) %in% c("69")),]

#plot shannon diversity
metadata2 %>%
  filter(!is.na(shannon_entropy)) %>%
  ggplot(aes(x=Timepoint, y=shannon_entropy, color=Location)) +
  geom_boxplot() + geom_jitter(shape=16, size= 3, position = position_jitter(0.2)) + lims(y = c(0,35)) +
  xlab("Time Point") +
  ylab("Shannon Diversity") +
  theme_q2r() + # try other themes like theme_bw() or theme_classic()
  scale_color_viridis_d(name="Location") # use different color scale which is color blind friendly
#ggsave("Shannon_by_time.pdf", height=3, width=4, device="pdf") # save a PDF 3 inches by 4 inches

metadata3 %>%
  filter(!is.na(shannon_entropy)) %>%
  arrange(shannon_entropy) %>%
  mutate(loctime=factor(loctime, levels=c("NS-T0","NS-T14","ME-T0","ME-T14","NH-T0","NH-T14","Field-T0","Field-T14","MA-T0","MA-T14","SC-T0","SC-T14","FL-T0","FL-T14"))) %>%
  ggplot(aes(x=loctime, y=shannon_entropy, color=loctime)) +
  geom_boxplot() + geom_jitter(shape=16, size= 2, position = position_jitter(0.2))  +
  xlab("Time Point") +
  ylab("Shannon Diversity") +
  theme_q2r() + # try other themes like theme_bw() or theme_classic()
  scale_color_viridis_d(name="Location") # use different color scale which is color blind friendly

#I like this one - box plot
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

metadata3 %>%
  filter(!is.na(shannon_entropy)) %>%
  arrange(shannon_entropy) %>%
  mutate(ns_time=factor(ns_time, levels=c("north_T0","north_T14","control_T0","control_T14","south_T0","south_T14"))) %>%
  ggplot(aes(x=ns_time, y=shannon_entropy, color=Timepoint)) +
  geom_boxplot() + geom_jitter(shape=16, size= 2, position = position_jitter(0.2))  +
  xlab("North vs South") +
  ylab("Shannon Diversity") +
  theme_bw() + # try other themes like theme_bw() or theme_classic()
  scale_color_manual(values=c("#669BBC", "#003049"))

metadata3 %>%
  filter(!is.na(shannon_entropy)) %>%
  arrange(shannon_entropy) %>%
  mutate(ns_time=factor(ns_time, levels=c("north_T0","north_T14","control_T0","control_T14","south_T0","south_T14"))) %>%
  ggplot(aes(x=ns_time, y=shannon_entropy, color=Timepoint)) +
  geom_boxplot() + geom_jitter(shape=16, size= 2, position = position_jitter(0.2))  +
  xlab("North vs South") +
  ylab("Shannon Diversity") +
  theme_bw() + # try other themes like theme_bw() or theme_classic()
  scale_color_manual(values=c("#669BBC", "#003049"))

#multipanneled barplot for each location divided by timepoint - I like this one
metadata3 %>%
  filter(!is.na(shannon_entropy)) %>%
  arrange(shannon_entropy) %>%
  mutate(Location=factor(Location, levels=c("NS","ME","NH","Field","MA","SC","FL"))) %>%
  ggplot(aes(x=Timepoint, y=shannon_entropy, fill=Timepoint)) +
  stat_summary(geom="bar", fun.data=mean_se, color="black") + #here black is the outline for the bars
  geom_jitter(shape=21, width=0.2, height=0) +
  coord_cartesian(ylim=c(4,8)) + # adjust y-axis
  facet_grid(~Location) + # create a panel for each body site
  xlab("Timepoint") +
  ylab("Shannon Diversity") +
  theme_q2r() +
  scale_fill_manual(values=c("#669BBC","#003049")) + #specify custom colors
  theme(legend.position="none") #remove the legend as it isn't needed

#Tukeys:
#subset data
sub_T0<- subset(metadata3, metadata3$Timepoint  == "T0")
sub_T14<- subset(metadata3, metadata3$Timepoint  == "T14")

sub_NS<- subset(metadata3, metadata3$Location  == "NS")
sub_ME<- subset(metadata3, metadata3$Location  == "ME")
sub_NH<- subset(metadata3, metadata3$Location  == "NH")
sub_MA<- subset(metadata3, metadata3$Location  == "MA")
sub_SC<- subset(metadata3, metadata3$Location  == "SC")
sub_FL<- subset(metadata3, metadata3$Location  == "FL")

#tukey 
full_mod<-lm(shannon_entropy ~ loctime,metadata3)
anova(full_mod)
#          Df Sum Sq Mean Sq F value    Pr(>F)    
#loctime   13 16.357 1.25826  3.4569 0.0006392 ***
#Residuals 54 19.655 0.36398 
boxplot(shannon_entropy ~ loctime, data = metadata3,ylab="Richness")

full_mod<-lm(shannon_entropy ~ ns_time,metadata3)
anova(full_mod)
#ns_time    5 14.505 2.90107  8.3631 4.33e-06 ***
#Residuals 62 21.507 0.34689 
boxplot(shannon_entropy ~ ns_time, data = metadata3,ylab="Richness")

T0_mod<-lm(shannon_entropy ~ Location,sub_T0)
anova(T0_mod)
T14_mod<-lm(shannon_entropy ~ Location,sub_T14)
anova(T14_mod)

NS_mod<-lm(shannon_entropy ~ Timepoint,sub_NS)
anova(NS_mod)
ME_mod<-lm(shannon_entropy ~ Timepoint,sub_ME)
anova(ME_mod)
NH_mod<-lm(shannon_entropy ~ Timepoint,sub_NH)
anova(NH_mod)
MA_mod<-lm(shannon_entropy ~ Timepoint,sub_MA)
anova(MA_mod)
SC_mod<-lm(shannon_entropy ~ Timepoint,sub_SC)
anova(SC_mod)
FL_mod<-lm(shannon_entropy ~ Timepoint,sub_FL)
anova(FL_mod)

library(agricolae)
full_compar<-HSD.test(full_mod, "loctime")

full_compar<-HSD.test(full_mod, "ns_time")

T0_compar<-HSD.test(T0_mod, "Location")
T14_compar<-HSD.test(T14_mod, "Location")

NS_compar<-HSD.test(NS_mod, "Timepoint")
ME_compar<-HSD.test(ME_mod, "Timepoint")
NH_compar<-HSD.test(NH_mod, "Timepoint")
MA_compar<-HSD.test(MA_mod, "Timepoint")
SC_compar<-HSD.test(SC_mod, "Timepoint")
FL_compar<-HSD.test(FL_mod, "Timepoint")



#### faith pd 

metadata<-read_q2metadata("metadata-file-2022_v2.tsv")
faith<-read_qza("faith_pd_vector.qza")
faith$data$SampleID
#rename column names (currently V1 and V2)
names(faith$data)[names(faith$data)=="V1"] <- "SampleID" 
names(faith$data)[names(faith$data)=="V2"] <- "faith_pd" 

#faith<-faith$data %>% rownames_to_column("sample-id") # this moves the sample names to a new column that matches the metadata and allows them to be merged
#rename sample id column in shannon file
#names(faith)[1]<-"SampleID"
#see how many samples shared between shannon data and metadata
gplots::venn(list(metadata=metadata$SampleID, faith=faith$data$SampleID))

#remove samples that were exluded: 
#NH-T14-FIELD-B5 (40)
#NH-T14-FIELD-B1 (36)
#Blank-2-meso (72)

#remove from metadata 
metadata2<-metadata[!(row.names(metadata) %in% c("36","40","72")),]
gplots::venn(list(metadata=metadata2$SampleID, faith=faith$data$SampleID))

#add shannon data to metadata file
metadata2<-
  metadata2 %>% 
  left_join(faith$data,by="SampleID")
head(metadata2)

#remove blank
metadata3<-metadata2[!(row.names(metadata2) %in% c("69")),]

#I like this one - box plot
metadata3 %>%
  filter(!is.na(faith_pd)) %>%
  arrange(faith_pd) %>%
  mutate(Location=factor(Location, levels=c("NS","ME","NH","Field","MA","SC","FL"))) %>%
  ggplot(aes(x=Location, y=faith_pd, color=Timepoint)) +
  geom_boxplot() + geom_jitter(shape=16, size= 2, position = position_jitter(0.2))  +
  xlab("Location") +
  ylab("Faith PD") +
  theme_bw() + # try other themes like theme_bw() or theme_classic()
  scale_color_manual(values=c("#669BBC", "#003049"))

metadata3 %>%
  filter(!is.na(faith_pd)) %>%
  arrange(faith_pd) %>%
  mutate(Location=factor(Location, levels=c("NS","ME","NH","Field","MA","SC","FL"))) %>%
  ggplot(aes(x=Timepoint, y=faith_pd, fill=Timepoint)) +
  geom_boxplot() + geom_jitter(shape=16, size= 2, position = position_jitter(0.2))  +
  xlab("Location") +
  ylab("Faith PD") +
  facet_grid(~Location) +
  theme_bw() + # try other themes like theme_bw() or theme_classic()
  scale_fill_manual(values=c("#83C5BE","#006D77"))



#multipanneled barplot for each location divided by timepoint - I like this one
metadata3 %>%
  filter(!is.na(faith_pd)) %>%
  arrange(faith_pd) %>%
  mutate(Location=factor(Location, levels=c("NS","ME","NH","Field","MA","SC","FL"))) %>%
  ggplot(aes(x=Timepoint, y=faith_pd, fill=Timepoint)) +
  stat_summary(geom="bar", fun.data=mean_se, color="black") + #here black is the outline for the bars
  geom_jitter(shape=21, width=0.2, height=0) +
  coord_cartesian(ylim=c(50,250)) + # adjust y-axis
  facet_grid(~Location) + # create a panel for each body site
  xlab("Timepoint") +
  ylab("Faith PD") +
  theme_q2r() +
  scale_fill_manual(values=c("#83C5BE","#006D77")) + #specify custom colors
  theme(legend.position="none") #remove the legend as it isn't needed



#Tukeys:
#subset data
sub_T0<- subset(metadata3, metadata3$Timepoint  == "T0")
sub_T14<- subset(metadata3, metadata3$Timepoint  == "T14")

sub_NS<- subset(metadata3, metadata3$Location  == "NS")
sub_ME<- subset(metadata3, metadata3$Location  == "ME")
sub_NH<- subset(metadata3, metadata3$Location  == "NH")
sub_MA<- subset(metadata3, metadata3$Location  == "MA")
sub_SC<- subset(metadata3, metadata3$Location  == "SC")
sub_FL<- subset(metadata3, metadata3$Location  == "FL")

#tukey 
full_mod<-lm(faith_pd ~ loctime,metadata3)
anova(full_mod)

T0_mod<-lm(faith_pd ~ Location,sub_T0)
anova(T0_mod)
T14_mod<-lm(faith_pd ~ Location,sub_T14)
anova(T14_mod)

NS_mod<-lm(faith_pd ~ Timepoint,sub_NS)
anova(NS_mod)
ME_mod<-lm(faith_pd ~ Timepoint,sub_ME)
anova(ME_mod)
NH_mod<-lm(faith_pd ~ Timepoint,sub_NH)
anova(NH_mod)
MA_mod<-lm(faith_pd ~ Timepoint,sub_MA)
anova(MA_mod)
SC_mod<-lm(faith_pd ~ Timepoint,sub_SC)
anova(SC_mod)
FL_mod<-lm(faith_pd ~ Timepoint,sub_FL)
anova(FL_mod)

library(agricolae)
full_compar<-HSD.test(full_mod, "loctime")

T0_compar<-HSD.test(T0_mod, "Location")
T14_compar<-HSD.test(T14_mod, "Location")

NS_compar<-HSD.test(NS_mod, "Timepoint")
ME_compar<-HSD.test(ME_mod, "Timepoint")
NH_compar<-HSD.test(NH_mod, "Timepoint")
MA_compar<-HSD.test(MA_mod, "Timepoint")
SC_compar<-HSD.test(SC_mod, "Timepoint")
FL_compar<-HSD.test(FL_mod, "Timepoint")




#### Evenness

metadata<-read_q2metadata("metadata-file-2022_v2.tsv")
even<-read_qza("evenness_vector.qza")

even<-even$data %>% rownames_to_column("sample-id") # this moves the sample names to a new column that matches the metadata and allows them to be merged
#rename sample id column in shannon file
names(even)[1]<-"SampleID"
#see how many samples shared between shannon data and metadata
gplots::venn(list(metadata=metadata$SampleID, even=even$SampleID))


#remove samples that were exluded: 
#NH-T14-FIELD-B5 (40)
#NH-T14-FIELD-B1 (36)
#Blank-2-meso (72)

#remove from metadata 
metadata2<-metadata[!(row.names(metadata) %in% c("36","40","72")),]
gplots::venn(list(metadata=metadata2$SampleID, even=even$SampleID))

#add shannon data to metadata file
metadata2<-
  metadata2 %>% 
  left_join(even,by="SampleID")
head(metadata2)

#remove blank
metadata3<-metadata2[!(row.names(metadata2) %in% c("69")),]

#I like this one - box plot
metadata3 %>%
  filter(!is.na(pielou_evenness)) %>%
  arrange(pielou_evenness) %>%
  mutate(Location=factor(Location, levels=c("NS","ME","NH","Field","MA","SC","FL"))) %>%
  ggplot(aes(x=Location, y=pielou_evenness, color=Timepoint)) +
  geom_boxplot() + geom_jitter(shape=16, size= 2, position = position_jitter(0.2))  +
  xlab("Location") +
  ylab("Pielou Evenness") +
  theme_bw() + # try other themes like theme_bw() or theme_classic()
  scale_color_manual(values=c("#669BBC", "#003049"))

metadata3 %>%
  filter(!is.na(pielou_evenness)) %>%
  arrange(pielou_evenness) %>%
  mutate(Location=factor(Location, levels=c("NS","ME","NH","Field","MA","SC","FL"))) %>%
  ggplot(aes(x=Timepoint, y=pielou_evenness, color=Timepoint)) +
  geom_boxplot() + geom_jitter(shape=16, size= 2, position = position_jitter(0.2))  +
  xlab("Location") +
  ylab("Pielou Evenness") +
  facet_grid(~Location) +
  theme_bw() + # try other themes like theme_bw() or theme_classic()
  scale_color_manual(values=c("#669BBC", "#003049"))

#multipanneled barplot for each location divided by timepoint - I like this one
metadata3 %>%
  filter(!is.na(pielou_evenness)) %>%
  arrange(pielou_evenness) %>%
  mutate(Location=factor(Location, levels=c("NS","ME","NH","Field","MA","SC","FL"))) %>%
  ggplot(aes(x=Timepoint, y=pielou_evenness, fill=Timepoint)) +
  stat_summary(geom="bar", fun.data=mean_se, color="black") + #here black is the outline for the bars
  geom_jitter(shape=21, width=0.2, height=0) +
  coord_cartesian(ylim=c(0.55,0.75)) + # adjust y-axis
  facet_grid(~Location) + # create a panel for each body site
  xlab("Timepoint") +
  ylab("Pielou Evenness") +
  theme_q2r() +
  scale_fill_manual(values=c("#669BBC","#003049")) + #specify custom colors
  theme(legend.position="none") #remove the legend as it isn't needed

#Tukeys:
#subset data
sub_T0<- subset(metadata3, metadata3$Timepoint  == "T0")
sub_T14<- subset(metadata3, metadata3$Timepoint  == "T14")

sub_NS<- subset(metadata3, metadata3$Location  == "NS")
sub_ME<- subset(metadata3, metadata3$Location  == "ME")
sub_NH<- subset(metadata3, metadata3$Location  == "NH")
sub_MA<- subset(metadata3, metadata3$Location  == "MA")
sub_SC<- subset(metadata3, metadata3$Location  == "SC")
sub_FL<- subset(metadata3, metadata3$Location  == "FL")

#tukey 
full_mod<-lm(pielou_evenness ~ loctime,metadata3)
anova(full_mod)

T0_mod<-lm(pielou_evenness ~ Location,sub_T0)
anova(T0_mod)
T14_mod<-lm(pielou_evenness ~ Location,sub_T14)
anova(T14_mod)

NS_mod<-lm(pielou_evenness ~ Timepoint,sub_NS)
anova(NS_mod)
ME_mod<-lm(pielou_evenness ~ Timepoint,sub_ME)
anova(ME_mod)
NH_mod<-lm(pielou_evenness ~ Timepoint,sub_NH)
anova(NH_mod)
MA_mod<-lm(pielou_evenness ~ Timepoint,sub_MA)
anova(MA_mod)
SC_mod<-lm(pielou_evenness ~ Timepoint,sub_SC)
anova(SC_mod)
FL_mod<-lm(pielou_evenness ~ Timepoint,sub_FL)
anova(FL_mod)

library(agricolae)
full_compar<-HSD.test(full_mod, "loctime")

T0_compar<-HSD.test(T0_mod, "Location")
T14_compar<-HSD.test(T14_mod, "Location")

NS_compar<-HSD.test(NS_mod, "Timepoint")
ME_compar<-HSD.test(ME_mod, "Timepoint")
NH_compar<-HSD.test(NH_mod, "Timepoint")
MA_compar<-HSD.test(MA_mod, "Timepoint")
SC_compar<-HSD.test(SC_mod, "Timepoint")
FL_compar<-HSD.test(FL_mod, "Timepoint")





### Beta Diversity - PCoA ###
library(tidyverse)
library(qiime2R)

## Unifrac ##

#read in files
metadata<-read_q2metadata("metadata-file-2022_v2.tsv")
uwunifrac<-read_qza("unweighted_unifrac_pcoa_results.qza")
shannon<-read_qza("shannon_vector.qza")$data %>% rownames_to_column("SampleID") 


uwunifrac$data$Vectors %>%
  select(SampleID, PC1, PC2) %>%
  left_join(metadata,by="SampleID") %>%
  left_join(shannon,by="SampleID") %>%
  ggplot(aes(x=PC1, y=PC2, color=Location, shape=Timepoint, size=shannon_entropy)) +
  geom_point(alpha=0.5) + #alpha controls transparency and helps when points are overlapping
  theme_q2r() +
  scale_shape_manual(values=c(4,16,1), name="Timepoint") + #see http://www.sthda.com/sthda/RDoc/figure/graphs/r-plot-pch-symbols-points-in-r.png for numeric shape codes
  scale_size_continuous(name="Shannon Diversity") +
  scale_color_manual(values= c("#001219", "#005f73", "#0a9396","#e9d8a6", "#ee9b00","#bb3e03", "#ae2012", "#9b2226"))



uwunifrac$data$Vectors %>%
  select(SampleID, PC1, PC2) %>%
  left_join(metadata,by="SampleID") %>%
  left_join(shannon,by="SampleID") %>%
  arrange(shannon_entropy) %>%
  mutate(Location=factor(Location, levels=c("NS","ME","NH","Field","MA","SC","FL"))) %>%
  ggplot(aes(x=PC1, y=PC2, color=Timepoint, shape=Location, size=shannon_entropy)) +
  geom_point(alpha=0.9) + #alpha controls transparency and helps when points are overlapping
  theme_q2r() +
  scale_shape_manual(values=c(4,16,1,15,17,18,2,5), name="Location") + #see http://www.sthda.com/sthda/RDoc/figure/graphs/r-plot-pch-symbols-points-in-r.png for numeric shape codes
  scale_size_continuous(name="Shannon Diversity") +
  scale_color_manual(values= c("#e9d8a6","#94d2bd","#003049"))



## Jaccard ##

#read in files
metadata<-read_q2metadata("metadata-file-2022_v2.tsv")
jaccard<-read_qza("jaccard_pcoa_results.qza")
shannon<-read_qza("shannon_vector.qza")$data %>% rownames_to_column("SampleID") 


jaccard$data$Vectors %>%
  select(SampleID, PC1, PC2) %>%
  left_join(metadata,by="SampleID") %>%
  left_join(shannon,by="SampleID") %>%
  ggplot(aes(x=PC1, y=PC2, color=Location, shape=Timepoint, size=shannon_entropy)) +
  geom_point(alpha=0.5) + #alpha controls transparency and helps when points are overlapping
  theme_q2r() +
  scale_shape_manual(values=c(4,16,1), name="Timepoint") + #see http://www.sthda.com/sthda/RDoc/figure/graphs/r-plot-pch-symbols-points-in-r.png for numeric shape codes
  scale_size_continuous(name="Shannon Diversity") +
  scale_fill_manual(values= c("#001219", "#005f73", "#0a9396","#e9d8a6", "#ee9b00","#bb3e03", "#ae2012", "#9b2226"))


jaccard$data$Vectors %>%
  select(SampleID, PC1, PC2) %>%
  left_join(metadata,by="SampleID") %>%
  left_join(shannon,by="SampleID") %>%
  arrange(shannon_entropy) %>%
  mutate(Location=factor(Location, levels=c("NS","ME","NH","Field","MA","SC","FL"))) %>%
  ggplot(aes(x=PC1, y=PC2, color=Timepoint, shape=Location, size=shannon_entropy)) +
  geom_point(alpha=0.8) + #alpha controls transparency and helps when points are overlapping
  theme_q2r() +
  scale_shape_manual(values=c(4,16,1,15,17,18,2,5), name="Location") + #see http://www.sthda.com/sthda/RDoc/figure/graphs/r-plot-pch-symbols-points-in-r.png for numeric shape codes
  scale_size_continuous(name="Shannon Diversity") +
  scale_color_manual(values= c("#e9d8a6","#94d2bd","#003049"))



## bray curtis ##

#read in files
metadata<-read_q2metadata("metadata-file-2022_v2.tsv")
bray<-read_qza("bray_curtis_pcoa_results.qza")
shannon<-read_qza("shannon_vector.qza")$data %>% rownames_to_column("SampleID") 


bray$data$Vectors %>%
  select(SampleID, PC1, PC2) %>%
  left_join(metadata,by="SampleID") %>%
  left_join(shannon,by="SampleID") %>%
  arrange(shannon_entropy) %>%
  mutate(Location=factor(Location, levels=c("NS","ME","NH","Field","MA","SC","FL"))) %>%
  ggplot(aes(x=PC1, y=PC2, color=Timepoint, shape=Location, size=shannon_entropy)) +
  geom_point(alpha=0.8) + #alpha controls transparency and helps when points are overlapping
  theme_q2r() +
  scale_shape_manual(values=c(4,16,1,15,17,18,2,5), name="Location") + #see http://www.sthda.com/sthda/RDoc/figure/graphs/r-plot-pch-symbols-points-in-r.png for numeric shape codes
  scale_size_continuous(name="Shannon Diversity") +
  scale_color_manual(values= c("#e9d8a6","#94d2bd","#003049"))




### Heatmap of top 30 most abudnant taxonomic features ###

library(tidyverse)
library(qiime2R)

#read in data
metadata<-read_q2metadata("metadata-file-2022_v2.tsv")
SVs<-read_qza("table.qza")$data
taxonomy<-read_qza("taxonomy_v2.qza")$data

#subset data frame to remove all unassigned
taxonomy2<- subset(taxonomy, taxonomy$Taxon !="Unassigned")

SVs<-apply(SVs, 2, function(x) x/sum(x)*100) #convert to percent



SVsToPlot<-  
  data.frame(MeanAbundance=rowMeans(SVs)) %>% #find the average abundance of a SV
  rownames_to_column("Feature.ID") %>%
  arrange(desc(MeanAbundance)) %>%
  top_n(50, MeanAbundance) %>%
  pull(Feature.ID) #extract only the names from the table


SVsToPlot2<-  
  data.frame(MeanAbundance=rowMeans(SVs)) %>% #find the average abundance of a SV
  rownames_to_column("Feature.ID") %>%
  pull(Feature.ID)

write.csv(SVs, file= "SVs.csv")

#Locations
SVs %>%
  as.data.frame() %>%
  rownames_to_column("Feature.ID") %>%
  gather(-Feature.ID, key="SampleID", value="Abundance") %>%
  mutate(Feature.ID=if_else(Feature.ID %in% SVsToPlot,  Feature.ID, "Remainder")) %>% #flag features to be collapsed
  group_by(SampleID, Feature.ID) %>%
  summarize(Abundance=sum(Abundance), by=NULL) %>%
  left_join(metadata) %>%
  mutate(NormAbundance=log10(Abundance+0.01)) %>% # do a log10 transformation after adding a 0.01% pseudocount. Could also add 1 read before transformation to percent
  left_join(taxonomy) %>%
  mutate(Feature=paste(Feature.ID, Taxon)) %>%
  mutate(Feature=gsub("[kpcofgs]__", "", Feature)) %>% # trim out leading text from taxonomy string
  arrange(Feature) %>%
  mutate(Location=factor(Location, levels=c("NS","ME","NH","Field","MA","SC","FL","blank"))) %>%
  ggplot(aes(x=SampleID, y=Feature, fill=NormAbundance)) +
  geom_tile() +
  facet_grid(~Location, scales="free_x") +
  theme_q2r() +
  theme(axis.text.x=element_text(angle=45, hjust=1)) +
  scale_fill_viridis_c(name="log10(% Abundance)")

SVs %>%
  as.data.frame() %>%
  rownames_to_column("Feature.ID") %>%
  gather(-Feature.ID, key="SampleID", value="Abundance") %>%
  mutate(Feature.ID=if_else(Feature.ID %in% SVsToPlot,  Feature.ID, "Remainder")) %>% #flag features to be collapsed
  group_by(SampleID, Feature.ID) %>%
  summarize(Abundance=sum(Abundance), by=NULL) %>%
  left_join(metadata) %>%
  mutate(NormAbundance=log10(Abundance+0.01)) %>% # do a log10 transformation after adding a 0.01% pseudocount. Could also add 1 read before transformation to percent
  left_join(taxonomy) %>%
  mutate(Feature=paste( Taxon)) %>%
  mutate(Feature=gsub("[kpcofgs]__", "", Feature)) %>% # trim out leading text from taxonomy string
  arrange(Feature) %>%
  mutate(Location=factor(Location, levels=c("NS","ME","NH","Field","MA","SC","FL","blank"))) %>%
  ggplot(aes(x=SampleID, y=Feature, fill=NormAbundance)) +
  geom_tile() +
  facet_grid(~Location, scales="free_x") +
  theme_q2r() +
  theme(axis.text.x=element_text(angle=45, hjust=1)) +
  scale_fill_viridis_c(name="log10(% Abundance)")


parse_taxonomy(taxonomy)
write.csv(parse_taxonomy(taxonomy), file= "taxonomy.file2.csv")
#timepoints
SVs %>%
  as.data.frame() %>%
  rownames_to_column("Feature.ID") %>%
  gather(-Feature.ID, key="SampleID", value="Abundance") %>%
  mutate(Feature.ID=if_else(Feature.ID %in% SVsToPlot,  Feature.ID, "Remainder")) %>% #flag features to be collapsed
  group_by(SampleID, Feature.ID) %>%
  summarize(Abundance=sum(Abundance)) %>%
  left_join(metadata) %>%
  mutate(NormAbundance=log10(Abundance+0.01)) %>% # do a log10 transformation after adding a 0.01% pseudocount. Could also add 1 read before transformation to percent
  left_join(taxonomy) %>%
  mutate(Feature=paste(Feature.ID, Taxon)) %>%
  mutate(Feature=gsub("[kpcofgs]__", "", Feature)) %>% # trim out leading text from taxonomy string
  arrange(Feature) %>%
  mutate(Timepoint=factor(Timepoint, levels=c("T0","T14","blank"))) %>%
  ggplot(aes(x=SampleID, y=Feature, fill=NormAbundance)) +
  geom_tile() +
  facet_grid(~Timepoint, scales="free_x") +
  theme_q2r() +
  theme(axis.text.x=element_text(angle=45, hjust=1)) +
  scale_fill_viridis_c(name="log10(% Abundance)")



## Phylum Level (Level 2)

#import data
phy_16s_barplot <- read.csv("~/Library/CloudStorage/GoogleDrive-sbirch1@charlotte.edu/My Drive/Postdoc_work/Mesocosom_exps/Mesocosom_Summer_2022/16s seq/QIIME Output/R/phylum_stacked_barplot_dat.csv", stringsAsFactors=TRUE)


#With Legend
phy_16s_barplot %>%
  arrange(Count) %>%
  mutate(Location=factor(Location, levels=c("FL-T14","FL-T0","SC-T14","SC-T0","MA-T14","MA-T0","Field-T14","Field-T0","NH-T14","NH-T0","ME-T14","ME-T0","NS-T14","NS-T0"))) %>%
  ggplot( aes(x=Count, y=Location, fill = Taxa)) + geom_bar(stat= "identity", color = "black") +theme_minimal() 

phy_16s_barplot %>%
  arrange(Count) %>%
  mutate(Location=factor(Location, levels=c("FL-T14","FL-T0","SC-T14","SC-T0","MA-T14","MA-T0","Field-T14","Field-T0","NH-T14","NH-T0","ME-T14","ME-T0","NS-T14","NS-T0"))) %>%
  ggplot( aes(x=Location, y=Count, fill = Taxa)) + geom_bar(stat= "identity", color = "black") +theme_minimal() 


#no legend
#Loc xaxis
phy_16s_barplot %>%
  arrange(Count) %>%
  mutate(Location=factor(Location, levels=c("FL-T14","FL-T0","SC-T14","SC-T0","MA-T14","MA-T0","Field-T14","Field-T0","NH-T14","NH-T0","ME-T14","ME-T0","NS-T14","NS-T0"))) %>%
  ggplot( aes(x=Location, y=Count, fill = Taxa)) + geom_bar(stat= "identity", show.legend = F,color = "ivory2") +theme_minimal() 

phy_16s_barplot %>%
  arrange(Count) %>%
  mutate(Location=factor(Location, levels=c("FL-T14","FL-T0","SC-T14","SC-T0","MA-T14","MA-T0","Field-T14","Field-T0","NH-T14","NH-T0","ME-T14","ME-T0","NS-T14","NS-T0"))) %>%
  ggplot( aes(x=Location, y=Count, fill = Taxa)) + geom_bar(stat= "identity", show.legend = F, position="fill",color = "ivory3") +theme_minimal() 


#Loc yaxis
phy_16s_barplot %>%
  arrange(Count) %>%
  mutate(Location=factor(Location, levels=c("FL-T14","FL-T0","SC-T14","SC-T0","MA-T14","MA-T0","Field-T14","Field-T0","NH-T14","NH-T0","ME-T14","ME-T0","NS-T14","NS-T0"))) %>%
  ggplot( aes(x=Count, y=Location, fill = Taxa)) + geom_bar(stat= "identity", show.legend = F) +theme_minimal() 


#relative percentage (add:  position="fill" in geom_bar)
phy_16s_barplot %>%
  arrange(Count) %>%
  mutate(Location=factor(Location, levels=c("FL-T14","FL-T0","SC-T14","SC-T0","MA-T14","MA-T0","Field-T14","Field-T0","NH-T14","NH-T0","ME-T14","ME-T0","NS-T14","NS-T0"))) %>%
  ggplot( aes(x=Count, y=Location, fill = Taxa)) + geom_bar(stat= "identity", show.legend = F, position="fill",color = "ivory3") +theme_minimal() 

library(viridis)




## Venn diagrams ## 

library(ggplot2)
library(ggVennDiagram)
library(RVenn)
library(ggvenn)
if (!require(devtools)) install.packages("devtools")
devtools::install_github("yanlinlin82/ggvenn")



#### GENUS #### 
#import data
genus_venn_dat <- read.csv("~/Library/CloudStorage/GoogleDrive-sbirch1@charlotte.edu/My Drive/Postdoc_work/Mesocosom_exps/Mesocosom_Summer_2022/16s seq/QIIME Output/R/genus_venn_dat.csv", stringsAsFactors=TRUE)

#change rownames to taxa
rownames(genus_venn_dat)<-genus_venn_dat$taxa
## Had to delete the rows with no class info "g__" --> there were 443

#turn each column into a set
set_NS.0 <- genus_venn_dat %>% filter(NH.NS.T0 >= 1) %>% rownames()
set_NS.14 <- genus_venn_dat %>% filter(NH.NS.T14 >= 1) %>% rownames()
set_field.0 <- genus_venn_dat %>% filter(NH.Field.T0 >= 1) %>% rownames()
set_ME.0 <- genus_venn_dat %>% filter(NH.ME.T0 >= 1) %>% rownames()
set_ME.14 <- genus_venn_dat %>% filter(NH.ME.T14 >= 1) %>% rownames()
set_NH.0 <- genus_venn_dat %>% filter(NH.NH.T0 >= 1) %>% rownames()
set_NH.14 <- genus_venn_dat %>% filter(NH.NH.T14 >= 1) %>% rownames()
set_MA.0 <- genus_venn_dat %>% filter(NH.MA.T0 >= 1) %>% rownames()
set_MA.14 <- genus_venn_dat %>% filter(NH.MA.T14 >= 1) %>% rownames()
set_SC.0 <- genus_venn_dat %>% filter(NH.SC.T0 >= 1) %>% rownames()
set_SC.14 <- genus_venn_dat %>% filter(NH.SC.T14 >= 1) %>% rownames()
set_FL.0 <- genus_venn_dat %>% filter(NH.FL.T0 >= 1) %>% rownames()
set_FL.14 <- genus_venn_dat %>% filter(NH.FL.T14 >= 1) %>% rownames()

#this works - make a list of the comparison and then graph with ggVenn

#NS
NS_field <- list(T0 = set_NS.0, T14 = set_NS.14, ctrl_Field = set_field.0 )
ggVennDiagram(NS_field,category.names=c("NS T0", "NS T14", "Field")) + scale_fill_gradient(low= "#3D405B",high = "#E07A5F") + ggtitle("Comparison of Genera present in NS T0 vs NS T14 vs Field")

#ME
ME_field <- list(T0 = set_ME.0, T14 = set_ME.14, ctrl_Field = set_field.0 )
ggVennDiagram(ME_field,category.names=c("ME T0", "ME T14", "Field")) + scale_fill_gradient(low= "#3D405B",high = "#E07A5F") + ggtitle("Comparison of Genera present in ME T0 vs ME T14 vs Field")

#NH
NH_field <- list(T0 = set_NH.0, T14 = set_NH.14, ctrl_Field = set_field.0 )
ggVennDiagram(NH_field,category.names=c("NH T0", "NH T14", "Field")) + scale_fill_gradient(low= "#3D405B",high = "#E07A5F") + ggtitle("Comparison of Genera present in NH T0 vs NH T14 vs Field")

#MA
MA_field <- list(T0 = set_MA.0, T14 = set_MA.14, ctrl_Field = set_field.0 )
ggVennDiagram(MA_field,category.names=c("MA T0", "MA T14", "Field")) + scale_fill_gradient(low= "#3D405B",high = "#E07A5F") + ggtitle("Comparison of Genera present in MA T0 vs MA T14 vs Field")

#SC
SC_field <- list(T0 = set_SC.0, T14 = set_SC.14, ctrl_Field = set_field.0 )
ggVennDiagram(SC_field,category.names=c("SC T0", "SC T14", "Field")) + scale_fill_gradient(low= "#3D405B",high = "#E07A5F") + ggtitle("Comparison of Genera present in SC T0 vs SC T14 vs Field")

#FL
FL_field <- list(T0 = set_FL.0, T14 = set_FL.14, ctrl_Field = set_field.0 )
ggVennDiagram(FL_field,category.names=c("FL T0", "FL T14", "Field")) + scale_fill_gradient(low= "#3D405B",high = "#E07A5F") + ggtitle("Comparison of Genera present in FL T0 vs FL T14 vs Field")


#turn into a ven object
NS_venn<-Venn(NS_field)
ME_venn<-Venn(ME_field)
NH_venn<-Venn(NH_field)
MA_venn<-Venn(MA_field)
SC_venn<-Venn(SC_field)
FL_venn<-Venn(FL_field)

str(NS_field)
str(NS_venn)

#use just list for ggvenn the venn object for set map - works 
ggvenn(NS_field) 
setmap(NS_venn, title = "Presence/Absence of Genera in Field x NS_T0 x NS_T14")
setmap(ME_venn, title = "Presence/Absence of Genera in Field x ME_T0 x ME_T14")
setmap(NH_venn, title = "Presence/Absence of Genera in Field x NH_T0 x NH_T14")
setmap(MA_venn, title = "Presence/Absence of Genera in Field x MA_T0 x MA_T14")
setmap(SC_venn, title = "Presence/Absence of Genera in Field x SC_T0 x SC_T14")
setmap(FL_venn, title = "Presence/Absence of Genera in Field x FL_T0 x FL_T14")


## make a list of all T0 dat and T14 dat
All_T14_dat = list(NS_T14 = set_NS.14, ME_T14 = set_ME.14, NH_T14 = set_NH.14, MA_T14 = set_MA.14, SC_T14 = set_SC.14, FL_T14 = set_FL.14, Field = set_field.0)
All_T0_dat = list(NS_T0 = set_NS.0, ME_T0 = set_ME.0, NH_T0 = set_NH.0, MA_T0 = set_MA.0, SC_T0 = set_SC.0, FL_T0 = set_FL.0, Field = set_field.0)

#make venn object
All_14_venn<-Venn(All_T14_dat) #For all T14 locations
All_T0_venn<-Venn(All_T0_dat) #For all T0 locations + Field

#make heatmap of all T0 and T14s
setmap(All_T0_venn, title= "Presence/Absence of Genera across all T0 groups")
setmap(All_14_venn, title= "Presence/Absence of Genera across all T14 groups")

#overlap of all sets
tot_overlap_T0s<-overlap(All_T0_venn) #tot overlap for all T0 locations (43 genera)
tot_overlap_T14s<-overlap(All_14_venn) #tot overlap for all T14 locations (108 genera)
tot_overlap_NS<-overlap(NS_venn)
tot_overlap_ME<-overlap(ME_venn)
tot_overlap_NH<-overlap(NH_venn)
tot_overlap_MA<-overlap(MA_venn)
tot_overlap_SC<-overlap(SC_venn)
tot_overlap_FL<-overlap(FL_venn)

#get total pairwise overlap of each sample
NS_T0vT14_tot <-overlap(NS_venn, c("T0","T14"))
NS_T0vField_tot<-overlap(NS_venn,c("T0","ctrl_Field"))
NS_T14vField_tot<-overlap(NS_venn,c("T14","ctrl_Field"))

ME_T0vT14_tot <-overlap(ME_venn, c("T0","T14"))
ME_T0vField_tot<-overlap(ME_venn,c("T0","ctrl_Field"))
ME_T14vField_tot<-overlap(ME_venn,c("T14","ctrl_Field"))

NH_T0vT14_tot <-overlap(NH_venn, c("T0","T14"))
NH_T0vField_tot<-overlap(NH_venn,c("T0","ctrl_Field"))
NH_T14vField_tot<-overlap(NH_venn,c("T14","ctrl_Field"))

MA_T0vT14_tot <-overlap(MA_venn, c("T0","T14"))
MA_T0vField_tot<-overlap(MA_venn,c("T0","ctrl_Field"))
MA_T14vField_tot<-overlap(MA_venn,c("T14","ctrl_Field"))

SC_T0vT14_tot <-overlap(SC_venn, c("T0","T14"))
SC_T0vField_tot<-overlap(SC_venn,c("T0","ctrl_Field"))
SC_T14vField_tot<-overlap(SC_venn,c("T14","ctrl_Field"))

FL_T0vT14_tot <-overlap(FL_venn, c("T0","T14"))
FL_T0vField_tot<-overlap(FL_venn,c("T0","ctrl_Field"))
FL_T14vField_tot<-overlap(FL_venn,c("T14","ctrl_Field"))


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
str(NS_list_venn)

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

#make a list of field and overlaps with T14
Field_T14_list = list(NS_Field_T14 = NS_T14vField_actual, ME_Field_T14 = ME_T0vT14_actual, NH_Field_T14 = NH_T0vField_actual, MA_Field_T14 = MA_T14vField_actual, SC_Field_T14 = SC_T14vField_actual, FL_Field_T14 = FL_T14vField_actual)
Field_T14_venn <- Venn(Field_T14_list)
Field_T14_overlap<-overlap(Field_T14_venn)
setmap(Field_T14_venn, title= "Presence/Absence of Genera Field vs T14s")

#north (60 tot overlap)
north_14s_list = list(NS_Field_T14 = NS_T14vField_actual, ME_Field_T14 = ME_T14vField_actual, NH_Field_T14 = NH_T14vField_actual)
north_venn<- Venn(north_14s_list)
north_overlap<-overlap(north_venn)
ggvenn(north_venn)
ggVennDiagram(north_14s_list,category.names=c("NS T14-Field", "ME T14-Field", "NH T14-Field")) + scale_fill_gradient(low= "#3D405B",high = "#E07A5F") + ggtitle("Comparison of Genera present in Field vs North 14s")
setmap(north_venn, title= "Presence/Absence of Genera Field vs North T14s")

#south (48 tot overlap)
south_14s_list = list(MA_Field_T14 = MA_T14vField_actual, SC_Field_T14 = SC_T14vField_actual, FL_Field_T14 = FL_T14vField_actual)
south_venn<- Venn(south_14s_list)
south_overlap<-overlap(south_venn)
ggvenn(south_venn)
ggVennDiagram(south_14s_list,category.names=c("MA T14-Field", "SC T14-Field", "FL T14-Field")) + scale_fill_gradient(low= "#3D405B",high = "#E07A5F") + ggtitle("Comparison of Genera present in Field vs South 14s")
setmap(south_venn, title= "Presence/Absence of Genera Field vs South T14s")


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
length(write_list_NS$Total_overlap) #60
length(write_list_NS$T0vT14_overlap) #31
length(write_list_NS$T0vField_overlap) #10
length(write_list_NS$T14vField_overlap) #112

# list to adjust <- longest list
length(write_list_NS$Total_overlap) <- length(write_list_NS$T14vField_overlap)
length(write_list_NS$T0vT14_overlap) <- length(write_list_NS$T14vField_overlap)
length(write_list_NS$T0vField_overlap) <- length(write_list_NS$T14vField_overlap)

#ME
length(write_list_ME$Total_overlap) #64
length(write_list_ME$T0vT14_overlap) #19
length(write_list_ME$T0vField_overlap) #8
length(write_list_ME$T14vField_overlap) #116

# list to adjust <- longest list
length(write_list_ME$Total_overlap) <- length(write_list_ME$T14vField_overlap)
length(write_list_ME$T0vT14_overlap) <- length(write_list_ME$T14vField_overlap)
length(write_list_ME$T0vField_overlap) <- length(write_list_ME$T14vField_overlap)

#NH
length(write_list_NH$Total_overlap) #62
length(write_list_NH$T0vT14_overlap) #28
length(write_list_NH$T0vField_overlap) #14
length(write_list_NH$T14vField_overlap) #112

# list to adjust <- longest list
length(write_list_NH$Total_overlap) <- length(write_list_NH$T14vField_overlap)
length(write_list_NH$T0vT14_overlap) <- length(write_list_NH$T14vField_overlap)
length(write_list_NH$T0vField_overlap) <- length(write_list_NH$T14vField_overlap)

#MA
length(write_list_MA$Total_overlap) #80
length(write_list_MA$T0vT14_overlap) #42
length(write_list_MA$T0vField_overlap) #40
length(write_list_MA$T14vField_overlap) #84

# list to adjust <- longest list
length(write_list_MA$Total_overlap) <- length(write_list_MA$T14vField_overlap)
length(write_list_MA$T0vT14_overlap) <- length(write_list_MA$T14vField_overlap)
length(write_list_MA$T0vField_overlap) <- length(write_list_MA$T14vField_overlap)

#SC
length(write_list_SC$Total_overlap) #80
length(write_list_SC$T0vT14_overlap) #53
length(write_list_SC$T0vField_overlap) #13
length(write_list_SC$T14vField_overlap) #99

# list to adjust <- longest list
length(write_list_SC$Total_overlap) <- length(write_list_SC$T14vField_overlap)
length(write_list_SC$T0vT14_overlap) <- length(write_list_SC$T14vField_overlap)
length(write_list_SC$T0vField_overlap) <- length(write_list_SC$T14vField_overlap)

#FL
length(write_list_FL$Total_overlap) #90
length(write_list_FL$T0vT14_overlap) #53
length(write_list_FL$T0vField_overlap) #12
length(write_list_FL$T14vField_overlap) #138

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

