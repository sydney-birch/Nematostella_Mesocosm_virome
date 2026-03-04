##Nematostella Heatmaps for Meso 2022 56 Immune genes 

## This is For Field, FL, MA, NH, NS, SC T0 & T14 (3 reps in each)

## The salmon output is in mapping dir 
## made lib to stages
## also import csv for homohydra map: Immune_genes_56.csv

library(tximport); library(readr); library(edgeR); library(dplyr); library(tidyr); library(gplots); library(limma)

setwd("~/Library/CloudStorage/GoogleDrive-sbirch1@charlotte.edu/My Drive/Postdoc_work/Mesocosom_exps/Mesocosom_Summer_2022/Transcriptomics/EdgeR_MAPPED_4-29-24/Host_immune_heatmap")
system('ls')

dir <- getwd()
list.files()


#library info file
devstages<-read.table("libraries_to_stages_2.txt",header=F,row.names=1)
dev1<-rownames(devstages)[which(devstages$V2==1)]; dev1files <- file.path(dir, "../mapping",dev1, "quant.sf"); names(dev1files)<-dev1 #NS 0
dev2 <-rownames(devstages)[which(devstages$V2==2)]; dev2files <- file.path(dir, "../mapping",dev2, "quant.sf"); names(dev2files)<-dev2 # NS 14
dev3 <-rownames(devstages)[which(devstages$V2==3)]; dev3files <- file.path(dir, "../mapping",dev3, "quant.sf"); names(dev3files)<-dev3 #ME 0
dev4 <-rownames(devstages)[which(devstages$V2==4)]; dev4files <- file.path(dir, "../mapping",dev4, "quant.sf"); names(dev4files)<-dev4 #ME 14
dev5<-rownames(devstages)[which(devstages$V2==5)]; dev5files <- file.path(dir, "../mapping",dev5, "quant.sf"); names(dev5files)<-dev5 #NH 0
dev6 <-rownames(devstages)[which(devstages$V2==6)]; dev6files <- file.path(dir, "../mapping",dev6, "quant.sf"); names(dev6files)<-dev6 #NH 14
dev7 <-rownames(devstages)[which(devstages$V2==7)]; dev7files <- file.path(dir, "../mapping",dev7, "quant.sf"); names(dev7files)<-dev7 #Field 0
dev8<-rownames(devstages)[which(devstages$V2==8)]; dev8files <- file.path(dir, "../mapping",dev8, "quant.sf"); names(dev8files)<-dev8 #Field 14
dev9 <-rownames(devstages)[which(devstages$V2==9)]; dev9files <- file.path(dir, "../mapping",dev9, "quant.sf"); names(dev9files)<-dev9 #MA 0
dev10 <-rownames(devstages)[which(devstages$V2==10)]; dev10files <- file.path(dir, "../mapping",dev10, "quant.sf"); names(dev10files)<-dev10 # MA 14
dev11 <-rownames(devstages)[which(devstages$V2==11)]; dev11files <- file.path(dir, "../mapping",dev11, "quant.sf"); names(dev11files)<-dev11 #SC 0
dev12<-rownames(devstages)[which(devstages$V2==12)]; dev12files <- file.path(dir, "../mapping",dev12, "quant.sf"); names(dev12files)<-dev12 #SC 14
dev13 <-rownames(devstages)[which(devstages$V2==13)]; dev13files <- file.path(dir, "../mapping",dev13, "quant.sf"); names(dev13files)<-dev13 #Fl 0
dev14 <-rownames(devstages)[which(devstages$V2==14)]; dev14files <- file.path(dir, "../mapping",dev14, "quant.sf"); names(dev14files)<-dev14 #FL 14

#Normalize Gene counts
txi.salmon<- tximport(c(dev1files,dev2files, dev3files, dev4files,dev5files,dev6files,dev7files,dev8files,dev9files, dev10files, dev11files,dev12files,dev13files,dev14files), type = "salmon", txOut=T)
cts <- txi.salmon$counts
head(cts)
b<-DGEList(counts=cts, group=c(1,1,1,2,2,2,3,3,3,4,4,4,5,5,5,6,6,6,7,7,7,8,8,8,9,9,9,10,10,10,11,11,11,12,12,12,13,13,13,14,14,14))
normcounts<-cpm(b, normalized.lib.sizes=T) #normalized library counts for all genes

#Set up Gene ID map: 
homohydra.map <- read.csv("~/Library/CloudStorage/GoogleDrive-sbirch1@charlotte.edu/My Drive/Postdoc_work/Mesocosom_exps/Mesocosom_Summer_2022/Transcriptomics/EdgeR_MAPPED_4-29-24/Host_immune_heatmap/Immune_genes_56.csv")

str(homohydra.map) #should be data frame
length(homohydra.map) #4 col with 54 enteries
length(unique(homohydra.map$UK_Seq_ID)) #54
length(unique(homohydra.map$Gene_Name)) #54

library(tidyr)
library(dplyr)

#select just TF genes
keep<-unique(homohydra.map$UK_Seq_ID); length(keep) #1 col with 54 enteries 
TFcts <-normcounts[keep,]
nrow(TFcts) #25
sort(homohydra.map$Gene_Name)  #doublecheck that IDs are correctly sampled - should be gene symbol #25 
sort(rownames(TFcts)) #25

library(gplots)
library(viridis)

#construct heatmap of genes in map across devs
#using row as color break 
rgb.palette <- colorRampPalette(c( "#5DADE2","black","#F4D03F"),space = "Lab")

#using row as color break (normal)
heatmap.2(TFcts,  scale="row", labRow = homohydra.map$UK_Seq_ID,    Colv=F,  trace="none",  dendrogram="row",  key=F,  col=rgb.palette(120),  density.info=NULL,  margins=c(5, 11),  lmat=rbind(4:3, 2:1),  lhei=c(1, 30),  lwid=c(1, 1),cexRow = 0.80);
graphics.off()

# take average (normalized) counts for each dev stage
tcts<-t(TFcts); dev<-factor(c(1,1,1,2,2,2,3,3,3,4,4,4,5,5,5,6,6,6,7,7,7,8,8,8,9,9,9,10,10,10,11,11,11,12,12,12,13,13,13,14,14,14));dat<-cbind(dev, tcts)
dev.means<-by(dat, dev, function(x) colMeans(x[, -1]))
dev.means2<-as.data.frame(t(sapply(dev.means, I)))

#heatmaps using row as color break (normal)
heatmap.2(t(dev.means2), labRow = homohydra.map$Gene_Name, scale="row", Colv=F, trace="none", dendrogram="row", key=TRUE, col=rgb.palette(120), density.info=NULL,  margins=c(2.5, 10.75),  lmat=rbind(4:3, 2:1),  lhei=c(1, 30),  lwid=c(1, 1), cexRow = 0.80);
graphics.off() 
#heatmap using row as color break - legend
heatmap.2(t(dev.means2), labRow = homohydra.map$Gene_Name, scale="row", Colv=F, trace="none", dendrogram="row", key=TRUE, col=rgb.palette(120), density.info=NULL);
graphics.off() 

#viridis colors. - gives legend
heatmap.2(t(dev.means2), labRow = homohydra.map$Gene_Name, scale="row", Colv=F, trace="none", dendrogram="row", key=TRUE, col=viridis, density.info=NULL);
graphics.off() 
#viridis colors - actual heatmap
heatmap.2(t(dev.means2), labRow = homohydra.map$Gene_Name, scale="row", Colv=F, trace="none", dendrogram="row", key=TRUE, col=viridis, density.info=NULL,  margins=c(2.5, 10.75),  lmat=rbind(4:3, 2:1),  lhei=c(1, 30),  lwid=c(1, 1), cexRow = 0.80);
graphics.off() 



#### Run on Timepoint seperation: 

library(tximport); library(readr); library(edgeR); library(dplyr); library(tidyr); library(gplots); library(limma)

setwd("~/Library/CloudStorage/GoogleDrive-sbirch1@charlotte.edu/My Drive/Postdoc_work/Mesocosom_exps/Mesocosom_Summer_2022/Transcriptomics/EdgeR_MAPPED_4-29-24/Host_immune_heatmap")
system('ls')

dir <- getwd()
list.files()


#library info file
devstages<-read.table("libraries_to_stages_3.txt",header=F,row.names=1)
dev1<-rownames(devstages)[which(devstages$V2==1)]; dev1files <- file.path(dir, "../mapping",dev1, "quant.sf"); names(dev1files)<-dev1 #NS 0
dev2 <-rownames(devstages)[which(devstages$V2==2)]; dev2files <- file.path(dir, "../mapping",dev2, "quant.sf"); names(dev2files)<-dev2 # NS 14
dev3 <-rownames(devstages)[which(devstages$V2==3)]; dev3files <- file.path(dir, "../mapping",dev3, "quant.sf"); names(dev3files)<-dev3 #ME 0
dev4 <-rownames(devstages)[which(devstages$V2==4)]; dev4files <- file.path(dir, "../mapping",dev4, "quant.sf"); names(dev4files)<-dev4 #ME 14
dev5<-rownames(devstages)[which(devstages$V2==5)]; dev5files <- file.path(dir, "../mapping",dev5, "quant.sf"); names(dev5files)<-dev5 #NH 0
dev6 <-rownames(devstages)[which(devstages$V2==6)]; dev6files <- file.path(dir, "../mapping",dev6, "quant.sf"); names(dev6files)<-dev6 #NH 14
dev7 <-rownames(devstages)[which(devstages$V2==7)]; dev7files <- file.path(dir, "../mapping",dev7, "quant.sf"); names(dev7files)<-dev7 #Field 0
dev8<-rownames(devstages)[which(devstages$V2==8)]; dev8files <- file.path(dir, "../mapping",dev8, "quant.sf"); names(dev8files)<-dev8 #Field 14
dev9 <-rownames(devstages)[which(devstages$V2==9)]; dev9files <- file.path(dir, "../mapping",dev9, "quant.sf"); names(dev9files)<-dev9 #MA 0
dev10 <-rownames(devstages)[which(devstages$V2==10)]; dev10files <- file.path(dir, "../mapping",dev10, "quant.sf"); names(dev10files)<-dev10 # MA 14
dev11 <-rownames(devstages)[which(devstages$V2==11)]; dev11files <- file.path(dir, "../mapping",dev11, "quant.sf"); names(dev11files)<-dev11 #SC 0
dev12<-rownames(devstages)[which(devstages$V2==12)]; dev12files <- file.path(dir, "../mapping",dev12, "quant.sf"); names(dev12files)<-dev12 #SC 14
dev13 <-rownames(devstages)[which(devstages$V2==13)]; dev13files <- file.path(dir, "../mapping",dev13, "quant.sf"); names(dev13files)<-dev13 #Fl 0
dev14 <-rownames(devstages)[which(devstages$V2==14)]; dev14files <- file.path(dir, "../mapping",dev14, "quant.sf"); names(dev14files)<-dev14 #FL 14

#Normalize Gene counts
txi.salmon<- tximport(c(dev1files,dev2files, dev3files, dev4files,dev5files,dev6files,dev7files,dev8files,dev9files, dev10files, dev11files,dev12files,dev13files,dev14files), type = "salmon", txOut=T)
cts <- txi.salmon$counts
head(cts)
b<-DGEList(counts=cts, group=c(1,1,1,2,2,2,3,3,3,4,4,4,5,5,5,6,6,6,7,7,7,8,8,8,9,9,9,10,10,10,11,11,11,12,12,12,13,13,13,14,14,14))
normcounts<-cpm(b, normalized.lib.sizes=T) #normalized library counts for all genes

#Set up Gene ID map: 
homohydra.map <- read.csv("~/Library/CloudStorage/GoogleDrive-sbirch1@charlotte.edu/My Drive/Postdoc_work/Mesocosom_exps/Mesocosom_Summer_2022/Transcriptomics/EdgeR_MAPPED_4-29-24/Host_immune_heatmap/Immune_genes_56.csv")

str(homohydra.map) #should be data frame
length(homohydra.map) #4 col with 54 enteries
length(unique(homohydra.map$UK_Seq_ID)) #54


library(tidyr)
library(dplyr)

#select just TF genes
keep<-unique(homohydra.map$UK_Seq_ID); length(keep) #1 col with 54 enteries 
TFcts <-normcounts[keep,]
nrow(TFcts) #25
sort(homohydra.map$Gene_Name)  #doublecheck that IDs are correctly sampled - should be gene symbol #25 
sort(rownames(TFcts)) #25

library(gplots)
library(viridis)

#construct heatmap of genes in map across devs
#using row as color break 
rgb.palette <- colorRampPalette(c( "#5DADE2","black","#F4D03F"),space = "Lab")

#using row as color break (normal)
heatmap.2(TFcts,  scale="row", labRow = homohydra.map$Gene_Name,    Colv=F,  trace="none",  dendrogram="row",  key=F,  col=rgb.palette(120),  density.info=NULL,  margins=c(5, 11),  lmat=rbind(4:3, 2:1),  lhei=c(1, 30),  lwid=c(1, 1),cexRow = 0.80);
graphics.off() 

# take average (normalized) counts for each dev stage
tcts<-t(TFcts); dev<-factor(c(1,1,1,2,2,2,3,3,3,4,4,4,5,5,5,6,6,6,7,7,7,8,8,8,9,9,9,10,10,10,11,11,11,12,12,12,13,13,13,14,14,14));dat<-cbind(dev, tcts)
dev.means<-by(dat, dev, function(x) colMeans(x[, -1]))
dev.means2<-as.data.frame(t(sapply(dev.means, I)))

#heatmaps using row as color break (normal)
heatmap.2(t(dev.means2), labRow = homohydra.map$Gene_Name, scale="row", Colv=F, trace="none", dendrogram="row", key=TRUE, col=rgb.palette(120), density.info=NULL,  margins=c(2.5, 10.75),  lmat=rbind(4:3, 2:1),  lhei=c(1, 30),  lwid=c(1, 1), cexRow = 0.80);
graphics.off() 
#heatmap using row as color break - legend
heatmap.2(t(dev.means2), labRow = homohydra.map$Gene_Name, scale="row", Colv=F, trace="none", dendrogram="row", key=TRUE, col=rgb.palette(120), density.info=NULL);
graphics.off() 

#viridis colors. - gives legend
heatmap.2(t(dev.means2), labRow = homohydra.map$Gene_Name, scale="row", Colv=F, trace="none", dendrogram="row", key=TRUE, col=viridis, density.info=NULL);
graphics.off() 
#viridis colors - actual heatmap
heatmap.2(t(dev.means2), labRow = homohydra.map$Gene_Name, scale="row", Colv=F, trace="none", dendrogram="row", key=TRUE, col=viridis, density.info=NULL,  margins=c(2.5, 10.75),  lmat=rbind(4:3, 2:1),  lhei=c(1, 30),  lwid=c(1, 1), cexRow = 0.80);
graphics.off() 

#no dendro - viridis - reps - disable auto row - order alphabetical
heatmap.2(t(dev.means2),  scale="row", labRow = sort(homohydra.map$Gene_Name),  Rowv=F,  Colv=F,trace="none",  dendrogram="none",  key=F,  col=viridis,  density.info=NULL,  margins=c(5, 11),  lmat=rbind(4:3, 2:1),  lhei=c(1, 30),  lwid=c(1, 1),cexRow = 0.80);
graphics.off() 
#no dendro - viridis - reps - disable auto row - order alphabetical - key
heatmap.2(t(dev.means2),  scale="row", labRow = sort(homohydra.map$Gene_Name),  Rowv=F,  Colv=F,trace="none",  dendrogram="none",  key=TRUE,  col=viridis,  density.info=NULL);
graphics.off() 


#no dendro - viridis - reps - disable auto row - order alphabetical - no z scroe = absolute expression not by row
heatmap.2(t(dev.means2),  scale="none", labRow = sort(homohydra.map$Gene_Name),  Rowv=F,  Colv=F,trace="none",  dendrogram="none",  key=F,  col=viridis,  density.info=NULL,  margins=c(5, 11),  lmat=rbind(4:3, 2:1),  lhei=c(1, 30),  lwid=c(1, 1),cexRow = 0.80);
graphics.off() 
#no dendro - viridis - reps - disable auto row - order alphabetical - no z scroe = absolute expression not by row - key
heatmap.2(t(dev.means2),  scale="none", labRow = sort(homohydra.map$Gene_Name),  Rowv=F,  Colv=F,trace="none",  dendrogram="none",  key=TRUE,  col=viridis,  density.info=NULL);
graphics.off() 






