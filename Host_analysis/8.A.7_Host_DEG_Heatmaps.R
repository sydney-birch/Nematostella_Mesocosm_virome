## Nematostella (Host) DEG Heatmaps

  ## Previously ran edger and got sig degs, in mapped dir
  ## Run prep_for_R_string_deg.py on sig degs (10_R_prep_work/heatmaps-Mapped_edgeR/T0_v_T14_DEGs)

library(tximport); library(readr); library(edgeR); library(dplyr); library(tidyr); library(gplots); library(limma)

setwd("~/Library/CloudStorage/GoogleDrive-sbirch1@charlotte.edu/My Drive/Postdoc_work/Mesocosom_exps/Mesocosom_Summer_2022/Transcriptomics/EdgeR_MAPPED_4-29-24")
system('ls')

dir <- getwd()
list.files()


#library info file
devstages<-read.table("libraries_to_stages_2.txt",header=F,row.names=1)
dev1<-rownames(devstages)[which(devstages$V2==1)]; dev1files <- file.path(dir, "mapping",dev1, "quant.sf"); names(dev1files)<-dev1 #NS T0
dev2 <-rownames(devstages)[which(devstages$V2==2)]; dev2files <- file.path(dir, "mapping",dev2, "quant.sf"); names(dev2files)<-dev2 #NS T14
dev3 <-rownames(devstages)[which(devstages$V2==3)]; dev3files <- file.path(dir, "mapping",dev3, "quant.sf"); names(dev3files)<-dev3 #ME T0
dev4 <-rownames(devstages)[which(devstages$V2==4)]; dev4files <- file.path(dir, "mapping",dev4, "quant.sf"); names(dev4files)<-dev4 #ME T14
dev5<-rownames(devstages)[which(devstages$V2==5)]; dev5files <- file.path(dir, "mapping",dev5, "quant.sf"); names(dev5files)<-dev5 #NH T0
dev6 <-rownames(devstages)[which(devstages$V2==6)]; dev6files <- file.path(dir, "mapping",dev6, "quant.sf"); names(dev6files)<-dev6 #NH T14
dev7 <-rownames(devstages)[which(devstages$V2==7)]; dev7files <- file.path(dir, "mapping",dev7, "quant.sf"); names(dev7files)<-dev7 #Field T0
dev8<-rownames(devstages)[which(devstages$V2==8)]; dev8files <- file.path(dir, "mapping",dev8, "quant.sf"); names(dev8files)<-dev8 #Field T14
dev9 <-rownames(devstages)[which(devstages$V2==9)]; dev9files <- file.path(dir, "mapping",dev9, "quant.sf"); names(dev9files)<-dev9 #MA T0
dev10 <-rownames(devstages)[which(devstages$V2==10)]; dev10files <- file.path(dir, "mapping",dev10, "quant.sf"); names(dev10files)<-dev10 #MA T14
dev11 <-rownames(devstages)[which(devstages$V2==11)]; dev11files <- file.path(dir, "mapping",dev11, "quant.sf"); names(dev11files)<-dev11 #SC T0
dev12<-rownames(devstages)[which(devstages$V2==12)]; dev12files <- file.path(dir, "mapping",dev12, "quant.sf"); names(dev12files)<-dev12 #SC T14
dev13 <-rownames(devstages)[which(devstages$V2==13)]; dev13files <- file.path(dir, "mapping",dev13, "quant.sf"); names(dev13files)<-dev13 #FL T0
dev14 <-rownames(devstages)[which(devstages$V2==14)]; dev14files <- file.path(dir, "mapping",dev14, "quant.sf"); names(dev14files)<-dev14 #FL T14

#Running on each genotype
## FL T0 vs T14 

txi.salmon<- tximport(c(dev13files,dev14files), type = "salmon", txOut=T)
cts <- txi.salmon$counts
head(cts)
b<-DGEList(counts=cts, group=c(13,13,13,14,14,14))
normcounts<-cpm(b, normalized.lib.sizes=T) #normalized library counts for all genes


#make a matrix of significant degs from this data set - run prep_sig_DEGs_for_heatmap.py in terminal to get sig degs formated in this way
sig_degs<-c("XR_007306819.1","XM_048721243.1","XM_048724741.1","XM_048731718.1","XR_007306253.1","XM_001631917.3","XM_032382958.2","XM_001635407.3","XM_001632539.3","XM_032363830.2","XM_032385697.2","XM_048733319.1") #Reduced 

#make into a data frame and change name of column 
sig_degs_df<-data.frame(sig_degs)
str(sig_degs_df)
names(sig_degs_df)[1] <- c("Stella_seqid")

#make data characters not factor so they can be compared
sig_degs_df$Stella_seqid<- as.character(sig_degs_df$Stella_seqid)
str(sig_degs_df)

#homohydra.map<- read.table("Hs_detect_light_stim_genes_R_output.txt",header=T, stringsAsFactors=F) #Full
homohydra.map<- read.table("full_uk_transcriptome_headers.txt",header=T, stringsAsFactors=F) #Reduced
str(homohydra.map) #should be data frame

library(tidyr)
library(dplyr)

#join sig degs to homohydra.map
str(homohydra.map)
sig_set<-inner_join(homohydra.map[,1:2],sig_degs_df,by="Stella_seqid")
nrow(sig_set)
length(unique(sig_set$Stella_seqid))


#homohydra.map<-unique_sig_set
#keep<-unique(homohydra.map$Hydractinia_seqid); length(keep) #1 col with 32 enteries 
homohydra.map<-sig_set
keep<-sig_set$Stella_seqid

#select just TF genes
TFcts <-normcounts[keep,]
nrow(TFcts) #25
sort(rownames(TFcts)) #25

#construct heatmap of genes in map across devs
#using row as color break 
rgb.palette <- colorRampPalette(c( "#5DADE2","black","#F4D03F"),space = "Lab")
#specifying color break 
rgb.palette2 <- colorRampPalette(c( "#5DADE2","black","#F4D03F"))(n=299) #,space = "Lab"
col_breaks = c(seq(0,33,length=100), seq(33.01,66,length=100), seq(66.01,100,length=100))
library(viridis)
library(gplots)


#using row as color break
heatmap.2(TFcts,  scale="row", labRow = homohydra.map$header,  Colv=F,  trace="none",  dendrogram="row",  key=T,  col=rgb.palette(120),  density.info=NULL, margins=c(5, 11), lmat=rbind(4:3, 2:1),  lhei=c(1, 30),  lwid=c(1, 1),cexRow = 0.80);
graphics.off() 
par(mfrow = c(1,1))
heatmap.2(TFcts,  scale="row", labRow = homohydra.map$header, Colv=F,  trace="none", col=rgb.palette, ColSideColors=rep(c('#83C5BE','#006D77'), each=3),dendrogram="row")


# take average (normalized) counts for each dev stage
tcts<-t(TFcts); dev<-factor(c(13,13,13,14,14,14));dat<-cbind(dev, tcts)
dev.means<-by(dat, dev, function(x) colMeans(x[, -1]))
dev.means2<-as.data.frame(t(sapply(dev.means, I)))

#heatmaps using row as color break
heatmap.2(t(dev.means2), labRow = homohydra.map$header, scale="row", Colv=F, trace="none", dendrogram="row", key=T, col=rgb.palette(120), density.info=NULL,  margins=c(2.5, 10.75),  lmat=rbind(4:3, 2:1),  lhei=c(1, 30),  lwid=c(1, 1), cexRow = 0.80);
heatmap.2(t(dev.means2), cellnote= round(t(dev.means2), digits = 3), notecol = "white", labRow = homohydra.map$header, scale="row", Colv=F, trace="none", dendrogram="row", key=F, col=rgb.palette(120), density.info=NULL)

          
          

#Running on each genotype
## SC T0 vs T14 

txi.salmon<- tximport(c(dev11files,dev12files), type = "salmon", txOut=T)
cts <- txi.salmon$counts
head(cts)
b<-DGEList(counts=cts, group=c(11,11,11,12,12,12))
normcounts<-cpm(b, normalized.lib.sizes=T) #normalized library counts for all genes


#make a matrix of significant degs from this data set - run prep_sig_DEGs_for_heatmap.py in terminal to get sig degs formated in this way
sig_degs<-c("XR_007305401.1","XR_007305681.1","XM_048731672.1","XM_001634477.3","XR_004292988.2","XM_032381817.2","XM_032381798.2","XM_032370603.2","XM_032375039.2","XR_007307226.1","XM_048726115.1","XR_007306627.1","XM_048727377.1","XM_048723243.1","XM_048731381.1","XM_001637966.3","XM_048731023.1","XM_048731718.1","XM_048734697.1","XM_048731326.1","XR_007306626.1","XM_001638417.3","XM_048724741.1","XM_048720927.1","XM_001636428.3","XM_032378265.2","XM_048734254.1","XM_048731952.1","XR_004292730.1","XM_048725419.1","XM_032371824.2","XM_048725205.1","XM_048730561.1","XM_032367576.2","XR_007307117.1","XR_007307495.1","XM_001631087.3","XM_032366754.2","XM_001633659.3","XM_032382501.2","XM_048728941.1","XM_032387003.2","XM_001626957.3") #Reduced 

#make into a data frame and change name of column 
sig_degs_df<-data.frame(sig_degs)
str(sig_degs_df)
names(sig_degs_df)[1] <- c("Stella_seqid")

#make data characters not factor so they can be compared
sig_degs_df$Stella_seqid<- as.character(sig_degs_df$Stella_seqid)
str(sig_degs_df)

#homohydra.map<- read.table("Hs_detect_light_stim_genes_R_output.txt",header=T, stringsAsFactors=F) #Full
homohydra.map<- read.table("full_uk_transcriptome_headers.txt",header=T, stringsAsFactors=F) #Reduced
str(homohydra.map) #should be data frame

library(tidyr)
library(dplyr)

#join sig degs to homohydra.map
str(homohydra.map)
sig_set<-inner_join(homohydra.map[,1:2],sig_degs_df,by="Stella_seqid")
nrow(sig_set)
length(unique(sig_set$Stella_seqid))

homohydra.map<-sig_set
keep<-sig_set$Stella_seqid

#select just TF genes
TFcts <-normcounts[keep,]
nrow(TFcts) #25
sort(rownames(TFcts)) #25

#construct heatmap of genes in map across devs
#using row as color break 
rgb.palette <- colorRampPalette(c( "#5DADE2","black","#F4D03F"),space = "Lab")
#specifying color break 
rgb.palette2 <- colorRampPalette(c( "#5DADE2","black","#F4D03F"))(n=299) #,space = "Lab"
col_breaks = c(seq(0,33,length=100), seq(33.01,66,length=100), seq(66.01,100,length=100))

library(gplots)

#using row as color break
heatmap.2(TFcts,  scale="row", labRow = homohydra.map$header,  Colv=F,  trace="none",  dendrogram="row",  key=T,  col=rgb.palette(120),  density.info=NULL, margins=c(5, 11), lmat=rbind(4:3, 2:1),  lhei=c(1, 30),  lwid=c(1, 1),cexRow = 0.80);
graphics.off() 
par(mfrow = c(1,1))
heatmap.2(TFcts,  scale="row", labRow = homohydra.map$header, Colv=F,  trace="none", col=rgb.palette, ColSideColors=rep(c('#83C5BE','#006D77'), each=3),dendrogram="row")



# take average (normalized) counts for each dev stage
tcts<-t(TFcts); dev<-factor(c(11,11,11,12,12,12));dat<-cbind(dev, tcts)
dev.means<-by(dat, dev, function(x) colMeans(x[, -1]))
dev.means2<-as.data.frame(t(sapply(dev.means, I)))

#heatmaps using row as color break
heatmap.2(t(dev.means2), labRow = homohydra.map$header, scale="row", Colv=F, trace="none", dendrogram="row", key=T, col=rgb.palette(120), density.info=NULL,  margins=c(2.5, 10.75),  lmat=rbind(4:3, 2:1),  lhei=c(1, 30),  lwid=c(1, 1), cexRow = 0.80);
heatmap.2(t(dev.means2), cellnote= round(t(dev.means2), digits = 3), notecol = "white", labRow = homohydra.map$header, scale="row", Colv=F, trace="none", dendrogram="row", key=F, col=rgb.palette(120), density.info=NULL)



#Running on each genotype
## MA T0 vs T14 

txi.salmon<- tximport(c(dev9files,dev10files), type = "salmon", txOut=T)
cts <- txi.salmon$counts
head(cts)
b<-DGEList(counts=cts, group=c(9,9,9,10,10,10))
normcounts<-cpm(b, normalized.lib.sizes=T) #normalized library counts for all genes


#make a matrix of significant degs from this data set - run prep_sig_DEGs_for_heatmap.py in terminal to get sig degs formated in this way
sig_degs<-c("XM_001632539.3","XR_007306819.1","XR_004292988.2","XM_032385756.2","XM_001628020.3","XM_001635041.3","XM_001625470.3","XM_001636932.3","XM_001631917.3","XM_032378581.2","XM_032362107.2","XM_048724741.1","XM_001637875.3","XM_048721533.1","XM_048730053.1","XM_032362443.2","XM_001630092.3","XM_048721527.1","XR_004294420.2","XM_048733907.1","XM_032364263.2","XM_032365359.2","XM_048730855.1","XM_032379762.2","XM_001631945.2") #Reduced 

#make into a data frame and change name of column 
sig_degs_df<-data.frame(sig_degs)
str(sig_degs_df)
names(sig_degs_df)[1] <- c("Stella_seqid")

#make data characters not factor so they can be compared
sig_degs_df$Stella_seqid<- as.character(sig_degs_df$Stella_seqid)
str(sig_degs_df)

#homohydra.map<- read.table("Hs_detect_light_stim_genes_R_output.txt",header=T, stringsAsFactors=F) #Full
homohydra.map<- read.table("full_uk_transcriptome_headers.txt",header=T, stringsAsFactors=F) #Reduced
str(homohydra.map) #should be data frame

library(tidyr)
library(dplyr)

#join sig degs to homohydra.map
str(homohydra.map)
sig_set<-inner_join(homohydra.map[,1:2],sig_degs_df,by="Stella_seqid")
nrow(sig_set)
length(unique(sig_set$Stella_seqid))


homohydra.map<-sig_set
keep<-sig_set$Stella_seqid

#select just TF genes
TFcts <-normcounts[keep,]
nrow(TFcts) #25
sort(rownames(TFcts)) #25

#construct heatmap of genes in map across devs
#using row as color break 
rgb.palette <- colorRampPalette(c( "#5DADE2","black","#F4D03F"),space = "Lab")
#specifying color break 
#rgb.palette2 <- colorRampPalette(c( "#5DADE2","black","#F4D03F"))(n=299) #,space = "Lab"
#col_breaks = c(seq(0,33,length=100), seq(33.01,66,length=100), seq(66.01,100,length=100))

library(gplots)

#using row as color break
graphics.off() 
par(mfrow = c(1,1))
heatmap.2(TFcts,  scale="row", labRow = homohydra.map$header,  Colv=F,  trace="none",  dendrogram="row",  key=T,  col=rgb.palette(120),  density.info=NULL, margins=c(5, 11), lmat=rbind(4:3, 2:1),  lhei=c(1, 30),  lwid=c(1, 1),cexRow = 0.80);
graphics.off() 
par(mfrow = c(1,1))
heatmap.2(TFcts,  scale="row", labRow = homohydra.map$header, Colv=F,  trace="none", col=rgb.palette, ColSideColors=rep(c('#83C5BE','#006D77'), each=3),dendrogram="row")
graphics.off() 
par(mfrow = c(1,1))


# take average (normalized) counts for each dev stage
tcts<-t(TFcts); dev<-factor(c(9,9,9,10,10,10));dat<-cbind(dev, tcts)
dev.means<-by(dat, dev, function(x) colMeans(x[, -1]))
dev.means2<-as.data.frame(t(sapply(dev.means, I)))

#heatmaps using row as color break
heatmap.2(t(dev.means2), labRow = homohydra.map$header, scale="row", Colv=F, trace="none", dendrogram="row", key=T, col=rgb.palette(120), density.info=NULL,  margins=c(2.5, 10.75),  lmat=rbind(4:3, 2:1),  lhei=c(1, 30),  lwid=c(1, 1), cexRow = 0.80);
heatmap.2(t(dev.means2), cellnote= round(t(dev.means2), digits = 3), notecol = "white", labRow = homohydra.map$header, scale="row", Colv=F, trace="none", dendrogram="row", key=F, col=rgb.palette(120), density.info=NULL)


          

#Running on each genotype
## Field T0 vs T14 

txi.salmon<- tximport(c(dev7files,dev8files), type = "salmon", txOut=T)
cts <- txi.salmon$counts
head(cts)
b<-DGEList(counts=cts, group=c(7,7,7,8,8,8))
normcounts<-cpm(b, normalized.lib.sizes=T) #normalized library counts for all genes


#make a matrix of significant degs from this data set - run prep_sig_DEGs_for_heatmap.py in terminal to get sig degs formated in this way
sig_degs<-c("XM_001629965.3","XM_032366283.2","XM_032384803.2","XM_048734625.1","XM_048732439.1","XM_032381131.2","XM_032385682.2","XM_032377463.2","XM_032371824.2","XM_048725001.1","XM_048724214.1","XM_048726145.1","XM_048733179.1","XM_001634132.3","XM_048726383.1","XM_001632982.3","XM_048719951.1","XM_001629432.3","XM_048723130.1","XM_032366282.2","XM_032386219.2","XM_032364280.2","XM_032381223.2","XR_007306796.1","XM_001636455.3") #Reduced 

#make into a data frame and change name of column 
sig_degs_df<-data.frame(sig_degs)
str(sig_degs_df)
names(sig_degs_df)[1] <- c("Stella_seqid")

#make data characters not factor so they can be compared
sig_degs_df$Stella_seqid<- as.character(sig_degs_df$Stella_seqid)
str(sig_degs_df)

#homohydra.map<- read.table("Hs_detect_light_stim_genes_R_output.txt",header=T, stringsAsFactors=F) #Full
homohydra.map<- read.table("full_uk_transcriptome_headers.txt",header=T, stringsAsFactors=F) #Reduced
str(homohydra.map) #should be data frame

library(tidyr)
library(dplyr)

#join sig degs to homohydra.map
str(homohydra.map)
sig_set<-inner_join(homohydra.map[,1:2],sig_degs_df,by="Stella_seqid")
nrow(sig_set)
length(unique(sig_set$Stella_seqid))

homohydra.map<-sig_set
keep<-sig_set$Stella_seqid

#select just TF genes
TFcts <-normcounts[keep,]
nrow(TFcts) #25
sort(rownames(TFcts)) #25

#construct heatmap of genes in map across devs
#using row as color break 
rgb.palette <- colorRampPalette(c( "#5DADE2","black","#F4D03F"),space = "Lab")
#specifying color break 
#rgb.palette2 <- colorRampPalette(c( "#5DADE2","black","#F4D03F"))(n=299) #,space = "Lab"
#col_breaks = c(seq(0,33,length=100), seq(33.01,66,length=100), seq(66.01,100,length=100))

library(gplots)

#using row as color break
graphics.off() 
par(mfrow = c(1,1))
heatmap.2(TFcts,  scale="row", labRow = homohydra.map$header,  Colv=F,  trace="none",  dendrogram="row",  key=T,  col=rgb.palette(120),  density.info=NULL, margins=c(5, 11), lmat=rbind(4:3, 2:1),  lhei=c(1, 30),  lwid=c(1, 1),cexRow = 0.80);
graphics.off() 
par(mfrow = c(1,1))
heatmap.2(TFcts,  scale="row", labRow = homohydra.map$header, Colv=F,  trace="none", col=rgb.palette, ColSideColors=rep(c('#83C5BE','#006D77'), each=3),dendrogram="row")
graphics.off() 
par(mfrow = c(1,1))


# take average (normalized) counts for each dev stage
tcts<-t(TFcts); dev<-factor(c(7,7,7,8,8,8));dat<-cbind(dev, tcts)
dev.means<-by(dat, dev, function(x) colMeans(x[, -1]))
dev.means2<-as.data.frame(t(sapply(dev.means, I)))

#heatmaps using row as color break
heatmap.2(t(dev.means2), labRow = homohydra.map$header, scale="row", Colv=F, trace="none", dendrogram="row", key=T, col=rgb.palette(120), density.info=NULL,  margins=c(2.5, 10.75),  lmat=rbind(4:3, 2:1),  lhei=c(1, 30),  lwid=c(1, 1), cexRow = 0.80);
heatmap.2(t(dev.means2), cellnote= round(t(dev.means2), digits = 3), notecol = "white", labRow = homohydra.map$header, scale="row", Colv=F, trace="none", dendrogram="row", key=F, col=rgb.palette(120), density.info=NULL)

          
          
          
          
#Running on each genotype
## NH T0 vs T14 

txi.salmon<- tximport(c(dev5files,dev6files), type = "salmon", txOut=T)
cts <- txi.salmon$counts
head(cts)
b<-DGEList(counts=cts, group=c(5,5,5,6,6,6))
normcounts<-cpm(b, normalized.lib.sizes=T) #normalized library counts for all genes


#make a matrix of significant degs from this data set - run prep_sig_DEGs_for_heatmap.py in terminal to get sig degs formated in this way
sig_degs<-c("XM_032369564.2","XM_048731718.1","XR_004295687.2","XM_001628020.3","XM_032382924.2","XM_032381836.2","XM_048733772.1","XM_032364762.2","XM_048725094.1","XM_001636034.3","XM_048725419.1","XM_001641466.3","XM_001627599.3","XM_048723501.1","XM_032384424.2","XM_001632539.3","XM_001624935.3","XM_048730410.1","XM_032371878.2","XM_032373961.2","XM_048729070.1","XM_001622599.3","XM_001629293.3","XM_048727740.1","XM_001639572.3","XM_032372917.2","XM_032386556.2","XM_048723557.1","XM_032385144.2","XM_032375534.2","XM_001636932.3","XM_032385537.2","XM_048726238.1","XM_032372536.2","XM_048729334.1","XM_048720799.1","XM_032379332.2","XM_048734688.1","XM_032381469.2","XM_032365761.2","XM_048726145.1","XM_048726018.1","XM_048728999.1") #Reduced 

#make into a data frame and change name of column 
sig_degs_df<-data.frame(sig_degs)
str(sig_degs_df)
names(sig_degs_df)[1] <- c("Stella_seqid")

#make data characters not factor so they can be compared
sig_degs_df$Stella_seqid<- as.character(sig_degs_df$Stella_seqid)
str(sig_degs_df)

#homohydra.map<- read.table("Hs_detect_light_stim_genes_R_output.txt",header=T, stringsAsFactors=F) #Full
homohydra.map<- read.table("full_uk_transcriptome_headers.txt",header=T, stringsAsFactors=F) #Reduced
str(homohydra.map) #should be data frame

library(tidyr)
library(dplyr)

#join sig degs to homohydra.map
str(homohydra.map)
sig_set<-inner_join(homohydra.map[,1:2],sig_degs_df,by="Stella_seqid")
nrow(sig_set)
length(unique(sig_set$Stella_seqid))

homohydra.map<-sig_set
keep<-sig_set$Stella_seqid

#select just TF genes
TFcts <-normcounts[keep,]
nrow(TFcts) #25
sort(rownames(TFcts)) #25

#construct heatmap of genes in map across devs
#using row as color break 
rgb.palette <- colorRampPalette(c( "#5DADE2","black","#F4D03F"),space = "Lab")
#specifying color break 
#rgb.palette2 <- colorRampPalette(c( "#5DADE2","black","#F4D03F"))(n=299) #,space = "Lab"
#col_breaks = c(seq(0,33,length=100), seq(33.01,66,length=100), seq(66.01,100,length=100))

library(gplots)

#using row as color break
graphics.off() 
par(mfrow = c(1,1))
heatmap.2(TFcts,  scale="row", labRow = homohydra.map$header,  Colv=F,  trace="none",  dendrogram="row",  key=T,  col=rgb.palette(120),  density.info=NULL, margins=c(5, 11), lmat=rbind(4:3, 2:1),  lhei=c(1, 30),  lwid=c(1, 1),cexRow = 0.80);
graphics.off() 
par(mfrow = c(1,1))
heatmap.2(TFcts,  scale="row", labRow = homohydra.map$header, Colv=F,  trace="none", col=rgb.palette, ColSideColors=rep(c('#83C5BE','#006D77'), each=3),dendrogram="row")
graphics.off() 
par(mfrow = c(1,1))


# take average (normalized) counts for each dev stage
tcts<-t(TFcts); dev<-factor(c(5,5,5,6,6,6));dat<-cbind(dev, tcts)
dev.means<-by(dat, dev, function(x) colMeans(x[, -1]))
dev.means2<-as.data.frame(t(sapply(dev.means, I)))

#heatmaps using row as color break
heatmap.2(t(dev.means2), labRow = homohydra.map$header, scale="row", Colv=F, trace="none", dendrogram="row", key=T, col=rgb.palette(120), density.info=NULL,  margins=c(2.5, 10.75),  lmat=rbind(4:3, 2:1),  lhei=c(1, 30),  lwid=c(1, 1), cexRow = 0.80);
heatmap.2(t(dev.means2), cellnote= round(t(dev.means2), digits = 3), notecol = "white", labRow = homohydra.map$header, scale="row", Colv=F, trace="none", dendrogram="row", key=F, col=rgb.palette(120), density.info=NULL)

          
          
          
          
          
#Running on each genotype
## ME T0 vs T14 

txi.salmon<- tximport(c(dev3files,dev4files), type = "salmon", txOut=T)
cts <- txi.salmon$counts
head(cts)
b<-DGEList(counts=cts, group=c(3,3,3,4,4,4))
normcounts<-cpm(b, normalized.lib.sizes=T) #normalized library counts for all genes


#make a matrix of significant degs from this data set - run prep_sig_DEGs_for_heatmap.py in terminal to get sig degs formated in this way
sig_degs<-c("XR_007306877.1","XM_048729070.1","XR_007306918.1","XM_001631917.3","XM_001636932.3","XM_048722657.1","XM_001624482.3","XM_001623270.3","XM_001631422.3","XR_007306253.1","XM_048733772.1","XM_048734432.1","XM_032382924.2","XR_007312385.1","XM_001626671.3","XM_032371825.2","XM_048721249.1","XR_007305260.1","XM_048732914.1","XR_007306916.1","XM_032377265.2","XM_001636392.3","XM_048731960.1","XM_032375382.2","XM_048730246.1","XM_001634075.3","XM_001627064.3","XM_001629293.3","XM_048720915.1","XM_048724741.1","XM_048730410.1","XM_032365879.2","XM_048723806.1","XM_048730683.1","XM_048728999.1","XR_007312194.1","XR_007314268.1","XM_001626679.3","XM_048725109.1","XM_048723718.1","XM_001632982.3","XM_001635407.3","XM_032386726.2","XR_007312590.1","XM_032377385.2","XM_048726078.1","XM_032364009.2","XM_032386007.2","XM_048725092.1","XM_048721953.1","XM_001638554.3","XR_007312192.1","XM_048721179.1","XM_001640639.3","XM_048734604.1","XM_001638165.3","XM_032385144.2","XM_032367448.2","XM_001626241.3","XM_001628740.3","XM_032387399.2","XM_048724100.1","XM_001641929.3","XM_032377993.2","XM_032363372.2","XM_001633287.3","XM_032386927.2","XM_032383745.2","XR_007307164.1","XM_048733738.1","XM_001631860.3","XM_048723498.1","XM_048721957.1","XM_048729371.1","XM_048720580.1","XM_048722692.1","XM_001636654.3","XM_048730578.1","XM_001626957.3","XM_048725787.1","XM_048729734.1","XM_048726050.1","XM_001627279.3","XM_001627402.3","XM_048727207.1","XM_032364772.2","XM_032376148.2","XM_048732427.1","XM_001634488.3","XM_048725227.1","XM_048720696.1","XM_048725404.1","XM_048723705.1","XM_032369445.2","XM_032365207.2","XM_032385117.2","XM_048722917.1","XM_032364263.2","XM_001629238.3","XM_048725366.1","XM_048723494.1","XM_032370494.2","XM_032380541.2","XM_001630669.3","XM_048728876.1","XM_048732896.1","XR_007307142.1","XM_048730601.1","XM_032386970.2","XM_048719709.1","XM_048721527.1","XM_032373921.2","XM_001636405.3","XR_007305305.1","XM_032379630.2","XM_001636514.3","XM_048725113.1","XM_048723488.1","XM_048721268.1","XM_032371278.2","XM_048721809.1","XM_048721789.1","XM_048731565.1","XM_048727188.1","XM_001629396.3","XR_007306627.1","XM_032367550.2","XM_048723496.1","XM_032367112.2","XM_048723504.1","XM_048728499.1","XM_001636071.3","XM_048730906.1","XM_001635071.3","XM_048721466.1","XM_048728827.1","XM_048733363.1","XM_048724444.1","XM_048732509.1","XM_048732374.1","XM_001629503.3","XM_032378599.2","XR_004292921.2","XM_001626090.3","XR_004297034.2","XM_032373987.2","XM_032366219.2","XM_032380321.2","XM_001630315.3","XM_032363258.2","XM_048729040.1","XM_032371555.2","XM_032365218.2","XM_001637845.3","XM_001636301.3","XM_048720698.1","XM_048729222.1","XM_001636575.3","XM_032384667.2","XM_001636598.3","XM_048724558.1","XM_032376542.2","XM_048721132.1","XM_048731707.1","XM_048721333.1","XM_032366595.2","XR_004293447.2","XM_048723128.1","XM_048726619.1","XM_001638302.3","XM_048723524.1","XM_032381054.2","XM_048732713.1","XM_032364294.2","XM_001630398.3","XM_001626638.3","XM_048721497.1","XM_032376480.2","XM_001638489.3","XM_048722775.1","XR_007305189.1","XM_048734583.1","XM_001619797.3","XM_032380348.2","XM_001637591.3","XM_001626761.3","XM_048729738.1","XM_032379655.2","XM_032370869.2","XR_007305698.1","XM_048721139.1","XM_048721936.1","XM_001638134.3","XM_001624294.3","XM_032376256.2","XM_048723918.1","XR_007307557.1","XM_032363460.2","XM_048733814.1","XM_001633333.3","XM_048733990.1","XR_007306626.1","XM_048732164.1","XM_001627740.3","XM_048719665.1","XM_001647469.3","XM_032374527.2","XM_001629388.3","XM_048723801.1","XR_004296830.2","XM_001633335.3","XM_001623710.2","XM_032364360.2","XM_048731936.1","XM_032363960.2","XM_032381139.2","XM_001623467.3","XM_048729302.1","XM_032370334.2","XM_032381347.2","XM_032366159.2","XR_007306629.1","XM_001631991.3","XM_001637837.3","XM_001634064.3","XM_048724736.1","XM_048734532.1","XM_001636086.3","XM_032366603.2","XR_007306415.1","XM_032384093.2","XM_001639580.3","XM_048728792.1","XM_048726713.1","XM_048727178.1","XM_048725767.1","XM_048725562.1","XM_048729564.1","XM_001627083.3","XM_048733448.1","XR_007313731.1","XM_032383015.2","XM_001638034.3","XR_007312046.1","XM_048731663.1","XR_004294420.2","XM_001634477.3","XM_048733739.1","XM_001626176.3","XM_032364427.2","XM_001640381.3","XM_032371550.2","XR_004296070.2","XM_032367265.2","XM_048726077.1","XM_032375167.2","XM_048723244.1","XM_048733318.1","XM_001633843.3","XM_032379762.2","XM_032379332.2","XM_032379062.2","XM_048733547.1","XM_032365683.2","XM_048720802.1","XM_048729559.1","XM_032366265.2","XM_032385174.2","XM_001630073.3","XM_001629087.2","XM_001628266.3","XM_001628128.3","XM_032386257.2","XM_032366858.2","XM_001641749.3","XM_001647486.3","XM_048722456.1","XM_032387288.2","XM_032363159.2","XM_032382712.2","XM_048724725.1","XM_048720920.1","XM_048723441.1","XM_032367590.2","XM_048731280.1","XM_048722130.1","XM_048726173.1","XM_001628583.3","XM_032382123.2","XM_048725998.1","XM_001633796.3","XM_048721326.1","XM_048722628.1","XM_001636034.3","XM_032375357.2","XM_032370047.2","XM_001634542.3","XM_032384885.2","XM_048729010.1","XM_032363570.2","XM_048721305.1","XM_048720024.1","XM_048721253.1","XM_048728389.1","XM_048728988.1","XM_001619813.3","XM_001629597.3","XM_048719979.1","XM_048732822.1","XM_032379789.2","XM_032363538.2","XM_048721531.1","XM_048725785.1","XM_032370281.2","XM_032379211.2","XM_048733407.1","XM_001640912.3","XM_001633518.3","XM_048729633.1","XM_001625714.3","XM_048723760.1","XR_004297331.2","XM_032378741.2","XM_048729791.1","XR_007314261.1","XM_048732700.1","XM_048729736.1","XM_001633659.3","XM_001631518.3","XM_032371824.2","XR_007306989.1","XM_032377681.2","XM_048729205.1","XM_032382727.2","XR_007308206.1","XM_048723477.1","XM_048728488.1","XM_048721863.1","XM_048728498.1","XM_001630411.3","XM_032364762.2","XM_032383259.2","XR_007314072.1","XM_001635222.3","XM_001640140.3","XM_048732033.1","XM_032379829.2","XM_032382011.2","XM_032378053.2","XM_032376414.2","XM_032383977.2","XM_048731564.1","XM_001625226.3","XM_032375101.2","XM_001636428.3","XM_032381131.2","XM_001635264.3","XM_048720611.1","XM_001637279.3","XM_032381281.2","XM_001627857.3","XM_048720649.1","XM_032375534.2","XM_048723167.1","XM_032378782.2","XM_032379717.2","XM_032381920.2","XM_032364001.2","XM_048727062.1","XM_032379946.2","XM_001640047.3","XM_001622639.3","XM_048720311.1","XM_032384120.2","XM_001637865.3","XM_001641466.3","XM_048721893.1") #Reduced 

#make into a data frame and change name of column 
sig_degs_df<-data.frame(sig_degs)
str(sig_degs_df)
names(sig_degs_df)[1] <- c("Stella_seqid")

sig_degs_df <- read.csv("~/Library/CloudStorage/GoogleDrive-sbirch1@charlotte.edu/My Drive/Postdoc_work/Mesocosom_exps/Mesocosom_Summer_2022/Transcriptomics/EdgeR_MAPPED_4-29-24/deg heatmaps/ME_degs.txt-R_string copy.csv", sep="")
str(sig_degs_df)
#make data characters not factor so they can be compared
sig_degs_df$Stella_seqid<- as.character(sig_degs_df$Stella_seqid)
str(sig_degs_df)

#homohydra.map<- read.table("Hs_detect_light_stim_genes_R_output.txt",header=T, stringsAsFactors=F) #Full
homohydra.map<- read.table("full_uk_transcriptome_headers.txt",header=T, stringsAsFactors=F) #Reduced
str(homohydra.map) #should be data frame

library(tidyr)
library(dplyr)

#join sig degs to homohydra.map
str(homohydra.map)
sig_set<-inner_join(homohydra.map[,1:2],sig_degs_df,by="Stella_seqid")
nrow(sig_set)
length(unique(sig_set$Stella_seqid))

homohydra.map<-sig_set
keep<-sig_set$Stella_seqid

#select just TF genes
TFcts <-normcounts[keep,]
nrow(TFcts) #25
sort(rownames(TFcts)) #25

#construct heatmap of genes in map across devs
#using row as color break 
rgb.palette <- colorRampPalette(c( "#5DADE2","black","#F4D03F"),space = "Lab")
#specifying color break 
#rgb.palette2 <- colorRampPalette(c( "#5DADE2","black","#F4D03F"))(n=299) #,space = "Lab"
#col_breaks = c(seq(0,33,length=100), seq(33.01,66,length=100), seq(66.01,100,length=100))

library(gplots)

#using row as color break
graphics.off() 
par(mfrow = c(1,1))
heatmap.2(TFcts,  scale="row", labRow = homohydra.map$header,  Colv=F,  trace="none",  dendrogram="row",  key=T,  col=rgb.palette(120),  density.info=NULL, margins=c(5, 11), lmat=rbind(4:3, 2:1),  lhei=c(1, 30),  lwid=c(1, 1),cexRow = 0.80);
graphics.off() 
par(mfrow = c(1,1))
heatmap.2(TFcts,  scale="row", labRow = homohydra.map$header, Colv=F,  trace="none", col=rgb.palette, ColSideColors=rep(c('#83C5BE','#006D77'), each=3),dendrogram="row")
graphics.off() 
par(mfrow = c(1,1))


# take average (normalized) counts for each dev stage
tcts<-t(TFcts); dev<-factor(c(3,3,3,4,4,4));dat<-cbind(dev, tcts)
dev.means<-by(dat, dev, function(x) colMeans(x[, -1]))
dev.means2<-as.data.frame(t(sapply(dev.means, I)))

#heatmaps using row as color break
heatmap.2(t(dev.means2), labRow = homohydra.map$header, scale="row", Colv=F, trace="none", dendrogram="row", key=T, col=rgb.palette(120), density.info=NULL,  margins=c(2.5, 10.75),  lmat=rbind(4:3, 2:1),  lhei=c(1, 30),  lwid=c(1, 1), cexRow = 0.80);
heatmap.2(t(dev.means2), cellnote= round(t(dev.means2), digits = 3), notecol = "white", labRow = homohydra.map$header, scale="row", Colv=F, trace="none", dendrogram="row", key=F, col=rgb.palette(120), density.info=NULL)


          
          
          
          
          
#Running on each genotype
## NS T0 vs T14 

txi.salmon<- tximport(c(dev1files,dev2files), type = "salmon", txOut=T)
cts <- txi.salmon$counts
head(cts)
b<-DGEList(counts=cts, group=c(1,1,1,2,2,2))
normcounts<-cpm(b, normalized.lib.sizes=T) #normalized library counts for all genes


#make a matrix of significant degs from this data set - run prep_sig_DEGs_for_heatmap.py in terminal to get sig degs formated in this way
sig_degs<-c("XM_048722159.1","XM_048722160.1","XM_048725419.1","XM_032371448.2","XM_032365057.2","XM_048729797.1","XM_001631917.3","XM_032381798.2","XM_048733772.1","XR_007306253.1","XM_032384595.2","XM_048733739.1","XM_048721139.1","XM_001631422.3","XR_004296111.2","XM_032364702.2","XM_048722170.1","XM_032380597.2","XR_007314237.1","XM_001626916.3","XM_032381817.2","XR_007307279.1","XM_001636932.3","XM_048731924.1","XM_001628740.3","XM_048730410.1","XM_032381836.2","XM_032377385.2","XM_032375894.2","XM_048725998.1","XM_001636779.3","XM_048725113.1","XR_007312608.1","XM_048724181.1","XM_032373961.2","XM_048720646.1","XM_032364677.2","XM_048724444.1","XM_001639370.3","XM_001637574.3","XM_032385159.2","XM_032363883.2","XM_048732559.1","XR_007306819.1","XM_032384412.2","XM_001624583.3","XM_048720444.1","XM_048722169.1","XM_048734432.1","XM_032386927.2","XR_007306442.1","XM_032378637.2","XM_032362243.2","XM_032365054.2","XM_032367285.2","XM_032385537.2","XM_048727755.1","XM_048725092.1","XM_048720966.1","XM_032377986.2","XM_032366183.2","XR_007307164.1","XM_001633601.3","XM_032367765.2","XM_048721466.1","XM_032381131.2","XM_048726696.1","XM_032385144.2","XM_001623270.3","XM_048732914.1","XM_032387420.2","XM_001636816.3","XM_032364675.2","XM_048720927.1","XR_007306627.1","XM_048732544.1","XM_048725404.1","XM_048729734.1","XM_001630068.3","XM_048733318.1","XR_007313759.1","XM_032383008.2","XM_032386693.2","XM_001641929.3","XM_048729371.1","XM_048732270.1","XM_048727213.1","XM_001625226.3","XM_048721268.1","XM_048723557.1","XM_048723628.1","XM_048727409.1","XM_032372897.2","XM_048725085.1","XM_048719902.1","XM_048720462.1","XM_032373970.2","XM_032369016.2","XM_001622287.3","XM_048724548.1","XM_048734086.1","XM_048720799.1","XM_032380321.2","XM_001641807.3","XM_048722917.1","XM_032366022.2","XM_048720920.1","XR_007308206.1","XM_001619797.3","XM_048724741.1","XM_001626671.3","XM_032366283.2","XM_048723498.1","XM_048723494.1","XM_001630411.3","XM_048733738.1","XM_032384093.2","XM_032386726.2","XM_001640910.3","XM_032380381.2","XM_048732713.1","XM_032382924.2","XM_048724186.1","XM_032366282.2","XM_032375406.2","XM_032369459.2","XM_048721305.1","XM_048725109.1","XM_032381933.2","XM_032383013.2","XM_001628020.3","XM_032363196.2","XM_032373875.2","XM_001633335.3","XM_032373529.2","XR_007307206.1","XM_032385462.2","XM_048734670.1","XM_048728026.1","XM_032375362.2","XM_032380975.2","XM_001626761.3","XR_004295687.2","XM_048731319.1","XM_048721809.1","XM_048730053.1","XM_048728555.1","XM_001634477.3","XM_032363960.2","XM_001630315.3","XM_048730433.1","XM_001627064.3","XM_048729796.1","XM_048734697.1","XM_048724737.1","XM_001630535.3","XM_048732936.1","XM_048723705.1","XM_048731204.1","XM_048724165.1","XM_032377758.1","XM_048725000.1","XM_001635555.3","XM_048733179.1","XR_007312590.1","XR_004295015.2","XM_048726312.1","XM_032376510.2","XM_048728572.1","XM_048726741.1","XM_032380122.2","XM_001627599.3","XM_048726510.1","XM_001627279.3","XM_048725767.1","XM_032371555.2","XM_032383978.2","XM_032370287.2","XR_007312399.1","XM_048726927.1","XR_007305249.1","XR_007312408.1","XM_032370281.2","XM_001641530.3","XM_048721957.1","XM_048733312.1","XM_048728976.1","XM_001633033.3","XM_048728792.1","XM_048723326.1","XM_001636071.3","XR_007306849.1","XM_048725562.1","XR_007312194.1","XM_048733990.1","XM_032365107.2","XM_048728422.1","XM_001633659.3","XR_007307210.1","XM_032378741.2","XM_048728827.1","XM_032384129.2","XM_048721533.1","XM_048722692.1","XM_032373921.2","XR_007312578.1","XM_048726161.1","XM_048723373.1","XM_032370601.2","XM_048731784.1","XM_048725725.1","XM_001626638.3","XM_048725647.1","XM_032364294.2","XM_048728999.1","XM_032364145.2","XM_048722628.1","XM_048726054.1","XR_004293910.2","XM_048722503.1","XM_032372059.2","XM_032364263.2","XM_032363258.2","XM_032363587.2","XM_032377074.2","XM_048720696.1","XM_001633762.3","XM_032381139.2","XM_001630669.3","XM_001639511.3","XM_001625696.3","XM_048731718.1","XM_001641466.3","XM_048725306.1","XM_048726448.1","XM_032363157.2","XM_001638917.3","XR_004292988.2","XR_007312561.1","XM_048723473.1","XM_048731735.1","XM_048723718.1","XM_032375470.2","XM_001636575.3","XM_032364888.2","XM_032361894.2","XM_001632149.3","XM_032376414.2","XM_048720321.1","XM_048721243.1","XR_004295126.2","XM_032366219.2","XM_032385756.2","XM_048728488.1","XM_001638165.3","XM_001641640.3","XM_032386904.2","XM_032373182.2","XM_032379789.2","XR_007306602.1","XM_001638034.3","XM_048731707.1","XM_048728359.1","XM_048726055.1","XM_048733807.1","XM_048733697.1","XM_048723918.1","XR_004290963.2","XM_048730727.1","XM_048727207.1","XM_048723007.1","XM_001627857.3","XM_001626679.3","XM_048731800.1","XM_048732164.1","XM_048726662.1","XM_048723116.1","XM_032364009.2","XM_001636283.3","XM_001640140.3","XM_032386394.2","XM_048724215.1","XM_001627285.3","XM_032367463.2","XM_001636482.3","XM_001627402.3","XM_032369557.2","XR_004292730.1","XM_001627998.3","XM_001629087.2","XM_048727934.1","XM_032363359.2","XM_032383991.2","XM_048724214.1","XM_048720477.1","XM_032375246.2","XM_048729736.1","XM_001638478.3","XM_032366541.2","XM_048720915.1","XM_032381347.2","XM_032382290.2","XM_032362676.2","XM_001633290.3","XM_032377004.2","XM_048726076.1","XM_032380490.2","XM_001637837.3","XM_032380618.2","XM_048724447.1","XM_048720624.1","XM_048734600.1","XM_048720800.1","XM_001631950.2","XM_048722756.1","XM_032382789.2","XM_001623354.3","XM_048722050.1","XM_048726232.1","XM_048721227.1","XM_032367579.2","XM_032384803.2","XM_048721953.1","XM_032380616.2","XM_048721497.1","XM_001641062.3","XM_001625668.3","XM_001636497.3","XM_032387003.2","XM_048727187.1","XM_048719783.1","XM_048719790.1","XM_032372999.2","XM_048726450.1","XR_007314261.1","XM_048729302.1","XM_032369021.2","XM_048726260.1","XM_048728648.1","XM_032377749.2","XM_048721326.1","XM_001640210.2","XR_007306801.1","XM_032381213.2","XM_048727529.1","XR_004291135.2","XM_001632539.3","XM_048723524.1","XM_001635334.3","XM_001641219.2","XM_001628684.3","XM_048728498.1","XR_007305690.1","XM_048725787.1","XM_048728540.1","XM_048722237.1","XM_032378053.2","XM_048724818.1","XM_001638061.3","XM_032384128.2","XM_048724703.1","XM_032365846.2","XM_032363302.2","XR_004295937.2","XM_032381281.2","XM_048732762.1","XM_048734755.1","XM_001639649.3","XM_001636428.3","XM_032370154.2","XM_048729552.1","XM_032370147.2","XM_048722347.1","XM_048732113.1","XM_001628583.3","XM_048732380.1","XM_032378265.2","XM_032362917.2","XM_048724198.1","XM_001632693.3","XM_048726144.1","XM_001628860.3","XM_048726358.1","XM_048724464.1","XM_001631878.3","XM_001629168.3","XM_032374003.2","XM_032377522.2","XM_032367284.2","XM_032377228.2","XM_001639626.3","XM_001630918.3","XM_048722252.1","XM_048720827.1","XM_032376622.2","XM_048727482.1","XM_032371552.2","XM_001638268.3","XM_048722397.1","XM_048730480.1","XM_001631769.3","XM_001636928.3","XM_048729222.1","XM_048728223.1","XM_032383959.2","XM_048733073.1","XM_001634109.3","XM_048732210.1","XM_001631435.3","XM_048733279.1","XM_048722583.1","XM_048725061.1","XM_032363026.2","XM_001629963.3","XM_048723773.1","XM_032375039.2","XM_001625853.3","XM_048733814.1","XM_048725563.1","XR_004296397.2","XM_032376256.2","XM_048726568.1","XM_048722623.1","XM_048727081.1","XM_032365481.2","XM_048732578.1","XM_032366562.2","XR_007311891.1","XM_032370334.2","XM_048723500.1","XM_032385660.2","XR_007306415.1","XM_001623803.3","XM_001638417.3","XM_048726694.1","XM_048721789.1","XR_007307075.1","XM_048721603.1","XM_032373987.2","XM_048734208.1","XM_048727208.1","XM_001639625.3","XM_048726770.1","XM_032373225.2","XM_001628994.2","XM_032364761.2","XM_001641857.3","XM_048721249.1","XM_032379717.2","XM_001640639.3","XM_032370869.2","XM_048727819.1","XM_048731921.1","XM_001635222.3","XM_001636241.3","XM_001631958.3","XM_048728308.1","XM_032369575.2","XM_048734161.1","XM_001632792.3","XM_032385249.2","XM_032383259.2","XM_001623828.3","XM_048727495.1","XM_032370803.2","XM_048729355.1","XM_048723737.1","XM_048729556.1","XM_032384471.2") #Reduced 

#make into a data frame and change name of column 
sig_degs_df<-data.frame(sig_degs)
str(sig_degs_df)
names(sig_degs_df)[1] <- c("Stella_seqid")

sig_degs_df <- read.csv("~/Library/CloudStorage/GoogleDrive-sbirch1@charlotte.edu/My Drive/Postdoc_work/Mesocosom_exps/Mesocosom_Summer_2022/Transcriptomics/EdgeR_MAPPED_4-29-24/deg heatmaps/NS_degs.txt-R_string copy.csv", sep="")
str(sig_degs_df)
#make data characters not factor so they can be compared
sig_degs_df$Stella_seqid<- as.character(sig_degs_df$Stella_seqid)
str(sig_degs_df)

#homohydra.map<- read.table("Hs_detect_light_stim_genes_R_output.txt",header=T, stringsAsFactors=F) #Full
homohydra.map<- read.table("full_uk_transcriptome_headers.txt",header=T, stringsAsFactors=F) #Reduced
str(homohydra.map) #should be data frame

library(tidyr)
library(dplyr)

#join sig degs to homohydra.map
str(homohydra.map)
sig_set<-inner_join(homohydra.map[,1:2],sig_degs_df,by="Stella_seqid")
nrow(sig_set)
length(unique(sig_set$Stella_seqid))


homohydra.map<-sig_set
keep<-sig_set$Stella_seqid

#select just TF genes
TFcts <-normcounts[keep,]
nrow(TFcts) #25
sort(rownames(TFcts)) #25

#construct heatmap of genes in map across devs
#using row as color break 
rgb.palette <- colorRampPalette(c( "#5DADE2","black","#F4D03F"),space = "Lab")
#specifying color break 
#rgb.palette2 <- colorRampPalette(c( "#5DADE2","black","#F4D03F"))(n=299) #,space = "Lab"
#col_breaks = c(seq(0,33,length=100), seq(33.01,66,length=100), seq(66.01,100,length=100))

library(gplots)

#pdf("BestTFhits all and averaged.pdf",height=11,width=8)  
#using row as color break
graphics.off() 
par(mfrow = c(1,1))
heatmap.2(TFcts,  scale="row", labRow = homohydra.map$header,  Colv=F,  trace="none",  dendrogram="row",  key=T,  col=rgb.palette(120),  density.info=NULL, margins=c(5, 11), lmat=rbind(4:3, 2:1),  lhei=c(1, 30),  lwid=c(1, 1),cexRow = 0.80);
graphics.off() 
par(mfrow = c(1,1))
heatmap.2(TFcts,  scale="row", labRow = homohydra.map$header, Colv=F,  trace="none", col=rgb.palette, ColSideColors=rep(c('#83C5BE','#006D77'), each=3),dendrogram="row")
graphics.off() 
par(mfrow = c(1,1))


# take average (normalized) counts for each dev stage
tcts<-t(TFcts); dev<-factor(c(1,1,1,2,2,2));dat<-cbind(dev, tcts)
dev.means<-by(dat, dev, function(x) colMeans(x[, -1]))
dev.means2<-as.data.frame(t(sapply(dev.means, I)))

#heatmaps using row as color break
heatmap.2(t(dev.means2), labRow = homohydra.map$header, scale="row", Colv=F, trace="none", dendrogram="row", key=T, col=rgb.palette(120), density.info=NULL,  margins=c(2.5, 10.75),  lmat=rbind(4:3, 2:1),  lhei=c(1, 30),  lwid=c(1, 1), cexRow = 0.80);
heatmap.2(t(dev.means2), cellnote= round(t(dev.means2), digits = 3), notecol = "white", labRow = homohydra.map$header, scale="row", Colv=F, trace="none", dendrogram="row", key=F, col=rgb.palette(120), density.info=NULL)
          
          
