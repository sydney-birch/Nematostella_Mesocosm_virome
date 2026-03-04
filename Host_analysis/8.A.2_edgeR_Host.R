##EdgeR - Nematostella MAPPED READS- 7 locations, 2 time points, 3 reps 

setwd("~/Library/CloudStorage/GoogleDrive-sbirch1@charlotte.edu/My Drive/Postdoc_work/Mesocosom_exps/Mesocosom_Summer_2022/Transcriptomics/EdgeR_MAPPED_4-29-24")
#source("https://bioconductor.org/biocLite.R"); biocLite("tximport"); install.packages("readr"); biocLite("edgeR")    #<-----uncomment & run if any not installed yet

library(tximport); library(readr); library(edgeR); library(limma)

#in your cwd make a text file named libraries_to_stages.txt - this should have 2 columns - the first is the name
#of your sample and the second (separated by a space) is the group number it belongs to (so the stage)

#Make a dir in your cwd named mapping - place all of your output dirs for each sample that was outputted by salmon in this dir
dir <- getwd()
list.files()

#library info file - based on location and time point (14)
devstages<-read.table("libraries_to_stages_2.txt",header=F,row.names=1)
dev1<-rownames(devstages)[which(devstages$V2==1)]; dev1files <- file.path(dir, "mapping",dev1, "quant.sf"); names(dev1files)<-dev1
dev2 <-rownames(devstages)[which(devstages$V2==2)]; dev2files <- file.path(dir, "mapping",dev2, "quant.sf"); names(dev2files)<-dev2
dev3 <-rownames(devstages)[which(devstages$V2==3)]; dev3files <- file.path(dir, "mapping",dev3, "quant.sf"); names(dev3files)<-dev3
dev4 <-rownames(devstages)[which(devstages$V2==4)]; dev4files <- file.path(dir, "mapping",dev4, "quant.sf"); names(dev4files)<-dev4
dev5<-rownames(devstages)[which(devstages$V2==5)]; dev5files <- file.path(dir, "mapping",dev5, "quant.sf"); names(dev5files)<-dev5
dev6 <-rownames(devstages)[which(devstages$V2==6)]; dev6files <- file.path(dir, "mapping",dev6, "quant.sf"); names(dev6files)<-dev6
dev7 <-rownames(devstages)[which(devstages$V2==7)]; dev7files <- file.path(dir, "mapping",dev7, "quant.sf"); names(dev7files)<-dev7
dev8<-rownames(devstages)[which(devstages$V2==8)]; dev8files <- file.path(dir, "mapping",dev8, "quant.sf"); names(dev8files)<-dev8
dev9 <-rownames(devstages)[which(devstages$V2==9)]; dev9files <- file.path(dir, "mapping",dev9, "quant.sf"); names(dev9files)<-dev9
dev10 <-rownames(devstages)[which(devstages$V2==10)]; dev10files <- file.path(dir, "mapping",dev10, "quant.sf"); names(dev10files)<-dev10
dev11 <-rownames(devstages)[which(devstages$V2==11)]; dev11files <- file.path(dir, "mapping",dev11, "quant.sf"); names(dev11files)<-dev11
dev12<-rownames(devstages)[which(devstages$V2==12)]; dev12files <- file.path(dir, "mapping",dev12, "quant.sf"); names(dev12files)<-dev12
dev13 <-rownames(devstages)[which(devstages$V2==13)]; dev13files <- file.path(dir, "mapping",dev13, "quant.sf"); names(dev13files)<-dev13
dev14 <-rownames(devstages)[which(devstages$V2==14)]; dev14files <- file.path(dir, "mapping",dev14, "quant.sf"); names(dev14files)<-dev14

##pick colors for each library type
stage1col<-"#DE3C00"
stage2col <-"#FFA07A"
stage3col <-"#08548F"
stage4col <-"aquamarine3"
stage5col<-"#FFC125"
stage6col <-"#8B4500"
stage7col <-"#8B668B"
stage8col<-"#FFB6C1"
stage9col <-"#556B2F"
stage10col <-"#A2CD5A"
stage11col <-"#CDC1C5"
stage12col<-"#8B8386"
stage13col <-"#FF82AB"
stage14col <-"#8B475D"


################################
### Multidimensional scaling ###
################################

runEdgeRmds<-function(salmonquantfiles, groups, colors){
  ##  read in files with tximport
  txi.salmon<- tximport(salmonquantfiles, type = "salmon", txOut=T, dropInfReps=TRUE)
  cts <- txi.salmon$counts
  print(colSums(cts))
  
  keep <- rowSums(cpm(cts)>10.0) >= length(groups)
  cts<- cts[keep,]
  dim(cts)
  print(colSums(cts))
  
  group <- groups
  y <- DGEList(counts=cts ,group=group)
  y <- calcNormFactors(y)
  y <-estimateCommonDisp(y)
  y <- estimateTagwiseDisp(y, prior.n=1.79) # TagwiseDisp n-value should be close to: 50/(#samples - #groups) = 50/(42-14) = 50/28 =1.79
  
  
  plotMDS.DGEList(y , main = "MDS Plot for Count Data", labels = colnames(y),col=colors)
  legend("topright",legend=paste("Location ",levels(as.factor(groups))),text.col=colors[seq(1,length(colors),3)]) #change the last number to how many reps you have
  
}
#pdf(file="edgeRplots.pdf",width=8, height=11)
cols<-c(stage1col, stage2col, stage3col, stage4col,stage5col,stage6col,stage7col,stage8col, stage9col, stage10col, stage11col,stage12col,stage13col,stage14col)
runEdgeRmds(c(dev1files,dev2files, dev3files, dev4files,dev5files,dev6files,dev7files,dev8files,dev9files, dev10files, dev11files,dev12files,dev13files,dev14files), c(1,1,1,2,2,2,3,3,3,4,4,4,5,5,5,6,6,6,7,7,7,8,8,8,9,9,9,10,10,10,11,11,11,12,12,12,13,13,13,14,14,14),rep(cols,each=3)) #change the last number to how many groups you have




##############################################
### pairwise comparisons between libraries ###
##############################################

runEdgeRpwise<-function(salmonquantfiles, groups,colors){
  ##  read in files with tximport
  txi.salmon<- tximport(salmonquantfiles, type = "salmon", txOut=T, dropInfReps=TRUE)
  cts <- txi.salmon$counts
  print(colSums(cts))
  
  keep <- rowSums(cpm(cts)>10.0) >= length(groups)
  cts<- cts[keep,]
  dim(cts)
  print(colSums(cts))
  
  y <- DGEList(counts=cts ,group=groups)
  y <- calcNormFactors(y)
  
  y<-estimateCommonDisp(y)
  # TagwiseDisp n-value should be close to: 50/(#samples - #groups) = 50/(36-3) = 50/33 =1.515152
  y <- estimateTagwiseDisp(y)#prior.n=1.79
  
  group<-levels(as.factor(groups))		
  et<-exactTest(y, pair=c(group[1],group[2]))
  tab<-summary(de <- decideTestsDGE(et, p=0.01, adjust="BH"))
  n<-tab[1]+tab[3]
  detags <- rownames(y)[as.logical(de)]
  
  
  plotSmear(et, de.tags=detags, main="DGE Exact Test")
  abline(h = c(-2, 2), col = "blue")
  abline(h = c(-4, 4), col = "blue")
  
  plotMDS.DGEList(y , main = "MDS Plot for Count Data", labels = colnames(y), col=colors)
  textcol<-colors[seq(1,length(colors),3)]
  legend("center",legend=paste("Stage ",group),text.col= textcol,bty="n")
  
  plotBCV(y, main="BCV plot")
  
  meanVarPlot <- plotMeanVar(estimateCommonDisp(y) , show.raw.vars=TRUE,
                             show.tagwise.vars=TRUE,
                             show.binned.common.disp.vars=FALSE,
                             show.ave.raw.vars=FALSE , NBline=TRUE,
                             nbins=100,
                             pch=16,
                             xlab="Mean Expresion (Log10)",
                             ylab="Variance (Log10)",
                             main="Mean-Variance Plot")
  
  #positive FC: higher in group2 than 1 
  return(topTags(et, n=n))
}

options(max.print = 1000000000)

### T0 vs T14 ####

#NS
#evaluate stage1 vs stage2 :: positive logFC values mean high in stage 2 than stage 1
par(mfrow=c(2,2),oma=c(1,1,2,0))
s1vs2<-runEdgeRpwise(c(dev1files,dev2files), c(1,1,1,2,2,2),rep(c(stage1col, stage2col),each=3))
title("Stage NS T0 vs NS T14",outer=T)
nrow(s1vs2)
write.table(s1vs2, "NS_T0vT14_DEGs_edgeR_output.txt", sep="\t")

#ME
#evaluate stage3 vs stage4
par(mfrow=c(2,2),oma=c(1,1,2,0))
s3vs4<-runEdgeRpwise(c(dev3files,dev4files), c(3,3,3,4,4,4), rep(c(stage3col, stage4col),each=3))
title("Stage ME T0 vs ME T14",outer=T)
nrow(s3vs4)
write.table(s3vs4, "ME_T0vT14_DEGs_edgeR_output.txt", sep="\t")

#NH
#evaluate stage5 vs stage6
par(mfrow=c(2,2),oma=c(1,1,2,0))
s5vs6<-runEdgeRpwise(c(dev5files,dev6files), c(5,5,5,6,6,6), rep(c(stage5col, stage6col),each=3))
title("Stage NH T0 vs NH T14",outer=T)
nrow(s5vs6)
write.table(s5vs6, "NH_T0vT14_DEGs_edgeR_output.txt", sep="\t")

#Field
#evaluate stage7 vs stage8 :: positive logFC values mean high in stage 2 than stage 1
par(mfrow=c(2,2),oma=c(1,1,2,0))
s7vs8<-runEdgeRpwise(c(dev7files,dev8files), c(7,7,7,8,8,8),rep(c(stage7col, stage8col),each=3))
title("Stage Field T0 vs Field T14",outer=T)
nrow(s7vs8)
write.table(s7vs8, "FIELD_T0vT14_DEGs_edgeR_output.txt", sep="\t")

#MA
#evaluate stage9 vs stage10 :: positive logFC values mean high in stage 2 than stage 1
par(mfrow=c(2,2),oma=c(1,1,2,0))
s9vs10<-runEdgeRpwise(c(dev9files,dev10files), c(9,9,9,10,10,10),rep(c(stage9col, stage10col),each=3))
title("Stage MA T0 vs MA T14",outer=T)
nrow(s9vs10)
write.table(s9vs10, "MA_T0vT14_DEGs_edgeR_output.txt", sep="\t")

#SC
#evaluate stage11 vs stage12 :: positive logFC values mean high in stage 2 than stage 1
par(mfrow=c(2,2),oma=c(1,1,2,0))
s11vs12<-runEdgeRpwise(c(dev11files,dev12files), c(11,11,11,12,12,12),rep(c(stage11col, stage12col),each=3))
title("Stage SC T0 vs SC T14",outer=T)
nrow(s11vs12)
write.table(s11vs12, "SC_T0vT14_DEGs_edgeR_output.txt", sep="\t")

#FL
#evaluate stage13 vs stage14 :: positive logFC values mean high in stage 2 than stage 1
par(mfrow=c(2,2),oma=c(1,1,2,0))
s13vs14<-runEdgeRpwise(c(dev13files,dev14files), c(13,13,13,14,14,14),rep(c(stage13col, stage14col),each=3))
title("Stage FL T0 vs FL T14",outer=T)
nrow(s13vs14)
write.table(s13vs14, "FL_T0vT14_DEGs_edgeR_output.txt", sep="\t")

#Barplot
par(mfrow=c(1,1))
barplot(c(nrow(s1vs2), nrow(s3vs4),nrow(s5vs6),nrow(s7vs8),nrow(s9vs10),nrow(s11vs12),nrow(s13vs14)),names.arg=c("NS", "ME", "NH","Field", "MA", "SC", "FL"), main="Number of differentially expressed transcripts (T0 vs T14)")

