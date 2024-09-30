
rm(list = ls())
set.seed(0)

library(ggplot2)
library(cowplot)
library(scales)
library(tidyr)
library(grid)
library(gridExtra)
library(factoextra)
library(dplyr)
library(data.table)
library(tidyverse)

source("/Users/juanmacias/Library/CloudStorage/Box-Box/LawsonLab/JuanMacias/Misc_Code/CS_Aesthetic.R")
source("/Users/juanmacias/Documents/GitHub/General_Code/R_Scripts/General/UsefulPlottingFunctions.R")

setwd("/Users/juanmacias/Library/CloudStorage/Box-Box/LawsonLab/JuanMacias/Postdoc/HPRC/HiC/CompartmentAnalysis_10kb/")

# Import data
AllCompartmentData<-read.delim(file = "All_ReSmoothed_eigenvector_10000_30.bed", header = FALSE)
colnames(AllCompartmentData)<-c("Chr","Start","Stop","Score","Sample","Genome")

# Change assembly names
AllCompartmentData$Genome[grep(pattern = ".mat", x = AllCompartmentData$Genome)]<-"Maternal"
AllCompartmentData$Genome[grep(pattern = ".pat", x = AllCompartmentData$Genome)]<-"Paternal"
AllCompartmentData$Score[AllCompartmentData$Score=="."]<-NA

# Replace . with NA
AllCompartmentData$Genome[AllCompartmentData$Genome=="GRCh38"]<-"hg38"

# Make numeric
AllCompartmentData$Score<-as.numeric(AllCompartmentData$Score)

# Make wide
AllCompartmentData.Wide<-spread(AllCompartmentData, Genome, Score)

# Remove any bins that do not have a value for all for genomes
AllCompartmentData.Wide<-subset(AllCompartmentData.Wide, !is.na(CHM13) & !is.na(hg38) & !is.na(Maternal) & !is.na(Paternal) )

# Filter out ChrX
AllCompartmentData.Wide<-subset(AllCompartmentData.Wide, Chr !="chrX")

# Order chromosome
AllCompartmentData.Wide$Chr<-factor(x = AllCompartmentData.Wide$Chr, levels = c("chr1","chr2","chr3","chr4","chr5","chr6","chr7","chr8","chr9","chr10","chr11","chr12","chr13","chr14","chr15","chr16","chr17","chr18","chr19","chr20","chr21","chr22"), ordered = TRUE)

# Plot HG00621 Maternal vs hg38
ExampleScatter_Chr18<-ggplot(subset(AllCompartmentData.Wide, Sample=="HG00621" & Chr=="chr18"), aes(x=hg38,y=Maternal))+
  geom_abline(slope = 1,intercept = 0,color="darkgrey")+
  geom_hline(yintercept = 0, color="orange")+
  geom_vline(xintercept = 0, color="orange")+
  geom_point(size=0.1)+
  CS.THEME+
  ggtitle("HG00621: Maternal vs hg38")+
  xlab("hg38 compartment score")+
  ylab("Maternal compartment score")+
  ggtitle("chr18")

ggsave2("ExampleScatter_Chr18.png", ExampleScatter_Chr18, units = "in", height = 3.5, width = 3.5)

ExampleScatter_Chr21<-ggplot(subset(AllCompartmentData.Wide, Sample=="HG00621" & Chr=="chr21"), aes(x=hg38,y=Maternal))+
  geom_abline(slope = 1,intercept = 0,color="darkgrey")+
  geom_hline(yintercept = 0, color="orange")+
  geom_vline(xintercept = 0, color="orange")+
  geom_point(size=0.1)+
  CS.THEME+
  ggtitle("HG00621: Maternal vs hg38")+
  xlab("hg38 compartment score")+
  ylab("Maternal compartment score")+
  ggtitle("chr21")

ggsave2("ExampleScatter_Chr21.png", ExampleScatter_Chr21, units = "in", height = 3.5, width = 3.5)

ExampleScatter_Chr22<-ggplot(subset(AllCompartmentData.Wide, Sample=="HG00621" & Chr=="chr22"), aes(x=hg38,y=Maternal))+
  geom_abline(slope = 1,intercept = 0,color="darkgrey")+
  geom_hline(yintercept = 0, color="orange")+
  geom_vline(xintercept = 0, color="orange")+
  geom_point(size=0.1)+
  CS.THEME+
  ggtitle("HG00621: Maternal vs hg38")+
  xlab("hg38 compartment score")+
  ylab("Maternal compartment score")+
  ggtitle("chr22")

ggsave2("ExampleScatter_Chr22.png", ExampleScatter_Chr22, units = "in", height = 3.5, width = 3.5)

ExampleScatter_Chr9<-ggplot(subset(AllCompartmentData.Wide, Sample=="HG00621" & Chr=="chr9"), aes(x=hg38,y=Maternal))+
  geom_abline(slope = 1,intercept = 0,color="darkgrey")+
  geom_hline(yintercept = 0, color="orange")+
  geom_vline(xintercept = 0, color="orange")+
  geom_point(size=0.1)+
  CS.THEME+
  ggtitle("HG00621: Maternal vs hg38")+
  xlab("hg38 compartment score")+
  ylab("Maternal compartment score")+
  ggtitle("chr9")

ggsave2("ExampleScatter_Chr9.png", ExampleScatter_Chr9, units = "in", height = 3.5, width = 3.5)




# Plot split scatter plots
CHM13.SplitByChr<-ggplot(AllCompartmentData.Wide, aes(x=hg38,y=CHM13))+
  geom_abline(slope = 1,intercept = 0,color="darkgrey")+
  geom_hline(yintercept = 0, color="orange")+
  geom_vline(xintercept = 0, color="orange")+
  geom_point(size=0.1)+
  facet_grid(Sample~Chr)+
  CS.THEME+
  ggtitle("CHM13 vs hg38")

Maternal.SplitByChr<-ggplot(AllCompartmentData.Wide, aes(x=hg38,y=Maternal))+
  geom_abline(slope = 1,intercept = 0,color="darkgrey")+
  geom_hline(yintercept = 0, color="orange")+
  geom_vline(xintercept = 0, color="orange")+
  geom_point(size=0.1)+
  facet_grid(Sample~Chr)+
  CS.THEME+
  ggtitle("Maternal vs hg38")

Paternal.SplitByChr<-ggplot(AllCompartmentData.Wide, aes(x=hg38,y=Paternal))+
  geom_abline(slope = 1,intercept = 0,color="darkgrey")+
  geom_hline(yintercept = 0, color="orange")+
  geom_vline(xintercept = 0, color="orange")+
  geom_point(size=0.1)+
  facet_grid(Sample~Chr)+
  CS.THEME+
  ggtitle("Paternal vs hg38")

SplitByChrPanel<-plot_grid(CHM13.SplitByChr,Maternal.SplitByChr,Paternal.SplitByChr,labels = "AUTO",nrow = 3)
ggsave2("SplitByChrPanel.png", SplitByChrPanel, units = "in", height = 27, width = 33)


# Calculate R^2s for continuous values
ListOfChromosomes<-unique(AllCompartmentData.Wide$Chr)

# Define function
calcRsquared<-function(DATAFRAME,SAMPLE,CHROMOSOME,REF1,REF2){
  SET<-subset(DATAFRAME,Chr==CHROMOSOME & Sample==SAMPLE)
  REF1SET<-data.frame(Ref1=SET[,colnames(SET) == REF1])
  REF2SET<-data.frame(Ref2=SET[,colnames(SET) == REF2])
  COMBINED<-cbind(REF1SET,REF2SET)
  Rsqd<-summary(lm(Ref1~Ref2, COMBINED))[8]$r.squared
  
  data.frame(
    Sample=SAMPLE,
    Chromosome=CHROMOSOME,
    Query=REF1,
    Reference=REF2,
    Rsqd=Rsqd
  )
  
}

ListOfRuns<-unique(AllCompartmentData.Wide[,c(1,4)])

ContinuousRsqd.dataframe<-rbind(
  lapply(1:nrow(ListOfRuns), function(INDEX)  calcRsquared(AllCompartmentData.Wide,ListOfRuns[INDEX,2],ListOfRuns[INDEX,1],"CHM13","hg38") ) %>% bind_rows(),
  lapply(1:nrow(ListOfRuns), function(INDEX)  calcRsquared(AllCompartmentData.Wide,ListOfRuns[INDEX,2],ListOfRuns[INDEX,1],"Maternal","hg38") ) %>% bind_rows(),
  lapply(1:nrow(ListOfRuns), function(INDEX)  calcRsquared(AllCompartmentData.Wide,ListOfRuns[INDEX,2],ListOfRuns[INDEX,1],"Paternal","hg38") ) %>% bind_rows()
)

ContinuousRsqdPlot<-ggplot(ContinuousRsqd.dataframe,aes(x=Query,y=Rsqd))+
  geom_violin()+
  geom_jitter(size=1,width = 0.1)+
  facet_grid(.~Sample)+
  theme(axis.text.x = element_text(angle = 45,hjust = 1))+
  CS.THEME+
  ylab(expression(paste("R"^2)))+
  scale_y_continuous(limits = c(0,1))


ggsave2(filename = "ContinuousRsqdPlot.png",plot = ContinuousRsqdPlot, units = "in", width = 8, height = 3)


# Calculate R^2s binarized

AllCompartmentData.Wide.Orig<-AllCompartmentData.Wide

AllCompartmentData.Wide$CHM13[AllCompartmentData.Wide$CHM13>0]<-1
AllCompartmentData.Wide$CHM13[AllCompartmentData.Wide$CHM13<0]<--1
AllCompartmentData.Wide$hg38[AllCompartmentData.Wide$hg38>0]<-1
AllCompartmentData.Wide$hg38[AllCompartmentData.Wide$hg38<0]<--1
AllCompartmentData.Wide$Maternal[AllCompartmentData.Wide$Maternal>0]<-1
AllCompartmentData.Wide$Maternal[AllCompartmentData.Wide$Maternal<0]<--1
AllCompartmentData.Wide$Paternal[AllCompartmentData.Wide$Paternal>0]<-1
AllCompartmentData.Wide$Paternal[AllCompartmentData.Wide$Paternal<0]<--1

BinarizedRsqd.dataframe<-rbind(
  lapply(1:nrow(ListOfRuns), function(INDEX)  calcRsquared(AllCompartmentData.Wide,ListOfRuns[INDEX,2],ListOfRuns[INDEX,1],"CHM13","hg38") ) %>% bind_rows(),
  lapply(1:nrow(ListOfRuns), function(INDEX)  calcRsquared(AllCompartmentData.Wide,ListOfRuns[INDEX,2],ListOfRuns[INDEX,1],"Maternal","hg38") ) %>% bind_rows(),
  lapply(1:nrow(ListOfRuns), function(INDEX)  calcRsquared(AllCompartmentData.Wide,ListOfRuns[INDEX,2],ListOfRuns[INDEX,1],"Paternal","hg38") ) %>% bind_rows()
)

BinarizedRsqdPlot<-ggplot(BinarizedRsqd.dataframe,aes(x=Query,y=Rsqd))+
  geom_violin()+
  geom_jitter(size=1,width = 0.1)+
  facet_grid(.~Sample)+
  theme(axis.text.x = element_text(angle = 45,hjust = 1))+
  CS.THEME+
  ylab(expression(paste("R"^2)))+
  scale_y_continuous(limits = c(0,1))

ggsave2(filename = "BinarizedRsqdPlot.png",plot = BinarizedRsqdPlot, units = "in", width = 8, height = 3)




# Calculate per bin difference

AllCompartmentData.Wide$CHM13.hg38<-((AllCompartmentData.Wide$CHM13+1)/(AllCompartmentData.Wide$hg38+1))-1
AllCompartmentData.Wide$Maternal.hg38<-((AllCompartmentData.Wide$Maternal+1)/(AllCompartmentData.Wide$hg38+1))-1
AllCompartmentData.Wide$Paternal.hg38<-((AllCompartmentData.Wide$Paternal+1)/(AllCompartmentData.Wide$hg38+1))-1

PositionDifference.CHM13<-ggplot(AllCompartmentData.Wide,aes(x=Start/1000000,y=CHM13.hg38*100))+
  geom_hline(yintercept = 0, color="red")+
  geom_line()+
  facet_grid(Sample~Chr)+
  xlab("Position (Mb)")+
  CS.THEME+
  ylab("% Changed")

PositionDifference.Maternal<-ggplot(AllCompartmentData.Wide,aes(x=Start/1000000,y=Maternal.hg38*100))+
  geom_hline(yintercept = 0, color="red")+
  geom_line()+
  facet_grid(Sample~Chr)+
  xlab("Position (Mb)")+
  CS.THEME+
  ylab("% Changed")

PositionDifference.Paternal<-ggplot(AllCompartmentData.Wide,aes(x=Start/1000000,y=Paternal.hg38*100))+
  geom_hline(yintercept = 0, color="red")+
  geom_line()+
  facet_grid(Sample~Chr)+
  xlab("Position (Mb)")+
  CS.THEME+
  ylab("% Changed")

PositioDiferencePanel<-plot_grid(PositionDifference.CHM13,PositionDifference.Maternal,PositionDifference.Paternal,nrow = 3, labels = "AUTO")


ggsave2("PositioDiferencePanel.png", PositioDiferencePanel, units = "in", height = 27, width = 40)


# I should re-make the above as circos plots


