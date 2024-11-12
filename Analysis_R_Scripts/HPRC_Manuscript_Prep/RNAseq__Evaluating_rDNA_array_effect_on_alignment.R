


rm(list = ls())
set.seed(0)

library(ggplot2)
require(cowplot)
library(scales)
library(gridExtra)
library(grid)
library(dplyr)
library(tidyr)
library(lsr)
library(magick)
library(pdftools)
library(tidyverse)

source("/Users/juanmacias/Documents/GitHub/JuanMacias_General_Code/R_Scripts/General/UsefulPlottingFunctions.R")
source("/Users/juanmacias/Library/CloudStorage/Box-Box/LawsonLab/JuanMacias/Misc_Code/CS_Aesthetic.R")

InputDirectory<-"/Users/juanmacias/Library/CloudStorage/Box-Box/LawsonLab/JuanMacias/Postdoc/HPRC/RNAseq/rDNA_Array/"
OutputDirectory<-"/Users/juanmacias/Library/CloudStorage/Box-Box/LawsonLab/JuanMacias/Postdoc/Manuscripts/BenchmarkingPaper/Supplementary_Figure_Files/RNAseq/"
setwd(InputDirectory)

### Import data
AlignmentResults<-read.delim(file = "AlignerResults.txt", header = TRUE, sep = "\t")

ExpansionData<-subset(AlignmentResults, is.na(uniquely.mapped_Subjunc) )
AlignmentResults<-subset(AlignmentResults, !is.na(uniquely.mapped_Subjunc) )[,-c(15:18)]

AlignmentResults.Long<-gather(AlignmentResults, Class, Percent, uniquely.mapped_STAR:unmapped_BBMap)
AlignmentResults.Long$Placement<-unlist(lapply(AlignmentResults.Long$Class, function(CLASS) str_split(CLASS,"_")[[1]][1] ))
AlignmentResults.Long$Aligner<-unlist(lapply(AlignmentResults.Long$Class, function(CLASS) str_split(CLASS,"_")[[1]][2] ))

AlignmentResults.Long$Placement[AlignmentResults.Long$Placement=="uniquely.mapped"]<-"UniquelyMapped"
AlignmentResults.Long$Placement[AlignmentResults.Long$Placement=="multi.mapped"]<-"Multimapped"
AlignmentResults.Long$Placement[AlignmentResults.Long$Placement=="unmapped"]<-"Unmapped"


setwd(OutputDirectory)


### Generate marginals plot
MarginalsPlot<-ggplot(AlignmentResults.Long,aes(x=Aligner,y=100*Percent,fill=Placement))+
  geom_bar(stat = 'identity')+
  facet_grid(Genome~Sample*Replicate)+
  ylab("% of mapped reads")+
  theme(axis.text.x = element_text(angle = 30,hjust=1))+
  scale_fill_viridis(discrete = TRUE, option = "cividis")+
  theme(legend.position = "bottom")

ggsave2(filename = "MarginalsPlot.tmp.png", plot = plot_grid(MarginalsPlot,labels = "A"), units = "in", width = 12, height = 4)



### Generate effect of array expansion
ExpansionData<-ExpansionData[,-c(6:7,9:13)]

colnames(ExpansionData)[5]<-"UniquelyMapped"
colnames(ExpansionData)[6]<-"Multimapped"
colnames(ExpansionData)[8]<-"Too.many.loci"
colnames(ExpansionData)[9]<-"Too.many.mismatches"
colnames(ExpansionData)[10]<-"Too.short"
colnames(ExpansionData)[11]<-"Other"

ExpansionData$Genome[ExpansionData$Genome=="hg38.expanded"]<-"hg38\nexpanded"

UniquelyMapped<-ggplot(ExpansionData, aes(x=Genome,y=UniquelyMapped*100,fill=Replicate))+
  geom_bar(stat="identity", position = position_dodge())+
  CS.THEME+
  scale_fill_manual(values = c(CS.ColorPalette[11,1],CS.ColorPalette[12,1]))+
  #theme(axis.text.x = element_text(angle = 30,hjust=1))+
  ylab("%Uniquely mapped\n")+
  facet_grid(~Sample)+
  scale_y_continuous(limits = c(0,75), expand = c(0,0))

MultiMapped<-ggplot(ExpansionData, aes(x=Genome,y=Multimapped*100,fill=Replicate))+
  geom_bar(stat="identity", position = position_dodge())+
  CS.THEME+
  scale_fill_manual(values = c(CS.ColorPalette[11,1],CS.ColorPalette[12,1]))+
  #theme(axis.text.x = element_text(angle = 30,hjust=1))+
  ylab("%Multi-mapped mapped\n")+
  facet_grid(~Sample)+
  scale_y_continuous(limits = c(0,75), expand = c(0,0))

TooManyLoci<-ggplot(ExpansionData, aes(x=Genome,y=Too.many.loci*100,fill=Replicate))+
  geom_bar(stat="identity", position = position_dodge())+
  CS.THEME+
  scale_fill_manual(values = c(CS.ColorPalette[11,1],CS.ColorPalette[12,1]))+
  #theme(axis.text.x = element_text(angle = 30,hjust=1))+
  ylab("%Unmapped:\nToo many loci")+
  facet_grid(~Sample)+
  scale_y_continuous(limits = c(0,75), expand = c(0,0))

TooManyMismatches<-ggplot(ExpansionData, aes(x=Genome,y=Too.many.mismatches*100,fill=Replicate))+
  geom_bar(stat="identity", position = position_dodge())+
  CS.THEME+
  scale_fill_manual(values = c(CS.ColorPalette[11,1],CS.ColorPalette[12,1]))+
  #theme(axis.text.x = element_text(angle = 30,hjust=1))+
  ylab("%Unmapped:\nToo many mismatches")+
  facet_grid(~Sample)+
  scale_y_continuous(limits = c(0,75), expand = c(0,0))

Tooshort<-ggplot(ExpansionData, aes(x=Genome,y=Too.short*100,fill=Replicate))+
  geom_bar(stat="identity", position = position_dodge())+
  CS.THEME+
  scale_fill_manual(values = c(CS.ColorPalette[11,1],CS.ColorPalette[12,1]))+
  #theme(axis.text.x = element_text(angle = 30,hjust=1))+
  ylab("%Unmapped:\nToo short")+
  facet_grid(~Sample)+
  scale_y_continuous(limits = c(0,75), expand = c(0,0))

Other<-ggplot(ExpansionData, aes(x=Genome,y=Other*100,fill=Replicate))+
  geom_bar(stat="identity", position = position_dodge())+
  CS.THEME+
  scale_fill_manual(values = c(CS.ColorPalette[11,1],CS.ColorPalette[12,1]))+
  #theme(axis.text.x = element_text(angle = 30,hjust=1))+
  ylab("%Unmapped:\nOther")+
  facet_grid(~Sample)+
  scale_y_continuous(limits = c(0,75), expand = c(0,0))

ExpansionData_Panel<-plot_grid(UniquelyMapped,MultiMapped,TooManyLoci,TooManyMismatches,Tooshort,Other,labels = c("B","C","D","E","F","G"))

ggsave2(filename = "ExpansionData_Panel.tmp.png", plot = ExpansionData_Panel, units = "in", width = 15, height = 6)



### Array length

ArraySize<-data.frame(
  Genome=c("hg38","CHM13","HG00621 Maternal","HG00621 Paternal","HG00741 Maternal","HG00741 Paternal","HG01952 Maternal","HG01952 Paternal","HG01987 Maternal","HG01987 Paternal","HG03516 Maternal","HG03516 Paternal"),
  AlignmentLength=c(161,2500,239,70,423,227,165,179,185,133,423,0)
)

ArraySize<-ArraySize[order(ArraySize$AlignmentLength),]

PlotArrayLength<-ggplot(ArraySize,aes(y=AlignmentLength,x=Genome))+
  geom_bar(stat = "identity")+
  coord_flip()+
  CS.THEME+
  scale_y_continuous(limits = c(0,2600), expand = c(0,0))+
  ylab("Alignment length (Kb)")


### Array examples




setwd(InputDirectory)

### Plot parts
#CHM13_Alignment <- ggdraw() + draw_image("CHM13_vs_all_FP671120_alignment_sequence.png",scale = 0.9)
#HG01952Maternal_Alignment <- ggdraw() + draw_image("HG01952mat_vs_all_FP671120_alignment_sequence.png",scale = 0.9)
#HG00621Maternal_Alignment <- ggdraw() + draw_image("HG00621mat_vs_all_FP671120_alignment_sequence.png",scale = 0.9)

CHM13_Alignment <- ggdraw() + draw_image("FP671120_FP236383_vs_CHM13.svg",scale = 0.9)
HG01952Maternal_Alignment <- ggdraw() + draw_image("FP671120_FP236383_v_HG01952_maternal.svg",scale = 0.9)
HG00621Maternal_Alignment <- ggdraw() + draw_image("FP671120_FP236383_v_HG00621_maternal.svg",scale = 0.9)


setwd(OutputDirectory)

ArrayLengthPanel<-plot_grid(PlotArrayLength,CHM13_Alignment,HG01952Maternal_Alignment,HG00621Maternal_Alignment,labels = c("H","I","J","K") )

ggsave2(filename = "ArrayLengthPanel.tmp.png",plot = ArrayLengthPanel, units = "in", height = 8, width = 20)


PlotTop <- ggdraw() + draw_image("MarginalsPlot.tmp.png",scale = 1)
PlotMiddle <- ggdraw() + draw_image("ExpansionData_Panel.tmp.png",scale = 1)
PlotBottom <- ggdraw() + draw_image("ArrayLengthPanel.tmp.png",scale = 1)

rDNA_Array_Panel<-plot_grid(PlotTop,PlotMiddle,PlotBottom,ncol=1)
ggsave2(filename = "rDNA_Array_Panel.png", plot = rDNA_Array_Panel, units = "in", height = 12, width = 10)

file.remove("MarginalsPlot.tmp.png")
file.remove("ExpansionData_Panel.tmp.png")
file.remove("ArrayLengthPanel.tmp.png")

