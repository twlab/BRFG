
rm(list = ls())
set.seed(0)

library(ggplot2)
require(cowplot)
library(scales)
library(gridExtra)
library(grid)
library(dplyr)
library(tidyr)

source("/Users/juanmacias/Library/CloudStorage/Box-Box/LawsonLab/JuanMacias/Misc_Code/CS_Aesthetic.R")

setwd("/Users/juanmacias/Library/CloudStorage/Box-Box/LawsonLab/JuanMacias/Postdoc/HPRC/WGBS/")

PlacementData<-read.delim("/Users/juanmacias/Library/CloudStorage/Box-Box/LawsonLab/JuanMacias/Postdoc/HPRC/WGBS/Combined_Counts_placement.txt",header=FALSE,sep = "\t")
colnames(PlacementData)<-c("Counts","Code","File")


PlacementData$File<-as.character(PlacementData$File)

### Any NA entries?
sum(is.na(PlacementData[,1]))
### Filter them out...
PlacementData<-PlacementData[!is.na(PlacementData[,1]),]

PlacementData$Sample<-matrix(unlist(strsplit(PlacementData$File,split = "_")),nrow=length(PlacementData$File),byrow = TRUE)[,3]
PlacementData$Replicate<-matrix(unlist(strsplit(PlacementData$File,split = "_")),nrow=length(PlacementData$File),byrow = TRUE)[,4]
PlacementData$Query.Genome<-matrix(unlist(strsplit(PlacementData$File,split = "_")),nrow=length(PlacementData$File),byrow = TRUE)[,2]
PlacementData$Reference.Genome<-"hg38"

table(PlacementData$Query.Genome)

PlacementData$Query<-NA

PlacementData$Query[PlacementData$Code==1]<-"Unmapped"
PlacementData$Reference[PlacementData$Code==1]<-"Unmapped"

PlacementData$Query[PlacementData$Code==2]<-"Uniquely mapped"
PlacementData$Reference[PlacementData$Code==2]<-"Unmapped"

PlacementData$Query[PlacementData$Code==3]<-"multimapped"
PlacementData$Reference[PlacementData$Code==3]<-"Unmapped"

PlacementData$Query[PlacementData$Code==4]<-"Unmapped"
PlacementData$Reference[PlacementData$Code==4]<-"Uniquely mapped"

PlacementData$Query[PlacementData$Code==5]<-"Uniquely mapped"
PlacementData$Reference[PlacementData$Code==5]<-"Uniquely mapped"

PlacementData$Query[PlacementData$Code==6]<-"multimapped"
PlacementData$Reference[PlacementData$Code==6]<-"Uniquely mapped"

PlacementData$Query[PlacementData$Code==7]<-"Unmapped"
PlacementData$Reference[PlacementData$Code==7]<-"multimapped"

PlacementData$Query[PlacementData$Code==8]<-"Uniquely mapped"
PlacementData$Reference[PlacementData$Code==8]<-"multimapped"

PlacementData$Query[PlacementData$Code==9]<-"multimapped"
PlacementData$Reference[PlacementData$Code==9]<-"multimapped"


PlacementData$Query<-factor(PlacementData$Query,levels=c("Unmapped","Uniquely mapped","multimapped"))
PlacementData$Reference<-factor(PlacementData$Reference,levels=c("Unmapped","Uniquely mapped","multimapped"))




PlacementData$Percent<-NA

for (FILE in unique(PlacementData$File)) {
  PlacementData$Percent[PlacementData$File==FILE]<-PlacementData$Counts[PlacementData$File==FILE]/sum(PlacementData$Counts[PlacementData$File==FILE])
}

PlacementData$Query.Genome[grep("mat",PlacementData$Query.Genome)]<-"Maternal"
PlacementData$Query.Genome[grep("pat",PlacementData$Query.Genome)]<-"Paternal"



## add replicate labels
# Make sure these are correct!
PlacementData$Library<-PlacementData$Replicate

PlacementData[PlacementData$Replicate =="lib3" & PlacementData$Sample =="HG00621",]$Replicate<-"B"
PlacementData[PlacementData$Replicate =="lib2" & PlacementData$Sample =="HG00621",]$Replicate<-"A"

PlacementData[PlacementData$Replicate =="lib2" & PlacementData$Sample =="HG00741",]$Replicate<-"B"
PlacementData[PlacementData$Replicate =="lib1" & PlacementData$Sample =="HG00741",]$Replicate<-"A"

PlacementData[PlacementData$Replicate =="lib3" & PlacementData$Sample =="HG01952",]$Replicate<-"B"
PlacementData[PlacementData$Replicate =="lib1" & PlacementData$Sample =="HG01952",]$Replicate<-"A"

PlacementData[PlacementData$Replicate =="lib3" & PlacementData$Sample =="HG01978",]$Replicate<-"B"
PlacementData[PlacementData$Replicate =="lib1" & PlacementData$Sample =="HG01978",]$Replicate<-"A"

PlacementData[PlacementData$Replicate =="lib2" & PlacementData$Sample =="HG03516",]$Replicate<-"B"
PlacementData[PlacementData$Replicate =="lib1" & PlacementData$Sample =="HG03516",]$Replicate<-"A"


table(PlacementData$Sample)



subset(PlacementData,Sample=="HG00621" & Replicate=="A" & Query.Genome=="CHM13" & Reference.Genome=="hg38")




PlacementData$Agreement<-NA
PlacementData$Agreement[PlacementData$Query == PlacementData$Reference]<-"Agree"
PlacementData$Agreement[PlacementData$Query != PlacementData$Reference]<-"Disagree"

CodePlot.1<-ggplot( subset(PlacementData, Code==1),aes(x=Query.Genome,y=100*Percent,fill=Replicate))+
  CS.THEME+
  geom_bar(stat="identity", width=.75, position = "dodge")+
  facet_grid(~Sample)+
  CS.THEME+
  theme(axis.text.x = element_text(angle = 30,hjust=1))+
  scale_fill_manual(values = c(CS.ColorPalette[11,1],CS.ColorPalette[12,1]))+
  ylab("Percent of reads")+
  ggtitle("Reference: Unmapped, Query: Unmapped")

CodePlot.2<-ggplot( subset(PlacementData, Code==2),aes(x=Query.Genome,y=100*Percent,fill=Replicate))+
  CS.THEME+
  geom_bar(stat="identity", width=.75, position = "dodge")+
  facet_grid(~Sample)+
  CS.THEME+
  theme(axis.text.x = element_text(angle = 30,hjust=1))+
  scale_fill_manual(values = c(CS.ColorPalette[11,1],CS.ColorPalette[12,1]))+
  ylab("Percent of reads")+
  ggtitle("Reference: Unmapped, Query: Uniquely mapped")

CodePlot.3<-ggplot( subset(PlacementData, Code==3),aes(x=Query.Genome,y=100*Percent,fill=Replicate))+
  CS.THEME+
  geom_bar(stat="identity", width=.75, position = "dodge")+
  facet_grid(~Sample)+
  CS.THEME+
  theme(axis.text.x = element_text(angle = 30,hjust=1))+
  scale_fill_manual(values = c(CS.ColorPalette[11,1],CS.ColorPalette[12,1]))+
  ylab("Percent of reads")+
  ggtitle("Reference: Unmapped, Query: Multimapped")

CodePlot.4<-ggplot( subset(PlacementData, Code==4),aes(x=Query.Genome,y=100*Percent,fill=Replicate))+
  CS.THEME+
  geom_bar(stat="identity", width=.75, position = "dodge")+
  facet_grid(~Sample)+
  CS.THEME+
  theme(axis.text.x = element_text(angle = 30,hjust=1))+
  scale_fill_manual(values = c(CS.ColorPalette[11,1],CS.ColorPalette[12,1]))+
  ylab("Percent of reads")+
  ggtitle("Reference: Uniquely mapped, Query: Unmapped")

CodePlot.5<-ggplot( subset(PlacementData, Code==5),aes(x=Query.Genome,y=100*Percent,fill=Replicate))+
  CS.THEME+
  geom_bar(stat="identity", width=.75, position = "dodge")+
  facet_grid(~Sample)+
  CS.THEME+
  theme(axis.text.x = element_text(angle = 30,hjust=1))+
  scale_fill_manual(values = c(CS.ColorPalette[11,1],CS.ColorPalette[12,1]))+
  ylab("Percent of reads")+
  ggtitle("Reference: Uniquely mapped, Query: Uniquely mapped")

CodePlot.6<-ggplot( subset(PlacementData, Code==6),aes(x=Query.Genome,y=100*Percent,fill=Replicate))+
  CS.THEME+
  geom_bar(stat="identity", width=.75, position = "dodge")+
  facet_grid(~Sample)+
  CS.THEME+
  theme(axis.text.x = element_text(angle = 30,hjust=1))+
  scale_fill_manual(values = c(CS.ColorPalette[11,1],CS.ColorPalette[12,1]))+
  ylab("Percent of reads")+
  ggtitle("Reference: Uniquely mapped, Query: Multimapped")

CodePlot.7<-ggplot( subset(PlacementData, Code==7),aes(x=Query.Genome,y=100*Percent,fill=Replicate))+
  CS.THEME+
  geom_bar(stat="identity", width=.75, position = "dodge")+
  facet_grid(~Sample)+
  CS.THEME+
  theme(axis.text.x = element_text(angle = 30,hjust=1))+
  scale_fill_manual(values = c(CS.ColorPalette[11,1],CS.ColorPalette[12,1]))+
  ylab("Percent of reads")+
  ggtitle("Reference: Multimapped, Query: Unmapped")

CodePlot.8<-ggplot( subset(PlacementData, Code==8),aes(x=Query.Genome,y=100*Percent,fill=Replicate))+
  CS.THEME+
  geom_bar(stat="identity", width=.75, position = "dodge")+
  facet_grid(~Sample)+
  CS.THEME+
  theme(axis.text.x = element_text(angle = 30,hjust=1))+
  scale_fill_manual(values = c(CS.ColorPalette[11,1],CS.ColorPalette[12,1]))+
  ylab("Percent of reads")+
  ggtitle("Reference: Multimapped, Query: Uniquely mapped")

CodePlot.9<-ggplot( subset(PlacementData, Code==9),aes(x=Query.Genome,y=100*Percent,fill=Replicate))+
  CS.THEME+
  geom_bar(stat="identity", width=.75, position = "dodge")+
  facet_grid(~Sample)+
  CS.THEME+
  theme(axis.text.x = element_text(angle = 30,hjust=1))+
  scale_fill_manual(values = c(CS.ColorPalette[11,1],CS.ColorPalette[12,1]))+
  ylab("Percent of reads")+
  ggtitle("Reference: Multimapped, Query: Multimapped")

FatePlacementBarplots.Panel<-plot_grid(CodePlot.1,CodePlot.2,CodePlot.3,CodePlot.4,CodePlot.5,CodePlot.6,CodePlot.7,CodePlot.8,CodePlot.9, ncol = 3, labels = c("S","T","U","V","W","X","Y","Z","AA"))
#FatePlacementBarplots.Panel<-plot_grid(CodePlot.1,CodePlot.2,CodePlot.4,CodePlot.5, ncol = 2, labels = "AUTO")


ggsave2("FatePlacementBarplots.Panel.png",plot = FatePlacementBarplots.Panel, units = "in", width = 27, height = 10)

calculateMarginal.Query<-function(QUERY.GENOME,SAMPLE,REPLICATE){
  data.frame(
    Sample=SAMPLE,
    Replicate=REPLICATE,
    Genome=QUERY.GENOME,
    Unmapped=sum(subset(PlacementData,Query.Genome==QUERY.GENOME & Sample ==SAMPLE & Replicate==REPLICATE & Query=="Unmapped")$Percent),
    UniquelyMapped=sum(subset(PlacementData,Query.Genome==QUERY.GENOME & Sample ==SAMPLE & Replicate==REPLICATE & Query=="Uniquely mapped")$Percent),
    Multimapped=sum(subset(PlacementData,Query.Genome==QUERY.GENOME & Sample ==SAMPLE & Replicate==REPLICATE & Query=="multimapped")$Percent)
  )
}

sum(
  subset(PlacementData, Query.Genome=="CHM13" & Sample =="HG00621" & Replicate=="A" & Query=="Unmapped")$Percent
  )

calculateMarginal.Reference<-function(REFERENCE.GENOME,SAMPLE,REPLICATE,QUERY.GENOME){
  data.frame(
    Sample=SAMPLE,
    Replicate=REPLICATE,
    Genome=REFERENCE.GENOME,
    Unmapped=sum(subset(PlacementData,Reference.Genome==REFERENCE.GENOME & Sample ==SAMPLE & Replicate==REPLICATE & Reference=="Unmapped" & Query.Genome==QUERY.GENOME)$Percent),
    UniquelyMapped=sum(subset(PlacementData,Reference.Genome==REFERENCE.GENOME & Sample ==SAMPLE & Replicate==REPLICATE & Reference=="Uniquely mapped" & Query.Genome==QUERY.GENOME)$Percent),
    Multimapped=sum(subset(PlacementData,Reference.Genome==REFERENCE.GENOME & Sample ==SAMPLE & Replicate==REPLICATE & Reference=="multimapped" & Query.Genome==QUERY.GENOME)$Percent)
  )
}





MarginalsDF<-lapply(unique(PlacementData$Query.Genome),
       function(QUERY.GENOME) lapply(
                                  unique(PlacementData$Sample),
                                  function(SAMPLE) lapply(
                                                      unique(PlacementData$Replicate),
                                                      function(REPLICATE) calculateMarginal.Query(QUERY.GENOME,SAMPLE,REPLICATE)
                                                        ) %>% bind_rows()
                              ) %>% bind_rows()
    ) %>% bind_rows()


MarginalsDF<-rbind(
  MarginalsDF,
  lapply(unique(PlacementData$Reference.Genome),
         function(REFERENCE.GENOME) lapply(
           unique(PlacementData$Sample),
           function(SAMPLE) lapply(
             unique(PlacementData$Replicate),
             function(REPLICATE) calculateMarginal.Reference(REFERENCE.GENOME,SAMPLE,REPLICATE,"CHM13")
           ) %>% bind_rows()
         ) %>% bind_rows()
  ) %>% bind_rows()
)





MarginalsDF.Long<-gather(MarginalsDF,Placement,Percent,Unmapped:Multimapped)

MarginalsDF.Long$Genome<-as.character(MarginalsDF.Long$Genome)


MarginalsPlot<-ggplot(MarginalsDF.Long,aes(x=Genome,y=100*Percent,fill=Placement))+
  geom_bar(stat = 'identity')+
  facet_grid(Replicate~Sample)+
  ylab("% of mapped reads")+
  theme(axis.text.x = element_text(angle = 30,hjust=1))+
  scale_fill_viridis(discrete = TRUE, option = "cividis")+
  theme(legend.position = "bottom")

ggsave2(filename = "MarginalsPlacementPlot.png",plot = MarginalsPlot,units = "in", height = 4, width = 10)


DP.CHM13<-ggplot(subset(PlacementData,Query.Genome=="CHM13"),aes(x=Query,y=Reference))+
  geom_point(aes(size = Percent))+
  CS.THEME+
  theme(axis.text.x = element_text(angle = 30,hjust=1))+
  facet_grid(Replicate~Sample)+
  ggtitle("CHM13")

DP.Paternal<-ggplot(subset(PlacementData,Query.Genome=="Paternal"),aes(x=Query,y=Reference))+
  geom_point(aes(size = Percent))+
  CS.THEME+
  theme(axis.text.x = element_text(angle = 30,hjust=1))+
  facet_grid(Replicate~Sample)+
  ggtitle("Paternal")

DP.Maternal<-ggplot(subset(PlacementData,Query.Genome=="Maternal"),aes(x=Query,y=Reference))+
  geom_point(aes(size = Percent))+
  CS.THEME+
  theme(axis.text.x = element_text(angle = 30,hjust=1))+
  facet_grid(Replicate~Sample)+
  ggtitle("Maternal")

DP.ALL<-plot_grid(DP.CHM13,DP.Paternal,DP.Maternal,nrow=3,labels = "AUTO")

ggsave2("DP.ALL.png",plot = DP.ALL, units = "in", width = 9, height = 11)





####


RatioData<-lapply(
  unique(PlacementData$File),
  function(FILE) data.frame(
    file=FILE,
    #PlacementRatio=sum(PlacementData$Counts[PlacementData$File==FILE & PlacementData$Agreement=="Disagree"])/sum(PlacementData$Counts[PlacementData$File==FILE & PlacementData$Agreement=="Agree"]),
    #Unique.Multi.Ratio=PlacementData$Percent[PlacementData$File==FILE & PlacementData$Code==6]/PlacementData$Percent[PlacementData$File==FILE & PlacementData$Code==8],
    Disagreement.Percent=100*sum(PlacementData$Counts[PlacementData$File==FILE & PlacementData$Agreement=="Disagree"])/sum(PlacementData$Counts[PlacementData$File==FILE])
  )
  ) %>% bind_rows()




RatioData$file<-as.character(RatioData$file)

RatioData$Sample<-unlist(lapply(RatioData$file,function(FILE) strsplit(FILE,split = "_")[[1]][[3]] ))
RatioData$Replicate<-unlist(lapply(RatioData$file,function(FILE) strsplit(FILE,split = "_")[[1]][[4]] ))
RatioData$Query<-unlist(lapply(RatioData$file,function(FILE) strsplit(FILE,split = "_")[[1]][[2]] ))

RatioData$Reference<-"hg38"
RatioData$Query[grep("mat",RatioData$Query)]<-"Maternal"
RatioData$Query[grep("pat",RatioData$Query)]<-"Paternal"


RatioData[RatioData$Replicate =="lib3" & RatioData$Sample =="HG00621",]$Replicate<-"B"
RatioData[RatioData$Replicate =="lib2" & RatioData$Sample =="HG00621",]$Replicate<-"A"

RatioData[RatioData$Replicate =="lib2" & RatioData$Sample =="HG00741",]$Replicate<-"B"
RatioData[RatioData$Replicate =="lib1" & RatioData$Sample =="HG00741",]$Replicate<-"A"

RatioData[RatioData$Replicate =="lib3" & RatioData$Sample =="HG01952",]$Replicate<-"B"
RatioData[RatioData$Replicate =="lib1" & RatioData$Sample =="HG01952",]$Replicate<-"A"

RatioData[RatioData$Replicate =="lib3" & RatioData$Sample =="HG01978",]$Replicate<-"B"
RatioData[RatioData$Replicate =="lib1" & RatioData$Sample =="HG01978",]$Replicate<-"A"

RatioData[RatioData$Replicate =="lib2" & RatioData$Sample =="HG03516",]$Replicate<-"B"
RatioData[RatioData$Replicate =="lib1" & RatioData$Sample =="HG03516",]$Replicate<-"A"

mean(RatioData[,2])
sd(RatioData[,2])

FateChangePlot<-ggplot(RatioData,aes(x=Query,y=Disagreement.Percent,fill=Replicate))+
  CS.THEME+
  geom_bar(stat="identity", width=.75, position = "dodge")+
  facet_grid(.~Sample)+
  scale_y_continuous(expand = c(0,0),limits = c(0,25))+
  ylab("%Reads changing fates")+
  scale_fill_manual(values = c(CS.ColorPalette[11,1],CS.ColorPalette[12,1]))+
  theme(legend.position = "bottom")+
  ggtitle("WGBS - changes in read fate")


ggsave2("FateChangePlot.png",plot = plot_grid(FateChangePlot,ncol=1,labels = "C"), units = "in",width = 19,height = 4)














