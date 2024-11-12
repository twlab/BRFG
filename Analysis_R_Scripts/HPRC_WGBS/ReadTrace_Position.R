


rm(list = ls())
set.seed(0)

library(ggplot2)
require(cowplot)
library(scales)
library(gridExtra)
library(grid)
library(dplyr)
library(tidyr)
library(networkD3)

source("/Users/juanmacias/Library/CloudStorage/Box-Box/LawsonLab/JuanMacias/Misc_Code/CS_Aesthetic.R")
setwd("/Users/juanmacias/Library/CloudStorage/Box-Box/LawsonLab/JuanMacias/Postdoc/HPRC/WGBS/ReadPositionTrace/")

importChromosomePlacementFile<-function(FILENUMBER){
  FILELIST<-dir("ProcessedFiles/")[grep("ChromosomePlaceMentTable.txt", dir("ProcessedFiles/"))]
  
  SAMPLE<-strsplit(FILELIST[FILENUMBER], split = "_")[[1]][6]
  
  REPLICATE<-strsplit(FILELIST[FILENUMBER], split = "_")[[1]][7]
  
  GENOME<-strsplit(FILELIST[FILENUMBER], split = "_")[[1]][8]
  REFGENOME<-strsplit(FILELIST[FILENUMBER], split = "_")[[1]][9]
  #print(SAMPLE)
  FILENAME<-paste( c("ProcessedFiles/Lifted_Code_5_Reads",GENOME,SAMPLE,REPLICATE,GENOME,REFGENOME,"Positions_ChromosomePlaceMentTable.txt"), collapse = "_")
  
  print(FILENAME)
  POSITIONCHANGEDF<-read.delim(FILENAME,header = FALSE)
  colnames(POSITIONCHANGEDF)<-c("Chromosome.Reference","Chromosome.Query","Chromosome.ReferenceLiftedToQuery")
  
  POSITIONCHANGEDF$Frequency<-unlist(lapply(
    trimws(POSITIONCHANGEDF$Chromosome.Reference, which = "left"),
    function(LINE) strsplit(LINE,split = " ")[[1]][1]
  ))
  
  POSITIONCHANGEDF$Chromosome.Reference<-unlist(lapply(
    trimws(POSITIONCHANGEDF$Chromosome.Reference, which = "left"),
    function(LINE) strsplit(LINE,split = " ")[[1]][2]
  ))
  
  
  POSITIONCHANGEDF$Sample<-SAMPLE
  POSITIONCHANGEDF$Replicate<-REPLICATE
  POSITIONCHANGEDF$Query.Genome<-GENOME
  POSITIONCHANGEDF$Reference.Genome<-REFGENOME

  return(POSITIONCHANGEDF)
}



PositionChromosomeDataframe<-lapply(
  1:30,
  function(FILENUMBER) importChromosomePlacementFile(FILENUMBER)
) %>% bind_rows()

str(PositionChromosomeDataframe)
PositionChromosomeDataframe$Frequency<-as.numeric(PositionChromosomeDataframe$Frequency)
str(PositionChromosomeDataframe)


### Filter out chrM entries
PositionChromosomeDataframe<-subset(PositionChromosomeDataframe, Chromosome.Reference!="chrM" & Chromosome.Query!="chrM" )



### Change names
PositionChromosomeDataframe$Query.Genome[grep("mat",PositionChromosomeDataframe$Query.Genome)]<-"Maternal"
PositionChromosomeDataframe$Query.Genome[grep("pat",PositionChromosomeDataframe$Query.Genome)]<-"Paternal"

### Change library names
LibraryConversionTable<-read.delim(file = "Library_Conversion_Table.txt",header = FALSE)
colnames(LibraryConversionTable)<-c("Sample","Replicate","Library")

GetLibraryReplicate<-function(SAMPLE,LIBRARY){
  return(subset(LibraryConversionTable, Sample==SAMPLE & Library==LIBRARY)$Replicate)
}

PositionChromosomeDataframe$Library<-PositionChromosomeDataframe$Replicate
PositionChromosomeDataframe$Replicate<-unlist(lapply(1:nrow(PositionChromosomeDataframe), function(INDEX) GetLibraryReplicate(PositionChromosomeDataframe$Sample[INDEX],PositionChromosomeDataframe$Replicate[INDEX]) ))






#merge(PositionChromosomeDataframe,LibraryConversionTable,by=c("Sample","Replicate"))




### Same chromosome before liftover
calcSameChromBeforeLiftPercent<-function(SAMPLE,REPLICATE,QUERY){
  TotalNumberOfEntries <- sum( subset(PositionChromosomeDataframe, Sample==SAMPLE & Replicate==REPLICATE & Query.Genome==QUERY)$Frequency )
  SameChromBeforeLift <- sum( subset(PositionChromosomeDataframe, Sample==SAMPLE & Replicate==REPLICATE & Query.Genome==QUERY & Chromosome.Reference==Chromosome.Query)$Frequency )
  
  Output<-data.frame(
    Sample=SAMPLE,
    Replicate=REPLICATE,
    Query=QUERY,
    SameChr=SameChromBeforeLift,
    Total=TotalNumberOfEntries
  )
  
  return(Output)
}

ChromosomeComparisonBeforelift.Dataframe<-lapply(
  unique(PositionChromosomeDataframe$Sample),
  function(SAMPLE) lapply(unique(PositionChromosomeDataframe$Replicate),
                          function(REPLICATE) lapply(unique(PositionChromosomeDataframe$Query.Genome),
                                                     function(QUERY) calcSameChromBeforeLiftPercent(SAMPLE,REPLICATE,QUERY)
      ) %>% bind_rows()
    ) %>% bind_rows()
  ) %>% bind_rows()


ChromosomeComparisonBeforelift.Dataframe$Percent<-ChromosomeComparisonBeforelift.Dataframe$SameChr/ChromosomeComparisonBeforelift.Dataframe$Total


Plot.SameChrBefore<-ggplot(ChromosomeComparisonBeforelift.Dataframe,aes(x=Query,y=Percent*100,fill=Replicate))+
  geom_boxplot(outlier.alpha = 0)+
  geom_jitter(position=position_jitterdodge(jitter.width = 0.25))+
  scale_fill_manual(values = c(CS.ColorPalette[11,1],CS.ColorPalette[12,1]))+
  CS.THEME+
  ylab("%Map to same chromosome")+
  ggtitle("Reads that map to\nsame chromosome before liftover")+
  theme(legend.position = "bottom")
  




### Can liftover

canLiftoverPercent<-function(SAMPLE,REPLICATE,QUERY){
  TotalNumberOfEntries <- sum( subset(PositionChromosomeDataframe, Sample==SAMPLE & Replicate==REPLICATE & Query.Genome==QUERY)$Frequency )
  CanLift <- sum( subset(PositionChromosomeDataframe, Sample==SAMPLE & Replicate==REPLICATE & Query.Genome==QUERY & Chromosome.ReferenceLiftedToQuery!="FAILEDLIFT")$Frequency )
  
  Output<-data.frame(
    Sample=SAMPLE,
    Replicate=REPLICATE,
    Query=QUERY,
    NoLiftover=CanLift,
    Total=TotalNumberOfEntries
  )
  
  return(Output)
}

CanLiftover.Dataframe<-lapply(
  unique(PositionChromosomeDataframe$Sample),
  function(SAMPLE) lapply(unique(PositionChromosomeDataframe$Replicate),
                          function(REPLICATE) lapply(unique(PositionChromosomeDataframe$Query.Genome),
                                                     function(QUERY) canLiftoverPercent(SAMPLE,REPLICATE,QUERY)
                          ) %>% bind_rows()
  ) %>% bind_rows()
) %>% bind_rows()

CanLiftover.Dataframe$Percent<-CanLiftover.Dataframe$NoLiftover/CanLiftover.Dataframe$Total

Plot.CanLiftover<-ggplot(CanLiftover.Dataframe,aes(x=Query,y=Percent*100,fill=Replicate))+
  geom_boxplot(outlier.alpha = 0)+
  geom_jitter(position=position_jitterdodge(jitter.width = 0.25))+
  scale_fill_manual(values = c(CS.ColorPalette[11,1],CS.ColorPalette[12,1]))+
  CS.THEME+
  ylab("%Reads that liftover")+
  ggtitle("Reads that do lift over\n")+
  theme(legend.position = "bottom")


### Same chromosome after liftover
calcSameChromAfterLiftPercent<-function(SAMPLE,REPLICATE,QUERY){
  TotalNumberOfEntries <- sum( subset(PositionChromosomeDataframe, Sample==SAMPLE & Replicate==REPLICATE & Query.Genome==QUERY)$Frequency )
  SameChromBeforeLift <- sum( subset(PositionChromosomeDataframe, Sample==SAMPLE & Replicate==REPLICATE & Query.Genome==QUERY & Chromosome.Reference==Chromosome.ReferenceLiftedToQuery)$Frequency )
  
  Output<-data.frame(
    Sample=SAMPLE,
    Replicate=REPLICATE,
    Query=QUERY,
    SameChr=SameChromBeforeLift,
    Total=TotalNumberOfEntries
  )
  
  return(Output)
}

ChromosomeComparisonAfterlift.Dataframe<-lapply(
  unique(PositionChromosomeDataframe$Sample),
  function(SAMPLE) lapply(unique(PositionChromosomeDataframe$Replicate),
                          function(REPLICATE) lapply(unique(PositionChromosomeDataframe$Query.Genome),
                                                     function(QUERY) calcSameChromAfterLiftPercent(SAMPLE,REPLICATE,QUERY)
                          ) %>% bind_rows()
  ) %>% bind_rows()
) %>% bind_rows()


ChromosomeComparisonAfterlift.Dataframe$Percent<-ChromosomeComparisonAfterlift.Dataframe$SameChr/ChromosomeComparisonAfterlift.Dataframe$Total


Plot.SameChrAfter<-ggplot(ChromosomeComparisonAfterlift.Dataframe,aes(x=Query,y=Percent*100,fill=Replicate))+
  geom_boxplot(outlier.alpha = 0)+
  geom_jitter(position=position_jitterdodge(jitter.width = 0.25))+
  scale_fill_manual(values = c(CS.ColorPalette[11,1],CS.ColorPalette[12,1]))+
  CS.THEME+
  ylab("%Map to same chromosome")+
  ggtitle("Reads that map to\nsame chromosome after liftover")+
  theme(legend.position = "bottom")







##### Position change
importPositionChangeFile<-function(FILENUMBER){
  FILELIST<-dir("ProcessedFiles/")[grep("PositionChangeTable.txt", dir("ProcessedFiles/"))]
  
  SAMPLE<-strsplit(FILELIST[FILENUMBER], split = "_")[[1]][6]
  
  REPLICATE<-strsplit(FILELIST[FILENUMBER], split = "_")[[1]][7]
  
  GENOME<-strsplit(FILELIST[FILENUMBER], split = "_")[[1]][8]
  REFGENOME<-strsplit(FILELIST[FILENUMBER], split = "_")[[1]][9]

  FILENAME<-paste( c("ProcessedFiles/Lifted_Code_5_Reads",GENOME,SAMPLE,REPLICATE,GENOME,REFGENOME,"Positions_PositionChangeTable.txt"), collapse = "_")
  print("Importing:")
  print(FILENAME)
  
  POSITIONCHANGEDF<-read.delim(FILENAME,header = FALSE)
  colnames(POSITIONCHANGEDF)<-c("PositionChangeValue")
  
  POSITIONCHANGEDF$Frequency<-unlist(lapply(
    trimws(POSITIONCHANGEDF$PositionChangeValue, which = "left"),
    function(LINE) strsplit(LINE,split = " ")[[1]][1]
  ))

  POSITIONCHANGEDF$PositionChangeValue<-unlist(lapply(
    trimws(POSITIONCHANGEDF$PositionChangeValue, which = "left"),
    function(LINE) strsplit(LINE,split = " ")[[1]][2]
  ))
  
  
  POSITIONCHANGEDF$Sample<-SAMPLE
  POSITIONCHANGEDF$Replicate<-REPLICATE
  POSITIONCHANGEDF$Query.Genome<-GENOME
  POSITIONCHANGEDF$Reference.Genome<-REFGENOME
  POSITIONCHANGEDF$Frequency<-as.numeric(POSITIONCHANGEDF$Frequency)
  POSITIONCHANGEDF$PositionChangeValue<-as.numeric(POSITIONCHANGEDF$PositionChangeValue)
  
  ### Change names
  POSITIONCHANGEDF$Query.Genome[grep("mat",POSITIONCHANGEDF$Query.Genome)]<-"Maternal"
  POSITIONCHANGEDF$Query.Genome[grep("pat",POSITIONCHANGEDF$Query.Genome)]<-"Paternal"
  
  ### Change library names
  POSITIONCHANGEDF$Library<-POSITIONCHANGEDF$Replicate
  POSITIONCHANGEDF$Replicate<-GetLibraryReplicate(POSITIONCHANGEDF$Sample[1],POSITIONCHANGEDF$Replicate[1])
  
  print("Caculating position change")
  TotalNumberOfEntries <- sum( POSITIONCHANGEDF$Frequency )
  SamePosition <- sum( subset(POSITIONCHANGEDF, PositionChangeValue=="0")$Frequency )
  
  Output<-data.frame(
    Sample=POSITIONCHANGEDF$Sample[1],
    Replicate=POSITIONCHANGEDF$Replicate[1],
    Query=POSITIONCHANGEDF$Query.Genome[1],
    SamePosition=SamePosition,
    Total=TotalNumberOfEntries
  )
  print("Done")
  return(Output)
  
}



PositionComparison.Dataframe<-lapply(
  1:30,
  function(FILENUMBER) importPositionChangeFile(FILENUMBER)
    ) %>% bind_rows()

head(PositionComparison.Dataframe)

PositionComparison.Dataframe$Percent<-PositionComparison.Dataframe$SamePosition/PositionComparison.Dataframe$Total


Plot.SamePosition<-ggplot(PositionComparison.Dataframe,aes(x=Query,y=Percent*100,fill=Replicate))+
  geom_boxplot(outlier.alpha = 0)+
  geom_jitter(position=position_jitterdodge(jitter.width = 0.25))+
  scale_fill_manual(values = c(CS.ColorPalette[11,1],CS.ColorPalette[12,1]))+
  CS.THEME+
  ylab("%Map to same position")+
  ggtitle("Reads that map to\nsame position after liftover")+
  theme(legend.position = "bottom")

mean(subset(PositionComparison.Dataframe,Query=="CHM13")$Percent)
sd(subset(PositionComparison.Dataframe,Query=="CHM13")$Percent)
mean(subset(PositionComparison.Dataframe,Query=="Maternal")$Percent)
sd(subset(PositionComparison.Dataframe,Query=="Maternal")$Percent)
mean(subset(PositionComparison.Dataframe,Query=="Paternal")$Percent)
sd(subset(PositionComparison.Dataframe,Query=="Paternal")$Percent)


Panel<-plot_grid(
  Plot.SameChrBefore,
  Plot.CanLiftover,
  Plot.SameChrAfter,
  Plot.SamePosition,
  ncol=4,
  labels = c("L","M","N","O")
)

ggsave2(filename = "Read_Position_Trace_Panel.png",plot = Panel,units = "in", height = 5, width = 18)


