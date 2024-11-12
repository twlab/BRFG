
rm(list=ls())
set.seed(0)

require("ggplot2")
require("cowplot")
require("viridis")
require("reshape2")
require("dplyr")
require("tidyverse")

setwd("/Users/juanmacias/Library/CloudStorage/Box-Box/LawsonLab/JuanMacias/Postdoc/Experimental/")
source("/Users/juanmacias/Library/CloudStorage/Box-Box/LawsonLab/JuanMacias/Misc_Code/CS_Aesthetic.R")

theme_set(theme_cowplot())


Links<-read.delim("bschr8.gfa.links", header=FALSE, sep="\t", stringsAsFactors=FALSE)
colnames(Links)<-c("Linetype","SourceSeg", "SourceSeg.Orient", "TargetSeg", "TargetSeg.Orient", "Overlap")

Segments<-read.delim("bschr8.gfa.segments", header=FALSE, sep="\t", stringsAsFactors=FALSE)
colnames(Segments)<-c("Linetype", "SegID", "SegSeq","StableSequence","Offset","Rank")

# remove unnecessary columns
Segments<-Segments[,-ncol(Segments)]
Segments<-Segments[,-1]
Links<-Links[,-6]
Links<-Links[,-1]

# convert character columns to numeric
Segments$SegID<-as.numeric(Segments$SegID)
Segments$IsRef<-ifelse(Segments$StableSequence!="",TRUE,FALSE)



Links$SourceSeg<-as.numeric(Links$SourceSeg)
Links$TargetSeg<-as.numeric(Links$TargetSeg)
Links$SourceIsRef<-Segments$IsRef[match(Links$SourceSeg, Segments$SegID)]
Links$TargetIsRef<-Segments$IsRef[match(Links$TargetSeg, Segments$SegID)]


Segments$Offset<-unlist(lapply(Segments$Offset, function(X) as.numeric(gsub("SO:i:","",X))))
Segments$Offset<-as.numeric(Segments$Offset)

Segments$SegSeqLength<-Segments$Offset+nchar(Segments$SegSeq)


# Add score
Links$Score<-abs(Links$TargetSeg-Links$SourceSeg)
Links[Links$Score>10,]$Score<-11

Links$Score<-as.character(Links$Score)
Links$Score[Links$Score=="11"]<-">10"
Links$Score<-factor(Links$Score, levels=c("1","2","3","4","5","6","7","8","9","10",">10"), ordered=TRUE)

Links$SourceSegPosition<-Segments$SegSeqLength[match(Links$SourceSeg, Segments$SegID)]

Links$TargetSegSeqLength<-nchar(Segments$SegSeq[match(Links$TargetSeg, Segments$SegID)])

Links<-Links[order(Links$TargetSegSeqLength),]

# Set x scale to 0-11, wit a tick every 1, no subticks
GraphStructurePlot<-ggplot(subset(Links, SourceIsRef==TRUE & TargetIsRef==FALSE), aes(x=Score, y=SourceSegPosition/1000000, color=TargetSegSeqLength)) +
  geom_jitter(size=0.2, width=0.4,alpha=1)+
  scale_color_viridis(option="magma", begin = 0.85, end = 0.1, trans="log10", name="Targe segment\nsequence length (bp)")+
  theme(legend.position="bottom",
        panel.grid.minor = element_blank(),
        panel.grid.minor.x = element_blank()#, 
        #panel.grid.minor.y = element_blank()
        )+
  #scale_x_discrete(breaks = seq(1, 11, by = 1), expand = c(0,0))+
  CS.THEME+
  scale_y_continuous(limits = c(0,ceiling(max(subset(Links, SourceIsRef==TRUE)$SourceSegPosition/1000000))), expand = c(0,0))+
  ylab("Position (Mb)")+
  xlab("Change in segmentID")+
  theme(legend.text = element_text(size=9))

# save plot
ggsave("GraphStructurePlot.png", GraphStructurePlot, width=5, height=7, units="in", dpi=300)







