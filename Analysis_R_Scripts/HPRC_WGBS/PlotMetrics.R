
rm(list = ls())
set.seed(0)

library(ggplot2)
library(cowplot)
library(scales)

source("/Users/juanmacias/Library/CloudStorage/Box-Box/LawsonLab/JuanMacias/Misc_Code/CS_Aesthetic.R")
setwd("/Users/juanmacias/Library/CloudStorage/Box-Box/LawsonLab/JuanMacias/Postdoc/HPRC/WGBS/")


AlignmentMetrics<-read.delim("AlignmentMetrics.txt",header = TRUE)
AlignmentMetrics$ref[grep("_pat",AlignmentMetrics$ref)]<-"Paternal"
AlignmentMetrics$ref[grep("_mat",AlignmentMetrics$ref)]<-"Maternal"




Plot.A<-ggplot(AlignmentMetrics,aes(x=ref,y=total/1000000,fill=replicate))+
  geom_bar(stat="identity", position = position_dodge())+
  ylab("Total number\nof reads (M)")+
  CS.THEME+
  scale_y_continuous(expand = c(0,0),limits = c(0,300))+
  facet_grid(~sample)+
  theme(axis.text.x = element_text(angle = 30,hjust=1))+
  scale_fill_manual(values = c(CS.ColorPalette[10,1],CS.ColorPalette[11,1],CS.ColorPalette[12,1]))

Plot.B<-ggplot(AlignmentMetrics,aes(x=ref,y=after_trim/1000000,fill=replicate))+
  geom_bar(stat="identity", position = position_dodge())+
  ylab("Number of reads\nafter trim (M)")+
  CS.THEME+
  scale_y_continuous(expand = c(0,0),limits = c(0,300))+
  facet_grid(~sample)+
  theme(axis.text.x = element_text(angle = 30,hjust=1))+
  scale_fill_manual(values = c(CS.ColorPalette[10,1],CS.ColorPalette[11,1],CS.ColorPalette[12,1]))

Plot.C<-ggplot(AlignmentMetrics,aes(x=ref,y=mapped_perc,fill=replicate))+
  geom_bar(stat="identity", position = position_dodge())+
  ylab("%Reads mapped")+
  CS.THEME+
  scale_y_continuous(expand = c(0,0),limits = c(0,1))+
  facet_grid(~sample)+
  theme(axis.text.x = element_text(angle = 30,hjust=1))+
  scale_fill_manual(values = c(CS.ColorPalette[10,1],CS.ColorPalette[11,1],CS.ColorPalette[12,1]))

Plot.D<-ggplot(AlignmentMetrics,aes(x=ref,y=unique_perc,fill=replicate))+
  geom_bar(stat="identity", position = position_dodge())+
  ylab("%Uniquely mapped")+
  CS.THEME+
  scale_y_continuous(expand = c(0,0),limits = c(0,1))+
  facet_grid(~sample)+
  theme(axis.text.x = element_text(angle = 30,hjust=1))+
  scale_fill_manual(values = c(CS.ColorPalette[10,1],CS.ColorPalette[11,1],CS.ColorPalette[12,1]))

Plot.E<-ggplot(AlignmentMetrics,aes(x=ref,y=multi_perc,fill=replicate))+
  geom_bar(stat="identity", position = position_dodge())+
  ylab("%Multi-mapped")+
  CS.THEME+
  scale_y_continuous(expand = c(0,0),limits = c(0,1))+
  facet_grid(~sample)+
  theme(axis.text.x = element_text(angle = 30,hjust=1))+
  scale_fill_manual(values = c(CS.ColorPalette[10,1],CS.ColorPalette[11,1],CS.ColorPalette[12,1]))

Plot.F<-ggplot(AlignmentMetrics,aes(x=ref,y=dedup_perc,fill=replicate))+
  geom_bar(stat="identity", position = position_dodge())+
  ylab("%Uniquly mapped\nreads after de-duplication")+
  CS.THEME+
  scale_y_continuous(expand = c(0,0),limits = c(0,1))+
  facet_grid(~sample)+
  theme(axis.text.x = element_text(angle = 30,hjust=1))+
  scale_fill_manual(values = c(CS.ColorPalette[10,1],CS.ColorPalette[11,1],CS.ColorPalette[12,1]))



AlignmentMetricsPlot<-plot_grid(Plot.A,Plot.B,Plot.C,Plot.D,Plot.E,Plot.F, labels = "AUTO", ncol=2)


ggsave2(filename = "AlignmentMetricsPlot.png",plot = AlignmentMetricsPlot, units = "in", width = 17, height = 8)






Coveragemetrics<-read.delim("Genome_coverage_metrics.txt",header = TRUE)
Coveragemetrics$ref[grep("_pat",Coveragemetrics$ref)]<-"Paternal"
Coveragemetrics$ref[grep("_mat",Coveragemetrics$ref)]<-"Maternal"



Plot.G<-ggplot(Coveragemetrics,aes(x=ref,y=max_coverage,fill=replicate))+
  geom_bar(stat="identity", position = position_dodge())+
  ylab("Max coverage")+
  CS.THEME+
  scale_y_continuous(expand = c(0,0),limits = c(0,10))+
  facet_grid(~sample)+
  theme(axis.text.x = element_text(angle = 30,hjust=1))+
  scale_fill_manual(values = c(CS.ColorPalette[10,1],CS.ColorPalette[11,1],CS.ColorPalette[12,1]))

Plot.H<-ggplot(Coveragemetrics,aes(x=ref,y=c_cov,fill=replicate))+
  geom_bar(stat="identity", position = position_dodge())+
  ylab("C coverage")+
  CS.THEME+
  scale_y_continuous(expand = c(0,0),limits = c(0,6))+
  facet_grid(~sample)+
  theme(axis.text.x = element_text(angle = 30,hjust=1))+
  scale_fill_manual(values = c(CS.ColorPalette[10,1],CS.ColorPalette[11,1],CS.ColorPalette[12,1]))

Plot.I<-ggplot(Coveragemetrics,aes(x=ref,y=cg_cov,fill=replicate))+
  geom_bar(stat="identity", position = position_dodge())+
  ylab("CpG coverage")+
  CS.THEME+
  scale_y_continuous(expand = c(0,0),limits = c(0,6))+
  facet_grid(~sample)+
  theme(axis.text.x = element_text(angle = 30,hjust=1))+
  scale_fill_manual(values = c(CS.ColorPalette[10,1],CS.ColorPalette[11,1],CS.ColorPalette[12,1]))


CoverageMetricsPlot<-plot_grid(Plot.G,Plot.H,Plot.I, labels = "AUTO", ncol=1)

ggsave2(filename = "CoverageMetricsPlot.png",plot = CoverageMetricsPlot, units = "in", width = 17/2, height = 8)












ConversionMetrics<-read.delim("Conversion_Metrics.txt",header = TRUE)
ConversionMetrics$ref[grep("_pat",ConversionMetrics$ref)]<-"Paternal"
ConversionMetrics$ref[grep("_mat",ConversionMetrics$ref)]<-"Maternal"


Plot.J<-ggplot(ConversionMetrics,aes(x=ref,y=lambda_rate,fill=replicate))+
  geom_bar(stat="identity", position = position_dodge())+
  ylab("Lambda conversion\nrate")+
  CS.THEME+
  scale_y_continuous(expand = c(0,0),limits = c(0,1))+
  facet_grid(~sample)+
  theme(axis.text.x = element_text(angle = 30,hjust=1))+
  scale_fill_manual(values = c(CS.ColorPalette[10,1],CS.ColorPalette[11,1],CS.ColorPalette[12,1]))

Plot.K<-ggplot(ConversionMetrics,aes(x=ref,y=CA_met,fill=replicate))+
  geom_bar(stat="identity", position = position_dodge())+
  ylab("CA methylation")+
  CS.THEME+
  scale_y_continuous(expand = c(0,0),limits = c(0,0.1))+
  facet_grid(~sample)+
  theme(axis.text.x = element_text(angle = 30,hjust=1))+
  scale_fill_manual(values = c(CS.ColorPalette[10,1],CS.ColorPalette[11,1],CS.ColorPalette[12,1]))

Plot.L<-ggplot(ConversionMetrics,aes(x=ref,y=CC_met,fill=replicate))+
  geom_bar(stat="identity", position = position_dodge())+
  ylab("CC methylation")+
  CS.THEME+
  scale_y_continuous(expand = c(0,0),limits = c(0,0.1))+
  facet_grid(~sample)+
  theme(axis.text.x = element_text(angle = 30,hjust=1))+
  scale_fill_manual(values = c(CS.ColorPalette[10,1],CS.ColorPalette[11,1],CS.ColorPalette[12,1]))

Plot.M<-ggplot(ConversionMetrics,aes(x=ref,y=CT_met,fill=replicate))+
  geom_bar(stat="identity", position = position_dodge())+
  ylab("CT methylation")+
  CS.THEME+
  scale_y_continuous(expand = c(0,0),limits = c(0,0.1))+
  facet_grid(~sample)+
  theme(axis.text.x = element_text(angle = 30,hjust=1))+
  scale_fill_manual(values = c(CS.ColorPalette[10,1],CS.ColorPalette[11,1],CS.ColorPalette[12,1]))

Plot.N<-ggplot(ConversionMetrics,aes(x=ref,y=CH_conv_rate,fill=replicate))+
  geom_bar(stat="identity", position = position_dodge())+
  ylab("CH conversion rate")+
  CS.THEME+
  scale_y_continuous(expand = c(0,0),limits = c(0,1))+
  facet_grid(~sample)+
  theme(axis.text.x = element_text(angle = 30,hjust=1))+
  scale_fill_manual(values = c(CS.ColorPalette[10,1],CS.ColorPalette[11,1],CS.ColorPalette[12,1]))


ConversionMetricsPlot<-plot_grid(Plot.J,Plot.K,Plot.L,Plot.M,Plot.N,labels = "AUTO", ncol=2)

ggsave2(filename = "ConversionMetricsPlot.png",plot = ConversionMetricsPlot, units = "in", width = 17, height = 8)




AlignmentMetrics$barcode==Coveragemetrics$barcode
AlignmentMetrics$barcode==ConversionMetrics$barcode



library(factoextra)

PCA.Subset<-cbind(AlignmentMetrics[,c(12,14,16,20)],Coveragemetrics[,c(9,10,11)],ConversionMetrics[,c(13,14,15,16,17)])




#row.names(PCA.Subset)<-TARGET.2.DATA$file
# PCA.Subset[,5]<-PCA.Subset[,5]/max(PCA.Subset[,5])
# PCA.Subset[,6]<-PCA.Subset[,6]/max(PCA.Subset[,6])
# PCA.Subset[,7]<-PCA.Subset[,7]/max(PCA.Subset[,7])
# PCA.Subset[,8]<-PCA.Subset[,8]/max(PCA.Subset[,8])
# PCA.Subset[,9]<-PCA.Subset[,9]/max(PCA.Subset[,9])
# PCA.Subset[,10]<-PCA.Subset[,10]/max(PCA.Subset[,10])
# PCA.Subset[,11]<-PCA.Subset[,11]/max(PCA.Subset[,11])
# PCA.Subset[,12]<-PCA.Subset[,12]/max(PCA.Subset[,12])
# PCA.Subset[,13]<-PCA.Subset[,13]/max(PCA.Subset[,13])

# Perfor PCA
res.pca<-prcomp(PCA.Subset, scale=TRUE)
#Plot.PCs<-fviz_eig(res.pca)
#Plot.PCs
#ggsave2(filename = "QC_PCs.png",plot = Plot.PCs,units = "in",width = 6/2, height = 6/2)

# Get the partial eta squared values
library(tidyr)
PCA.VarObject<-cbind(res.pca$x,AlignmentMetrics[,c(4,5,1)])
colnames(PCA.VarObject)[13:15]<-c("Sample","Replicate","Genome")

manova_model<-manova(cbind(PC1,PC2,PC3,PC4,PC5,PC6,PC7,PC8,PC9,PC10,PC11,PC12)~Sample*Genome+Replicate,PCA.VarObject)
manova_model.Summary<-summary(manova_model)
aov_results <- summary.aov(manova_model)

NUM.PCs<-12
PC.PartialEtaSquareds<-lapply(1:NUM.PCs,
                              function(PC){aov_results[[PC]][[2]][1:5]/sum(aov_results[[PC]][[2]])}
) %>% do.call(rbind,.) %>% as.data.frame 

colnames(PC.PartialEtaSquareds)<-c("Sample","Genome","Sample:Genome","Replicate","Residual")

PC.Statistics<-cbind(as.factor(1:NUM.PCs),get_eig(res.pca),PC.PartialEtaSquareds)
colnames(PC.Statistics)<-c("PC","Eigenvalue","VariancePercent","CumulativeVariancePercent","PartialEtaSquared.Sample","PartialEtaSquared.Genome","PartialEtaSquared.Sample:Genome","PartialEtaSquared.Replicate","PartialEtaSquared.Residual")
rownames(PC.Statistics)<-NULL

PC.Statistics$VariancePercent.Sample<-PC.Statistics$VariancePercent*PC.Statistics$PartialEtaSquared.Sample
PC.Statistics$VariancePercent.Genome<-PC.Statistics$VariancePercent*PC.Statistics$PartialEtaSquared.Genome
PC.Statistics$VariancePercent.Sample.Genome<-PC.Statistics$VariancePercent*PC.Statistics$`PartialEtaSquared.Sample:Genome`
PC.Statistics$VariancePercent.Replicate<-PC.Statistics$VariancePercent*PC.Statistics$PartialEtaSquared.Replicate
PC.Statistics$VariancePercent.Residual<-PC.Statistics$VariancePercent*PC.Statistics$PartialEtaSquared.Residual

# Save the PC statistics
write.table(PC.Statistics, file = "PC.QC.Statistics.txt", sep = "\t", row.names = FALSE, quote = FALSE)

# Make long df to plot variance percents
PC.Statistics.VariancePercent.Long<-gather(PC.Statistics[,c(1,10:14)],key="VarianceType",value="VariancePercent",VariancePercent.Sample:VariancePercent.Residual)
# clean up names
PC.Statistics.VariancePercent.Long$VarianceType<-gsub("VariancePercent.","",PC.Statistics.VariancePercent.Long$VarianceType)
PC.Statistics.VariancePercent.Long$VarianceType<-gsub("Sample.Genome","Sample:Genome",PC.Statistics.VariancePercent.Long$VarianceType)

PC.Statistics.VariancePercent.Long$VarianceType<-factor(PC.Statistics.VariancePercent.Long$VarianceType, levels = c("Sample","Genome","Sample:Genome","Replicate","Residual"), ordered = TRUE)

Plot.PCs<-ggplot(PC.Statistics.VariancePercent.Long, aes(fill=VarianceType, y=VariancePercent, x=PC)) + 
  geom_bar(position="stack", stat="identity")+
  scale_fill_manual(values = c("#ff7f5d","#2d5879","#b2823f","#aaaaaa","#333434"))+
  #scale_fill_manual(values = c("#2d5879","#aaaaaa","#333434","#ff7f5d","#b2823f"))+
  ylab("Percent of variance")+
  xlab("Principal component")+
  ggtitle("WGBS - QC metrics\n%Variance captured by PCs partitioned by factor")+
  CS.THEME+
  theme(legend.position="bottom")+
  scale_y_continuous(limits = c(0,round(max(PC.Statistics$VariancePercent)*1.05)), expand = c(0,0))+
  theme(legend.key.size = unit(0.2, "cm"))

ggsave2(filename = "QC_PCs.png",plot = Plot.PCs,units = "in",width = 7, height = 6)

# Plot the PCA
groups.Sample <- as.factor(AlignmentMetrics$sample)
PCA.Plot.Sample<-fviz_pca_ind(res.pca,
             col.ind = groups.Sample,
             palette = viridis(5, end=0.9),
             repel = TRUE,
             addEllipses = TRUE,
             ellipse.type = "convex",
             title="WGBS - Samples"
)+
CS.THEME+
  theme(legend.position="bottom",legend.text = element_text(size=8,hjust = 0.7), legend.title = element_text(size=12),legend.key.size = unit(0.25, "cm"))


groups <- as.factor(AlignmentMetrics$ref)
PCA.Plot.Genome<-fviz_pca_ind(res.pca,
             col.ind = groups,
             palette = viridis(5,option = "magma", end = 0.8),
             repel = TRUE,
             title="WGBS - Genome"
)+
CS.THEME+
  theme(legend.position="bottom",legend.text = element_text(size=8,hjust = 0.7), legend.title = element_text(size=12),legend.key.size = unit(0.25, "cm"))


groups <- as.factor(AlignmentMetrics$replicate)
PCA.Plot.Replicate<-fviz_pca_ind(res.pca,
             col.ind = groups,
             palette = viridis(2,option = "cividis", end = 0.8),
             repel = TRUE,
             title="WGBS - Replicate"
)+
CS.THEME+
  theme(legend.position="bottom",legend.text = element_text(size=8,hjust = 0.7), legend.title = element_text(size=12),legend.key.size = unit(0.25, "cm"))

PCA.Panel<-plot_grid(PCA.Plot.Sample,PCA.Plot.Genome,PCA.Plot.Replicate,labels = c("G","H","I"),ncol=3)

ggsave2(filename = "QC.PCA.Panel.png",plot = PCA.Panel, units = "in", width = 17*0.75, height = 6*0.75)



PCA.Subset<-cbind(AlignmentMetrics[,c(1,4,5,12,14,16,20)],Coveragemetrics[,c(9,10,11)],ConversionMetrics[,c(13,14,15,16,17)])


summary(aov(mapped_perc~sample+ref,PCA.Subset))

### mapped_perc
100*summary(aov(mapped_perc~sample+ref,PCA.Subset))[[1]][,2]/sum(summary(aov(mapped_perc~sample+ref,PCA.Subset))[[1]][,2])

### max_coverage
100*summary(aov(max_coverage~sample+ref,PCA.Subset))[[1]][,2]/sum(summary(aov(max_coverage~sample+ref,PCA.Subset))[[1]][,2])





