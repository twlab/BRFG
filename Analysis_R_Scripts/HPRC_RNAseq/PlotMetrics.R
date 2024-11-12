
rm(list = ls())
set.seed(0)

library(dplyr)
library(tidyr)
library(ggplot2)
library(cowplot)
source("/Users/juanmacias/Library/CloudStorage/Box-Box/LawsonLab/JuanMacias/Misc_Code/CS_Aesthetic.R")
setwd("/Users/juanmacias/Library/CloudStorage/Box-Box/LawsonLab/JuanMacias/Postdoc/HPRC/RNAseq/")

AlignmentMetrics<-read.delim("AlignmentMetrics.txt",header = TRUE)

AlignmentMetrics$Replicate[AlignmentMetrics$Replicate==1]<-"A"
AlignmentMetrics$Replicate[AlignmentMetrics$Replicate==2]<-"B"


AlignmentMetrics$Annotation[AlignmentMetrics$Annotation=="liftoff"]<-"gencode"

AlignmentMetrics<-subset(AlignmentMetrics, Sample != "HG002")

colnames(AlignmentMetrics)




Number.of.input.reads.Plot<-ggplot(AlignmentMetrics,aes(x=Sample,y=Number.of.input.reads/1000000, fill=Replicate))+
  geom_bar(stat="identity", position = position_dodge())+
  CS.THEME+
  scale_fill_manual(values = c(CS.ColorPalette[11,1],CS.ColorPalette[12,1]))+
  theme(axis.text.x = element_text(angle = 30,hjust=1))+
  ylab("Number.of.input.reads (M)")

Average.mapped.length.Plot<-ggplot(AlignmentMetrics,aes(x=Genome,y=Average.mapped.length, fill=Replicate))+
  geom_bar(stat="identity", position = position_dodge())+
  CS.THEME+
  scale_fill_manual(values = c(CS.ColorPalette[11,1],CS.ColorPalette[12,1]))+
  theme(axis.text.x = element_text(angle = 30,hjust=1))+
  ylab("Average.mapped.length")+
  facet_grid(Annotation~Sample)



mean(subset(AlignmentMetrics, Annotation=="gencode" & Genome!="CHM13")[,14]+subset(AlignmentMetrics, Annotation=="gencode" & Genome!="CHM13")[,9])

mean(subset(AlignmentMetrics, Annotation=="gencode" & Genome=="CHM13")[,14]+subset(AlignmentMetrics, Annotation=="gencode" & Genome=="CHM13")[,9])


### Placement class breakdown

Uniquely.mapped.reads.Plot<-ggplot(AlignmentMetrics,aes(x=Genome,y=Uniquely.mapped.reads.., fill=Replicate))+
  geom_bar(stat="identity", position = position_dodge())+
  CS.THEME+
  scale_fill_manual(values = c(CS.ColorPalette[11,1],CS.ColorPalette[12,1]))+
  theme(axis.text.x = element_text(angle = 30,hjust=1))+
  ylab("Uniquely.mapped.reads%")+
  facet_grid(Annotation~Sample)

pct.of.reads.mapped.to.multiple.loci.Plot<-ggplot(AlignmentMetrics,aes(x=Genome,y=X...of.reads.mapped.to.multiple.loci, fill=Replicate))+
  geom_bar(stat="identity", position = position_dodge())+
  CS.THEME+
  scale_fill_manual(values = c(CS.ColorPalette[11,1],CS.ColorPalette[12,1]))+
  theme(axis.text.x = element_text(angle = 30,hjust=1))+
  ylab("%.of.reads.mapped.to.multiple.loci")+
  facet_grid(Annotation~Sample)


pct.of.reads.mapped.to.too.many.loci.Plot<-ggplot(AlignmentMetrics,aes(x=Genome,y=X...of.reads.mapped.to.too.many.loci, fill=Replicate))+
  geom_bar(stat="identity", position = position_dodge())+
  CS.THEME+
  scale_fill_manual(values = c(CS.ColorPalette[11,1],CS.ColorPalette[12,1]))+
  theme(axis.text.x = element_text(angle = 30,hjust=1))+
  ylab("%of.reads.mapped.to.too.many.loci")+
  facet_grid(Annotation~Sample)


pct.of.reads.unmapped..other.loci.Plot<-ggplot(AlignmentMetrics,aes(x=Genome,y=X...of.reads.unmapped..other, fill=Replicate))+
  geom_bar(stat="identity", position = position_dodge())+
  CS.THEME+
  scale_fill_manual(values = c(CS.ColorPalette[11,1],CS.ColorPalette[12,1]))+
  theme(axis.text.x = element_text(angle = 30,hjust=1))+
  ylab("%of.reads.unmapped..other")+
  facet_grid(Annotation~Sample)

pct.of.reads.unmapped..too.many.mismatches.Plot<-ggplot(AlignmentMetrics,aes(x=Genome,y=X...of.reads.unmapped..too.many.mismatches, fill=Replicate))+
  geom_bar(stat="identity", position = position_dodge())+
  CS.THEME+
  scale_fill_manual(values = c(CS.ColorPalette[11,1],CS.ColorPalette[12,1]))+
  theme(axis.text.x = element_text(angle = 30,hjust=1))+
  ylab("%of.reads.unmapped..too.many.mismatches")+
  facet_grid(Annotation~Sample)

pct.of.reads.unmapped..too.short.Plot<-ggplot(AlignmentMetrics,aes(x=Genome,y=X...of.reads.unmapped..too.short, fill=Replicate))+
  geom_bar(stat="identity", position = position_dodge())+
  CS.THEME+
  scale_fill_manual(values = c(CS.ColorPalette[11,1],CS.ColorPalette[12,1]))+
  theme(axis.text.x = element_text(angle = 30,hjust=1))+
  ylab("%of.reads.unmapped..too.short")+
  facet_grid(Annotation~Sample)


plot_grid(Uniquely.mapped.reads.Plot,pct.of.reads.mapped.to.multiple.loci.Plot,labels = "AUTO")


UnmappedRateMetricsPanel<-plot_grid(pct.of.reads.mapped.to.too.many.loci.Plot,pct.of.reads.unmapped..too.many.mismatches.Plot,pct.of.reads.unmapped..too.short.Plot,pct.of.reads.unmapped..other.loci.Plot,labels = "AUTO")

ggsave(filename = "Figures/ShortReadUnmappedRateMetricsPanel.png",plot = UnmappedRateMetricsPanel,units = "in", height = 10, width = 20)





Number.of.splices..AT.AC.Plot<-ggplot(AlignmentMetrics,aes(x=Genome,y=Number.of.splices..AT.AC, fill=Replicate))+
  geom_bar(stat="identity", position = position_dodge())+
  CS.THEME+
  scale_fill_manual(values = c(CS.ColorPalette[11,1],CS.ColorPalette[12,1]))+
  theme(axis.text.x = element_text(angle = 30,hjust=1))+
  ylab("Number.of.splices..AT.AC")+
  facet_grid(Annotation~Sample)

Number.of.splices..Annotated..sjdb..Plot<-ggplot(AlignmentMetrics,aes(x=Genome,y=Number.of.splices..Annotated..sjdb., fill=Replicate))+
  geom_bar(stat="identity", position = position_dodge())+
  CS.THEME+
  scale_fill_manual(values = c(CS.ColorPalette[11,1],CS.ColorPalette[12,1]))+
  theme(axis.text.x = element_text(angle = 30,hjust=1))+
  ylab("Number.of.splices..Annotated..sjdb.")+
  facet_grid(Annotation~Sample)

Number.of.splices..GC.AG.Plot<-ggplot(AlignmentMetrics,aes(x=Genome,y=Number.of.splices..GC.AG, fill=Replicate))+
  geom_bar(stat="identity", position = position_dodge())+
  CS.THEME+
  scale_fill_manual(values = c(CS.ColorPalette[11,1],CS.ColorPalette[12,1]))+
  theme(axis.text.x = element_text(angle = 30,hjust=1))+
  ylab("Number.of.splices..GC.AG")+
  facet_grid(Annotation~Sample)


Number.of.splices..GT.AG.Plot<-ggplot(AlignmentMetrics,aes(x=Genome,y=Number.of.splices..GT.AG, fill=Replicate))+
  geom_bar(stat="identity", position = position_dodge())+
  CS.THEME+
  scale_fill_manual(values = c(CS.ColorPalette[11,1],CS.ColorPalette[12,1]))+
  theme(axis.text.x = element_text(angle = 30,hjust=1))+
  ylab("Number.of.splices..GT.AG")+
  facet_grid(Annotation~Sample)

Number.of.splices..Non.canonical.Plot<-ggplot(AlignmentMetrics,aes(x=Genome,y=Number.of.splices..Non.canonical, fill=Replicate))+
  geom_bar(stat="identity", position = position_dodge())+
  CS.THEME+
  scale_fill_manual(values = c(CS.ColorPalette[11,1],CS.ColorPalette[12,1]))+
  theme(axis.text.x = element_text(angle = 30,hjust=1))+
  ylab("Number.of.splices..Non.canonical")+
  facet_grid(Annotation~Sample)

Number.of.splices..Total.Plot<-ggplot(AlignmentMetrics,aes(x=Genome,y=Number.of.splices..Total, fill=Replicate))+
  geom_bar(stat="identity", position = position_dodge())+
  CS.THEME+
  scale_fill_manual(values = c(CS.ColorPalette[11,1],CS.ColorPalette[12,1]))+
  theme(axis.text.x = element_text(angle = 30,hjust=1))+
  ylab("Number.of.splices..Total")+
  facet_grid(Annotation~Sample)





Deletion.average.length.Plot<-ggplot(AlignmentMetrics,aes(x=Genome,y=Deletion.average.length, fill=Replicate))+
  geom_bar(stat="identity", position = position_dodge())+
  CS.THEME+
  scale_fill_manual(values = c(CS.ColorPalette[11,1],CS.ColorPalette[12,1]))+
  theme(axis.text.x = element_text(angle = 30,hjust=1))+
  ylab("Deletion.average.length")+
  facet_grid(Annotation~Sample)


Insertion.average.length.Plot<-ggplot(AlignmentMetrics,aes(x=Genome,y=Insertion.average.length, fill=Replicate))+
  geom_bar(stat="identity", position = position_dodge())+
  CS.THEME+
  scale_fill_manual(values = c(CS.ColorPalette[11,1],CS.ColorPalette[12,1]))+
  theme(axis.text.x = element_text(angle = 30,hjust=1))+
  ylab("Insertion.average.length")+
  facet_grid(Annotation~Sample)


PCA.Subset<-AlignmentMetrics[,c(8:32)]
PCA.Subset<-PCA.Subset[,!(colnames(PCA.Subset) %in% colnames(PCA.Subset)[colSums(PCA.Subset)==0])]

library(factoextra)

res.pca<-prcomp(PCA.Subset, scale=TRUE)
fviz_eig(res.pca)


groups.Sample <- as.factor(AlignmentMetrics$Sample)
PCA.Plot.Sample<-fviz_pca_ind(res.pca,
                              col.ind = groups.Sample,
                              palette = viridis(6, end=0.9),
                              repel = TRUE,
                              title="Samples - PCA"
)+
  CS.THEME+
  theme(legend.position="bottom")

PCA.Plot.Sample

groups <- as.factor(AlignmentMetrics$Genome)
PCA.Plot.Genome<-fviz_pca_ind(res.pca,
                              col.ind = groups,
                              palette = viridis(4,option = "magma", end = 0.8),
                              repel = TRUE,
                              title="Genome - PCA"
)+
  CS.THEME+
  theme(legend.position="bottom")

PCA.Plot.Genome

groups <- as.factor(AlignmentMetrics$Replicate)
PCA.Plot.Replicate<-fviz_pca_ind(res.pca,
                                 col.ind = groups,
                                 palette = viridis(2,option = "cividis", end = 0.8),
                                 repel = TRUE,
                                 title="Replicate - PCA"
)+
  CS.THEME+
  theme(legend.position="bottom")


groups <- as.factor(AlignmentMetrics$Annotation)
PCA.Plot.Annotation<-fviz_pca_ind(res.pca,
                                  col.ind = groups,
                                  palette = viridis(4,option = "cividis", end = 0.8),
                                  repel = TRUE,
                                  title="Annotation - PCA"
)+
  CS.THEME+
  theme(legend.position="bottom")


AligmentMetric.PCA.Panel<-plot_grid(PCA.Plot.Sample,PCA.Plot.Genome,PCA.Plot.Replicate, PCA.Plot.Annotation,labels = "AUTO",ncol=2)

ggsave2(filename = "Figures/AligmentMetric.PCA.Panel.png",plot = AligmentMetric.PCA.Panel, units = "in", width = 12, height = 12)




###### Picard

PicardMetrics<-read.delim("PicardMetrics.txt",header = TRUE)

PicardMetrics$Replicate[PicardMetrics$Replicate==1]<-"A"
PicardMetrics$Replicate[PicardMetrics$Replicate==2]<-"B"

PicardMetrics$Annotation[PicardMetrics$Annotation=="LiftOff"]<-"Gencode"

PicardMetrics<-subset(PicardMetrics, Sample != "HG002")



pf_aligned_bases.Plot<-ggplot(PicardMetrics,aes(x=Genome,y=pf_aligned_bases, fill=Replicate))+
  geom_bar(stat="identity", position = position_dodge())+
  CS.THEME+
  scale_fill_manual(values = c(CS.ColorPalette[11,1],CS.ColorPalette[12,1]))+
  theme(axis.text.x = element_text(angle = 30,hjust=1))+
  ylab("pf_aligned_bases")+
  facet_grid(Annotation~Sample)

pct_ribosomal_bases.Plot<-ggplot(PicardMetrics,aes(x=Genome,y=pct_ribosomal_bases, fill=Replicate))+
  geom_bar(stat="identity", position = position_dodge())+
  CS.THEME+
  scale_fill_manual(values = c(CS.ColorPalette[11,1],CS.ColorPalette[12,1]))+
  theme(axis.text.x = element_text(angle = 30,hjust=1))+
  ylab("pct_ribosomal_bases")+
  facet_grid(Annotation~Sample)

pct_coding_bases.Plot<-ggplot(PicardMetrics,aes(x=Genome,y=pct_coding_bases, fill=Replicate))+
  geom_bar(stat="identity", position = position_dodge())+
  CS.THEME+
  scale_fill_manual(values = c(CS.ColorPalette[11,1],CS.ColorPalette[12,1]))+
  theme(axis.text.x = element_text(angle = 30,hjust=1))+
  ylab("pct_coding_bases")+
  facet_grid(Annotation~Sample)

pct_utr_bases.Plot<-ggplot(PicardMetrics,aes(x=Genome,y=pct_utr_bases, fill=Replicate))+
  geom_bar(stat="identity", position = position_dodge())+
  CS.THEME+
  scale_fill_manual(values = c(CS.ColorPalette[11,1],CS.ColorPalette[12,1]))+
  theme(axis.text.x = element_text(angle = 30,hjust=1))+
  ylab("pct_utr_bases")+
  facet_grid(Annotation~Sample)

pct_intronic_bases.Plot<-ggplot(PicardMetrics,aes(x=Genome,y=pct_intronic_bases, fill=Replicate))+
  geom_bar(stat="identity", position = position_dodge())+
  CS.THEME+
  scale_fill_manual(values = c(CS.ColorPalette[11,1],CS.ColorPalette[12,1]))+
  theme(axis.text.x = element_text(angle = 30,hjust=1))+
  ylab("pct_intronic_bases")+
  facet_grid(Annotation~Sample)


pct_intergenic_bases.Plot<-ggplot(PicardMetrics,aes(x=Genome,y=pct_intergenic_bases, fill=Replicate))+
  geom_bar(stat="identity", position = position_dodge())+
  CS.THEME+
  scale_fill_manual(values = c(CS.ColorPalette[11,1],CS.ColorPalette[12,1]))+
  theme(axis.text.x = element_text(angle = 30,hjust=1))+
  ylab("pct_intergenic_bases")+
  facet_grid(Annotation~Sample)


pct_mrna_bases.Plot<-ggplot(PicardMetrics,aes(x=Genome,y=pct_mrna_bases, fill=Replicate))+
  geom_bar(stat="identity", position = position_dodge())+
  CS.THEME+
  scale_fill_manual(values = c(CS.ColorPalette[11,1],CS.ColorPalette[12,1]))+
  theme(axis.text.x = element_text(angle = 30,hjust=1))+
  ylab("pct_mrna_bases")+
  facet_grid(Annotation~Sample)

pct_usable_bases.Plot<-ggplot(PicardMetrics,aes(x=Genome,y=pct_usable_bases, fill=Replicate))+
  geom_bar(stat="identity", position = position_dodge())+
  CS.THEME+
  scale_fill_manual(values = c(CS.ColorPalette[11,1],CS.ColorPalette[12,1]))+
  theme(axis.text.x = element_text(angle = 30,hjust=1))+
  ylab("pct_usable_bases")+
  facet_grid(Annotation~Sample)


PicardTypesPanel<-plot_grid(
  pct_ribosomal_bases.Plot,
  pct_coding_bases.Plot,
  pct_utr_bases.Plot,
  pct_intronic_bases.Plot,
  pct_intergenic_bases.Plot,
  pct_mrna_bases.Plot,
  pct_usable_bases.Plot,
  labels="AUTO"
)

ggsave(filename = "Figures/ShortReadPicardTypesPanel.png",plot = PicardTypesPanel, units = "in", width = 30, height = 20)



median_cv_coverage.Plot<-ggplot(PicardMetrics,aes(x=Genome,y=median_cv_coverage, fill=Replicate))+
  geom_bar(stat="identity", position = position_dodge())+
  CS.THEME+
  scale_fill_manual(values = c(CS.ColorPalette[11,1],CS.ColorPalette[12,1]))+
  theme(axis.text.x = element_text(angle = 30,hjust=1))+
  ylab("median_cv_coverage")+
  facet_grid(Annotation~Sample)

median_5prime_bias.Plot<-ggplot(PicardMetrics,aes(x=Genome,y=median_5prime_bias, fill=Replicate))+
  geom_bar(stat="identity", position = position_dodge())+
  CS.THEME+
  scale_fill_manual(values = c(CS.ColorPalette[11,1],CS.ColorPalette[12,1]))+
  theme(axis.text.x = element_text(angle = 30,hjust=1))+
  ylab("median_5prime_bias")+
  facet_grid(Annotation~Sample)

median_3prime_bias.Plot<-ggplot(PicardMetrics,aes(x=Genome,y=median_3prime_bias, fill=Replicate))+
  geom_bar(stat="identity", position = position_dodge())+
  CS.THEME+
  scale_fill_manual(values = c(CS.ColorPalette[11,1],CS.ColorPalette[12,1]))+
  theme(axis.text.x = element_text(angle = 30,hjust=1))+
  ylab("median_3prime_bias")+
  facet_grid(Annotation~Sample)

median_5prime_to_3prime_bias.Plot<-ggplot(PicardMetrics,aes(x=Genome,y=median_5prime_to_3prime_bias, fill=Replicate))+
  geom_bar(stat="identity", position = position_dodge())+
  CS.THEME+
  scale_fill_manual(values = c(CS.ColorPalette[11,1],CS.ColorPalette[12,1]))+
  theme(axis.text.x = element_text(angle = 30,hjust=1))+
  ylab("median_5prime_to_3prime_bias")+
  facet_grid(Annotation~Sample)

CoverageBiasPicardPanel<-plot_grid(
  median_cv_coverage.Plot,
  median_5prime_bias.Plot,
  median_3prime_bias.Plot,
  median_5prime_to_3prime_bias.Plot,
  labels = "AUTO"
)

ggsave(filename = "Figures/ShortReadCoverageBiasPicardPanel.png",plot = CoverageBiasPicardPanel, units = "in", height = 12, width = 20)



Picard.PCA.Subset<-PicardMetrics[,c(13:22)]
Picard.PCA.Subset<-Picard.PCA.Subset[,!(colnames(Picard.PCA.Subset) %in% colnames(Picard.PCA.Subset)[colSums(Picard.PCA.Subset)==0])]

library(factoextra)

Picard.res.pca<-prcomp(Picard.PCA.Subset, scale=TRUE)
fviz_eig(Picard.res.pca)


groups.Sample <- as.factor(PicardMetrics$Sample)
Picard.PCA.Plot.Sample<-fviz_pca_ind(Picard.res.pca,
                              col.ind = groups.Sample,
                              palette = viridis(6, end=0.9),
                              repel = TRUE,
                              title="Samples - PCA"
)+
  CS.THEME+
  theme(legend.position="bottom")

Picard.PCA.Plot.Sample

groups <- as.factor(PicardMetrics$Genome)
Picard.PCA.Plot.Genome<-fviz_pca_ind(Picard.res.pca,
                              col.ind = groups,
                              palette = viridis(4,option = "magma", end = 0.8),
                              repel = TRUE,
                              title="Genome - PCA"
)+
  CS.THEME+
  theme(legend.position="bottom")

Picard.PCA.Plot.Genome

groups <- as.factor(PicardMetrics$Replicate)
Picard.PCA.Plot.Replicate<-fviz_pca_ind(Picard.res.pca,
                                 col.ind = groups,
                                 palette = viridis(2,option = "cividis", end = 0.8),
                                 repel = TRUE,
                                 title="Replicate - PCA"
)+
  CS.THEME+
  theme(legend.position="bottom")

Picard.PCA.Plot.Replicate

groups <- as.factor(PicardMetrics$Annotation)
Picard.PCA.Plot.Annotation<-fviz_pca_ind(Picard.res.pca,
                                  col.ind = groups,
                                  palette = viridis(4,option = "cividis", end = 0.8),
                                  repel = TRUE,
                                  title="Annotation - PCA"
)+
  CS.THEME+
  theme(legend.position="bottom")

Picard.PCA.Plot.Annotation

PicardMetrics.PCA.Panel<-plot_grid(Picard.PCA.Plot.Sample,Picard.PCA.Plot.Genome,Picard.PCA.Plot.Replicate, Picard.PCA.Plot.Annotation,labels = "AUTO",ncol=2)

ggsave2(filename = "Figures/PicardMetrics.PCA.Panel.png",plot = PicardMetrics.PCA.Panel, units = "in", width = 12, height = 12)




#### Combined QC and alignment metrics

subForPCA_Picard<-PicardMetrics[,c(2,3,4,5,13:22)]

subForPCA_Align<-AlignmentMetrics[,c(2,3,5,4,8:32)]

# Fix annotation names
subForPCA_Align$Annotation[subForPCA_Align$Annotation=="gencode"]<-"Gencode"

# Standardize reference names
subForPCA_Align$Genome[subForPCA_Align$Genome=="GRCh38"]<-"hg38"

# subset focused conditions
subForPCA_Picard<-subset(subForPCA_Picard, Sample != "HG002" & Annotation=="Gencode")
subForPCA_Align<-subset(subForPCA_Align, Sample != "HG002" & Annotation=="Gencode")

# remove annotations column
subForPCA_Align<-subForPCA_Align[,-4]
subForPCA_Picard<-subForPCA_Picard[,-4]

# order both objects

subForPCA_Picard<-subForPCA_Picard[order(paste(subForPCA_Picard$Sample,subForPCA_Picard$Genome,subForPCA_Picard$Replicate,sep = " ")),]
subForPCA_Align<-subForPCA_Align[order(paste(subForPCA_Align$Sample,subForPCA_Align$Genome,subForPCA_Align$Replicate,sep = " ")),]

# If properly ordered should all three return 0
sum(subForPCA_Align$Sample != subForPCA_Picard$Sample)
sum(subForPCA_Align$Genome != subForPCA_Picard$Genome)
sum(subForPCA_Align$Replicate != subForPCA_Picard$Replicate)



# Combine
PCA_Combined<-bind_cols(subForPCA_Picard,subForPCA_Align[,-1:-3])

PCA_Combined_Groups<-PCA_Combined[,1:3]
PCA_Combined<-PCA_Combined[,-1:-3]

PCA_Combined<-PCA_Combined[,!(colnames(PCA_Combined) %in% colnames(PCA_Combined)[colSums(PCA_Combined)==0])]


1:ncol(PCA_Combined)

CombinationPairs<-combinat::combn(1:ncol(PCA_Combined),2)

runCor<-function(DATA,COL1,COL2){
  rho=cor(x = DATA[,COL1], y = DATA[,COL2])
  return(
    data.frame(
      rho=rho,
      Col1=COL1,
      Col2=COL2
    )
  )
}


CorrelationData<-lapply(1:ncol(CombinationPairs), function(X) runCor(PCA_Combined,CombinationPairs[1,X],CombinationPairs[2,X]) ) %>% bind_rows()

CorrelationData$ColName1<-colnames(PCA_Combined)[CorrelationData$Col1]
CorrelationData$ColName2<-colnames(PCA_Combined)[CorrelationData$Col2]

subset(CorrelationData, abs(rho)>0.4 )

### I should filter these to remove columns which are highly correlate with other columns

ggplot(CorrelationData, aes(x=ColName1,y=ColName2,fill=abs(rho) ))+
  geom_tile()+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1))

library(factoextra)

dim(PCA_Combined)
Combined.res.pca<-prcomp(PCA_Combined, scale=TRUE)

#Plot.PCs<-fviz_eig(Combined.res.pca)
#Plot.PCs

#ggsave2(filename = "QC_PCs.png",plot = Plot.PCs,units = "in",width = 6/2, height = 6/2)

# Get the partial eta squared values
library(tidyr)
PCA.VarObject<-cbind(Combined.res.pca$x,PCA_Combined_Groups[,1:3])

NUM.PCs<-19 # Must decide how many PCs to include and set formula below appropriately
manova_model<-manova(cbind(PC1,PC2,PC3,PC4,PC5,PC6,PC7,PC8,PC9,PC10,PC11,PC12,PC13,PC14,PC15,PC16,PC17,PC18,PC19)~Sample*Genome+Replicate, PCA.VarObject)
manova_model.Summary<-summary(manova_model)
aov_results <- summary.aov(manova_model)

PC.PartialEtaSquareds<-lapply(1:NUM.PCs,
                              function(PC){aov_results[[PC]][[2]][1:5]/sum(aov_results[[PC]][[2]])}
) %>% do.call(rbind,.) %>% as.data.frame 

colnames(PC.PartialEtaSquareds)<-c("Sample","Genome","Sample:Genome","Replicate","Residual")

PC.Statistics<-cbind(as.factor(1:NUM.PCs),get_eig(res.pca)[1:NUM.PCs,],PC.PartialEtaSquareds)
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
  ggtitle("RNA-seq - QC metrics\n%Variance captured by PCs partitioned by factor")+
  CS.THEME+
  theme(legend.position="bottom")+
  scale_y_continuous(limits = c(0,round(max(PC.Statistics$VariancePercent)*1.05)), expand = c(0,0))+
  theme(legend.key.size = unit(0.2, "cm"))

ggsave2(filename = "QC_PCs.png",plot = Plot.PCs,units = "in",width = 7, height = 6)


# Plot PCAs
groups.Sample <- as.factor(PCA_Combined_Groups$Sample)
Combined.PCA.Plot.Sample<-fviz_pca_ind(Combined.res.pca,
                                     col.ind = groups.Sample,
                                     palette = viridis(6, end=0.9),
                                     repel = TRUE,
                                     title="RNA-seq - Samples"
)+
  CS.THEME+
  theme(legend.position="bottom",legend.text = element_text(size=8,hjust = 0.7), legend.title = element_text(size=12),legend.key.size = unit(0.25, "cm"))

Combined.PCA.Plot.Sample

groups <- as.factor(PCA_Combined_Groups$Genome)
Combined.PCA.Plot.Genome<-fviz_pca_ind(Combined.res.pca,
                                     col.ind = groups,
                                     palette = viridis(4,option = "magma", end = 0.8),
                                     repel = TRUE,
                                     title="RNA-seq - Genome"
)+
  CS.THEME+
  theme(legend.position="bottom",legend.text = element_text(size=8,hjust = 0.7), legend.title = element_text(size=12),legend.key.size = unit(0.25, "cm"))

Combined.PCA.Plot.Genome

groups <- as.factor(PCA_Combined_Groups$Replicate)
Combined.PCA.Plot.Replicate<-fviz_pca_ind(Combined.res.pca,
                                        col.ind = groups,
                                        palette = viridis(2,option = "cividis", end = 0.8),
                                        repel = TRUE,
                                        title="RNA-seq - Replicate"
)+
  CS.THEME+
  theme(legend.position="bottom",legend.text = element_text(size=8,hjust = 0.7), legend.title = element_text(size=12),legend.key.size = unit(0.25, "cm"))

Combined.PCA.Plot.Replicate

CombinedMetrics.PCA.Panel<-plot_grid(Combined.PCA.Plot.Sample,Combined.PCA.Plot.Genome,Combined.PCA.Plot.Replicate,labels = c("D","E","F"),ncol=3)

ggsave2(filename = "Figures/CombinedMetrics.PCA.Panel.png",plot = CombinedMetrics.PCA.Panel, units = "in", width = 17*0.75, height = 6*0.75)



