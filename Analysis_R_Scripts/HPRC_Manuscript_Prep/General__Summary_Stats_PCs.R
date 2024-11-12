
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

source("/Users/juanmacias/Library/CloudStorage/Box-Box/LawsonLab/JuanMacias/Misc_Code/CS_Aesthetic.R")
setwd("/Users/juanmacias/Library/CloudStorage/Box-Box/LawsonLab/JuanMacias/Postdoc/Manuscripts/BenchmarkingPaper/Supplementary_Tables_Files/")

### Load data QC/Alignment
#### QC metrics
ATACseq.QC<-read.delim("/Users/juanmacias/Library/CloudStorage/Box-Box/LawsonLab/JuanMacias/Postdoc/HPRC/ATACseq/AutomaticTesting/Latest/PC.Statistics.txt", header = TRUE, stringsAsFactors = FALSE)
RNAseq.QC<-read.delim("/Users/juanmacias/Library/CloudStorage/Box-Box/LawsonLab/JuanMacias/Postdoc/HPRC/RNAseq/PC.QC.Statistics.txt", header = TRUE, stringsAsFactors = FALSE)
WGBS.QC<-read.delim("/Users/juanmacias/Library/CloudStorage/Box-Box/LawsonLab/JuanMacias/Postdoc/HPRC/WGBS/PC.QC.Statistics.txt", header = TRUE, stringsAsFactors = FALSE)
HiC.QC<-read.delim("/Users/juanmacias/Library/CloudStorage/Box-Box/LawsonLab/JuanMacias/Postdoc/HPRC/HiC/PC.Statistics.txt", header = TRUE, stringsAsFactors = FALSE)
#### Functional outputs
ATACseq<-read.delim("/Users/juanmacias/Library/CloudStorage/Box-Box/LawsonLab/JuanMacias/Postdoc/HPRC/ATACseq/AutomaticTesting/Latest/PeakComparison/PC.Statistics.txt", header = TRUE, stringsAsFactors = FALSE)
RNAseq<-read.delim("/Users/juanmacias/Library/CloudStorage/Box-Box/LawsonLab/JuanMacias/Postdoc/HPRC/RNAseq/PC.NoPangenome.Statistics.txt", header = TRUE, stringsAsFactors = FALSE)
WGBS<-read.delim("/Users/juanmacias/Library/CloudStorage/Box-Box/LawsonLab/JuanMacias/Postdoc/MethylGrapher/PC.Statistics.txt", header = TRUE, stringsAsFactors = FALSE)


### ATACseq % metric var captured by all PCs
sum(ATACseq.QC[,3])
### ATACseq % of captured var attributable to sample
sum(ATACseq.QC[,10])
### ATACseq % of captured var attributable to sample:genome
sum(ATACseq.QC[,11])
### ATACseq % of captured var attributable to genome
sum(ATACseq.QC[,12])

### RNAseq % metric var captured by all PCs
sum(RNAseq.QC[,3])
### RNAseq % of captured var attributable to sample
sum(RNAseq.QC[,10])
### RNAseq % of captured var attributable to sample:genome
sum(RNAseq.QC[,11])
### RNAseq % of captured var attributable to genome
sum(RNAseq.QC[,12])

### WGBS % metric var captured by all PCs
sum(WGBS.QC[,3])
### WGBS % of captured var attributable to sample
sum(WGBS.QC[,10])
### WGBS % of captured var attributable to sample:genome
sum(WGBS.QC[,11])
### WGBS % of captured var attributable to genome
sum(WGBS.QC[,12])

### HiC % metric var captured by all PCs
sum(HiC.QC[,3])
### HiC % of captured var attributable to sample
sum(HiC.QC[,8])
### HiC % of captured var attributable to sample:genome
sum(HiC.QC[,9])
### HiC % of captured var attributable to genome
sum(HiC.QC[,10])


### ATACseq % function value var captured by all PCs
sum(ATACseq[,3])
### ATACseq % of captured var attributable to sample
sum(ATACseq[,10])
### ATACseq % of captured var attributable to sample:genome
sum(ATACseq[,11])
### ATACseq % of captured var attributable to genome
sum(ATACseq[,12])

### RNAseq % function value var captured by all PCs
sum(RNAseq[,3])
### RNAseq % of captured var attributable to sample
sum(RNAseq[,10])
### RNAseq % of captured var attributable to sample:genome
sum(RNAseq[,11])
### RNAseq % of captured var attributable to genome
sum(RNAseq[,12])

### WGBS % function value var captured by all PCs
sum(WGBS[,3])
### WGBS % of captured var attributable to sample
sum(WGBS[,10])
### WGBS % of captured var attributable to sample:genome
sum(WGBS[,11])
### WGBS % of captured var attributable to genome
sum(WGBS[,12])
