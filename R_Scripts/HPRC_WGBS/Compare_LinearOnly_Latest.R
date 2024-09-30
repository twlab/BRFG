
rm(list = ls())
set.seed(0)

library(ggplot2)
library(cowplot)
library(dplyr)
library(tidyr)
library(viridis)
source("/Users/juanmacias/Library/CloudStorage/Box-Box/LawsonLab/JuanMacias/Misc_Code/CS_Aesthetic.R")
library(factoextra)

setwd("/Users/juanmacias/Library/CloudStorage/Box-Box/LawsonLab/JuanMacias/Postdoc/MethylGrapher/")

### Import and bundle methylC files

# hg38
hg38.HG00621.A<-read.table("BenchmarkingSamplesOutputsLatest/Binned_AllIntersectSet_hg38.HG00621_lib2.CG.methylC.gz")[,c(1,2,4)]
colnames(hg38.HG00621.A)<-c("Chromosome","Start","Methylation")
hg38.HG00621.A$Sample<-"HG00621"
hg38.HG00621.A$Replicate<-"A"
hg38.HG00621.A$Method<-"hg38"

hg38.HG00621.B<-read.table("BenchmarkingSamplesOutputsLatest/Binned_AllIntersectSet_hg38.HG00621_lib3.CG.methylC.gz")[,c(1,2,4)]
colnames(hg38.HG00621.B)<-c("Chromosome","Start","Methylation")
hg38.HG00621.B$Sample<-"HG00621"
hg38.HG00621.B$Replicate<-"B"
hg38.HG00621.B$Method<-"hg38"

hg38.HG00741.A<-read.table("BenchmarkingSamplesOutputsLatest/Binned_AllIntersectSet_HG00741_BRep1_CG.graph.Precomputed.methylC.gz")[,c(1,2,4)]
colnames(hg38.HG00741.A)<-c("Chromosome","Start","Methylation")
hg38.HG00741.A$Sample<-"HG00741"
hg38.HG00741.A$Replicate<-"A"
hg38.HG00741.A$Method<-"hg38"

hg38.HG00741.B<-read.table("BenchmarkingSamplesOutputsLatest/Binned_AllIntersectSet_HG00741_BRep2_CG.graph.Precomputed.methylC.gz")[,c(1,2,4)]
colnames(hg38.HG00741.B)<-c("Chromosome","Start","Methylation")
hg38.HG00741.B$Sample<-"HG00741"
hg38.HG00741.B$Replicate<-"B"
hg38.HG00741.B$Method<-"hg38"

hg38.HG01952.A<-read.table("BenchmarkingSamplesOutputsLatest/Binned_AllIntersectSet_HG01952_BRep1_CG.graph.Precomputed.methylC.gz")[,c(1,2,4)]
colnames(hg38.HG01952.A)<-c("Chromosome","Start","Methylation")
hg38.HG01952.A$Sample<-"HG01952"
hg38.HG01952.A$Replicate<-"A"
hg38.HG01952.A$Method<-"hg38"

hg38.HG01952.B<-read.table("BenchmarkingSamplesOutputsLatest/Binned_AllIntersectSet_HG01952_BRep2_CG.graph.Precomputed.methylC.gz")[,c(1,2,4)]
colnames(hg38.HG01952.B)<-c("Chromosome","Start","Methylation")
hg38.HG01952.B$Sample<-"HG01952"
hg38.HG01952.B$Replicate<-"B"
hg38.HG01952.B$Method<-"hg38"

hg38.HG01978.A<-read.table("BenchmarkingSamplesOutputsLatest/Binned_AllIntersectSet_HG01978_BRep1_CG.graph.Precomputed.methylC.gz")[,c(1,2,4)]
colnames(hg38.HG01978.A)<-c("Chromosome","Start","Methylation")
hg38.HG01978.A$Sample<-"HG01978"
hg38.HG01978.A$Replicate<-"A"
hg38.HG01978.A$Method<-"hg38"

hg38.HG01978.B<-read.table("BenchmarkingSamplesOutputsLatest/Binned_AllIntersectSet_HG01978_BRep2_CG.graph.Precomputed.methylC.gz")[,c(1,2,4)]
colnames(hg38.HG01978.B)<-c("Chromosome","Start","Methylation")
hg38.HG01978.B$Sample<-"HG01978"
hg38.HG01978.B$Replicate<-"B"
hg38.HG01978.B$Method<-"hg38"

hg38.HG03516.A<-read.table("BenchmarkingSamplesOutputsLatest/Binned_AllIntersectSet_HG03516_BRep1_CG.graph.Precomputed.methylC.gz")[,c(1,2,4)]
colnames(hg38.HG03516.A)<-c("Chromosome","Start","Methylation")
hg38.HG03516.A$Sample<-"HG03516"
hg38.HG03516.A$Replicate<-"A"
hg38.HG03516.A$Method<-"hg38"

hg38.HG03516.B<-read.table("BenchmarkingSamplesOutputsLatest/Binned_AllIntersectSet_HG03516_BRep2_CG.graph.Precomputed.methylC.gz")[,c(1,2,4)]
colnames(hg38.HG03516.B)<-c("Chromosome","Start","Methylation")
hg38.HG03516.B$Sample<-"HG03516"
hg38.HG03516.B$Replicate<-"B"
hg38.HG03516.B$Method<-"hg38"

hg38Combined<-rbind(
  hg38.HG00621.A,
  hg38.HG00741.A,
  hg38.HG01952.A,
  hg38.HG01978.A,
  hg38.HG03516.A,
  hg38.HG00621.B,
  hg38.HG00741.B,
  hg38.HG01952.B,
  hg38.HG01978.B,
  hg38.HG03516.B
)

rm(
  hg38.HG00621.A,
  hg38.HG00741.A,
  hg38.HG01952.A,
  hg38.HG01978.A,
  hg38.HG03516.A,
  hg38.HG00621.B,
  hg38.HG00741.B,
  hg38.HG01952.B,
  hg38.HG01978.B,
  hg38.HG03516.B
)

table(hg38Combined[,c(4:6)])


# CHM13
CHM13.HG00621.A<-read.table("BenchmarkingSamplesOutputsLatest/Binned_AllIntersectSet_HG00621_BRep1_CG.graph.Precomputed.methylC.gz")[,c(1,2,4)]
colnames(CHM13.HG00621.A)<-c("Chromosome","Start","Methylation")
CHM13.HG00621.A$Sample<-"HG00621"
CHM13.HG00621.A$Replicate<-"A"
CHM13.HG00621.A$Method<-"CHM13"

CHM13.HG00621.B<-read.table("BenchmarkingSamplesOutputsLatest/Binned_AllIntersectSet_HG00621_BRep2_CG.graph.Precomputed.methylC.gz")[,c(1,2,4)]
colnames(CHM13.HG00621.B)<-c("Chromosome","Start","Methylation")
CHM13.HG00621.B$Sample<-"HG00621"
CHM13.HG00621.B$Replicate<-"B"
CHM13.HG00621.B$Method<-"CHM13"

CHM13.HG00741.A<-read.table("BenchmarkingSamplesOutputsLatest/Binned_AllIntersectSet_HG00741_BRep1_CG.graph.Precomputed.methylC.gz")[,c(1,2,4)]
colnames(CHM13.HG00741.A)<-c("Chromosome","Start","Methylation")
CHM13.HG00741.A$Sample<-"HG00741"
CHM13.HG00741.A$Replicate<-"A"
CHM13.HG00741.A$Method<-"CHM13"

CHM13.HG00741.B<-read.table("BenchmarkingSamplesOutputsLatest/Binned_AllIntersectSet_HG00741_BRep2_CG.graph.Precomputed.methylC.gz")[,c(1,2,4)]
colnames(CHM13.HG00741.B)<-c("Chromosome","Start","Methylation")
CHM13.HG00741.B$Sample<-"HG00741"
CHM13.HG00741.B$Replicate<-"B"
CHM13.HG00741.B$Method<-"CHM13"

CHM13.HG01952.A<-read.table("BenchmarkingSamplesOutputsLatest/Binned_AllIntersectSet_HG01952_BRep1_CG.graph.Precomputed.methylC.gz")[,c(1,2,4)]
colnames(CHM13.HG01952.A)<-c("Chromosome","Start","Methylation")
CHM13.HG01952.A$Sample<-"HG01952"
CHM13.HG01952.A$Replicate<-"A"
CHM13.HG01952.A$Method<-"CHM13"

CHM13.HG01952.B<-read.table("BenchmarkingSamplesOutputsLatest/Binned_AllIntersectSet_HG01952_BRep2_CG.graph.Precomputed.methylC.gz")[,c(1,2,4)]
colnames(CHM13.HG01952.B)<-c("Chromosome","Start","Methylation")
CHM13.HG01952.B$Sample<-"HG01952"
CHM13.HG01952.B$Replicate<-"B"
CHM13.HG01952.B$Method<-"CHM13"

CHM13.HG01978.A<-read.table("BenchmarkingSamplesOutputsLatest/Binned_AllIntersectSet_HG01978_BRep1_CG.graph.Precomputed.methylC.gz")[,c(1,2,4)]
colnames(CHM13.HG01978.A)<-c("Chromosome","Start","Methylation")
CHM13.HG01978.A$Sample<-"HG01978"
CHM13.HG01978.A$Replicate<-"A"
CHM13.HG01978.A$Method<-"CHM13"

CHM13.HG01978.B<-read.table("BenchmarkingSamplesOutputsLatest/Binned_AllIntersectSet_HG01978_BRep2_CG.graph.Precomputed.methylC.gz")[,c(1,2,4)]
colnames(CHM13.HG01978.B)<-c("Chromosome","Start","Methylation")
CHM13.HG01978.B$Sample<-"HG01978"
CHM13.HG01978.B$Replicate<-"B"
CHM13.HG01978.B$Method<-"CHM13"

CHM13.HG03516.A<-read.table("BenchmarkingSamplesOutputsLatest/Binned_AllIntersectSet_HG03516_BRep1_CG.graph.Precomputed.methylC.gz")[,c(1,2,4)]
colnames(CHM13.HG03516.A)<-c("Chromosome","Start","Methylation")
CHM13.HG03516.A$Sample<-"HG03516"
CHM13.HG03516.A$Replicate<-"A"
CHM13.HG03516.A$Method<-"CHM13"

CHM13.HG03516.B<-read.table("BenchmarkingSamplesOutputsLatest/Binned_AllIntersectSet_HG03516_BRep2_CG.graph.Precomputed.methylC.gz")[,c(1,2,4)]
colnames(CHM13.HG03516.B)<-c("Chromosome","Start","Methylation")
CHM13.HG03516.B$Sample<-"HG03516"
CHM13.HG03516.B$Replicate<-"B"
CHM13.HG03516.B$Method<-"CHM13"

CHM13Combined<-rbind(
  CHM13.HG00621.A,
  CHM13.HG00741.A,
  CHM13.HG01952.A,
  CHM13.HG01978.A,
  CHM13.HG03516.A,
  CHM13.HG00621.B,
  CHM13.HG00741.B,
  CHM13.HG01952.B,
  CHM13.HG01978.B,
  CHM13.HG03516.B
)

rm(
  CHM13.HG00621.A,
  CHM13.HG00741.A,
  CHM13.HG01952.A,
  CHM13.HG01978.A,
  CHM13.HG03516.A,
  CHM13.HG00621.B,
  CHM13.HG00741.B,
  CHM13.HG01952.B,
  CHM13.HG01978.B,
  CHM13.HG03516.B
)

table(CHM13Combined[,c(4:6)])


# mat
mat.HG00621.A<-read.table("BenchmarkingSamplesOutputsLatest/Binned_AllIntersectSet_HG00621_BRep1_CG.graph.Precomputed.methylC.gz")[,c(1,2,4)]
colnames(mat.HG00621.A)<-c("Chromosome","Start","Methylation")
mat.HG00621.A$Sample<-"HG00621"
mat.HG00621.A$Replicate<-"A"
mat.HG00621.A$Method<-"mat"

mat.HG00621.B<-read.table("BenchmarkingSamplesOutputsLatest/Binned_AllIntersectSet_HG00621_BRep2_CG.graph.Precomputed.methylC.gz")[,c(1,2,4)]
colnames(mat.HG00621.B)<-c("Chromosome","Start","Methylation")
mat.HG00621.B$Sample<-"HG00621"
mat.HG00621.B$Replicate<-"B"
mat.HG00621.B$Method<-"mat"

mat.HG00741.A<-read.table("BenchmarkingSamplesOutputsLatest/Binned_AllIntersectSet_HG00741_BRep1_CG.graph.Precomputed.methylC.gz")[,c(1,2,4)]
colnames(mat.HG00741.A)<-c("Chromosome","Start","Methylation")
mat.HG00741.A$Sample<-"HG00741"
mat.HG00741.A$Replicate<-"A"
mat.HG00741.A$Method<-"mat"

mat.HG00741.B<-read.table("BenchmarkingSamplesOutputsLatest/Binned_AllIntersectSet_HG00741_BRep2_CG.graph.Precomputed.methylC.gz")[,c(1,2,4)]
colnames(mat.HG00741.B)<-c("Chromosome","Start","Methylation")
mat.HG00741.B$Sample<-"HG00741"
mat.HG00741.B$Replicate<-"B"
mat.HG00741.B$Method<-"mat"

mat.HG01952.A<-read.table("BenchmarkingSamplesOutputsLatest/Binned_AllIntersectSet_HG01952_BRep1_CG.graph.Precomputed.methylC.gz")[,c(1,2,4)]
colnames(mat.HG01952.A)<-c("Chromosome","Start","Methylation")
mat.HG01952.A$Sample<-"HG01952"
mat.HG01952.A$Replicate<-"A"
mat.HG01952.A$Method<-"mat"

mat.HG01952.B<-read.table("BenchmarkingSamplesOutputsLatest/Binned_AllIntersectSet_HG01952_BRep2_CG.graph.Precomputed.methylC.gz")[,c(1,2,4)]
colnames(mat.HG01952.B)<-c("Chromosome","Start","Methylation")
mat.HG01952.B$Sample<-"HG01952"
mat.HG01952.B$Replicate<-"B"
mat.HG01952.B$Method<-"mat"

mat.HG01978.A<-read.table("BenchmarkingSamplesOutputsLatest/Binned_AllIntersectSet_HG01978_BRep1_CG.graph.Precomputed.methylC.gz")[,c(1,2,4)]
colnames(mat.HG01978.A)<-c("Chromosome","Start","Methylation")
mat.HG01978.A$Sample<-"HG01978"
mat.HG01978.A$Replicate<-"A"
mat.HG01978.A$Method<-"mat"

mat.HG01978.B<-read.table("BenchmarkingSamplesOutputsLatest/Binned_AllIntersectSet_HG01978_BRep2_CG.graph.Precomputed.methylC.gz")[,c(1,2,4)]
colnames(mat.HG01978.B)<-c("Chromosome","Start","Methylation")
mat.HG01978.B$Sample<-"HG01978"
mat.HG01978.B$Replicate<-"B"
mat.HG01978.B$Method<-"mat"

mat.HG03516.A<-read.table("BenchmarkingSamplesOutputsLatest/Binned_AllIntersectSet_HG03516_BRep1_CG.graph.Precomputed.methylC.gz")[,c(1,2,4)]
colnames(mat.HG03516.A)<-c("Chromosome","Start","Methylation")
mat.HG03516.A$Sample<-"HG03516"
mat.HG03516.A$Replicate<-"A"
mat.HG03516.A$Method<-"mat"

mat.HG03516.B<-read.table("BenchmarkingSamplesOutputsLatest/Binned_AllIntersectSet_HG03516_BRep2_CG.graph.Precomputed.methylC.gz")[,c(1,2,4)]
colnames(mat.HG03516.B)<-c("Chromosome","Start","Methylation")
mat.HG03516.B$Sample<-"HG03516"
mat.HG03516.B$Replicate<-"B"
mat.HG03516.B$Method<-"mat"

matCombined<-rbind(
  mat.HG00621.A,
  mat.HG00741.A,
  mat.HG01952.A,
  mat.HG01978.A,
  mat.HG03516.A,
  mat.HG00621.B,
  mat.HG00741.B,
  mat.HG01952.B,
  mat.HG01978.B,
  mat.HG03516.B
)

rm(
  mat.HG00621.A,
  mat.HG00741.A,
  mat.HG01952.A,
  mat.HG01978.A,
  mat.HG03516.A,
  mat.HG00621.B,
  mat.HG00741.B,
  mat.HG01952.B,
  mat.HG01978.B,
  mat.HG03516.B
)

table(matCombined[,c(4:6)])


# pat
pat.HG00621.A<-read.table("BenchmarkingSamplesOutputsLatest/Binned_AllIntersectSet_HG00621_BRep1_CG.graph.Precomputed.methylC.gz")[,c(1,2,4)]
colnames(pat.HG00621.A)<-c("Chromosome","Start","Methylation")
pat.HG00621.A$Sample<-"HG00621"
pat.HG00621.A$Replicate<-"A"
pat.HG00621.A$Method<-"pat"

pat.HG00621.B<-read.table("BenchmarkingSamplesOutputsLatest/Binned_AllIntersectSet_HG00621_BRep2_CG.graph.Precomputed.methylC.gz")[,c(1,2,4)]
colnames(pat.HG00621.B)<-c("Chromosome","Start","Methylation")
pat.HG00621.B$Sample<-"HG00621"
pat.HG00621.B$Replicate<-"B"
pat.HG00621.B$Method<-"pat"

pat.HG00741.A<-read.table("BenchmarkingSamplesOutputsLatest/Binned_AllIntersectSet_HG00741_BRep1_CG.graph.Precomputed.methylC.gz")[,c(1,2,4)]
colnames(pat.HG00741.A)<-c("Chromosome","Start","Methylation")
pat.HG00741.A$Sample<-"HG00741"
pat.HG00741.A$Replicate<-"A"
pat.HG00741.A$Method<-"pat"

pat.HG00741.B<-read.table("BenchmarkingSamplesOutputsLatest/Binned_AllIntersectSet_HG00741_BRep2_CG.graph.Precomputed.methylC.gz")[,c(1,2,4)]
colnames(pat.HG00741.B)<-c("Chromosome","Start","Methylation")
pat.HG00741.B$Sample<-"HG00741"
pat.HG00741.B$Replicate<-"B"
pat.HG00741.B$Method<-"pat"

pat.HG01952.A<-read.table("BenchmarkingSamplesOutputsLatest/Binned_AllIntersectSet_HG01952_BRep1_CG.graph.Precomputed.methylC.gz")[,c(1,2,4)]
colnames(pat.HG01952.A)<-c("Chromosome","Start","Methylation")
pat.HG01952.A$Sample<-"HG01952"
pat.HG01952.A$Replicate<-"A"
pat.HG01952.A$Method<-"pat"

pat.HG01952.B<-read.table("BenchmarkingSamplesOutputsLatest/Binned_AllIntersectSet_HG01952_BRep2_CG.graph.Precomputed.methylC.gz")[,c(1,2,4)]
colnames(pat.HG01952.B)<-c("Chromosome","Start","Methylation")
pat.HG01952.B$Sample<-"HG01952"
pat.HG01952.B$Replicate<-"B"
pat.HG01952.B$Method<-"pat"

pat.HG01978.A<-read.table("BenchmarkingSamplesOutputsLatest/Binned_AllIntersectSet_HG01978_BRep1_CG.graph.Precomputed.methylC.gz")[,c(1,2,4)]
colnames(pat.HG01978.A)<-c("Chromosome","Start","Methylation")
pat.HG01978.A$Sample<-"HG01978"
pat.HG01978.A$Replicate<-"A"
pat.HG01978.A$Method<-"pat"

pat.HG01978.B<-read.table("BenchmarkingSamplesOutputsLatest/Binned_AllIntersectSet_HG01978_BRep2_CG.graph.Precomputed.methylC.gz")[,c(1,2,4)]
colnames(pat.HG01978.B)<-c("Chromosome","Start","Methylation")
pat.HG01978.B$Sample<-"HG01978"
pat.HG01978.B$Replicate<-"B"
pat.HG01978.B$Method<-"pat"

pat.HG03516.A<-read.table("BenchmarkingSamplesOutputsLatest/Binned_AllIntersectSet_HG03516_BRep1_CG.graph.Precomputed.methylC.gz")[,c(1,2,4)]
colnames(pat.HG03516.A)<-c("Chromosome","Start","Methylation")
pat.HG03516.A$Sample<-"HG03516"
pat.HG03516.A$Replicate<-"A"
pat.HG03516.A$Method<-"pat"

pat.HG03516.B<-read.table("BenchmarkingSamplesOutputsLatest/Binned_AllIntersectSet_HG03516_BRep2_CG.graph.Precomputed.methylC.gz")[,c(1,2,4)]
colnames(pat.HG03516.B)<-c("Chromosome","Start","Methylation")
pat.HG03516.B$Sample<-"HG03516"
pat.HG03516.B$Replicate<-"B"
pat.HG03516.B$Method<-"pat"

patCombined<-rbind(
  pat.HG00621.A,
  pat.HG00741.A,
  pat.HG01952.A,
  pat.HG01978.A,
  pat.HG03516.A,
  pat.HG00621.B,
  pat.HG00741.B,
  pat.HG01952.B,
  pat.HG01978.B,
  pat.HG03516.B
)

rm(
  pat.HG00621.A,
  pat.HG00741.A,
  pat.HG01952.A,
  pat.HG01978.A,
  pat.HG03516.A,
  pat.HG00621.B,
  pat.HG00741.B,
  pat.HG01952.B,
  pat.HG01978.B,
  pat.HG03516.B
)

table(patCombined[,c(4:6)])

print("Test for any differences, 0 means they match correctly")

sum(!(table(matCombined[,c(4:6)])==table(patCombined[,c(4:6)])))
sum(!(table(matCombined[,c(4:6)])==table(CHM13Combined[,c(4:6)])))
sum(!(table(matCombined[,c(4:6)])==table(hg38Combined[,c(4:6)])))




table(hg38Combined[,c(4:6)])
table(CHM13Combined[,c(4:6)])
table(patCombined[,c(4:6)])
table(matCombined[,c(4:6)])



AllCombined<-rbind(hg38Combined,CHM13Combined,patCombined,matCombined)

rm(hg38Combined,CHM13Combined,patCombined,matCombined)

table(AllCombined[,c(4:6)])

AllCombined$Key<-paste(AllCombined$Chromosome,AllCombined$Start, sep="_")

# Filter out empty bins
AllCombined<-subset(AllCombined, Methylation != ".")

# Ensure value types
AllCombined$Methylation<-as.numeric(AllCombined$Methylation)

# Spread
AllCombined.Wide<-spread(AllCombined[,c(3:7)],Key,Methylation)

## Remove Long formatted file...
rm(AllCombined)


### Ignore hg38
#AllCombined.Wide<-AllCombined.Wide[AllCombined.Wide$Method!="hg38",]

### Rename
#AllCombined.Wide$Method[AllCombined.Wide$Method=="Bismark"]<-"hg38"
AllCombined.Wide$Method[AllCombined.Wide$Method=="mat"]<-"Maternal"
AllCombined.Wide$Method[AllCombined.Wide$Method=="pat"]<-"Paternal"

### Make PCA and plot
PCA.Subset<-AllCombined.Wide[,c(4:ncol(AllCombined.Wide))]


print("Filter to exclude NAs")
print("Sites before filtering...")
print(dim(PCA.Subset))

PCA.Subset<-PCA.Subset[,unlist(apply(PCA.Subset, 2, function(X) sum(!is.na(X))==nrow(PCA.Subset) ))]

print("Sites after filtering...")
print(dim(PCA.Subset))

print("Filter to exclude zero variance sites")
print("Sites before filtering...")
print(dim(PCA.Subset))

PCA.Subset<-PCA.Subset[,unlist(apply(PCA.Subset, 2, function(X) sd(X)!=0 ))]
print("Sites after filtering...")
print(dim(PCA.Subset))


## Run PCA
res.pca<-prcomp(PCA.Subset, scale=TRUE)

groups.Sample <- as.factor(AllCombined.Wide$Sample)
groups.Genome <- as.factor(AllCombined.Wide$Method)
groups.Replicate <- as.factor(AllCombined.Wide$Replicate)


# Get the partial eta squared values
library(tidyr)
PCA.VarObject<-cbind(res.pca$x,data.frame(Sample=groups.Sample),data.frame(Replicate=groups.Replicate),data.frame(Genome=groups.Genome))

manova_model<-manova(cbind(PC1,PC2,PC3,PC4,PC5)~Sample*Genome+Replicate,PCA.VarObject)
manova_model.Summary<-summary(manova_model)
aov_results <- summary.aov(manova_model)

NUM.PCs<-5
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
write.table(PC.Statistics, file = "PC.Statistics.txt", sep = "\t", row.names = FALSE, quote = FALSE)

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
  ggtitle("WGBS - Methylation\n%Variance captured by PCs partitioned by factor")+
  CS.THEME+
  theme(legend.position="bottom")+
  scale_y_continuous(limits = c(0,round(max(PC.Statistics$VariancePercent)*1.05)), expand = c(0,0))+
  theme(legend.key.size = unit(0.2, "cm"))

ggsave2(filename = "Methylation_PCs.png",plot = Plot.PCs,units = "in",width = 9, height = 6)




#Plot.PCs<-fviz_eig(res.pca)

#ggsave2(filename = "AllRuns_PCs.png",plot = Plot.PCs,units = "in",width = 6/2, height = 6/2)



PCA.Plot.Sample<-fviz_pca_ind(res.pca,
                              col.ind = groups.Sample,
                              palette = viridis(5, end=0.9),
                              repel = TRUE,
                              addEllipses = TRUE,
                              ellipse.type = "convex",
                              title="Methylation - Samples"
)+
  CS.THEME+
  theme(legend.position="right")

#### Edit in pangenome label


PCA.Plot.Genome<-fviz_pca_ind(res.pca,
                              col.ind = groups.Genome,
                              palette = viridis(5,option = "magma", end = 0.8),
                              repel = TRUE,
                              title="Methylation - Genome"
)+
  CS.THEME+
  theme(legend.position="right")


PCA.Plot.Replicate<-fviz_pca_ind(res.pca,
                                 col.ind = groups.Replicate,
                                 palette = viridis(2,option = "cividis", end = 0.8),
                                 repel = TRUE,
                                 title="Methylation - Replicate"
)+
  CS.THEME+
  theme(legend.position="right")

PCA.Panel<-plot_grid(PCA.Plot.Sample,PCA.Plot.Genome,PCA.Plot.Replicate,labels = c("G","H","I"),ncol=3)

ggsave2(filename = "AllRuns_OnReferenceCpGMethylation.PCA.Panel.png",plot = PCA.Panel, units = "in", width = 17*1, height = 6*0.75)
