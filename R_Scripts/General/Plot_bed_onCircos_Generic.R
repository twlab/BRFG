
rm(list = ls())
set.seed(123)
# Load necessary libraries
library(circlize)
library(readr)

# Capture command line arguments
args <- commandArgs(trailingOnly = TRUE)

# Assuming the first argument is the path to the BED file and the second is the assembly name
#bed_file_path <- args[1]
#assembly_name <- args[2]

bed_file_path <- "/Users/juanmacias/Library/CloudStorage/Box-Box/LawsonLab/JuanMacias/Postdoc/Assembly_Benchmarking/ContigAlignments/HG00621.HG00621.hap1.maternal.noCS.resorted.bed.gz"
assembly_name <- "HG00621.maternal"



# Path to the chromosome lengths file
chromosome_lengths_path <- "/Users/juanmacias/Library/CloudStorage/Box-Box/LawsonLab/JuanMacias/Postdoc/Assembly_Benchmarking/ChromosomeSizes/GenomeLengths.csv"

# Read chromosome lengths and select the relevant column based on the assembly name
read_chromosome_lengths <- function(lengths_file_path, assembly_name) {
  chromosome_lengths_df <- read_csv(lengths_file_path)
  if (!assembly_name %in% colnames(chromosome_lengths_df)) {
    stop("Assembly name not found in the chromosome lengths file.")
  }
  chromosome_lengths <- data.frame(
    chrom = chromosome_lengths_df[[1]], # Assuming the first column is chromosome names
    length = chromosome_lengths_df[[assembly_name]]
  )
  return(chromosome_lengths)
}

# Initialize circos plot with custom chromosomes using the specific assembly
chromosome_lengths <- read_chromosome_lengths(chromosome_lengths_path, assembly_name)

# Exclude chromosomes with NA lengths
valid_chromosomes <- !is.na(chromosome_lengths$length)
chromosome_lengths <- chromosome_lengths[valid_chromosomes, ]




# Convert chromosome lengths to a matrix with 2 columns for xlim
xlim_matrix <- matrix(c(rep(0, nrow(chromosome_lengths)), chromosome_lengths$length), ncol = 2, byrow = FALSE)
rownames(xlim_matrix) <- chromosome_lengths$chrom

# Initialize circos plot
circos.initialize(factors = chromosome_lengths$chrom, xlim = xlim_matrix)

as.list(setNames(chromosome_lengths$length, chromosome_lengths$chrom))

# Ensure to define the function read_compressed_bed_and_prepare_for_circlize here
# Function: read_compressed_bed_and_prepare_for_circlize(bed_file_path) { ... }


bed_list<-read.delim("/Users/juanmacias/Library/CloudStorage/Box-Box/LawsonLab/JuanMacias/Postdoc/Assembly_Benchmarking/ContigAlignments/HG00621.HG00621.hap1.maternal.noCS.resorted.bed", header = FALSE, sep = "\t", col.names = c("chr", "start", "end", "name", "score", "strand"))


# Prepare data
#bed_list <- read_compressed_bed_and_prepare_for_circlize(bed_file_path)

as.data.frame(bed_list)


names(bed_list)




bed_list 


#bed = generateRandomBed(nr = 200)

circos.genomicTrackPlotRegion(as.data.frame(bed_list))
circos.genomicRect(as.data.frame(bed_list)[,2:3],value = as.data.frame(bed_list)[,5])


circos.genomicRect(as.data.frame(bed_list$chr1)[,2:3],value = as.data.frame(bed_list$chr1)[,5])


circos.genomicTrackPlotRegion(as.data.frame(bed_list$chr2))
circos.genomicRect(as.data.frame(bed_list$chr2)[,2:3],value = as.data.frame(bed_list$chr2)[,5])



# Plot each chromosome segment from BED
for(chrom in names(bed_list)) {
  for(i in 1:nrow(bed_list[[chrom]])) {
    circos.genomicRect(region = as.data.frame(bed_list[[chrom]])[,1:2])
  }
}

# Customize the plot if needed
# circos.trackPlotRegion(...)
# circos.text(...)
# circos.link(...)

# Clear circos plot after plotting to reset parameters
circos.clear()