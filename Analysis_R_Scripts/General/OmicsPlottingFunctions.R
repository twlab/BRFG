


library(ggplot2)
require(cowplot)
library(scales)
library(gridExtra)
library(grid)
library(dplyr)
library(tidyr)
library(zoo)

source("/Users/juanmacias/Library/CloudStorage/Box-Box/LawsonLab/JuanMacias/Misc_Code/CS_Aesthetic.R")

# Import Chromosome Sizes
ChromosomeSize<-read.delim("/Users/juanmacias/Library/CloudStorage/Box-Box/LawsonLab/JuanMacias/Postdoc/HPRC/HPRC_General/GenomeLengths.csv",sep = ",")

# In this R script I will define and describe R functions useful for plotting omics results

## manhattan like plot
### With optional highlights
### With optional vertical and horizontal bars
### With option for custom chromosome sizes
### Plotting per chromosome

ChromoAtlasScatter<-function(
    OBJ,
    CHROMOSOME,
    HIGHLIGHT=NULL,
    VERTICALLINES=NULL, 
    HORIZONTALLINES=NULL,
    CHROMSIZES=NA,
    Xlabel=NA,
    Ylabel=NA,
    PointSize=1,
    PositionScalingFactor=1,
    ValueScalingFactor=1,
    FitGAM=FALSE,
    GAMLineColor="black",
    FILLVARIABLE=FALSE,
    PointAlpha=1,
    FillMax=NULL,
    FillMin=NULL
    ){
  
  # Check whether custom chromosomes will be used. If not use hg38 main chromosomes.
  if(is.na(CHROMSIZES)){
    
    CHROMSIZES<-data.frame(
      Chr=ChromosomeSize$Chromsome,
      Length=ChromosomeSize$hg38
    )
    
  }
  
  # Check if CHROMOSOME exists in CHROMSIZES, if not return error
  if(CHROMOSOME %in% CHROMSIZES[,1]){
    
    if(CHROMOSOME %in% OBJ[,1]){
      
      if(FILLVARIABLE){
        
        # Rename OBJ columns
        colnames(OBJ)<-c("Chromosome","Position","Value","ColorValue")
        OBJ$Position<-as.numeric(as.character(OBJ$Position))
        OBJ$Value<-as.numeric(as.character(OBJ$Value))
        OBJ$ColorValue<-as.numeric(as.character(OBJ$ColorValue))
        
      } else {
        
        # Rename OBJ columns
        colnames(OBJ)<-c("Chromosome","Position","Value")
        OBJ$Position<-as.numeric(as.character(OBJ$Position))
        OBJ$Value<-as.numeric(as.character(OBJ$Value))
        
      }
      
      # Define CHROMLENGTH
      CHROMLENGTH<-CHROMSIZES[CHROMSIZES[,1]==CHROMOSOME,2]
      
      # Subset down to selected CHROMOSOME
      OBJ<-subset(OBJ, Chromosome==CHROMOSOME)
      
      if(FILLVARIABLE){
        
        if(!is.null(FillMax) & !is.null(FillMin)){
          
          PLOT<-ggplot(OBJ, aes(x=Position/PositionScalingFactor,y=Value/ValueScalingFactor,color=ColorValue))+
            geom_point(size=PointSize, stroke=0, alpha=PointAlpha)+
            CS.THEME+
            scale_x_continuous(limits = c(0,CHROMLENGTH/PositionScalingFactor), expand = c(0,0))+
            facet_grid(~Chromosome)+
            scale_color_viridis(begin = 0.9, end = 0.1, option = "cividis", limits=c(FillMin,FillMax))
          
        } else {
          
          PLOT<-ggplot(OBJ, aes(x=Position/PositionScalingFactor,y=Value/ValueScalingFactor,color=ColorValue))+
            geom_point(size=PointSize, stroke=0, alpha=PointAlpha)+
            CS.THEME+
            scale_x_continuous(limits = c(0,CHROMLENGTH/PositionScalingFactor), expand = c(0,0))+
            facet_grid(~Chromosome)+
            scale_color_viridis(begin = 0.9, end = 0.1, option = "cividis")
          
        }
        
      } else {
        
        PLOT<-ggplot(OBJ, aes(x=Position/PositionScalingFactor,y=Value/ValueScalingFactor))+
          geom_point(size=PointSize, stroke=0, color="grey50", alpha=PointAlpha)+
          CS.THEME+
          scale_x_continuous(limits = c(0,CHROMLENGTH/PositionScalingFactor), expand = c(0,0))+
          facet_grid(~Chromosome)
        
      }
      
      if(!is.na(Xlabel)){
        PLOT<-PLOT+
          xlab(Xlabel)
      }
      
      if(!is.na(Ylabel)){
        PLOT<-PLOT+
          ylab(Ylabel)
      }
      
      if(!is.null(HORIZONTALLINES)){
        PLOT<-PLOT+
          geom_hline(yintercept = HORIZONTALLINES/ValueScalingFactor, color="red", linetype="dashed", size=2)
        
        # Find min Value
        MIN.YVALUE<-min(c(0, HORIZONTALLINES/ValueScalingFactor,OBJ$Value/ValueScalingFactor), na.rm = TRUE)
        
        if(MIN.YVALUE<0){
          MIN.YVALUE<-MIN.YVALUE*1.1
        }
        
        # Find max Value
        MAX.YVALUE<-max(c(0, HORIZONTALLINES/ValueScalingFactor,OBJ$Value/ValueScalingFactor), na.rm = TRUE)
        
        # Adjust y plotting margin
        PLOT<-PLOT+
          scale_y_continuous(limits = c(MIN.YVALUE,MAX.YVALUE*1.1), expand = c(0,0))
        
      } else {
        # Find min Value
        MIN.YVALUE<-min(0, OBJ$Value/ValueScalingFactor, na.rm = TRUE)
        
        if(MIN.YVALUE<0){
          MIN.YVALUE<-MIN.YVALUE*1.1
        }
        
        # Find max Value
        MAX.YVALUE<-max(0, OBJ$Value/ValueScalingFactor, na.rm = TRUE)
        
        # Adjust y plotting margin
        PLOT<-PLOT+
          scale_y_continuous(limits = c(MIN.YVALUE,MAX.YVALUE*1.1), expand = c(0,0))
        
      }
      
      if(!is.null(VERTICALLINES)){
        PLOT<-PLOT+
          geom_vline(xintercept = VERTICALLINES/PositionScalingFactor)
      }
      
      if(FitGAM){
        PLOT<-PLOT+
          geom_smooth(method = "gam", color=GAMLineColor, se=FALSE, formula = y ~ s(x, bs = "cs") )
      }
      
      if(!is.null(HIGHLIGHT) & CHROMOSOME %in% HIGHLIGHT[,1]){
        colnames(HIGHLIGHT)<-c("Chr","Start","End")
        HIGHLIGHT<-subset(HIGHLIGHT,Chr==CHROMOSOME)
        
        PLOT<-PLOT+
          geom_rect(data = HIGHLIGHT, inherit.aes = FALSE, aes(xmin = Start/PositionScalingFactor, xmax = End/PositionScalingFactor, ymin = -Inf, ymax = Inf), fill = "#55575725", alpha = 0.3)
      }
      
      return(PLOT)
    } else {
      print("Warning: Chromsome is valid, but data object lacks any points to plot on it...")
      PLOT<-ggplot()
    }

  } else {
    
    print("Error: CHROMOSOME missing from CHROMSIZES")
    return(NULL)
  }
}
