
plotLayout<-function(LAYOUT.OBJ){
  DF<-LAYOUT.OBJ
  # convert to data.frame and gather
  DF <- as.data.frame(DF)
  DF$row <- rownames(DF)
  DF <- DF %>% gather(key='column', value='value', -row)
  
  DF$row<-as.numeric(DF$row)
  
  DF$column<-factor(DF$column,unique(DF$column),ordered = TRUE)
  
  ggplot(DF,aes(x=column, y=row, fill=value)) + 
    geom_raster() +
    scale_fill_viridis_c()+
    theme(axis.text.y = element_blank())
}
