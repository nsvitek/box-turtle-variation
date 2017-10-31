#function to determine how many PCs are nonzero. 
#input: 
#     "an", or results from the home-baked anderson function. 
nPC<-function(an){
  result<-which(an$percent==0.00) %>% min-1
  return(result)
}

### draw the error ellipse around landmark constellations. Updated version of 'errbubble' in sensitivity_utils.R
#Input: x & y vals for ellipse center, the values returned from errbvals()
#   can also set rgb() vals for col and border
#Output: a drawing layer of error ellipses
plot.errell<-function(shape,errcalc,color.border="gray",color.bg="gray",alpha.border=0.1,alpha.bg=0.5){
  rgbbor<-rgb(t(col2rgb(color.border)/255),alpha=alpha.border)
  rgbcol<-rgb(t(col2rgb(color.bg)/255),alpha=alpha.bg) #change color of residuals
    for (i in 1:nrow(errcalc)){
    draw.ellipse(shape[i,1],shape[i,2],a=errcalc[i,1],
                 b=errcalc[i,2],col=rgbcol,border=rgbbor)	
  }
}


#function to subset metadata from the complete dataset
#input: 
#      all.data: 2D data frame of all data, both shape coordinates and metadata
#      shape.data: 2D data frame of only shape coordinates
#      filter: the which() statement of which specimens to include
#      metadata.columns: vector of which columns are metadata, not shape coordinates
#output: a 2D data frame of the metadata for the specimens in the dataset of interest. 
makemetadata<-function(alldata,shapedata,filter,metadata_columns){
  if(length(which(is.na(shapedata[filter,1])))>=1){
    metadata<-alldata[filter,metadata_columns] %>% 
      .[-which(is.na(shapedata[filter,1])),]
  } else {
    metadata<-alldata[filter,metadata_columns]
  }
  metadata<-droplevels(metadata)
}

### make the mean shape of a group of shapes
#input: matrix where rows = specimens, cols= points; # dimensions in dataset (2,3)
#output: n,m,k,meanshape
#default m=2 (2D)
mshp<-function(x,m=2){
  n<-nrow(x)
  k<-ncol(x)/m
  meanshape<-t(matrix(colMeans(x),m,k))
  result<-list(n,m,k,meanshape)
  names(result)<-c("n","m","k","meanshape")
  return(result)
}

# #function to plot shape, with or without standard deviation
# # shapenames = vector of shapes. Or list of shapes?
# function(shapes,colors,plot.bubbles=TRUE){
#   xlims<-range(shapes) %>% +c(-0.015,0.015)
#   ylims<-range(shapes) %>% +c(-0.015,0.015)
#   plot(shapes[[1]],asp=1,pch=21,bg=colors[1],axes=FALSE,xlab="",ylab="",
#      xlim=xlims, ylim=ylims,cex=0.5)
#   for (i in 1:length(shapes)){ ####NEEDS TO CREATE OBJECTS THAT WON'T BE OVERRITTEN
#     rgbcol<-rgb(t(col2rgb(colors[1])/255),alpha=.5)
#     rgbbor<-rgb(t(col2rgb(colors[1])/255),alpha=.1)
#     eb_car<-errbubble(x=carolina_mshp$meanshape[,1],y=carolina_mshp$meanshape[,2],sd_carolina)
#   }
# rgbcol<-rgb(t(col2rgb(colors[2])/255),alpha=.5)
# rgbbor<-rgb(t(col2rgb(colors[2])/255),alpha=.1)
# eb_coa<-errbubble(x=coahuila_mshp$meanshape[,1],y=coahuila_mshp$meanshape[,2],sd_coahuila)
# lines(carolina_mshp$meanshape[conx,],col=spp_col[1],lwd=2)
# lines(coahuila_mshp$meanshape[conx,],col=spp_col[2],lwd=2)
# points(carolina_mshp$meanshape,asp=1,pch=21,bg=spp_col[1],cex=0.5)
# points(coahuila_mshp$meanshape,asp=1,pch=21,bg=spp_col[2],cex=0.5)
# }