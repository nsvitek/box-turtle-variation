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

#Adam Rountrey's rotateAMatrix function, from UC Berkeley's Morphometrics Short Course
rotateAMatrix = function(mat,degrees,rotCenX,rotCenY){
  radns = degrees/360*2*pi
  rotmat = matrix(c(cos(radns),sin(radns),-sin(radns),cos(radns)),2,2)
  tempmat = matrix(NA,dim(mat)[1],dim(mat)[2])
  tempmat[,1] = mat[,1] - rotCenX
  tempmat[,2] = mat[,2] - rotCenY
  tempmat = as.matrix(tempmat)%*%rotmat
  tempmat[,1] = tempmat[,1] + rotCenX
  tempmat[,2] = tempmat[,2] + rotCenY
  tempmat
}