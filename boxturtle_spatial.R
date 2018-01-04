# developing spatial eigenvector analyses (SEVM) or Moran Eigenvector (MEM) filters
# to apply to multivariate models of morphometric data.
library(sp) #making spatial object, coords()
library(spdep) #neighbors, all the x2y functions
library(geosphere) #distm

# Principal Components Analysis
PCAssp<-prcomp(ssp_set,scale.=FALSE) #explore
an_ssp<-anderson(PCAssp$sdev)

# # make data spatial
ssp.space<-cbind(ssp_metadata,PCAssp$x[,c(1:PCc)])
coordinates(ssp.space)<-c("longitude","latitude")
# class(ssp_metadata)
# First, find a way to calculate Moran's I and view correlograms ------

# #calculates global Moran's I
near.six<-knearneigh(ssp.space,k=16,RANN=F)
neighbors<-knn2nb(near.six)
plot(neighbors,coordinates(ssp.space))
spaceweights<-nb2listw(neighbors)
moran.mc(ssp_metadata$lat_cs,spaceweights,nsim=99)
moran.mc(PCAssp$x[,1],spaceweights,nsim=99)
# #make correlogram
moran.corr<-sp.correlogram(neighbors,ssp_metadata$lat_cs,order=10,method="I",zero.policy=TRUE)
moran.corr<-sp.correlogram(neighbors,PCAssp$x[,1],order=10,method="I",zero.policy=TRUE)

plot(moran.corr)
#results are qualitatively similar to SAM with 10 default distance classes but not quantitatively
#next neighbor test
# neighbors2<-dnearneigh(cbind(ssp_metadata$longitude,ssp_metadata$latitude),d1=0,d2=282,longlat=TRUE)

# ?nblag #might be worth trying: https://cran.r-project.org/web/packages/spdep/vignettes/nb.pdf

# Next, find a way to create spatial eigenvectors -------
# compute distance matrix
locality_dist<-distm(ssp.space@coords,fun=distVincentyEllipsoid) #output in meters
locality_dist2<-locality_dist
locality_dist2[which(locality_dist2>=283*1000)]<-0 #approximation of threshold from sam

# PCNM axes of the dist. matrix (from 'vegan' package)
sevm.test1<-pcnm(locality_dist)
pcnm.axes <- sevm.test1$vectors
(sevm.test1$values/sum(sevm.test1$values)*100)[1:10] #% variance explained? does that work with negative eigenvalues?
#note that eigenvalues in sevm.test1 don't match SAM-based eigenvalues

str(sevm.test1)
# Can you re-create results from SAM? -------

sam.text<-read.table("../test.txt",sep="\t",header = TRUE)
str(sam.text)

plot(ssp.space@coords,pch=21,bg=terrain.colors(100)[cut(sevm.test1$vectors[,89],100)]) #plot an axis
plot(ssp.space@coords,pch=21,bg=terrain.colors(100)[cut(sam.text[,1],100)]) #plot an axis
#looks like pcnm provides same answers as SAM, except sometimes axes flipped.
#note that the threshold for sevm.test1 is same as truncation distance in SAM.



plot(ssp.space@coords,pch=21,bg=terrain.colors(100)[cut(sevm.test1$vectors[,89],12)]) #plot an axis
plot(ssp.space@coords,pch=21,bg=filterr[cut(sevm.test1$vectors[,1],12)],alpha=0.3,cex=2) #plot an axis
gg.ssp<-cbind(ssp.space@coords,sevm.test1$vectors) %>% as.data.frame
colnames(gg.ssp)
ggplot(data=gg.ssp, aes(x=longitude,y=latitude)) +
  geom_point(aes(fill=PCNM1),size=4,alpha=0.3,shape=21)+
  scale_fill_gradient(low="white",high="blue")+
  theme_minimal()

colorgunblue<-colorRampPalette(c(rgb(1,1,1,0.01),rgb(0,0,1,0.5)),alpha=TRUE)
colorgunred<-colorRampPalette(c(rgb(1,1,1,0.01),rgb(1,0,0,0.5)),alpha=TRUE)
colorgungreen<-colorRampPalette(c(rgb(1,1,1,0.01),rgb(0,1,0,0.5)),alpha=TRUE)
plot(ssp.space@coords,pch=21,bg=colorgungreen(12)[cut(sevm.test1$vectors[,1],12)],cex=2)
points(ssp.space@coords,pch=21,bg=colorgunred(12)[cut(sevm.test1$vectors[,2],12)],cex=2,add=TRUE)
points(ssp.space@coords,pch=21,bg=colorgunblue(12)[cut(sevm.test1$vectors[,3],12)],cex=2,add=TRUE)
#3 is too many. 

plot(ssp.space@coords,pch=21,bg=colorgunblue(12)[cut(sevm.test1$vectors[,3],12)],cex=2)
points(ssp.space@coords,pch=21,bg=colorgunred(12)[cut(sevm.test1$vectors[,2],12)],cex=2,add=TRUE)
#kind of cool: spatial filters 2 and 3 might account for spotty significance of triunguis and bauri
#note that major/gulf coast may also be mixed in.
gg.ssp$PCNM1
# Pick filters ------
library(packfor)
# install.packages("packfor", repos="http://R-Forge.R-project.org")

forward.sel(PCAssp$x[,c(1:PCc)], sevm.test1$vectors, K = nrow(sevm.test1$vectors) - 1, R2thresh = 0.99, adjR2thresh = 0.99,nperm = 999, 
            R2more = 0.001, alpha = 0.05, Xscale = TRUE, Ycenter = TRUE, Yscale
            = FALSE)

#try selecting only 3 variables?
fs<-forward.sel(PCAssp$x[,c(1:PCc)], sevm.test1$vectors, 3, R2thresh = 0.99, adjR2thresh = 0.99,nperm = 999, 
            R2more = 0.001, alpha = 0.05, Xscale = TRUE, Ycenter = TRUE, Yscale
            = FALSE) #not very interesting.
str(fs)
fs$variables

# Apply to multivariate models ------

#try again according to reviewer + Collyer recommendations
ssp.gpa<-arrayspecs(ssp_set,p=ncol(ssp_set)/2,k=2) %>% gpagen
ssp2.gdf<-geomorph.data.frame(ssp.gpa,ssp=ssp_metadata$ssp,
                              PCNM1=sevm.test1$vectors[,1],PCNM2=sevm.test1$vectors[,2],
                              PCNM3=sevm.test1$vectors[,3],PCNM89=sevm.test1$vectors[,89],
                              PCNM10=sevm.test1$vectors[,10],PCNM13=sevm.test1$vectors[,13],
                              PCNM23=sevm.test1$vectors[,23],PCNM33=sevm.test1$vectors[,33],
                              PCNM6=sevm.test1$vectors[,6],PCNM4=sevm.test1$vectors[,4])
ssp2.gdf$Csize<-ssp_metadata[,cs_metadata_col]

procD.lm(coords~PCNM1+PCNM2+PCNM3+PCNM89+PCNM10+
           PCNM13+PCNM23+PCNM4+log(Csize)+ssp,iter=499,data=ssp2.gdf)

advanced.procD.lm(coords~(PCNM1+PCNM2+PCNM3+PCNM89+PCNM10+PCNM13+PCNM23+PCNM4+log(Csize)),
                  ~(PCNM1+PCNM2+PCNM3+PCNM89+PCNM10+PCNM13+PCNM23+PCNM4+log(Csize))+ssp,iter=499,data=ssp2.gdf) 

# advanced.procD.lm(coords~(PCNM1+PCNM2+PCNM3+PCNM89+PCNM10+PCNM13+PCNM23+PCNM4+log(Csize)),
#                   ~(PCNM1+PCNM2+PCNM3+PCNM89+PCNM10+PCNM13+PCNM23+PCNM4+log(Csize))+ssp,iter=499,
#                   groups=~ssp,slope=~(PCNM1+PCNM2+PCNM3+PCNM89+PCNM10+PCNM13+PCNM23+PCNM4+log(Csize)),data=ssp2.gdf) 
# #does not work. Can't make a slope for more than 1 variable plus log(cs)

advanced.procD.lm(coords~(PCNM1+log(Csize)),  ~(PCNM1+log(Csize))+ssp, iter=499,
                  groups=~ssp,slope=~(PCNM1+log(Csize)),data=ssp2.gdf) #some ssp still significant

advanced.procD.lm(coords~(PCNM2+log(Csize)),  ~(PCNM2+log(Csize))+ssp, iter=499,
                  groups=~ssp,slope=~(PCNM2+log(Csize)),data=ssp2.gdf) #ssp not significant

advanced.procD.lm(coords~(PCNM3+log(Csize)),  ~(PCNM3+log(Csize))+ssp, iter=499,
                  groups=~ssp,slope=~(PCNM3+log(Csize)),data=ssp2.gdf) #ssp not significant



#original
advanced.procD.lm(coords~log(Csize),~log(Csize)+ssp,
                  groups=~ssp,slope=~log(Csize),iter=499,data=ssp.gdf) 
#result: dorsal: differences between triunguis-major and triunguis-carolina group
#       lateral: differences between bauri-carolina, bauri-triunguis, maybe major-triunguis?
#     posterior: nothing significant


advanced.procD.lm(coords~log(Csize)+ssp,~log(Csize)*ssp,
                  groups=~ssp,slope=~log(Csize),angle.type="deg",iter=499,data=ssp.gdf) 
#dorsal: only major/triunguis significant
#lateral: same as model above, plus slopes are different (means can't use a single allometric model?)
#posterior: no reason to check this one