# SSP: PCA ------
# Principal Components Analysis (look at data)
PCAssp<-prcomp(ssp_set,scale.=FALSE) #explore
an_ssp<-anderson(PCAssp$sdev)

# commented out because doesn't need revision. Code used for figures.
# cairo_pdf(paste("ssp_PCA_",view,".pdf",sep=""),width = 3.6, height = 3.6,family ="Arial")
# par(mar=c(3.9,3.6,.5,.5))
# plot(PCAssp$x[,1:2],cex=1,
#      bg=ssp_col[ssp_metadata$ssp],pch=(20+as.numeric(ssp_metadata$ssp)),
#      xlab="",ylab="")
# mtext(paste("PC 1 (",round(an_ssp$percent[1],3),"%)",sep=""),side=1,line=2,cex=1.2)
# mtext(paste("PC 2 (",round(an_ssp$percent[2],3),"%)",sep=""),side=2,line=2,cex=1.2)
# # legend('bottomleft',cex=1,pch=c(21:24),pt.bg=ssp_col,legend=levels(ssp_metadata$ssp))
# dev.off()
# embed_fonts(paste("ssp_PCA_",view,".pdf",sep=""))

# commented out because doesn't need revision. Code used for figures.
# cairo_pdf("ssp_PCA_legend.pdf",width = 3.6, height = 3.6,family ="Arial")
# par(mar=c(.5,.5,.5,.5))
# plot(PCAssp$x[,1:2],type="n",axes=F,xlab="",ylab="")
# legend('center',cex=2,pch=c(21:24),pt.bg=ssp_col,legend=levels(ssp_metadata$ssp))
# dev.off()
# embed_fonts("ssp_PCA_legend.pdf")

# SSP: k-means -------
classifier<-as.character(ssp_metadata$ssp)
cluster.data<-PCAssp$x[,1:PCc]
mc<-4
BIC<-mclustBIC(cluster.data,G=1:mc)
plot(BIC)
one.cluster<-max(BIC[1,])
more.clusters<-max(BIC[c(2:4),],na.rm=TRUE) #BIC score of best model with more than one cluster
how.much<-one.cluster-more.clusters
best.k<-((which(BIC==max(BIC,na.rm=TRUE)) /mc)-floor((which(BIC==max(BIC,na.rm=TRUE)) /mc)))*mc
kmeans.report<-rbind(one.cluster,more.clusters,how.much,best.k)
row.names(kmeans.report)<-c("BIC.best.one.cluster",
                            "BIC.best.multiple.clusters",
                            "deltaBIC","best.k")
write.csv(kmeans.report,paste("kmeans_results_",view,".csv",sep=""))
# summary(BIC)
mod1<-Mclust(cluster.data,x=BIC) #only one cluster, no need to model
summary(mod1)
table(classifier, mod1$classification) #not super helpful with one cluster
# mod1dr<-MclustDR(mod1,normalized=FALSE) #futher code is example of what to do
# summary(mod1dr) #if you have multiple clusters. Not helpful now.
# plot(mod1dr,what="pairs")
# plot(mod1dr, what = "boundaries", ngrid = 200)


# SSP: spatial 1: Morans I ------
# developing spatial eigenvector analyses (SEVM) or Moran Eigenvector (MEM) filters
# to apply to multivariate models of morphometric data.
# # make data spatial
ssp.space<-cbind(ssp_metadata,PCAssp$x[,c(1:PCc)])
coordinates(ssp.space)<-c("longitude","latitude")

# # #calculates global Moran's I
# near.sixteen<-knearneigh(ssp.space,k=16,RANN=F)
# neighbors<-knn2nb(near.sixteen)
# plot(neighbors,coordinates(ssp.space))
# spaceweights<-nb2listw(neighbors)
# moran.mc(ssp_metadata$lat_cs,spaceweights,nsim=99)
# moran.mc(PCAssp$x[,1],spaceweights,nsim=99)
# # #make correlogram
# moran.corr1<-sp.correlogram(neighbors,ssp_metadata$lat_cs,order=10,method="I",zero.policy=TRUE)
# moran.corr2<-sp.correlogram(neighbors,PCAssp$x[,1],order=10,method="I",zero.policy=TRUE)
# plot(moran.corr2)
# #results are qualitatively similar to SAM with 10 default distance classes but not quantitatively

# SSP: spatial 2: SEVM -------
# compute distance matrix
locality_dist<-distm(ssp.space@coords,fun=distVincentyEllipsoid) #output in meters

# PCNM axes of the distance matrix
sevm<-pcnm(locality_dist)
#note that eigenvalues in sevm.test1 don't match SAM-based eigenvalues
#also, see Dray's note about negative eigenvalues. 

# # Can you re-create results from SAM? -------
# sam.text<-read.table("../test.txt",sep="\t",header = TRUE)
# plot(ssp.space@coords,pch=21,bg=terrain.colors(100)[cut(sevm.test1$vectors[,89],100)]) #plot an axis
# plot(ssp.space@coords,pch=21,bg=terrain.colors(100)[cut(sam.text[,1],100)]) #plot an axis
# #looks like pcnm provides same answers as SAM, except sometimes axes flipped.
# #note that the threshold for sevm.test1 is same as truncation distance in SAM.
# Pick filters ------
# forward.sel(PCAssp$x[,c(1:PCc)], sevm.test1$vectors, K = nrow(sevm.test1$vectors) - 1, R2thresh = 0.99, adjR2thresh = 0.99,nperm = 999, 
#             R2more = 0.001, alpha = 0.05, Xscale = TRUE, Ycenter = TRUE, Yscale
#             = FALSE)

#try selecting only the top 3 variables. Keep the model simple
fs<-forward.sel(PCAssp$x[,c(1:PCc)], sevm$vectors, 3, R2thresh = 0.99, adjR2thresh = 0.99,nperm = 999,
                R2more = 0.001, alpha = 0.05, Xscale = TRUE, Ycenter = TRUE, Yscale
                = FALSE)
fs$variables

# # Visualize Filter Maps-------
# usa<-map_data('state')
# xlims<-range(ssp.space@bbox[1,])+c(-0.8,+0.8)
# ylims<-range(ssp.space@bbox[2,])+c(-0.8,+0.8)
# #
# # choices<-which(colnames(sevm$vectors) %in% fs$variables)
# ssp.space.plot<-cbind(ssp_metadata$longitude,ssp_metadata$latitude,sevm$vectors[,c(1,2,3)]) %>% as.data.frame
# colnames(ssp.space.plot)<-c("longitude","latitude","SE1","SE2","SE3")
# #can't figure out how to completely automate this plot in every view. Sorry.
# cairo_pdf("eigenmap1.pdf",width = 5.4, height = 4.4,family ="Arial",bg="transparent")
# ggplot(data=ssp.space.plot,aes(x=longitude,y=latitude))+
#   geom_polygon(data=usa,aes(x=long,y=lat,group=group),fill=NA,color="gray") +
#   geom_point(aes(fill=SE1),shape=21,size=2)+
#   scale_fill_gradientn(colours=map_colors) +
#   coord_cartesian(xlim=xlims,ylim=ylims)+
#   theme_classic() +
#   theme(axis.line=element_blank(),axis.text.x=element_blank(),
#         axis.text.y=element_blank(),axis.ticks=element_blank(),
#         axis.title.x=element_blank(),axis.title.y=element_blank(),
#         legend.position=c(0.9,0.3),text=element_text(size=10),
#         panel.background = element_rect(fill = "transparent",colour = NA), # or theme_blank()
#         panel.grid.minor = element_blank(),
#         panel.grid.major = element_blank(),
#         plot.background = element_rect(fill = "transparent",colour = NA))
# dev.off()
# embed_fonts("eigenmap1.pdf")
# 
# cairo_pdf("eigenmap2.pdf",width = 5.4, height = 4.4,family ="Arial",bg="transparent")
# ggplot(data=ssp.space.plot,aes(x=longitude,y=latitude))+
#   geom_polygon(data=usa,aes(x=long,y=lat,group=group),fill=NA,color="gray") +
#   geom_point(aes(fill=SE2),shape=21,size=2)+
#   scale_fill_gradientn(colours=map_colors) +
#   coord_cartesian(xlim=xlims,ylim=ylims)+
#   theme_classic() +
#   theme(axis.line=element_blank(),axis.text.x=element_blank(),
#         axis.text.y=element_blank(),axis.ticks=element_blank(),
#         axis.title.x=element_blank(),axis.title.y=element_blank(),
#         legend.position=c(0.9,0.3),text=element_text(size=10),
#         panel.background = element_rect(fill = "transparent",colour = NA), # or theme_blank()
#         panel.grid.minor = element_blank(),
#         panel.grid.major = element_blank(),
#         plot.background = element_rect(fill = "transparent",colour = NA))
# dev.off()
# embed_fonts("eigenmap2.pdf")
# 
# cairo_pdf("eigenmap3.pdf",width = 5.4, height = 4.4,family ="Arial",bg="transparent")
# ggplot(data=ssp.space.plot,aes(x=longitude,y=latitude))+
#   geom_polygon(data=usa,aes(x=long,y=lat,group=group),fill=NA,color="gray") +
#   geom_point(aes(fill=SE3),shape=21,size=2)+
#   scale_fill_gradientn(colours=map_colors) +
#   coord_cartesian(xlim=xlims,ylim=ylims)+
#   theme_classic() +
#   theme(axis.line=element_blank(),axis.text.x=element_blank(),
#         axis.text.y=element_blank(),axis.ticks=element_blank(),
#         axis.title.x=element_blank(),axis.title.y=element_blank(),
#         legend.position=c(0.9,0.3),text=element_text(size=10),
#         panel.background = element_rect(fill = "transparent",colour = NA), # or theme_blank()
#         panel.grid.minor = element_blank(),
#         panel.grid.major = element_blank(),
#         plot.background = element_rect(fill = "transparent",colour = NA))
# dev.off()
# embed_fonts("eigenmap3.pdf")

# SSP: stats4shape -----
#try again according to reviewer + Collyer recommendations

#first, make objects for Procrustes ANOVA
ssp.gpa<-arrayspecs(ssp_set,p=ncol(ssp_set)/2,k=2) %>% gpagen
ssp.gdf<-geomorph.data.frame(ssp.gpa,ssp=ssp_metadata$ssp,
                              PCNM1=sevm$vectors[,1],PCNM2=sevm$vectors[,2],
                              PCNM3=sevm$vectors[,3],PCNM5=sevm$vectors[,5])
ssp.gdf$Csize<-ssp_metadata[,cs_metadata_col]

#Make an aspatial model.
#Step 1 of Dr. Collyer's Little Guide
advanced.procD.lm(coords~log(Csize)+ssp,~log(Csize)*ssp,iter=999,data=ssp.gdf)

#Significant. Step 1b
advanced.procD.lm(coords~log(Csize)+ssp,~log(Csize)*ssp,
                  groups=~ssp,slope=~log(Csize),iter=999,data=ssp.gdf)
#signficant in lateral view in distance

#Not significant. Step 2. "Focus on the LS means If ANOVA does not return a significant result, groups are not different."
P.man.nospace<-advanced.procD.lm(coords~log(Csize),~log(Csize)+ssp,
                  groups=~ssp,slope=~log(Csize),iter=999,data=ssp.gdf)
P.man.nospace$P.slopes.dist
#results:
#dorsal: differences between triunguis-major and triunguis-carolina group
#lateral: differences between bauri-carolina, bauri-triunguis, maybe major-triunguis?
#posterior: nothing significant

#Incorporate spatial autocorrelation into the model
# top the EVs: posterior view: PCNM 1,2; lateral view: PCNM 1,2; dorsal view: PCNM 2,3
# go with 1 and 2, the two consistently important EVs? Or just 2? First separately.
# (sevm$values/sum(sevm$values)*100) #how to interpret with negative eigenvalues?

if(view=="lat"|view=="pos"){
  P.man.pair.space<-advanced.procD.lm(coords~(PCNM2+log(Csize)),  ~(PCNM2+log(Csize))+ssp, iter=999,
                    groups=~ssp,slope=~(PCNM2+log(Csize)),data=ssp.gdf)
  P.man.space<-procD.lm(coords~PCNM1+PCNM2+log(Csize), iter=999,data=ssp.gdf)

  #lateral: nothing
  #posterior: nothing, but there was nothing before.
} else if (view=="dor"){
  P.man.pair.space<-advanced.procD.lm(coords~(PCNM3+log(Csize)),  ~(PCNM3+log(Csize))+ssp, iter=999,
                                 groups=~ssp,slope=~(PCNM3+log(Csize)),data=ssp.gdf)
  P.man.space<-procD.lm(coords~PCNM2+PCNM3+log(Csize), iter=999,data=ssp.gdf)
  #dorsal: nothing.
}

P.man.pair.space$P.slopes.dist

P.man.results<-P.man.pair.space$P.slopes.dist
for (i in 1:ncol(P.man.results)){
  P.man.results[,i]<-p.adjust(P.man.results[,i],method="BH")
}
P.man.pair.nospace<-P.man.nospace$P.slopes.dist
for (i in 1:ncol(P.man.pair.nospace)){
  P.man.pair.nospace[,i]<-p.adjust(P.man.pair.nospace[,i],method="BH")
}

P.man.results[lower.tri(P.man.results)]<-P.man.pair.nospace[lower.tri(P.man.pair.nospace)]

write.csv(P.man.results,paste("Pman_pairwise_results_",view,".csv",sep=""))
write.csv(P.man.space$aov.table,paste("Pman_continuous_results_",view,".csv",sep=""))
# lower triangle is aspatial, upper triangle is spatial
# # try both
# advanced.procD.lm(coords~PCNM2+PCNM1+log(Csize),  ~PCNM2+PCNM1+log(Csize)+ssp, iter=999,
#                   groups=~ssp, slope=~(PCNM2+PCNM1+log(Csize)),data=ssp.gdf)
# #can't add in all 3 covariates into slope. Model still valid without all 3?

# # SSP: stats4size --------
# size.c<-ssp_metadata$carapace_length %>% log
# size.d<-ssp_metadata$dor_cs %>% log
# size.l<-ssp_metadata$lat_cs %>% log
# size.p<-ssp_metadata$pos_cs %>% log
# hist(size.c) #looks normal-ish
# shapiro.test(size.p)
# #c & p fail, d & l pass
#
# #first, look at spatial distribution of size
# usa<-map_data('state')
# xlims<-range(ssp.space@bbox[1,])+c(-0.8,+0.8)
# ylims<-range(ssp.space@bbox[2,])+c(-0.8,+0.8)
# ssp.size.plot<-cbind(ssp_metadata$longitude,ssp_metadata$latitude,size.d,size.l,size.p,size.c) %>% as.data.frame
# colnames(ssp.size.plot)<-c("longitude","latitude","dor","lat","pos","length")
# ggplot(data=ssp.size.plot,aes(x=longitude,y=latitude))+
#   geom_polygon(data=usa,aes(x=long,y=lat,group=group),fill=NA,color="gray") +
#   geom_point(aes(fill=length),shape=21,size=2)+
#   scale_fill_gradientn(colours=map_colors) +
#   coord_cartesian(xlim=xlims,ylim=ylims)+
#   theme_classic() +
#   theme(axis.line=element_blank(),axis.text.x=element_blank(),
#         axis.text.y=element_blank(),axis.ticks=element_blank(),
#         axis.title.x=element_blank(),axis.title.y=element_blank(),
#         legend.position=c(0.9,0.3),text=element_text(size=10),
#         panel.background = element_rect(fill = "transparent",colour = NA), # or theme_blank()
#         panel.grid.minor = element_blank(),
#         panel.grid.major = element_blank(),
#         plot.background = element_rect(fill = "transparent",colour = NA))
#
# summary(lm(ssp_metadata$pos_cs~sevm$vectors[,1]+sevm$vectors[,2]+sevm$vectors[,4]+ssp_metadata$ssp))
# summary(lm(ssp_metadata$lat_cs~sevm$vectors[,1]+sevm$vectors[,2]+sevm$vectors[,4]+ssp_metadata$ssp))
# summary(lm(ssp_metadata$dor_cs~sevm$vectors[,2]+sevm$vectors[,4]+ssp_metadata$ssp))
# summary(lm(ssp_metadata$carapace_length~sevm$vectors[,2]+sevm$vectors[,4]+ssp_metadata$ssp))

# SSP: Assignments -----
# # # Use CVAGen8 for CVA-based assignments tests with jackknife
# # # write tables for the two groups
# write.table(file=paste("cva_binary_",view,".txt",sep=""),cbind(ssp_all,ssp_all_metadata[,cs_metadata_col]),col.names=FALSE,row.names=FALSE,sep='\t')
# write.table(file=paste("cva_binary_grp_",view,".txt",sep=""),as.numeric(binary_all),col.names=FALSE,row.names=FALSE)
# write group ID's to compare to traditional ssp
# write.table(file=paste("cva_ssp_",view,"_grp.txt",sep=""),as.numeric(ssp_all_metadata$ssp),col.names=FALSE,row.names=FALSE)

# SSP: Plot Filter Shapes -----
# build model for EV2
model<-lm(as.matrix(ssp_set[,])~sevm$vectors[,2])
modm<-matrix(model$coefficients[2,],m,k,byrow=TRUE) #m in the linear model y=mx+b
# residuals<-model$residuals
modb<-matrix(model$coefficients[1,],m,k,byrow=TRUE) #b in the linear model y=mx+b

#when you remove predicted shape, what you are left are the residuals, yes?
testev2<-function(x){modm*x+modb} #make modern linear model
# testa(-.10) %>% plot(.,asp=1) #check

# figure out maximum and minimum filter value to plot predicted shapes
which.min<-which(sevm$vectors[,2]==min(sevm$vectors[,2]))
which.max<-which(sevm$vectors[,2]==max(sevm$vectors[,2]))
min_r<-testev2(sevm$vectors[which.min,2]) %>% rotateAMatrix(.,r1,0,0)
max_r<-testev2(sevm$vectors[which.max,2]) %>% rotateAMatrix(.,r1,0,0)

pdf(paste("mod_wireframe_EV2min_",view,".pdf",sep=""),width = 4.4, height = 4.4,useDingbats=FALSE)
par(mar=c(.01,.01,.01,.01))
plotRefToTarget(max_r,min_r,method="TPS",gridPars=gridPar(tar.pt.size=.5,grid.lwd=0.8))
points(min_r,pch=16)
dev.off()

pdf(paste("mod_wireframe_EV2max_",view,".pdf",sep=""),width = 4.4, height = 4.4,useDingbats=FALSE)
par(mar=c(.01,.01,.01,.01))
plotRefToTarget(min_r,max_r,method="TPS",gridPars=gridPar(tar.pt.size=.5,grid.lwd=0.8))
points(max_r,pch=15)
dev.off()

#repeat for EV3
model<-lm(as.matrix(ssp_set[,])~sevm$vectors[,3])
modm<-matrix(model$coefficients[2,],m,k,byrow=TRUE) #m in the linear model y=mx+b
# residuals<-model$residuals
modb<-matrix(model$coefficients[1,],m,k,byrow=TRUE) #b in the linear model y=mx+b

#when you remove predicted shape, what you are left are the residuals, yes?
testev3<-function(x){modm*x+modb} #make modern linear model
# testa(-.10) %>% plot(.,asp=1) #check

# figure out maximum and minimum filter value to plot predicted shapes
which.min<-which(sevm$vectors[,2]==min(sevm$vectors[,2]))
which.max<-which(sevm$vectors[,2]==max(sevm$vectors[,2]))
min_r<-testev3(sevm$vectors[which.min,2]) %>% rotateAMatrix(.,r1,0,0)
max_r<-testev3(sevm$vectors[which.max,2]) %>% rotateAMatrix(.,r1,0,0)

pdf(paste("mod_wireframe_EV3min_",view,".pdf",sep=""),width = 4.4, height = 4.4,useDingbats=FALSE)
par(mar=c(.01,.01,.01,.01))
plotRefToTarget(max_r,min_r,method="TPS",gridPars=gridPar(tar.pt.size=.5,grid.lwd=0.8))
points(min_r,pch=16)
dev.off()

pdf(paste("mod_wireframe_EV3max_",view,".pdf",sep=""),width = 4.4, height = 4.4,useDingbats=FALSE)
par(mar=c(.01,.01,.01,.01))
plotRefToTarget(min_r,max_r,method="TPS",gridPars=gridPar(tar.pt.size=.5,grid.lwd=0.8))
points(max_r,pch=15)
dev.off()
# # SSP: bgPCA --------
# # bin_residuals_mean<-colMeans(ssp_set)#compute grand mean
# # bin_residuals_centered<-as.matrix(ssp_set)-rep(1,nrow(ssp_set))%*%t(bin_residuals_mean) #Subtract that grandmean from each row
# #
# # #Calculate the group means
# # bin_means<-array(NA,dim=c(length(levels(binary)),ncol(ssp_set)))
# # for (i in 1:length(levels(binary))){
# #   bin_means[i,]<-colMeans(ssp_set[which(binary==levels(binary)[i]),])
# # }
# # B<-prcomp(bin_means)$rotation #eigenvectors
# # B2<-prcomp(bin_means)
# #
# # bgPCAbin<-bin_residuals_centered%*%B #Get the scores for all the individuals on the eigenvectors of the PCA of the means
# #
# # plot(bgPCAbin[,c(1,2)],cex=1.5,
# #      bg=ssp_col[binary],pch=21,xlab="",ylab="")
#
# # SAMdata --------
# #code was used to input data into SAM4.0 to check
# #whether spatial eigenvector modelling gave same results
# #in both SAM and vegan
# write1<-lapply(ssp_metadata,as.numeric) %>% as.data.frame
# write2<-cbind(write1,PCAssp$x[,c(1:PCc)])
# write2$major_growth_rings[which(is.na(write2$major_growth_rings))]<-9999
# colnames(write2)[which(colnames(write2)=="latitude")]<-"Latitude"
# colnames(write2)[which(colnames(write2)=="longitude")]<-"Longitude"
# write.csv(write2,paste("boxturtle_SAMdata_ssp_",view,".csv",sep=""),row.names=FALSE)
# getwd()