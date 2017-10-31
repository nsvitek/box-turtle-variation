#estimating measurement error, based on same way of measuring error in 
#submitted sensitivity analysis manuscript, which is itself based on 
#Lockwood et al. 2002 and Fruciano 2016

setwd(datadir)
setwd("../data-output/error")

#remove missing shapes
if(length(which(is.na(shape)))>=1){
  data_all<-data_all[-which(is.na(shape)),]
}
shape<-na.omit(shape)

#explore to get a sense of what you expect to find
PCA<-prcomp(shape,scale.=FALSE)
an<-anderson(PCA$sdev)

# Lockwood -----------------
# get list of which id's are ones with replicates
repnames<-which(data_all$err_rep==1) %>% data_all$id[.] %>% 
  unique %>% droplevels

repgroups<-rep(1,nrow(data_all)) #base vector of replicate groups, to be modified
for (i in 1:length(repnames)){repgroups[grep(repnames[i],data_all$id)]<-i+1}

errdist<-NULL
origdist<-as.vector(dist(PCA$x[which(repgroups==1),er_pcs],method="euclidean"))
for (j in 1:length(repnames)){
  if(j+1>=length(repnames)){break}
  errdist[[j]]<-c(as.vector(dist(PCA$x[which(repgroups==(j+1)),er_pcs],method="euclidean")))
  #problem here with errdist...what's going on?
  
} #calculate pairwise euclidean distances in PC space between non-replicate specimens

#look at distribution of distances
opar<-par()
jpeg(filename=paste("error_",view,"_hist.jpg",sep=""))
par(mfrow=c(1,2))
hist(origdist) #spot check: no crazy skewed/bimodal distributions
hist(unlist(errdist)) #spot check: no crazy skewed/bimodal distributions?
dev.off()
par<-opar

summary_groups<-matrix(c(mean(origdist),mean(unlist(errdist)),
                         sd(origdist),sd(unlist(errdist)),
                         min(origdist),min(unlist(errdist)),
                         max(origdist),max(unlist(errdist))),
                       nrow=4,byrow=TRUE)

colnames(summary_groups)<-c("between_unique_specimens","between_error_replicates")
summary_groups<-summary_groups %>% as.data.frame %>% 
  mutate(.,ratio=between_error_replicates/between_unique_specimens)
row.names(summary_groups)<-c("mean","standard deviation","minimum","maximum")

#write group-wise summary statistics to csv file
write.csv(summary_groups,file=paste(view,"_error_distances_summary-stats.csv",sep=""))

# # Solid graphing option, also from Kowalewski's tutorial:
jpeg(filename=paste("error_",view,"_hist_comparison.jpg",sep=""))
hist(origdist, xlab="distances between unique specimens", ylab="number of observations", 
     main="comparison of distances: \nbetween unique specimens vs. \nbetween error replicates\narrow = mean",
     col="black", breaks=100)
hist(origdist[origdist<=max(unlist(errdist))],breaks=round(summary_groups[4,3]*100,0),add=T, col="red",border="red")
arrows(mean(unlist(errdist)),1000,mean(unlist(errdist)),50,length=0.1, col="red", lwd=2, code=2)
points(mean(unlist(errdist)),1000,cex=1, pch=16, col="red")
dev.off()

# Fruciano -----------
testgdf<-geomorph.data.frame(coords=PCA$x,specimen=data_all$id)
errorANOVA<-procD.lm(coords~factor(specimen),data=testgdf,iter=999,RRPP=TRUE) %>% .$aov.table
repeatability<-((errorANOVA$MS[1]-errorANOVA$MS[2])/2.2)/(errorANOVA$MS[2]+((errorANOVA$MS[1]-errorANOVA$MS[2])/2.2))

jpeg(filename=paste("error_",view,"_pca.jpg",sep=""))
plot(PCA$x[,1:2],bg=c("black","red")[data_all$err_rep+1],pch=21)
title(paste("Repeatability =", round(repeatability,3)))
dev.off()

repPC<-find_repeatablePCs(PCA$x,data_all$id,2.2)

jpeg(filename=paste("error_",view,"_by_pc.jpg",sep=""))
repPC<-find_repeatablePCs(PCA$x,data_all$id,2)
title(paste("# repeatable PCs =", max(which(repPC>=0.9))))
dev.off()

plot(PCA$x[,1:2],bg=c("black","red")[data_all$err_rep+1],pch=21)
text(x=PCA$x[which(data_all$err_rep==1),1],y=PCA$x[which(data_all$err_rep==1),2],labels=data_all$id[which(data_all$err_rep==1)],pos=1)
# 5 in dor, 9 in lat, 13 in pos