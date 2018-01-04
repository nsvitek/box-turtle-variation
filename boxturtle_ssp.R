# SSP: PCA ------
# Principal Components Analysis (look at data)
PCAssp<-prcomp(ssp_set,scale.=FALSE) #explore
an_ssp<-anderson(PCAssp$sdev)

cairo_pdf(paste("ssp_PCA_",view,".pdf",sep=""),width = 3.6, height = 3.6,family ="Arial")
par(mar=c(3.9,3.6,.5,.5))
plot(PCAssp$x[,1:2],cex=1,
     bg=ssp_col[ssp_metadata$ssp],pch=(20+as.numeric(ssp_metadata$ssp)),
     xlab="",ylab="")
mtext(paste("PC 1 (",round(an_ssp$percent[1],3),"%)",sep=""),side=1,line=2,cex=1.2)
mtext(paste("PC 2 (",round(an_ssp$percent[2],3),"%)",sep=""),side=2,line=2,cex=1.2)
# legend('bottomleft',cex=1,pch=c(21:24),pt.bg=ssp_col,legend=levels(ssp_metadata$ssp))
dev.off()
embed_fonts(paste("ssp_PCA_",view,".pdf",sep=""))

# cairo_pdf("ssp_PCA_legend.pdf",width = 3.6, height = 3.6,family ="Arial")
# par(mar=c(.5,.5,.5,.5))
# plot(PCAssp$x[,1:2],type="n",axes=F,xlab="",ylab="")
# legend('center',cex=2,pch=c(21:24),pt.bg=ssp_col,legend=levels(ssp_metadata$ssp))
# dev.off()
# embed_fonts("ssp_PCA_legend.pdf")

# SSP: spatial 1: Morans I ------
# developing spatial eigenvector analyses (SEVM) or Moran Eigenvector (MEM) filters
# to apply to multivariate models of morphometric data.
# # make data spatial
ssp.space<-cbind(ssp_metadata,PCAssp$x[,c(1:PCc)])
coordinates(ssp.space)<-c("longitude","latitude")

# #calculates global Moran's I
near.sixteen<-knearneigh(ssp.space,k=16,RANN=F)
neighbors<-knn2nb(near.sixteen)
plot(neighbors,coordinates(ssp.space))
spaceweights<-nb2listw(neighbors)
moran.mc(ssp_metadata$lat_cs,spaceweights,nsim=99)
moran.mc(PCAssp$x[,1],spaceweights,nsim=99)
# #make correlogram
moran.corr1<-sp.correlogram(neighbors,ssp_metadata$lat_cs,order=10,method="I",zero.policy=TRUE)
moran.corr2<-sp.correlogram(neighbors,PCAssp$x[,1],order=10,method="I",zero.policy=TRUE)
plot(moran.corr2)
#results are qualitatively similar to SAM with 10 default distance classes but not quantitatively

# SSP: spatial 2: SEVM -------
# compute distance matrix
locality_dist<-distm(ssp.space@coords,fun=distVincentyEllipsoid) #output in meters
# locality_dist2<-locality_dist
# locality_dist2[which(locality_dist2>=283*1000)]<-0 #approximation of threshold from SAM

# PCNM axes of the distance matrix
sevm<-pcnm(locality_dist)
(sevm$values/sum(sevm$values)*100)[1:10] #% variance explained? does that work with negative eigenvalues?
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

# # Visualize -------
# usa<-map_data('state')
# xlims<-range(ssp.space@bbox[1,])+c(-0.8,+0.8)
# ylims<-range(ssp.space@bbox[2,])+c(-0.8,+0.8)
# 
# choices<-which(colnames(sevm$vectors) %in% fs$variables)
# ssp.space.plot<-cbind(ssp_metadata,sevm$vectors[,choices])
# #code for color choice from  Tim van Werkhoven for Paul Tol diverging color scheme
# map_gradient = NULL 
# for (x in seq(0,1,length.out=nrow(ssp_metadata))){
#   rcol <- 0.237 - 2.13*x + 26.92*x**2 - 65.5*x**3 + 63.5*x**4 - 22.36*x**5
#   gcol <- ((0.572 + 1.524*x - 1.811*x**2)/(1 - 0.291*x + 0.1574*x**2))**2
#   bcol <- 1/(1.579 - 4.03*x + 12.92*x**2 - 31.4*x**3 + 48.6*x**4 - 23.36*x**5)
#   map_gradient<-rbind(cols,c(rcol,gcol,bcol))
# }
# map_colors<-map_gradient %>% rgb(.)
# 
# #can't figure out how to completely automate this plot in every view. Sorry.
# cairo_pdf("eigenmap1.pdf",width = 5.4, height = 4.4,family ="Arial",bg="transparent")
# ggplot(data=ssp.space.plot,aes(x=longitude,y=latitude))+
#   geom_polygon(data=usa,aes(x=long,y=lat,group=group),fill=NA,color="gray") +
#   geom_point(aes(fill=PCNM1),shape=21,size=2)+
#   # scale_color_manual(values= ssp_col,name="Subspecies",
#   #                   labels=c("T. c. bauri","T. c. carolina","T. c. major","T. c. triunguis"))
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
#   geom_point(aes(fill=PCNM2),shape=21,size=2)+
#   # scale_color_manual(values= ssp_col,name="Subspecies",
#   #                   labels=c("T. c. bauri","T. c. carolina","T. c. major","T. c. triunguis"))
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
#   geom_point(aes(fill=PCNM3),shape=21,size=2)+
#   # scale_color_manual(values= ssp_col,name="Subspecies",
#   #                   labels=c("T. c. bauri","T. c. carolina","T. c. major","T. c. triunguis"))
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
# 
# 
# cairo_pdf("eigenmap5.pdf",width = 5.4, height = 4.4,family ="Arial",bg="transparent")
# ggplot(data=ssp.space.plot,aes(x=longitude,y=latitude))+
#   geom_polygon(data=usa,aes(x=long,y=lat,group=group),fill=NA,color="gray") +
#   geom_point(aes(fill=PCNM5),shape=21,size=2)+
#   # scale_color_manual(values= ssp_col,name="Subspecies",
#   #                   labels=c("T. c. bauri","T. c. carolina","T. c. major","T. c. triunguis"))
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
# embed_fonts("eigenmap5.pdf")

# SSP: stats4ssp -----
#try again according to reviewer + Collyer recommendations
# top the EVs: posterior view: PCNM 1,2,5; lateral view: PCNM 1,2,3; dorsal view: PCNM 2,3,1
# go with 1 and 2, the two consistently important EVs? Or just 2?

ssp.gpa<-arrayspecs(ssp_set,p=ncol(ssp_set)/2,k=2) %>% gpagen
ssp2.gdf<-geomorph.data.frame(ssp.gpa,ssp=ssp_metadata$ssp,
                              PCNM1=sevm.test1$vectors[,1],PCNM2=sevm.test1$vectors[,2],
                              PCNM3=sevm.test1$vectors[,3],PCNM5=sevm.test1$vectors[,5])
ssp2.gdf$Csize<-ssp_metadata[,cs_metadata_col]

advanced.procD.lm(coords~(PCNM2+log(Csize)),  ~(PCNM2+log(Csize))+ssp, iter=499,
                  groups=~ssp,slope=~(PCNM2+log(Csize)),data=ssp2.gdf) #ssp not significant

advanced.procD.lm(coords~(PCNM3+log(Csize)),  ~(PCNM3+log(Csize))+ssp, iter=499,
                  groups=~ssp,slope=~(PCNM3+log(Csize)),data=ssp2.gdf) #ssp not significant

advanced.procD.lm(coords~(PCNM1+log(Csize)),  ~(PCNM1+log(Csize))+ssp, iter=499,
                  groups=~ssp,slope=~(PCNM1+log(Csize)),data=ssp2.gdf) #some ssp still significant


advanced.procD.lm(coords~(PCNM2+PCNM1+log(Csize)),
                  ~(PCNM2+PCNM1+log(Csize))+ssp,iter=499,data=ssp2.gdf) 



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

#try again according to reviewer + Collyer recommendations

advanced.procD.lm(coords~log(Csize),~log(Csize)+ssp,
                  groups=~ssp,slope=~log(Csize),iter=999,data=ssp.gdf) 
#result: dorsal: differences between triunguis-major and triunguis-carolina group
#       lateral: differences between bauri-carolina, bauri-triunguis, maybe major-triunguis?
#     posterior: nothing significant


advanced.procD.lm(coords~log(Csize)+ssp,~log(Csize)*ssp,
                  groups=~ssp,slope=~log(Csize),angle.type="deg",iter=499,data=ssp.gdf) 
#dorsal: only major/triunguis significant
#lateral: same as model above, plus slopes are different (means can't use a single allometric model?)
#posterior: no reason to check this one

procD.lm(coords~log(Csize)*ssp,iter=499,data=ssp.gdf) 
# also, as a note, there's an interaction, but significant for slope -- reason for more than one allometric model? 
# not in dorsal view, i guess yes for bauri for lateral view, but the slopes were different before adding the interaction?

pm.report<-pMAN$aov.table[,c(1,4,5,7)]
for(i in 1:nrow(pm.report)){pm.report[i,c(4,8,12,16,20,24,28)]<-p.adjust(pm.report[i,c(4,8,12,16,20,24,28)], method = "bonferroni", n = (6+1)*3)} #"BH

write.csv(pm.report,paste("ssp_procMANOVA_pairwise_",view,".csv",sep=""))

# SSP: size? -----
s1<-aov(log(ssp_all_metadata$carapace_length)~ssp_all_metadata$ssp) %>% summary
s2<-aov(log(ssp_all_metadata$dor_cs)~ssp_all_metadata$ssp) %>% summary
s3<-aov(log(ssp_all_metadata$lat_cs)~ssp_all_metadata$ssp) %>% summary
s4<-aov(log(ssp_all_metadata$pos_cs)~ssp_all_metadata$ssp) %>% summary

stest<-rbind(s1[[1]][1,],s2[[1]][1,],s3[[1]][1,],s4[[1]][1,]) %>% 
  # cbind(.,rep(c("carapace_length","dorsal","lateral","posterior"),each=2))
  cbind(.,c("carapace_length","dorsal","lateral","posterior"))

stest[,5]<-p.adjust(stest[,5], method = "BH") #adjust for multiple testing
write.csv(stest,"ssp_size_anova.csv")

# #explore
# subset<-which(ssp_all_metadata$ssp!="major"&ssp_all_metadata$ssp!="triunguis")
# aov(log(ssp_all_metadata$carapace_length[subset])~ssp_all_metadata$ssp[subset]) %>% summary
# aov(log(ssp_all_metadata$dor_cs[subset])~ssp_all_metadata$ssp[subset]) %>% summary
# aov(log(ssp_all_metadata$lat_cs[subset])~ssp_all_metadata$ssp[subset]) %>% summary
# aov(log(ssp_all_metadata$pos_cs[subset])~ssp_all_metadata$ssp[subset]) %>% summary

#make plotting object for boxplot
all_cs<-c(ssp_all_metadata$dor_cs,ssp_all_metadata$lat_cs,ssp_all_metadata$pos_cs,ssp_all_metadata$carapace_length) %>% log %>% data.frame(.,ncol=1)
all_cs$group<-c(rep("dorsal CS",nrow(ssp_all_metadata)),rep("lateral CL",nrow(ssp_all_metadata)),rep("posterior CS",nrow(ssp_all_metadata)),rep("carapace length",nrow(ssp_all_metadata)))
all_cs$ssp<-rep(ssp_all_metadata$ssp,4)

cairo_pdf("ssp_size.pdf",width = 5.4, height = 4,family ="Arial")
ggplot(data=all_cs,aes(x=ssp)) +
  geom_boxplot(aes(y=.,fill=group)) +
  
  # geom_boxplot(binaxis="y",binwidth=0.02,stackdir="up",
  #              aes(y=log(dor_cs),fill=spp_col[1])) +
  # geom_boxplot(binaxis="y",binwidth=0.02,stackdir="center",
  #              aes(y=log(lat_cs),fill=spp_col[2])) +
  # geom_boxplot(binaxis="y",binwidth=0.02,stackdir="down",
  #              aes(y=log(pos_cs),fill=spp_col[3])) +
  scale_fill_manual(name="Measurement", labels=c("carapace length","dorsal CS",
                                                 "lateral CS",
                                                 "posterior CS"),
                    values=spp_col[c(1,2,3,4)]) +
  xlab("Subspecies") + 
  ylab("ln(Carapace Length)") +
  theme_classic()
dev.off()  
embed_fonts("ssp_size.pdf")

# SSP: Bayesian Clustering? -------
classifier<-as.character(ssp_metadata$ssp)
cluster.data<-PCAssp$x[,1:PCc]
BIC<-mclustBIC(cluster.data,G=1:4)
plot(BIC)
one.cluster<-max(BIC[1,])
more.clusters<-max(BIC[c(2:4),],na.rm=TRUE) #BIC score of best model with more than one cluster
how.much<-one.cluster-more.clusters
best.k<-floor(which(BIC==max(BIC,na.rm=TRUE)) /4 ) #4 used because up to 4 clusters tested
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

# SSP: Assignments -----
# # # Use CVAGen8 for CVA-based assignments tests with jackknife
# # # write tables for the two groups
# write.table(file=paste("cva_binary_",view,".txt",sep=""),cbind(ssp_all,ssp_all_metadata[,cs_metadata_col]),col.names=FALSE,row.names=FALSE,sep='\t')
# write.table(file=paste("cva_binary_grp_",view,".txt",sep=""),as.numeric(binary_all),col.names=FALSE,row.names=FALSE)
# write group ID's to compare to traditional ssp
# write.table(file=paste("cva_ssp_",view,"_grp.txt",sep=""),as.numeric(ssp_all_metadata$ssp),col.names=FALSE,row.names=FALSE)
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


