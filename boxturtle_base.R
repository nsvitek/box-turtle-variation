# Base script to pull together geometric morphometric analyses of Terrapene spp.

scriptsdir <- "C://cygwin/home/N.S/scripts"
# datadir <- "D:/Dropbox/Documents/research/turtles/Thesis/Vitek_YR_PublishTerrapene/data"
datadir <- "C:/Users/N.S/Dropbox/Documents/research/turtles/Thesis/Vitek_YR_PublishTerrapene/data"

# Load Dependencies ------------------------------------------------------------------
# load library dependencies
library(dplyr)
library(vegan)
library(geomorph)
library(adegenet) #for Monmonier
library(wesanderson) #for colors
library(RColorBrewer) #for colors
library(plotrix) #for draw.ellipse()
library(ggplot2)
library(gridExtra) #for ggplot
library(mclust) #k-means clustering
library(sp) #making spatial object, coords()
library(spdep) #neighbors, all the x2y functions
library(geosphere) #distm
library(packfor) #forward selection with multivariate response
library(maps)
library(ggrepel)
library(ggmap)
library(scales) #needed to make reverse log-10 y axis, as solved by stackexchange Brian Diggs
# install.packages("packfor", repos="http://R-Forge.R-project.org")
# install.packages("mclust")

# load custom dependencies
setwd(scriptsdir)
source("observer-free-morphotype-characterization/anderson.R") #code published with Vitek et al. 2017
source("observer-free-morphotype-characterization/figpiece3d.R") #code published with Vitek et al. 2017
source("observer-free-morphotype-characterization/sensitivity_utils.R") #code published with Vitek et al. 2017
source("observer-free-morphotype-characterization/find_repeatablePCs.R") #code published with Vitek et al. 2017
source("box-turtle-variation/boxturtle_utils.R")

# load Julien Claude's code:
source("Morphometrics_with_R.R")
# source("common.slope.test.R") #no longer necessary with revision

#settings to get non-Helvetica fonts to work with ggplot 
library(extrafont)
# Adjust the path to match your installation of Ghostscript
Sys.setenv(R_GSCMD = "C:/Program Files (x86)/gs/gs9.21/bin/gswin32c.exe")

# Load Data ------------------------------------------------------------------
setwd(datadir)
data_all<-read.csv("boxturtle_Rdata.csv", header=TRUE)
#Remove specimen erroneously landmarked
data_all<-data_all[-which(data_all$id=="CMS8901"),]
age_measures<-read.csv("boxturtlemod_age_measures.csv",header=TRUE)
data_uni<-filter(data_all,err_rep==0) #get rid of error replicates
metadata_columns<-c(1,58,147,214:231) #which columns in total dataset have metadata
k<-2 #number of dimensions
ssp_sd_id<-read.csv("ssp_sd_id.csv",header=TRUE) #sex assessments
# bp<-0.05/(choose(4,2)*3) #bonferroni-corrected p-value
opar<-par

# Colors ----------------------------------------------------------------
# spp_col<-wes_palette("Moonrise2")
fos_col<-c(brewer.pal(12,"Paired"),"white")
ssp_col<-brewer.pal(8,"Accent") %>% rev()
# view_col<-wes_palette("Darjeeling")
oss_col<-c("#FFFAD2","#F9BD7E","#ED875E","#AE1C3E") #are these tol colors?

# Set View --------------------------------------------------------------
viewoptions<-c("dor","lat","pos")
view<-viewoptions[2]
source(paste(scriptsdir,"box-turtle-variation/boxturtle_settings.R",sep="/"))

# Error -------------------------------------------------------------------
# ##### Evaluate human digitization error
# freeze<-ls() #snapshot of current environment
# source(paste(scriptsdir,"boxturtle_error_eval.R",sep="/"))
# rm(list = c(setdiff(ls(),freeze),"data_all")) #clean up environment
# #also remove data_all: will not be used further.
# 
# setwd("../data-output")

# linear: run once -----
# freeze<-ls() #snapshot of current environment
# source(paste(scriptsdir,"boxturtle_linear.R",sep="/")
# rm(list = c(setdiff(ls(),freeze),"data_all")) #clean up environment

# Subspecies and Geography ------
#run through first set of ssp analyses for each view.

source(paste(scriptsdir,"boxturtle_ssp.R",sep="/"))
print("Subspecies code sourced.")

# Allometry and Dimorphism -----

source(paste(scriptsdir,"boxturtle_asd.R",sep="/"))
print("Geographic/Dimorphism code sourced.")

# F+SSP: make ------
fssp_set<-rbind(ssp_set,fos_set)
fssp_metadata<-mutate(ssp_metadata,site2=site) %>% rbind(.,fos_metadata) %>% droplevels
PCAfssp<-prcomp(fssp_set,scale.=FALSE)
an_fssp<-anderson(PCAfssp$sdev)

#reorder sites alphabetically, modern at end.
fssp_metadata$site2<-factor(as.character(fssp_metadata$site2))
fssp_metadata$site2<-factor(fssp_metadata$site2,levels(fssp_metadata$site2)[c(1:9,11:12,10)])

fssp.gpa<-arrayspecs(fssp_set,p=ncol(fssp_set)/2,k=2) %>% gpagen
fssp.gdf<-geomorph.data.frame(fssp.gpa,fossil=fssp_metadata$fos_set)
fssp.gdf$Csize<-fssp_metadata[,cs_metadata_col]

# F+SSP: allometry test -------

#do two groups have difference in allometric slopes?
advanced.procD.lm(coords~log(Csize)+fossil,~log(Csize)*fossil,
                  groups=~fossil,slope=~log(Csize),angle.type="deg",iter=999,data=fssp.gdf) 
#posterior: yes, both slopes and lengths; lateral: yes, slopes are different (vector lengths at 0.08)
#dorsal, yes, groups and slopes are different.

# advanced.procD.lm(coords~log(Csize),~log(Csize)+fossil,
#                   groups=~fossil,slope=~log(Csize),angle.type="deg",iter=499,data=fssp2.gdf) 

# F+SSP: allometry model ------
# build allometric linear model for modern specimens to apply to fossil record

model<-lm(as.matrix(ssp_set[,])~log(ssp_metadata[,cs_metadata_col])) #allometric model with CS
# model<-lm(as.matrix(ssp_all[,])~log(ssp_all_metadata$carapace_length[])) #allometric model with CL
modm<-matrix(model$coefficients[2,],m,k,byrow=TRUE) #m in the linear model y=mx+b
residuals<-model$residuals
modb<-matrix(model$coefficients[1,],m,k,byrow=TRUE) #b in the linear model y=mx+b

#when you remove predicted shape, what you are left are the residuals, yes?
testa<-function(x){modm*x+modb} #make modern linear model

# testa(log(ssp_metadata[,cs_metadata_col][188])) %>% plot(.,asp=1) #check


# "correct" for modern allometric predictions
#test
lcs<-log(fssp_metadata[,cs_metadata_col])
predicted.shape<-array(NA,dim=c(m,k,nrow(fssp_set)))
for(i in 1:nrow(fssp_set)){predicted.shape[,,i]<-testa(lcs[i])}
residual.shape<-arrayspecs(fssp_set,m,k)-predicted.shape
mean.shape<-colMeans(fssp_set)%>% matrix(.,m,k,byrow=TRUE) %>% array(.,dim=c(m,k,nrow(fssp_set)))
remaining.shape<-(mean.shape+residual.shape)

#test
fssp_metadata[which(fssp_metadata$site=="melbourne"),cs_metadata_col] %>% log
testa(6.6) %>% plot(.,asp=1) #check
points(matrix(fssp_set[236,],nrow=m,ncol=k,byrow=TRUE),col="red")
points(remaining.shape[,,236],col="blue") #these three shapes should make sense. Do they?


# F+SSP: allometry modelling ------
# F+SSP: PCA-----
# PCAfos.allo<-prcomp(two.d.array(remaining.shape),scale.=FALSE)
PCAfssp.allo<-prcomp(two.d.array(remaining.shape),scale.=FALSE)
an_fssp.allo<-anderson(PCAfssp.allo$sdev)
plot(PCAfssp.allo$x[,1:2],bg=c(fos_col[1:11],"white")[fssp_metadata$site2],pch=c(21,22)[fssp_metadata$fos_set],
     cex=(fssp_metadata[,cs_metadata_col]/1000)*3)
# text(PCAfos.allo$x[,1],PCAfos.allo$x[,2],c(1:nrow(fos_set)))
legend('bottomleft',cex=.7,pch=c(21),pt.bg=fos_col,legend=levels(fssp_metadata$site2),pt.cex=2)

# F+SSP: stats-----
advanced.procD.lm(coords~log(Csize),~log(Csize)+site,
                  groups=~site,slope=~log(Csize),iter=499,data=fssp2.gdf) 
procD.lm(coords~log(Csize),RRPP=TRUE,iter=499,data=fssp2.gdf)


# F+SSP: disparity ------
# individual.disparity(PCAfssp$x[,c(1:7)]) #95% variance explained
disp_total<-individual.disparity(PCAfssp$x[,c(1:PCc)]) #all non-zero PCs
# individual.disparity(PCAfssp$x[,]) #all data. 1st very different from 2nd & 3rd

observed<-PCAfssp$x[which(fssp_metadata$ssp_set=="yes"),c(1:PCc)] %>% individual.disparity() 
replicates<-1000 #eventually, put at 10000?

# now, essentially, rarefaction for disparity: is the ssp dataset size large enough to capture standing disparity?
sampleN<-seq(10,nrow(ssp_set),by=10)   #start with 10, but maybe go to 5?

rare_disparity<-matrix(NA,nrow=length(sampleN),ncol=2) %>% as.data.frame
colnames(rare_disparity)<-c("mean","sd")

for(N in 1:length(sampleN)){
  distribution<-NULL #create holder for null resampled distribution
  for (i in 1:replicates) {						# start for-loop
    rm1<- sample(c(1:nrow(ssp_set)),size=sampleN[N], replace=TRUE) %>% #resample ssp specimens with replacement
      PCAfssp$x[which(fssp_metadata$ssp_set=="yes"),c(1:PCc)][.,] %>% #take PC scores of resamples
      individual.disparity  # calculate disparity for resamples
    distribution<- rbind(distribution,rm1)					# store all resampled disparities
  }									# end for-loop
  rare_disparity$mean[N]<-distribution %>% mean
  rare_disparity$sd[N]<-distribution %>% sd
  print(paste("sample size",sampleN[N],"complete."))
}
# rare_disparity<-mutate(rare_disparity,lo = mean-(sd*1.96), hi = mean+(sd*1.96)) #1.96 for two-tailed
rare_disparity<-mutate(rare_disparity,lo = mean, hi = mean+(sd*1.645)) #1.645 for one-tailed

# calculate changes to disparity and sample size with addition of different fossil sites
fos_disparity<-matrix(NA,nrow=length(levels(fos_metadata$site2)),ncol=2) %>% as.data.frame
colnames(fos_disparity)<-c("d","n")

for(i in 1:length(levels(fos_metadata$site2))){
  select_disp<-which(fssp_metadata$site2=="modern"|fssp_metadata$site2==levels(fos_metadata$site2)[i])
  fos_disparity$d[i]<-individual.disparity(PCAfssp$x[select_disp,c(1:PCc)]) #all non-zero PCs
  fos_disparity$n[i]<-length(which(fssp_metadata$site2==levels(fos_metadata$site2)[i]))+nrow(ssp_set)
}

# # put it all together in a plot
# cairo_pdf(paste(view,"_disparity.pdf",sep=""),width = 5.4, height = 2,family ="Arial")
# par(mar=c(3,3.3,.5,.5))
# plot(sampleN,rare_disparity$mean,ylim=c(min(rare_disparity$lo),max(rare_disparity$hi)),
#      xlim=c(min(sampleN),nrow(fssp_set)),xlab="",ylab="") #wasted space, mess with limits
# polygon(c(sampleN,rev(sampleN)), 
#         c(rare_disparity$lo,rev(rare_disparity$hi)),
#         col=adjustcolor("gray",alpha=0.2), border = NA)
# lines(sampleN,rare_disparity$mean,lwd=1,lty=1)	
# points(sampleN,rare_disparity$mean,pch=21,bg="black",cex=1)
# points(fos_disparity$n,fos_disparity$d,pch=21,bg=fos_col)
# # legend('bottomleft',cex=.8,pch=c(21),pt.bg=fos_col,legend=levels(fos_metadata$site2),pt.cex=2)
# points(nrow(fssp_set),disp_total)
# mtext("Sample Size",side=1,line=2,cex=1)
# mtext("Disparity",side=2,line=2,cex=1)
# dev.off()
# embed_fonts(paste(view,"_disparity.pdf",sep=""))

cairo_pdf(paste("disparity_",view,".pdf",sep=""),width = 5.4, height = 2,family ="Arial")
par(mar=c(3,3.3,.5,.5))
plot(sampleN,rare_disparity$mean,ylim=c(mean(rare_disparity$mean),mean(rare_disparity$mean)+0.0007),
     xlim=c(100,nrow(fssp_set)),xlab="",ylab="") #wasted space, mess with limits
polygon(c(sampleN,rev(sampleN)),
        c(rare_disparity$lo,rev(rare_disparity$hi)),
        col=adjustcolor("gray",alpha=0.2), border = NA)
lines(sampleN,rare_disparity$mean,lwd=1,lty=1)
points(sampleN,rare_disparity$mean,pch=21,bg="black",cex=1)
points(fos_disparity$n,fos_disparity$d,pch=21,bg=fos_col)
# legend('bottomleft',cex=.8,pch=c(21),pt.bg=fos_col,legend=levels(fos_metadata$site2),pt.cex=2)
points(nrow(fssp_set),disp_total,pch=22,bg="black")
mtext("Sample Size",side=1,line=2,cex=1)
mtext("Disparity",side=2,line=2,cex=1)
dev.off()
embed_fonts(paste("disparity_",view,".pdf",sep=""))

cairo_pdf("disparity_legend.pdf",width = 5.4, height = 7,family ="Arial")
plot(sampleN,rare_disparity$mean,type="n") #wasted space, mess with limits
legend('center',cex=.8,pch=c(21),pt.bg=fos_col,legend=levels(fos_metadata$site2),pt.cex=2)
dev.off()
embed_fonts("disparity_legend.pdf")



# # F: PCA ---------
# # Principal Components Analysis (look at data)
PCAfos<-prcomp(fos_set,scale.=FALSE)
# an<-anderson(PCAfos$sdev)
# 
# cairo_pdf(paste("fos_PCA_",view,".pdf",sep=""),width = 4.4, height = 4.4,family ="Arial")
# par(mar=c(3.9,3.6,.5,.5))
# plot(PCAfos$x[,c(1,2)],cex=1.5,
#      bg=fos_col[fos_metadata$site2],pch=21,
#      xlab="",ylab="")
# mtext(paste("PC 1 (",round(an$percent[1],3),"%)",sep=""),side=1,line=3,cex=2)
# mtext(paste("PC 2 (",round(an$percent[2],3),"%)",sep=""),side=2,line=2,cex=2)
# # legend('bottomleft',cex=.8,pch=c(21),pt.bg=fos_col,legend=levels(fos_metadata$site2),pt.cex=2)
# dev.off()
# 
# pdf("fos_shape_legend.pdf",useDingbats=FALSE) #make legend
# plot(PCAfos$x[,c(1,2)],type="n",axes=F,xlab="",ylab="")
# legend('center',cex=2,pch=c(21),pt.cex=3,pt.bg=fos_col,legend=levels(fos_metadata$site2))
# dev.off()

# fos_metadata[,cs_metadata_col]/(min(fos_metadata[,cs_metadata_col]))

#is this appropriate?
fos.gpa<-remaining.shape %>% gpagen
fos.gdf<-geomorph.data.frame(fos.gpa,site2=fos_metadata$site2)
fos.gdf$Csize<-fos_metadata[,cs_metadata_col]
advanced.procD.lm(coords~1,~site3,groups=~site3,iter=999,data=fos.gdf)


# F: k-means clustering -------
classifier<-as.character(fssp_metadata$site2)
X<-PCAfssp.allo$x[which(fssp_metadata$fos_set=="yes"),1:PCc]
mc<-5 #maximum cluster
BIC<-mclustBIC(X,G=1:mc)
plot(BIC)
one.cluster<-max(BIC[1,])
more.clusters<-max(BIC[c(2:4),],na.rm=TRUE) #BIC score of best model with more than one cluster
how.much<-one.cluster-more.clusters
best.k<-((which(BIC==max(BIC,na.rm=TRUE)) /mc)-floor((which(BIC==max(BIC,na.rm=TRUE)) /mc)))*mc #5 used because up to 5 clusters tested
kmeans.report<-rbind(one.cluster,more.clusters,how.much,best.k)
row.names(kmeans.report)<-c("BIC.best.one.cluster",
                            "BIC.best.multiple.clusters",
                            "deltaBIC","best.k")
write.csv(kmeans.report,paste("kmeans_results_",view,".csv",sep=""))

# procD.lm(coords~log(Csize)+site4,RRPP=TRUE,iter=499,data=fos.gdf)
# procD.lm(coords~log(Csize)+site5,RRPP=TRUE,iter=499,data=fos.gdf)
test.maybe<-advanced.procD.lm(coords~log(Csize),~log(Csize)+site3,
                              groups=~site3,slope=~log(Csize),iter=999,data=fos.gdf)
summary(test.maybe)

advanced.procD.lm(coords~log(Csize),~log(Csize)+site3,
                  groups=~site3,slope=~log(Csize),iter=999,data=fos.gdf)

advanced.procD.lm(coords~1,~site3,
                  groups=~site3,iter=999,data=fos.gdf)

#in lateral view Camelot2 is different if you take size into account, 
#if you don't take size into account, large ones are sort of/mostly different from much smaller ones, esp. Melbourne
#in posterior view, nothing is different with size into account, 
#if you don't take size into account, similar pattern to lateral -- differences seem to scale by size
#in dorsal view, with size, a random diff b/t ardis and camelot, ardis and haile, ardis and reddick
#without size, camelot2/friesanhahn, camelot2/haile, holocene/friesenhahn, melbourne/friesenhahn, seems random.

fssp.gpa<-arrayspecs(fssp_set,p=ncol(fssp_set)/2,k=2) %>% gpagen
# fssp2.gdf<-geomorph.data.frame(fssp.gpa,site=c(as.character(site5),rep("holocene",nrow(ssp_metadata))))
# fssp2.gdf<-geomorph.data.frame(fssp.gpa,site=c(as.character(site3),as.character(ssp_metadata$site)))
fssp2.gdf<-geomorph.data.frame(fssp.gpa,fossil=fssp_metadata$fos_set)
fssp2.gdf$Csize<-fssp_metadata[,cs_metadata_col]

#do two groups have difference in allometric slopes?
advanced.procD.lm(coords~log(Csize)+fossil,~log(Csize)*fossil,
                  groups=~fossil,slope=~log(Csize),angle.type="deg",iter=499,data=fssp2.gdf) 
#posterior: yes, both slopes and lengths; lateral: yes, slopes are different (vector lengths at 0.08)
#dorsal, yes, groups and slopes are different.

advanced.procD.lm(coords~log(Csize),~log(Csize)+fossil,
                  groups=~fossil,slope=~log(Csize),angle.type="deg",iter=499,data=fssp2.gdf) 


advanced.procD.lm(coords~log(Csize),~log(Csize)+site,
                  groups=~site,slope=~log(Csize),iter=499,data=fssp2.gdf) 
procD.lm(coords~log(Csize),RRPP=TRUE,iter=499,data=fssp2.gdf)

advanced.procD.lm(coords~1,~site,
                  groups=~site,iter=499,data=fssp2.gdf) 



# F: stats -----
#is this appropriate?
fssp.gpa<-remaining.shape %>% gpagen
fssp.gdf<-geomorph.data.frame(fssp.gpa,site2=fssp_metadata$site2,fossil=fssp_metadata$fos_set,
                              site3=c(rep("modern",nrow(ssp_set)),as.character(site3)))
# fos.gdf$Csize<-fos_metadata[,cs_metadata_col]
advanced.procD.lm(coords~1,~site2,groups=~site2,iter=999,data=fssp.gdf)


# # appendix -----
# #make lists of specimens in each dataset
# 
# #roundabout, but works
# raw_id<-age_measures %>% mutate(.,id=paste(Institution,Number,sep="")) %>% 
#   select(.,id) %>% arrange(.,id) #
# intermediate_id<-matrix("a",nrow=nrow(raw_id),ncol=2)
# for (row in 1:nrow(raw_id)){
#   intermediate_id[row,1]<-gsub("([A-Z/-]+)(.*)", "\\1",raw_id[row,], perl = TRUE)
#   intermediate_id[row,2]<-gsub("([A-Z/-]+)(.*)", "\\2",raw_id[row,], perl = TRUE)
# }
# colnames(intermediate_id)<-c("institution","number")
# intermediate_id[multi.mixedorder(intermediate_id[,1],intermediate_id[,2]),] %>% as.data.frame %>% 
#   mutate(.,id=paste(institution,number,sep=" ")) %>% select(.,id) %>% write.csv(.,"appendix_mat.csv")
# 
# #sd
# intermediate_id<-matrix(dim_metadata$id,nrow=length(dim_metadata$id),ncol=2)
# for (row in 1:length(dim_metadata$id)){
#   intermediate_id[row,1]<-gsub("([A-Z/-]+)(.*)", "\\1",dim_metadata$id[row], perl = TRUE)
#   intermediate_id[row,2]<-gsub("([A-Z/-]+)(.*)", "\\2",dim_metadata$id[row], perl = TRUE)
# }
# colnames(intermediate_id)<-c("institution","number")
# intermediate_id[multi.mixedorder(intermediate_id[,1],intermediate_id[,2]),] %>% as.data.frame %>% 
#   mutate(.,id=paste(institution,number,sep=" ")) %>% select(.,id) %>% write.csv(.,"appendix_dim.csv")
# 
# #ssp/geography
# intermediate_id<-matrix(ssp_metadata$id,nrow=length(ssp_metadata$id),ncol=2)
# for (row in 1:length(ssp_metadata$id)){
#   intermediate_id[row,1]<-gsub("([A-Z/-]+)(.*)", "\\1",ssp_metadata$id[row], perl = TRUE)
#   intermediate_id[row,2]<-gsub("([A-Z/-]+)(.*)", "\\2",ssp_metadata$id[row], perl = TRUE)
# }
# colnames(intermediate_id)<-c("institution","number")
# intermediate_id[multi.mixedorder(intermediate_id[,1],intermediate_id[,2]),] %>% as.data.frame %>% 
#   mutate(.,id=paste(institution,number,sep=" ")) %>% select(.,id) %>% write.csv(.,"appendix_ssp.csv")
# 
# #fossils
# intermediate_id<-matrix(fos_metadata$id,nrow=length(fos_metadata$id),ncol=2)
# for (row in 1:length(fos_metadata$id)){
#   intermediate_id[row,1]<-gsub("([A-Z/-]+)(.*)", "\\1",fos_metadata$id[row], perl = TRUE)
#   intermediate_id[row,2]<-gsub("([A-Z/-]+)(.*)", "\\2",fos_metadata$id[row], perl = TRUE)
# }
# colnames(intermediate_id)<-c("institution","number")
# intermediate_id[multi.mixedorder(intermediate_id[,1],intermediate_id[,2]),] %>% as.data.frame %>% 
#   mutate(.,id=paste(institution,number,sep=" ")) %>% select(.,id) %>% write.csv(.,"appendix_fos.csv")
# 
# 
# library(gtools)
# multi.mixedorder <- function(..., na.last = TRUE, decreasing = FALSE,stringsAsFactors=FALSE){
#   do.call(order, c(
#     lapply(list(...), function(l){
#       if(is.character(l)){
#         factor(l, levels=mixedsort(unique(l)))
#       } else {
#         l
#       }
#     }),
#     list(na.last = na.last, decreasing = decreasing)
#   ))
# }
# 
# # glyphs ---------
# # try out using glyphs to express centroid size
# #okay, so we essentially need a time series for each sequence that we can wrap into a star using polar coords
# #2*points on star + 1, making x_minor c(1,13) for a 6-pointed star
# # and the points need to scale [-1,1]
# # so the y_minor should go from (-measure/max(measure)) to (measure/max(measure))
# library(reshape)
# library(GGally)
# cs.mut<-log(fos_metadata[,cs_metadata_col]) %>% as.character %>% as.numeric %>% abs(.) #* -1
# PCAll<-prcomp(fos_set,scale.=FALSE)
# glyph.dat1<-cbind(PCAll$x[,1:2],cs.mut,fos_metadata$site2,fos_metadata$state) %>% as.data.frame
# colnames(glyph.dat1)<-c("PC1","PC2","cs","site2","dummy")
# range01 <- function(x){(x-min(x))/(max(x)-min(x))} #function from stack exchange to rescale [0,1]
# glyph.dat2<-mutate(glyph.dat1,stmin=1,stmax=(range01(cs))) #scale longitude so min is at 0 and max is at +-1
# glyph.dat3<-glyph.dat2[,c(1:ncol(glyph.dat2),
#                         rep(c(which(colnames(glyph.dat2)=="stmin"|colnames(glyph.dat2)=="stmax")),10),# number of reps controls shape, butin a weird polar way
#                         rep(c(which(colnames(glyph.dat2)=="stmin")),1))]
# glyph.dat4<-melt(glyph.dat3,id.vars=c("PC1","PC2","dummy","site2","cs"))
# # the Q in 2*Q+3 should match the integer in rep((1|2),Q) two lines above. May need to tweak Q to get shape.
# #in dor Q should be 20, in lat and pos should be 10
# glyph.dat5<-mutate(glyph.dat4,x=rep(seq(from=1,to=(2*10+3)*10,by=10),nrow(glyph.dat1))) 
# 
# glyph.shapes<-glyphs(data=glyph.dat5,x_major="PC1",y_minor="value",x_minor="x",y_major="PC2",
#                   polar=TRUE,height=0.005,width=0.005)
# 
# cairo_pdf(paste("fos_PCA2_",view,".pdf",sep=""),width = 4.4, height = 4.4,family ="Arial")
# ggplot(glyph.shapes, ggplot2::aes(gx, gy, group = gid,fill=factor(site2))) +
#   geom_path() +
#   geom_polygon()+
#   theme_bw() +
#   labs(x = "PC1", y = "PC2") +
#   scale_fill_manual(values=fos_col) +
#   coord_fixed(1.1)
# dev.off()
# delete me? -------
# site4<-site3
# levels(site4)<-c(levels(site3),"texas","florida")
# site4[which(site4=="friesenhahn cave"|site4=="ingleside")]<-"texas"
# site4[which(site4=="haile 8a"|site4=="reddick 1b")]<-"florida"
# site4<-droplevels(site4)
# site5<-site4
# site5[which(site5!="camelot2")]<-"holocene"
# site5<-droplevels(site5)

# # F: Haile -----
# # are Auffenberg's two Haile types significantly different? No. p > 0.1 all cases. 
# two.hailes<-as.character(fos_metadata$site)
# two.hailes[which(fos_metadata$id=="UF3148"|
#                    fos_metadata$id=="UF3150"|
#                    fos_metadata$id=="UF3143")]<-"haile2" #based on Auffenberg 1967
# subset<-which(two.hailes=="haile 8a"|two.hailes=="haile2")
# fos_metadata$id[subset]
# adonis(PCAfos$x[subset,c(1:PCc)]~two.hailes[subset],permutations=1000,method="euclidean")
# 
# # look at mean shapes vs each other and grand mean
# hailemean<-mshp(fos_set[which(fos_metadata$site=="haile 8a"),])$meanshape
# haile1mean<-mshp(fos_set[which(two.hailes=="haile 8a"),])$meanshape
# haile2mean<-mshp(fos_set[which(two.hailes=="haile2"),])$meanshape
# 
# #plot two comparative frameworks for proposed Haile morphs
# cairo_pdf(paste("haile_shape_",view,".pdf",sep=""),width = 7.6, height = 7.6)
# par(mfrow=c(1,2))
# par(mar=c(.5,.5,.5,.5))
# # could it be dimorphism? What are the odds of sampling only 2 males and 3 females or vice versa?
#this isn't an appropriate way to model--fix, eventually. 
# xlims<-c(min(c(hailemean[,1]-sd_female2[,1],hailemean[,1]-sd_male2[,1])),
#          max(c(hailemean[,1]+sd_female2[,1],hailemean[,1]+sd_male2[,1])))
# ylims<-c(min(c(hailemean[,2]-sd_female2[,2],hailemean[,2]-sd_male2[,2])),
#          max(c(hailemean[,2]+sd_female2[,2],hailemean[,2]+sd_male2[,2])))
# plot(hailemean,asp=1,axes=F,xlab="",ylab="",pch=21,bg="black")
# plot.errell(haile2mean,sd_male2,color.bg=spp_col[1],alpha.bg=0.5)
# plot.errell(haile1mean,sd_female2,color.bg=spp_col[2],alpha.bg=0.5)
# points(haile1mean[,1],haile1mean[,2],pch=21,bg="red",cex=0.6)
# points(haile2mean[,1],haile2mean[,2],pch=21,bg="blue",cex=0.6)
# 
# # not that different in size, use single allometric model
# temp<-fos_site_vs_ssp("haile 8a")
# xlims<-c(min(c(temp[,1,1]-sd_model2[,1,1],temp[,1,2]-sd_model2[,1,2])),
#          max(c(temp[,1,1]+sd_model2[,1,1],temp[,1,2]+sd_model2[,1,2])))
# ylims<-c(min(c(temp[,2,1]-sd_model2[,2,1],temp[,2,2]-sd_model2[,2,2])),
#          max(c(temp[,2,1]+sd_model2[,2,1],temp[,2,2]+sd_model2[,2,2])))
# # pdf(paste(view,levels(fos_metadata$site2)[i],"shape.pdf",sep="_"))
# #
# plot(temp[,,1],asp=1,axes=F,xlab="",ylab="",xlim=xlims,ylim=ylims,type="n")
# for (j in c(2,3)){
#   plot.errell(temp[,,j],sd_model2[,,j-1],color.bg=ssp_col[j-1],alpha.bg=0.6)
# }
# temppts2<-arrayspecs(fos_set, m, k)
# points(temppts2[,1,which(two.hailes=="haile 8a")],temppts2[,2,which(two.hailes=="haile 8a")],pch=21,bg="red",cex=0.6)
# points(temppts2[,1,which(two.hailes=="haile2")],temppts2[,2,which(two.hailes=="haile2")],pch=21,bg="blue",cex=0.6)
# # points(temp[,,1],pch=21,bg="black",cex=1)
# dev.off()
