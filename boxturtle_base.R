# Base script to pull together geometric morphometric analyses of Terrapene spp.

scriptsdir <- "C://cygwin/home/N.S/scripts"
# datadir <- "D:/Dropbox/Documents/research/turtles/Thesis/Vitek_YR_PublishTerrapene/data"
datadir <- "C:/Users/N.S/Dropbox/Documents/research/turtles/Thesis/Vitek_YR_PublishTerrapene/data"

# Load Dependencies ------------------------------------------------------------------
setwd(scriptsdir)
source("boxturtle_utils.R")
source("observer-free-morphotype-characterization/anderson.R")

# load dependencies
library(dplyr)
library(vegan)
library(geomorph)
library(adegenet) #for Monmonier
library(wesanderson) #for colors
library(RColorBrewer) #for colors
# library(scales) #for emulating ggplot colors
library(plotrix) #for draw.ellipse()
library(ggplot2)
library(gridExtra) #for ggplot
# install.packages("gridExtra")

setwd(scriptsdir)
source("observer-free-morphotype-characterization/anderson.R")
source("observer-free-morphotype-characterization/figpiece3d.R")
source("observer-free-morphotype-characterization/sensitivity_utils.R")
source("observer-free-morphotype-characterization/find_repeatablePCs.R")
source("boxturtle_utils.R")

#plus Julien Claude's code and common slopes test:
source("Morphometrics_with_R.R")
source("common.slope.test.R")

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
bp<-0.05/(choose(4,2)*3) #bonferroni-corrected p-value
opar<-par

# Colors ----------------------------------------------------------------
spp_col<-wes_palette("Moonrise2")
fos_col<-c(brewer.pal(12,"Paired"),"white")
ssp_col<-brewer.pal(8,"Accent") %>% rev()
view_col<-wes_palette("Darjeeling")
oss_col<-c("#FFFAD2","#F9BD7E","#ED875E","#AE1C3E")

# Set View --------------------------------------------------------------
viewoptions<-c("dor","lat","pos")
view<-viewoptions[3]
source(paste(scriptsdir,"boxturtle_settings.R",sep="/"))

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
binary_fssp<-c(binary,rep(3,nrow(fos_set)))

#reorder sites alphabetically, modern at end.
fssp_metadata$site2<-factor(as.character(fssp_metadata$site2))
fssp_metadata$site2<-factor(fssp_metadata$site2,levels(fssp_metadata$site2)[c(1:9,11:12,10)])

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


# F+SSP: stats-----
stats<-matrix(NA,nrow=length(levels(fos_metadata$site2)),ncol=4)

#test fossil sites against peninsular vs. more inland/continental box turtles
for (y in 1:length(levels(fos_metadata$site2))){
  subset<-list()
  subset[[1]]<-which(binary_fssp==1|fssp_metadata$site2==levels(fos_metadata$site2)[y]) #2=carolina
  subset[[2]]<-which(binary_fssp==2|fssp_metadata$site2==levels(fos_metadata$site2)[y]) #2=carolina
  for(z in 1:2){
    # make geomorph object for subset of critters to be tested
    fssp.gpa<-arrayspecs(fssp_set[subset[[z]],],p=ncol(fssp_set[subset[[z]],])/2,k=2) %>% gpagen
    fssp.gdf<-geomorph.data.frame(fssp.gpa,binary=binary_fssp[subset[[z]]],site=fssp_metadata$site2[subset[[z]]])
    fssp.gdf$Csize<-fssp_metadata[subset[[z]],cs_metadata_col]
    # test
    pMAN<-procD.lm(coords ~ log(Csize)*site,
                    iter = 999, data = fssp.gdf) #report R^2,F,p, df
    stats[y,z*2-1]<-pMAN$aov.table[2,5]
    stats[y,z*2]<-pMAN$aov.table[2,7]
    
  }
}

row.names(stats)<-levels(fos_metadata$site2)
colnames(stats)<-c("F_peninsular","p_peninsular","F_mainland","p_mainland")
stats[,c(2,4)]<-p.adjust(stats[,c(2,4)], method = "BH", n = length(stats[,c(2,4)])*3)
write.csv(stats,paste("fos_binary_procMANOVA_",view,".csv",sep=""))


# F+SSP: PCA-----

cairo_pdf(paste("fos_binary_PCA_",view,".pdf",sep=""),width = 4.4, height = 4.4,family ="Arial")
par(mar=c(3.9,3.6,.5,.5))
plot(PCAfssp$x[,1:2],pch=20+binary_fssp,bg=c(fos_col[1:11],"white")[fssp_metadata$site2],
     cex=log(fssp_metadata$carapace_length)/(mean(log(fssp_metadata$carapace_length))),xlab="",ylab="")
mtext(paste("PC 1 (",round(an_fssp$percent[1],3),"%)",sep=""),side=1,line=3,cex=2)
mtext(paste("PC 2 (",round(an_fssp$percent[2],3),"%)",sep=""),side=2,line=2,cex=2)
dev.off()
embed_fonts(paste("fos_binary_PCA_",view,".pdf",sep=""))

pdf("fos_ssp_legend.pdf",useDingbats=FALSE) #make legend
plot(PCAfssp$x[,1:2],type="n",axes=F,xlab="",ylab="")
legend('center',cex=1.5,pch=c(rep(23,11),21,22),pt.cex=3,pt.bg=c(fos_col[1:11],rep("white",2)),
       legend=c(levels(fssp_metadata$site2)[1:11],"modern peninsular","modern other"))
dev.off()
embed_fonts("fos_ssp_legend.pdf")


# F+SSP: PCA v2 -------
ternary<-binary_fssp
ternary[which(binary_fssp==1|
                fssp_metadata$site2=="camelot2"|
                fssp_metadata$site2=="devils den"|
                fssp_metadata$site2=="holocene"|
                fssp_metadata$site2=="vero")]<-"peninsular"
ternary[which(binary_fssp==2|
                fssp_metadata$site2=="melbourne")]<-"mainland"
ternary[which(fssp_metadata$site2=="ardis"|
              fssp_metadata$site2=="camelot"|
              fssp_metadata$site2=="friesenhahn cave"|
              fssp_metadata$site2=="haile 8a"|
              fssp_metadata$site2=="ingleside"|
              fssp_metadata$site2=="reddick 1b")]<-"neither"
ternary<-factor(ternary, levels=c("peninsular","mainland","neither"))

cairo_pdf(paste("fos_binary_PCAv2_",view,".pdf",sep=""),width = 4.4, height = 4.4,family ="Arial")
plot(PCAfssp$x[,1:2],pch=20+as.numeric(fssp_metadata$fos_set),bg=ssp_col[ternary],
     cex=log(fssp_metadata$carapace_length)/(mean(log(fssp_metadata$carapace_length))),xlab="",ylab="")
mtext(paste("PC 1 (",round(an_fssp$percent[1],3),"%)",sep=""),side=1,line=3,cex=2)
mtext(paste("PC 2 (",round(an_fssp$percent[2],3),"%)",sep=""),side=2,line=2,cex=2)
dev.off()
embed_fonts(paste("fos_binary_PCAv2_",view,".pdf",sep=""))

# F: stats -----
fos_pairs<-combn(levels(ternary),2)
stats2<-matrix(NA,nrow=ncol(fos_pairs),ncol=2)
rownames(stats2)<-paste(fos_pairs[1,],fos_pairs[2,],sep="_")
colnames(stats2)<-c("F","p")
for (y in 1:ncol(fos_pairs)){
  subset<-which(fssp_metadata$fos_set=="yes"&(ternary==fos_pairs[1,y]|ternary==fos_pairs[2,y]))
  subset<-which((ternary==fos_pairs[1,y]|ternary==fos_pairs[2,y])) 
  fssp.gpa<-arrayspecs(fssp_set[subset,],p=ncol(fssp_set[subset,])/2,k=2) %>% gpagen
  fssp.gdf<-geomorph.data.frame(fssp.gpa,binary=binary_fssp[subset],site=fssp_metadata$site2[subset])
  fssp.gdf$Csize<-fssp_metadata[subset,cs_metadata_col]
  # test
  pMAN<-procD.lm(coords ~ log(Csize)*site,
                 iter = 999, data = fssp.gdf) #report R^2,F,p, df
  stats2[y,1]<-pMAN$aov.table[2,5]
  stats2[y,2]<-pMAN$aov.table[2,7]
}
stats2[,2]<-p.adjust(stats2[,2], method = "BH", n = length(stats2[,2])*3)

write.csv(stats2,paste("fos_ternary_procMANOVA_",view,".csv",sep=""))

# F: allometry modelling ------
# custom function for making site-specific models for plotting
fos_site_vs_ssp<-function(searchterm){
  answer<-array(NA,dim=c(m,k,3))
  mlCS<-which(fos_metadata$site2==searchterm) %>% fos_metadata[.,cs_metadata_col] %>% mean %>% log
  answer[,,1]<-mshp(fos_set[which(fos_metadata$site2==searchterm),])$meanshape
  answer[,,2]<-testb(mlCS)
  answer[,,3]<-testn(mlCS)
  return(answer)
}

cairo_pdf(paste("fos_shape_",view,"2.pdf",sep=""),width = 7.6, height = 7.6) #/1.83
par(mfrow=c(3,4))
par(mar=c(.5,.5,.5,.5))
for(i in c(1:length(levels(fos_metadata$site2)))){
  temp<-fos_site_vs_ssp(levels(fos_metadata$site2)[i])
  xlims<-c(min(c(temp[,1,1]-sd_model2[,1,1],temp[,1,2]-sd_model2[,1,2])),
           max(c(temp[,1,1]+sd_model2[,1,1],temp[,1,2]+sd_model2[,1,2])))
  ylims<-c(min(c(temp[,2,1]-sd_model2[,2,1],temp[,2,2]-sd_model2[,2,2])),
           max(c(temp[,2,1]+sd_model2[,2,1],temp[,2,2]+sd_model2[,2,2])))
  # pdf(paste(view,levels(fos_metadata$site2)[i],"shape.pdf",sep="_"))
  # 
  plot(temp[,,1],asp=1,axes=F,xlab="",ylab="",xlim=xlims,ylim=ylims,type="n")
  for (j in c(2,3)){
    plot.errell(temp[,,j],sd_model2[,,j-1],color.bg=ssp_col[j-1],alpha.bg=0.6)
  }
  temppts<-fos_set[which(fos_metadata$site2==levels(fos_metadata$site2)[i]),]
  temppts2<-arrayspecs(temppts, m, k)
  # points(temppts2[,1,],temppts2[,2,],pch=21,bg=ssp_col[3],cex=0.5)
  # points(temp[,,1],pch=22,bg=ssp_col[3],cex=1)
  points(temppts2[,1,],temppts2[,2,],pch=21,bg="black",cex=0.3)
  points(temp[,,1],pch=21,bg="black",cex=1)
  # legend('bottomright',cex=.8,pch=c(21),pt.bg=fos_col,legend=levels(fos_metadata$site2),pt.cex=2)
  title(levels(fos_metadata$site2)[i])
  # dev.off()
}
dev.off()

samplesize<-fos_metadata %>% 
  group_by(site2) %>% 
  summarise(n=length(id))


fos_metadata$site2
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
# # F: PCA ---------
# # Principal Components Analysis (look at data)
# PCAfos<-prcomp(fos_set,scale.=FALSE)
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
# # F: bgPCA -----
# #bgPCA
# fos_residuals_mean<-colMeans(fos_set)#compute grand mean
# fos_residuals_centered<-as.matrix(fos_set)-rep(1,nrow(fos_set))%*%t(fos_residuals_mean) #Subtract that grandmean from each row
# 
# #Calculate the group means
# fos_means<-array(NA,dim=c(length(levels(fos_metadata$site2)),ncol(fos_set)))
# for (i in 1:length(levels(fos_metadata$site2))){
#   fos_means[i,]<-colMeans(fos_set[which(fos_metadata$site2==levels(fos_metadata$site2)[i]),])
# }
# B<-prcomp(fos_means)$rotation #eigenvectors
# B2<-prcomp(fos_means)
# 
# bgPCAfos<-fos_residuals_centered%*%B #Get the scores for all the individuals on the eigenvectors of the PCA of the means
# 
# cairo_pdf(paste("fos_bgPCA_",view,".pdf",sep=""),width = 4.4, height = 4.4,family ="Arial")
# par(mar=c(3.9,3.6,.5,.5))
# plot(bgPCAfos[,c(1,2)],cex=1.5,
#      bg=fos_col[fos_metadata$site2],pch=21,xlab="",ylab="")
# mtext("bgPC 1",side=1,line=3,cex=2)
# mtext("bgPC 2",side=2,line=2,cex=2)
# # points(B2$x[,c(1,2)],cex=1.5,bg=fos_col,pch=22)
# # legend('topleft',cex=.8,pch=c(21),pt.bg=fos_col,legend=levels(fos_metadata$site2),pt.cex=2)
# dev.off()

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


# # F: tree? -----
# library(phangorn)
# 
# an2<-anderson(B2$sdev)
# 
# D1<-dist(B2$x[,c(1:nPC(an2))], "euclidean")
# # D1<-dist(PCAfos$x[,c(1:PCc)], "euclidean") #better? c(1:PCc)
# # perform the clustering;
# cluster1<-hclust(D1, method="average") #other option is ward.D, ward.D2 or "average" for UPGMA
# # compute the distances along the branches
# D2<-cophenetic(cluster1) #for correlation testing, used in development
# # compute the correlation and write it to the screen
# cor(D1,D2) #<0.85 indicates significant distortion
# ## in all cases, correlation < 0.7
# # in tests with marsupial rep 180, "average" had highest correlation
# # draw the dendrogram
# # cluster1$labels<-as.character(fos_metadata$site2)
# cluster1$labels<-as.character(levels(fos_metadata$site2))
# plot(cluster1)
# 
# library(ape)
# # clusterpos<-cluster1
# # clusterlat<-cluster1
# # clusterdor<-cluster1
# test1<-consensus(c(as.phylo(clusterpos),as.phylo(clusterdor)))
# # ,as.phylo(clusterlat)
# plot(test1) #there is no consensus except ardis+camelot, even when using only >0.85 clusters
# 
