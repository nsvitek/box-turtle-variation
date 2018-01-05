# SD: PCA -----
# Initial exploration of data structure
PCAdim<-prcomp(dim_set,scale.=FALSE)
an_dim<-anderson(PCAdim$sdev)
cairo_pdf(paste("sd_PCA_",view,".pdf",sep=""),width = 1.8, height = 1.8,
          family ="Arial")
par(mar=c(3,3.3,.5,.5))
plot(PCAdim$x[,1:2],cex=dim_metadata[,cs_metadata_col]/(min(dim_metadata[,cs_metadata_col])),
     bg=ssp_col[dim_metadata$ssp],pch=(20+as.numeric(dim_metadata$sex)),
     xlab="",ylab="")
mtext(paste("PC 1 (",round(an_dim$percent[1],3),"%)",sep=""),side=1,line=2,cex=1)
mtext(paste("PC 2 (",round(an_dim$percent[2],3),"%)",sep=""),side=2,line=2,cex=1)
# legend('topleft',cex=1,pch=21,pt.bg=ssp_col,legend=levels(dim_metadata$ssp))
# legend('topright',cex=1,pch=c(21,22),pt.bg=ssp_col[1],legend=levels(dim_metadata$sex))
# title("Exploration: Sexual Dimorphism\nspecimens scaled by CS")
dev.off()
embed_fonts(paste("sd_PCA_",view,".pdf",sep=""))

# SD: statistics -----
#statistics: after accounting for size and geographic groups, is there sexual dimorphism?

#check, interaction beween the two covariates? Answer: No. 
advanced.procD.lm(coords ~ log(Csize) + ssp, ~log(Csize) * ssp,
                  groups=~ssp,slope=~log(Csize),iter = 999, data = dim.gdf) #report R^2,F,p, df

#no interactions here.
advanced.procD.lm(coords ~ (log(Csize) + ssp) + sex, ~(log(Csize) + ssp) * sex,
                  groups=~sex,slope=~log(Csize),iter = 999, data = dim.gdf) #report R^2,F,p, df


pMAN<-procD.lm(coords ~ (log(Csize) + ssp) * sex,iter = 999, data = dim.gdf)

for(i in 1:nrow(pMAN$aov.tab)){pMAN$aov.tab[i,7]<-p.adjust(pMAN$aov.tab[i,7], method = "BH", n = 3)} #correct for multiple tests, this time for 3 views
write.csv(pMAN$aov.table[,c(1,4,5,7)],paste("sd_procMANOVA_dimset_",view,".csv",sep=""))

# SD: CVA -----
# # Use CVAGen8 for CVA-based assignments tests with jackknife
# # write tables for known sex
# write.table(file=paste("cva_sd_dim_",view,".txt",sep=""),cbind(dim_set,dim_metadata[,cs_metadata_col]),col.names=FALSE,row.names=FALSE,sep='\t')
# write.table(file=paste("cva_sd_dim_",view,"_grp.txt",sep=""),as.numeric(dim_metadata$sex),col.names=FALSE,row.names=FALSE)
# # write group file for ssp dataset
# write.table(file=paste("cva_sd_ssp_",view,"_grp.txt",sep=""),as.numeric(ssp_all_metadata$sex),col.names=FALSE,row.names=FALSE) #grouping to check traditional ssp

# SD: model shapes -----
#plotting mean male and female shapes

#calculate mean shape for males
male_mshp<-mshp(dim_set[which(dim_metadata$sex=="male"),],m=2)

#calculate confidence intervals for male landmarks
sd_male2<-dim_set[which(dim_metadata$sex=="male"),] %>% 
  apply(.,2,function(x) sd(x)*1.96) %>% 
  matrix(.,ncol=2,byrow=TRUE)

#calculate mean shape for female
female_mshp<-mshp(dim_set[which(dim_metadata$sex=="female"),],m=2)

#calculate confidence intervals for female landmarks 
sd_female2<-dim_set[which(dim_metadata$sex=="female"),] %>% 
  apply(.,2,function(x) sd(x)*1.96) %>% 
  matrix(.,ncol=2,byrow=TRUE)

#calculate x and y limits
#plot the two mean species shape with standard deviation
# pdf(paste("sd_outline_",view,".pdf",sep=""),width = 4.4, height = 4.4,useDingbats=FALSE)
# par(mar=c(.5,.5,.5,.5))
# xlims<-c(min(c(female_mshp$meanshape[,1]-sd_female2[,1],male_mshp$meanshape[,1]-sd_male2[,1])),
#          max(c(female_mshp$meanshape[,1]+sd_female2[,1],male_mshp$meanshape[,1]+sd_male2[,1])))
# ylims<-c(min(c(female_mshp$meanshape[,2]-sd_female2[,2],male_mshp$meanshape[,2]-sd_male2[,2])),
#          max(c(female_mshp$meanshape[,2]+sd_female2[,2],male_mshp$meanshape[,2]+sd_male2[,2])))
# plot(male_mshp$meanshape,asp=1,pch=21,bg=spp_col[1],axes=FALSE,xlab="",ylab="",
#      xlim=xlims, ylim=ylims,cex=0.5,type="n")
# plot.errell(male_mshp$meanshape,sd_male2,color.bg=spp_col[1],alpha.bg=0.5)
# plot.errell(female_mshp$meanshape,sd_female2,color.bg=spp_col[2],alpha.bg=0.5)
# points(male_mshp$meanshape,asp=1,pch=21,bg=spp_col[1],cex=0.6)
# points(female_mshp$meanshape,asp=1,pch=22,bg=spp_col[2],cex=0.6)
# # title("SD, solid = allometrically corrected, red=female")
# dev.off()

# cairo_pdf("sd_shape_legend.pdf",width = 3.6, height = 3.6,family ="Arial") #make legend
# plot(male_mshp$meanshape,type="n",axes=F,xlab="",ylab="")
# legend('center',cex=2,pch=c(21,22),pt.cex=3,pt.bg=spp_col,legend=rev(levels(dim_metadata$sex)))
# dev.off()
# embed_fonts("sd_shape_legend.pdf")

# plot the way recommended by reviewers:
female_mshp_r<-rotateAMatrix(female_mshp$meanshape,r1,0,0)
male_mshp_r<-rotateAMatrix(male_mshp$meanshape,r1,0,0)

pdf(paste("sd_wireframe_male",view,".pdf",sep=""),width = 4.4, height = 4.4,useDingbats=FALSE)
par(mar=c(.01,.01,.01,.01))
plotRefToTarget(female_mshp_r,male_mshp_r,method="TPS",gridPars=gridPar(tar.pt.size=.5,grid.lwd=0.8))
points(male_mshp_r,pch=16)
dev.off()

pdf(paste("sd_wireframe_female",view,".pdf",sep=""),width = 4.4, height = 4.4,useDingbats=FALSE)
par(mar=c(.01,.01,.01,.01))
plotRefToTarget(male_mshp_r,female_mshp_r,method="TPS",gridPars=gridPar(tar.pt.size=.5,grid.lwd=0.8))
points(female_mshp_r,pch=15)
dev.off()