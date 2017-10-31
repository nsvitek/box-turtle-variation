# A: common slopes ------
#See Ch.11 of Zelditch et al workbook
# test for significant interaction a duplicate of sexual dimorphism model. Skipped.
# cs.MAN<-adonis(ssp_all~log(ssp_all_metadata[,cs_metadata_col])*binary_all,
#                permutations=1000,method="euclidean")$aov.tab[,c(1,4,5,6)]
# cs.MAN[,4]<-p.adjust(cs.MAN[,4], method = "BH", n = 3) #correct for multiple tests, this time for 3 views
# 
# 
# write.csv(cs.MAN,paste("binary_npMANOVA_cslopes_",view,".csv",sep=""))
# cs.pMAN<-procD.lm(coords ~ log(Csize) * binary,
#                 iter = 999, data = ssp.gdf)$aov.table[,c(1,4,5,7)]
# write.csv(cs.pMAN,paste("binary_procMANOVA_cslopes_",view,".csv",sep=""))

cslopes<-common.slope.test(log(ssp_all_metadata[,cs_metadata_col]),
                           ssp_all,binary_all,nperm=1000)
#Are both of the above are evidence that there is a different ontogenetic trajectory between the two?
ssp_cslopes<-advanced.procD.lm(coords ~ log(Csize) + binary, 
                               ~ log(Csize) * binary, 
                               groups = ~ binary, 
                               slope = ~ log(Csize), angle.type = "deg", iter = 999, data = ssp.gdf)

cs.report<-data.frame(cslopes.p=cslopes$p.value[1,2],
                      cslopes.diff=cslopes$obs.diff[1,2],
                      slope.angles=ssp_cslopes$slopes.angles[1,2],
                      Z.angles=ssp_cslopes$Z.angles[1,2],
                      P.angles=ssp_cslopes$P.angles[1,2]
                      ) 
cs.report[c(1,5)]<-p.adjust(cs.report[c(1,5)], method = "BH", n = 3) #correct for multiple tests, this time for 3 views
write.csv(cs.report,paste("binary_commonslopes_",view,".csv",sep=""))

# A: binary models ----
# build allometric linear models for the two spatial groups
modm<-modb<-sd_model1<-sd_model2<-array(NA,dim=c(m,k,length(levels(binary_all))))
residuals<-list()
for (i in 1:length(levels(binary_all))){
  subset<-which(binary_all==levels(binary_all)[i])
  model<-lm(as.matrix(ssp_all[subset,])~log(ssp_all_metadata[subset,cs_metadata_col])) #allometric model with CS
  # model<-lm(as.matrix(ssp_all[subset,])~log(ssp_all_metadata$carapace_length[subset])) #allometric model with CL
  modm[,,i]<-matrix(model$coefficients[2,],m,k,byrow=TRUE) #m in the linear model y=mx+b
  residuals[[i]]<-model$residuals
  modb[,,i]<-matrix(model$coefficients[1,],m,k,byrow=TRUE) #b in the linear model y=mx+b
  sd_model1[,,i]<-model$residuals %>%  # for confidence interval visualizations these may need to change, but note that it was also done in Popat et al. 2013
    apply(.,2,function(x) sd(x)*1) %>% matrix(.,ncol=2,byrow=TRUE) 
  sd_model2[,,i]<-model$residuals %>% 
    apply(.,2,function(x) sd(x)*1.96) %>% matrix(.,ncol=2,byrow=TRUE)
}
testb<-function(x){modm[,,1]*x+modb[,,1]} #make bauri linear model
testn<-function(x){modm[,,2]*x+modb[,,2]} #make 'all else' linear model

#find smallest, largest, and mean centroid size
standards<-ssp_all_metadata[,cs_metadata_col]  %>% summary #want 1st, 3th, 6th values
# standards<-ssp_all_metadata$carapace_length  %>% summary #want 1st, 3th, 6th values

cairo_pdf(paste("binary_shape_",view,".pdf",sep=""),width = 7.6, height = 2.6)
# postscript(paste("binary_shape_",view,".eps",sep=""),width = 5.4, height = 1.8)
par(mfrow=c(1,3))
par(mar=c(.5,.5,.5,.5))
for (z in c(1,4,6)){
  temp<-array(data=c(testb(log(standards[z])),testn(log(standards[z]))),dim=c(m,k,2))
  # pdf(paste("binary_shape_",view,z,".pdf",sep=""))
  # create x and y lims using sdmodel+temp?
  xlims<-c(min(c(temp[,1,1]-sd_model2[,1,1],temp[,1,2]-sd_model2[,1,2])),
  max(c(temp[,1,1]+sd_model2[,1,1],temp[,1,2]+sd_model2[,1,2])))
  ylims<-c(min(c(temp[,2,1]-sd_model2[,2,1],temp[,2,2]-sd_model2[,2,2])),
           max(c(temp[,2,1]+sd_model2[,2,1],temp[,2,2]+sd_model2[,2,2])))
  plot(temp[,,1],asp=1,axes=F,xlab="",ylab="",type="n",xlim=xlims,ylim=ylims)
  for(j in c(1,2)){plot.errell(temp[,,j],sd_model2[,,j],color.bg=spp_col[j],alpha.bg=0.5)}
  # for(j in c(1,2)){lines(temp[conx,,j],col="grey",lty=1+j)}
  # arrows(temp[,1,2],temp[,2,2],temp[,1,1],temp[,2,1],length=0.1,lwd=2,col="red")
  for(j in c(1,2)){points(temp[,,j],pch=20+j,bg=spp_col[j],cex=.6)}
  # dev.off()
}
dev.off()

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
nMAN<-adonis(PCAdim$x[,c(1:PCc)]~log(dim_metadata[,cs_metadata_col])*
               # binary_dim*
               dim_metadata$sex, permutations=1000, method="euclidean")
for(i in 1:nrow(nMAN$aov.tab)){nMAN$aov.tab[i,6]<-p.adjust(nMAN$aov.tab[i,6], method = "BH", n = 3)} #correct for multiple tests, this time for 3 views

write.csv(nMAN$aov.tab[,c(1,4,5,6)],paste("sd_npMANOVA_dimset_",view,".csv",sep=""))

pMAN<-procD.lm(coords ~ log(Csize)  * sex, #* binary
               iter = 999, data = dim.gdf) #report R^2,F,p, df
for(i in 1:nrow(pMAN$aov.tab)){pMAN$aov.tab[i,7]<-p.adjust(pMAN$aov.tab[i,7], method = "BH", n = 3)} #correct for multiple tests, this time for 3 views
write.csv(pMAN$aov.table[,c(1,4,5,7)],paste("sd_procMANOVA_dimset_",view,".csv",sep=""))

#statistics: again with larger dataset
nMAN<-adonis(PCAssp$x[,c(1:PCc)]~log(ssp_metadata[,cs_metadata_col])*
               # binary*
               ssp_metadata$sex, permutations=1000, method="euclidean")
for(i in 1:nrow(nMAN$aov.tab)){nMAN$aov.tab[i,6]<-p.adjust(nMAN$aov.tab[i,6], method = "BH", n = 3)} #correct for multiple tests, this time for 3 views

write.csv(nMAN$aov.tab[,c(1,4,5,6)],paste("sd_npMANOVA_sspset_",view,".csv",sep=""))

pMAN<-procD.lm(coords ~ log(Csize)  * sex, #* binary
               iter = 999, data = ssp.gdf) #report R^2,F,p, df
for(i in 1:nrow(pMAN$aov.tab)){pMAN$aov.tab[i,7]<-p.adjust(pMAN$aov.tab[i,7], method = "BH", n = 3)} #correct for multiple tests, this time for 3 views

write.csv(pMAN$aov.table[,c(1,4,5,7)],paste("sd_procMANOVA_sspset_",view,".csv",sep=""))

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
pdf(paste("sd_outline_",view,".pdf",sep=""),width = 4.4, height = 4.4,useDingbats=FALSE)
par(mar=c(.5,.5,.5,.5))
xlims<-c(min(c(female_mshp$meanshape[,1]-sd_female2[,1],male_mshp$meanshape[,1]-sd_male2[,1])),
         max(c(female_mshp$meanshape[,1]+sd_female2[,1],male_mshp$meanshape[,1]+sd_male2[,1])))
ylims<-c(min(c(female_mshp$meanshape[,2]-sd_female2[,2],male_mshp$meanshape[,2]-sd_male2[,2])),
         max(c(female_mshp$meanshape[,2]+sd_female2[,2],male_mshp$meanshape[,2]+sd_male2[,2])))
plot(male_mshp$meanshape,asp=1,pch=21,bg=spp_col[1],axes=FALSE,xlab="",ylab="",
     xlim=xlims, ylim=ylims,cex=0.5,type="n")
plot.errell(male_mshp$meanshape,sd_male2,color.bg=spp_col[1],alpha.bg=0.5)
plot.errell(female_mshp$meanshape,sd_female2,color.bg=spp_col[2],alpha.bg=0.5)
points(male_mshp$meanshape,asp=1,pch=21,bg=spp_col[1],cex=0.6)
points(female_mshp$meanshape,asp=1,pch=22,bg=spp_col[2],cex=0.6)
# title("SD, solid = allometrically corrected, red=female")
dev.off()

cairo_pdf("sd_shape_legend.pdf",width = 3.6, height = 3.6,family ="Arial") #make legend
plot(male_mshp$meanshape,type="n",axes=F,xlab="",ylab="")
legend('center',cex=2,pch=c(21,22),pt.cex=3,pt.bg=spp_col,legend=rev(levels(dim_metadata$sex)))
dev.off()
embed_fonts("sd_shape_legend.pdf")