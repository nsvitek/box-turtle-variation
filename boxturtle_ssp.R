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

# SSP: stats4ssp -----
#make pairwise combinations of subspecies
ssp.pair<-combn(levels(ssp_metadata$ssp),2)

# # nonparametrics(permutational) MANOVA with repeatible PCs
nMAN<-adonis(PCAssp$x[,c(1:PCc)]~log(ssp_metadata[,cs_metadata_col])*
         ssp_metadata$latitude*
         ssp_metadata$longitude*
         ssp_metadata$ssp, permutations=1000, method="euclidean")

#create base for pairwise analysis via for-loop
np.report<-nMAN$aov.tab[,c(1,4,5,6)]
for (i in 1:(ncol(ssp.pair))){
  ssp.subset<-which(ssp_metadata$ssp==ssp.pair[1,i]|
                      ssp_metadata$ssp==ssp.pair[2,i])
  nMANs<-adonis(PCAssp$x[ssp.subset,c(1:PCc)]~log(ssp_metadata[ssp.subset,cs_metadata_col])*
                  ssp_metadata$latitude[ssp.subset]*
                  ssp_metadata$longitude[ssp.subset]*
                  ssp_metadata$ssp[ssp.subset], permutations=1000, method="euclidean")
  np.report<-cbind(np.report,nMANs$aov.tab[,c(1,4,5,6)])
}

#correct for multiple tests, this time for 3 views for each pairwise + total data test [(6+1)*3]
for(i in 1:nrow(np.report)){np.report[i,c(4,8,12,16,20,24,28)]<-p.adjust(np.report[i,c(4,8,12,16,20,24,28)], method = "BH", n = (6+1)*3)} 

np.report<-rbind(rep(NA,ncol(np.report)),np.report) #make empty first row
np.report[1,1]<-"all"
for (i in 1:(ncol(ssp.pair))){
  np.report[1,(4*(i)-3)]<-paste(ssp.pair[1,i],"vs",ssp.pair[2,i],sep="_")
}
write.csv(np.report,paste("ssp_npMANOVA_pairwise_",view,".csv",sep=""))

# # Procrustes MANOVA for Goodall's F results with coordinate data
pMAN<-procD.lm(coords ~ log(Csize) * latitude * longitude * ssp,
               iter = 999, data = ssp.gdf) #report R^2,F,p, df

#create base for pairwise analysis via for-loop
pm.report<-pMAN$aov.table[,c(1,4,5,7)]
for (i in 1:(ncol(ssp.pair))){
  ssp.subset<-which(ssp_metadata$ssp==ssp.pair[1,i]|
                      ssp_metadata$ssp==ssp.pair[2,i])
  ssps.gpa<-arrayspecs(ssp_set[ssp.subset,],p=ncol(ssp_set)/2,k=2) %>% gpagen
  ssps.gdf<-geomorph.data.frame(ssps.gpa,binary=binary[ssp.subset],ssp=ssp_metadata$ssp[ssp.subset],
                                latitude=ssp_metadata$latitude[ssp.subset],
                                longitude=ssp_metadata$longitude[ssp.subset])
  ssps.gdf$Csize<-ssp_metadata[ssp.subset,cs_metadata_col]
  pMANs<-procD.lm(coords ~ log(Csize) * latitude * longitude * ssp,
                  iter = 999, data = ssps.gdf)
  pm.report<-cbind(pm.report,pMANs$aov.table[,c(1,4,5,7)])
}
for(i in 1:nrow(pm.report)){pm.report[i,c(4,8,12,16,20,24,28)]<-p.adjust(pm.report[i,c(4,8,12,16,20,24,28)], method = "BH", n = (6+1)*3)} 

pm.report<-rbind(rep(NA,ncol(pm.report)),pm.report) #make empty first row
pm.report[1,1]<-"all"
for (i in 1:(ncol(ssp.pair))){
  pm.report[1,(4*(i+1)-3)]<-paste(ssp.pair[1,i],"vs",ssp.pair[2,i],sep="_")
}
write.csv(pm.report,paste("ssp_procMANOVA_pairwise_",view,".csv",sep=""))

# SSP: Monmonier ------
#Monmonier analyses

xy1<-cbind(ssp_metadata$longitude,ssp_metadata$latitude)
d<-dist(PCAssp$x[,c(1:PCc)])
cn<-chooseCN(jitter(xy1),type=2,ask=FALSE) #type 1 connects FL peninsula to TX
# mon<-optimize.monmonier(xy1,d,cn,ntry=10) #look for threshold, preliminary results
mon<-optimize.monmonier(xy1,d,cn,ntry=10,threshold=thresh,return.best=FALSE)
nrun<-3
mon2<-monmonier(xy1, d, cn, threshold=thresh, bd.length=NULL, nrun=3,
          skip.local.diff=rep(mon,nrun), allowLoop=TRUE)

cairo_pdf(paste("ssp_Monmonier_",view,".pdf",sep=""),width = 5.4,height=3)
plot(mon2, add.arrow=FALSE,bwd=10,col=ssp_col[8])
points(mon2$xy,bg=ssp_col[ssp_metadata$ssp],pch=(20+as.numeric(ssp_metadata$ssp)))
dev.off()

# binary MANCOVA ----
#try again with post-Monmonier binary coding
ssp_metadata$binary<-ssp_metadata$ssp=="bauri"
nMAN<-adonis(PCAssp$x[,c(1:PCc)]~log(ssp_metadata[,cs_metadata_col])*
               ssp_metadata$latitude*
               ssp_metadata$longitude*
               ssp_metadata$binary, permutations=1000, method="euclidean")

#correct for multiple tests, this time for 3 views 
np.report<-nMAN$aov.tab[,c(1,4,5,6)]
for(i in 1:nrow(np.report)){np.report[i,4]<-p.adjust(np.report[i,4], method = "BH", n = 3)} 

write.csv(np.report,paste("binary_npMANOVA_pairwise_",view,".csv",sep=""))


pMAN<-procD.lm(coords ~ log(Csize) * latitude * longitude * binary,
               iter = 999, data = ssp.gdf) #report R^2,F,p, df
pm.report<-pMAN$aov.table[,c(1,4,5,7)]
for(i in 1:nrow(pm.report)){pm.report[i,4]<-p.adjust(np.report[i,4], method = "BH", n = 3)} 

write.csv(pm.report,paste("binary_procMANOVA_pairwise_",view,".csv",sep=""))

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


cairo_pdf("ssp_size.pdf",width = 5.4, height = 4,family ="Arial")
ggplot(data=ssp_all_metadata,aes(x=ssp,alpha=0.5)) +
  geom_violin(binaxis="y",binwidth=0.02,stackdir="center",
               aes(y=log(carapace_length),fill=spp_col[4])) +
  geom_violin(binaxis="y",binwidth=0.02,stackdir="up",
               aes(y=log(dor_cs),fill=spp_col[1])) +
  geom_violin(binaxis="y",binwidth=0.02,stackdir="center",
               aes(y=log(lat_cs),fill=spp_col[2])) +
  geom_violin(binaxis="y",binwidth=0.02,stackdir="down",
               aes(y=log(pos_cs),fill=spp_col[3])) +
  scale_fill_manual(name="Measurement", labels=c("carapace length","dorsal CS",
                                                 "lateral CS",
                                                 "posterior CS"),
                    values=spp_col[c(1,2,3,4)]) +
  xlab("Subspecies") + 
  ylab("ln(Carapace Length)") +
  theme_classic()
dev.off()  
embed_fonts("ssp_size.pdf")

# SSP: Assignments -----
# # # Use CVAGen8 for CVA-based assignments tests with jackknife
# # # write tables for the two groups
# write.table(file=paste("cva_binary_",view,".txt",sep=""),cbind(ssp_all,ssp_all_metadata[,cs_metadata_col]),col.names=FALSE,row.names=FALSE,sep='\t')
# write.table(file=paste("cva_binary_grp_",view,".txt",sep=""),as.numeric(binary_all),col.names=FALSE,row.names=FALSE)
# write group ID's to compare to traditional ssp
# write.table(file=paste("cva_ssp_",view,"_grp.txt",sep=""),as.numeric(ssp_all_metadata$ssp),col.names=FALSE,row.names=FALSE)