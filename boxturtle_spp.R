# Species: load -----------------------------------------------------------------
##### compare T. coahuila and T. carolina
##### pull specimens meant for species analyses
#set filter criteria
filterspp<-which(data_uni$spp=="coahuila"|data_uni$ssp_set=="yes")

#make dataset and metadataset, remove specimens missing for that view
if(length(which(is.na(shape[filterspp,1])))>=1){
  spp_metadata<-data_uni[filterspp,metadata_columns] %>% 
    .[-which(is.na(shape[filterspp,1])),]
} else {
  spp_metadata<-data_uni[filterspp,metadata_columns]
}
spp_set<-shape[filterspp,] %>% na.omit()

dim(spp_set)
dim(spp_metadata)

# SPP: explore -----
##### 2.1 Principal Components Analysis (look at data)
pdf(paste(view,"_spp_PCA.pdf",sep=""))
PCAspp<-prcomp(spp_set,scale.=FALSE)
an<-anderson(PCAspp$sdev)
plot(PCAspp$x[,1:2],cex=spp_metadata[,cs_metadata_col]/(min(spp_metadata[,cs_metadata_col])),
     bg=spp_col[spp_metadata$spp],pch=(20+as.numeric(spp_metadata$spp)),
     xlab=paste("PC 1 (",round(an$percent[1],3),"%)",sep=""),
     ylab=paste("PC 2 (",round(an$percent[2],3),"%)",sep=""))
legend('bottomleft',cex=1,pch=c(21,22),pt.bg=spp_col,legend=levels(spp_metadata$spp))
dev.off()

##### 2.1.1 CVA & Assignments tests
# #give up, doing it in CVAGen for mac. based on enough loop, use first 12 pc's for dor, lat, 7 for pos
# for (i in 1:length(an$percent)){
#   enough<-sum(an$percent[1:i])
#   if(enough>=95){
#     enough<-i
#     break
#   }
# }
# write.table(file=paste("cva_",view,"_spp.txt",sep=""),cbind(spp_set,spp_metadata[,cs_metadata_col]),col.names=FALSE,row.names=FALSE,sep='\t')
# write.table(file=paste("cva_",view,"_sppgrp.txt",sep=""),as.numeric(spp_metadata$spp),col.names=FALSE,row.names=FALSE)
##### 2.2 are the shapes significantly different? Do they look different?
test_sppshdif<-adonis(PCAspp$x~spp_metadata$spp,permutations=10000,method="euclidean")

#calculate mean shape for carolina
carolina_mshp<-mshp(spp_set[which(spp_metadata$spp=="carolina"),],m=2)

#calculate confidence intervals for carolina landmarks
sd_carolina<-spp_set[which(spp_metadata$spp=="carolina"),] %>% 
  apply(.,2,sd) %>% 
  matrix(.,ncol=2,byrow=TRUE)

#calculate mean shape for coahuila
coahuila_mshp<-mshp(spp_set[which(spp_metadata$spp=="coahuila"),],m=2)

#for future illustration purposes, specimen that is closest to the mean coahuila shape
stand<-colMeans(spp_set[which(spp_metadata$spp=="coahuila"),])
ard<-apply(spp_set,1,function(x) abs(x-stand)) %>% apply(.,2,sum)
standard<-spp_metadata$id[which(ard==min(ard))] #UF153944 for dorsal, posterior

#calculate confidence intervals for coahuila landmarks 
sd_coahuila<-spp_set[which(spp_metadata$spp=="coahuila"),] %>% 
  apply(.,2,sd) %>% 
  matrix(.,ncol=2,byrow=TRUE)

#calculate x and y limits
xlims<-range(c(carolina_mshp$meanshape[,1],coahuila_mshp$meanshape[,1])) %>% +c(-0.015,0.015)
ylims<-range(c(carolina_mshp$meanshape[,2],coahuila_mshp$meanshape[,2])) %>% +c(-0.015,0.015)

#plot the two mean species shape with standard deviation
pdf(paste(view,"_spp_outline.pdf",sep=""))
plot(carolina_mshp$meanshape,asp=1,pch=21,bg=spp_col[1],axes=FALSE,xlab="",ylab="",
     xlim=xlims, ylim=ylims,cex=0.5)
rgbcol<-rgb(t(col2rgb(spp_col[1])/255),alpha=.5)
rgbbor<-rgb(t(col2rgb(spp_col[1])/255),alpha=.1)
eb_car<-errbubble(x=carolina_mshp$meanshape[,1],y=carolina_mshp$meanshape[,2],sd_carolina)
rgbcol<-rgb(t(col2rgb(spp_col[2])/255),alpha=.5)
rgbbor<-rgb(t(col2rgb(spp_col[2])/255),alpha=.1)
eb_coa<-errbubble(x=coahuila_mshp$meanshape[,1],y=coahuila_mshp$meanshape[,2],sd_coahuila)
lines(carolina_mshp$meanshape[conx,],col=spp_col[1],lwd=2)
lines(coahuila_mshp$meanshape[conx,],col=spp_col[2],lwd=2)
points(carolina_mshp$meanshape,asp=1,pch=21,bg=spp_col[1],cex=0.5)
points(coahuila_mshp$meanshape,asp=1,pch=21,bg=spp_col[2],cex=0.5)
dev.off()

##### 4.7 can the [small] fossils be diagnosed as T. coahuila?

#make dataset and metadataset, remove specimens missing for that view
# q_metadata<-rbind(spp_metadata,fos_metadata)
# q_set<-rbind(spp_set,fos_set)

classify_metadata<-rbind(spp_metadata,smfos_metadata) %>% droplevels()
classify_set<-rbind(spp_set,smfos_set)

#Principal Components Analysis (look at data)
PCAclass<-prcomp(classify_set,scale.=FALSE)
an<-anderson(PCAclass$sdev)
plot(PCAclass$x[,1:2],cex=classify_metadata[,cs_metadata_col]/(min(classify_metadata[,cs_metadata_col])),
     bg=fos_col[classify_metadata$site],pch=(20+as.numeric(classify_metadata$spp)),
     xlab=paste("PC 1 (",round(an$percent[1],3),"%)",sep=""),
     ylab=paste("PC 2 (",round(an$percent[2],3),"%)",sep=""))
legend('bottomleft',cex=1,pch=c(21),pt.bg=fos_col,legend=levels(classify_metadata$site))

#look at size distribution
plot(classify_metadata[,cs_metadata_col],bg=fos_col[q_metadata$site],pch=(20+as.numeric(q_metadata$spp)),cex=1.5)
legend('bottomleft',cex=1,pch=c(21),pt.bg=fos_col,legend=levels(fos_metadata$site),pt.cex=2)
# SPP + F: Assigmnents -------------------------------------------------------------


# read in results of assignments tests, filter dataset, rinse, repeat. 
assnteststat<-read.csv(paste("cva_",view,"_smfospp_assnteststat.csv",sep=""),header=TRUE)
classify_meta<-cbind(classify_metadata,assnteststat)
# colnames(classify_meta)
# classify_meta[which(classify_meta$AssignedGroup==2),c(1,9,15,17,18,25)]

ssp_metadata2<-classify_metadata[which(classify_meta$Assigned==1 & classify_meta$ssp_set=="yes"),]
classify_metadata2<-rbind(ssp_metadata2,spp_metadata[which(spp_metadata$spp=="coahuila"),])
classify_set2<-shape[which(data_uni$id %in% classify_metadata2$id),]
dim(classify_metadata2)
dim(classify_set2)

write.table(file=paste("cva_",view,"_spp2.txt",sep=""),cbind(classify_set2,classify_metadata2[,cs_metadata_col]),col.names=FALSE,row.names=FALSE,sep='\t')
write.table(file=paste("cva_",view,"_spp2grp.txt",sep=""),as.numeric(classify_metadata2$spp),col.names=FALSE,row.names=FALSE)

assntest2stat<-read.csv(paste("cva_",view,"_smfospp2_assnteststat.csv",sep=""),header=TRUE)
classify_meta2<-rbind(classify_metadata2,smfos_metadata) %>% cbind(.,assntest2stat)

classify_meta2[which(classify_meta2$Assigned==2),c(1,9,15,17,18)]

ssp_metadata3<-classify_metadata2[which(classify_meta2$Assigned==1 & classify_meta2$ssp_set=="yes"),]
classify_metadata3<-rbind(ssp_metadata3,spp_metadata[which(spp_metadata$spp=="coahuila"),]) %>% droplevels()
classify_set3<-shape[which(data_uni$id %in% classify_metadata3$id),]
write.table(file=paste("cva_",view,"_spp3.txt",sep=""),cbind(classify_set3,classify_metadata3[,cs_metadata_col]),col.names=FALSE,row.names=FALSE,sep='\t')
write.table(file=paste("cva_",view,"_spp3grp.txt",sep=""),as.numeric(classify_metadata3$spp),col.names=FALSE,row.names=FALSE)

assntest3stat<-read.csv(paste("cva_",view,"_smfospp3_assnteststat.csv",sep=""),header=TRUE)
coords<-read.csv("cva_dor_smfospp3_cvascores.csv", header=TRUE)
classify_meta3<-rbind(classify_metadata3,smfos_metadata) %>% 
  cbind(.,assntest3stat,coords) %>% droplevels()

PCAclass3<-rbind(classify_set3,smfos_set) %>% prcomp(.,scale.=FALSE)
an<-anderson(PCAclass3$sdev)
#loop to figure out how many PCs explain 95% of shape variation
for (i in 1:length(an$percent)){
  enough<-sum(an$percent[1:i])
  if(enough>=95){
    enough<-i
    break
  }
}

#plot PCA
plot(PCAclass3$x[,1:2],cex=classify_meta3[,cs_metadata_col]/(min(classify_meta3[,cs_metadata_col])),
     bg=fos_col[classify_meta3$site],pch=(20+as.numeric(classify_meta3$spp)),
     xlab=paste("PC 1 (",round(an$percent[1],3),"%)",sep=""),
     ylab=paste("PC 2 (",round(an$percent[2],3),"%)",sep=""),
     main="PCA. points scaled by centroid size\nsquares = T. coahuila")
legend('bottomleft',cex=0.7,pch=c(21),pt.bg=fos_col,legend=levels(classify_meta3$site),pt.cex=1.5)

#plot CVA
plot(classify_meta3$CV1,classify_meta3$CV2,
     cex=classify_meta3[,cs_metadata_col]/(min(classify_meta3[,cs_metadata_col])),
     bg=fos_col[classify_meta3$site],pch=(20+as.numeric(classify_meta3$spp)),
     main="CVA. points scaled by centroid size\nsquares = T. coahuila")
legend('bottomleft',cex=0.7,pch=c(21),pt.bg=fos_col,legend=levels(classify_meta3$site),pt.cex=1.5)
# Small Fossil Shapes -----------------------------------------------------
#make dataset for small, odd fossils
smfos_metadata<-fos_metadata[which(fos_metadata$carapace_length<=150 &
                                     (fos_metadata$site=="camelot"|
                                        fos_metadata$site=="vero"|
                                        fos_metadata$site=="fort center"|
                                        fos_metadata$site=="grove's orange midden")),]
smfos_set<-fos_set[which(fos_metadata$id %in% smfos_metadata$id),]
dim(smfos_set)
dim(smfos_metadata)

cbind(as.character(smfos_metadata$id),as.character(smfos_metadata$site))
# # write "unknown" file for later CVA/assignments in CVAGen
# write.table(file=paste("cva_",view,"_fos.txt",sep=""),cbind(fos_set,fos_metadata[,cs_metadata_col]),col.names=FALSE,row.names=FALSE,sep='\t')
# write.table(file=paste("cva_",view,"_fosgrp.txt",sep=""),rep(7,nrow(fos_set),col.names=FALSE,row.names=FALSE),sep="\t")
# write.table(file=paste("cva_",view,"_smfos.txt",sep=""),cbind(smfos_set,smfos_metadata[,cs_metadata_col]),col.names=FALSE,row.names=FALSE,sep='\t')
# write.table(file=paste("cva_",view,"_smfosgrp.txt",sep=""),rep(7,nrow(smfos_set)),col.names=FALSE,row.names=FALSE,sep="\t")



#calculate mean shape. 
smfos_mshp<-mshp(smfos_set,m=2)

#calculate confidence intervals for small fossils landmarks
sd_smfos<-smfos_set %>% 
  apply(.,2,sd) %>% 
  matrix(.,ncol=2,byrow=TRUE)

#plot outline in comparison to two species
pdf(file=paste(view,"_spp_smfos_outline.pdf",sep=""))#,width=4,height=4,res=500,units="in")
plot(carolina_mshp$meanshape,asp=1,pch=21,bg=spp_col[1],axes=FALSE,xlab="",ylab="",
     xlim=xlims, ylim=ylims,cex=0.5)
rgbcol<-rgb(t(col2rgb(spp_col[1])/255),alpha=.5)
rgbbor<-rgb(t(col2rgb(spp_col[1])/255),alpha=.1)
eb_car<-errbubble(x=carolina_mshp$meanshape[,1],y=carolina_mshp$meanshape[,2],sd_carolina)
rgbcol<-rgb(t(col2rgb(spp_col[2])/255),alpha=.5)
rgbbor<-rgb(t(col2rgb(spp_col[2])/255),alpha=.1)
eb_coa<-errbubble(x=coahuila_mshp$meanshape[,1],y=coahuila_mshp$meanshape[,2],sd_coahuila)
rgbcol<-rgb(t(col2rgb(spp_col[3])/255),alpha=.5)
rgbbor<-rgb(t(col2rgb(spp_col[3])/255),alpha=.1)
eb_smfos<-errbubble(x=smfos_mshp$meanshape[,1],y=smfos_mshp$meanshape[,2],sd_smfos)
lines(carolina_mshp$meanshape[conx,],col=spp_col[1],lwd=2)
lines(coahuila_mshp$meanshape[conx,],col=spp_col[2],lwd=2)
lines(smfos_mshp$meanshape[conx,],col=spp_col[3],lwd=2)
points(carolina_mshp$meanshape,asp=1,pch=21,bg=spp_col[1],cex=0.5)
points(coahuila_mshp$meanshape,asp=1,pch=21,bg=spp_col[2],cex=0.5)
points(smfos_mshp$meanshape,asp=1,pch=21,bg=spp_col[3],cex=0.5)
legend('bottomright',cex=0.7,pch=21,pt.bg=spp_col,legend=c(levels(spp_metadata$spp),"small fossils"))
dev.off()


