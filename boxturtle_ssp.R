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
pMAN<-procD.lm(coords ~ log(Csize) * (latitude + longitude) * ssp,
               iter = 99, data = ssp.gdf) #report R^2,F,p, df

pMAN<-procD.lm(coords ~ log(Csize) * ssp, iter = 099, data = ssp.gdf)


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

#just for fun
advanced.procD.lm(coords~log(Csize)+ssp,~log(Csize)+latitude+ssp,
                  groups=~ssp,slope=~log(Csize),iter=499,data=ssp.gdf) #doesn't improve anything

advanced.procD.lm(coords~log(Csize),~log(Csize)+binary,
                  groups=~binary,slope=~log(Csize),iter=999,data=ssp.gdf) #no longer significant 

advanced.procD.lm(coords~1,~ssp,
                  groups=~ssp,iter=999,data=ssp.gdf) 



#create base for pairwise analysis via for-loop
pm.report<-pMAN$aov.table[,c(1,4,5,7)]
for (i in 1:(ncol(ssp.pair))){
  ssp.subset<-which(ssp_metadata$ssp==ssp.pair[1,i]|
                      ssp_metadata$ssp==ssp.pair[2,i])
  ssps.gpa<-arrayspecs(ssp_set[ssp.subset,],p=ncol(ssp_set)/2,k=2) %>% gpagen
  ssps.gdf<-geomorph.data.frame(ssps.gpa,ssp=ssp_metadata$ssp[ssp.subset],
                                latitude=ssp_metadata$latitude[ssp.subset],
                                longitude=ssp_metadata$longitude[ssp.subset])
  ssps.gdf$Csize<-ssp_metadata[ssp.subset,cs_metadata_col]
  pMANs<-procD.lm(coords ~ log(Csize) * ssp,
                  iter = 999, data = ssps.gdf)
  pm.report<-cbind(pm.report,pMANs$aov.table[,c(1,4,5,7)])
}
for(i in 1:nrow(pm.report)){pm.report[i,c(4,8,12,16,20,24,28)]<-p.adjust(pm.report[i,c(4,8,12,16,20,24,28)], method = "bonferroni", n = (6+1)*3)} #"BH

?p.adjust
pm.report<-rbind(rep(NA,ncol(pm.report)),pm.report) #make empty first row
pm.report[1,1]<-"all"
for (i in 1:(ncol(ssp.pair))){
  pm.report[1,(4*(i+1)-3)]<-paste(ssp.pair[1,i],"vs",ssp.pair[2,i],sep="_")
}
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

# SSP: 2BPLS?------------
#data must be in 3D array form, both sides
ssp_shape<-arrayspecs(ssp_set,p=m,k=2)
ssp_variable<-cbind(log(ssp_metadata[,cs_metadata_col]),ssp_metadata$ssp,
                    ssp_metadata$latitude,ssp_metadata$longitude) 

writeland.tps(ssp_shape,paste("ssp_",view,"_shape.TPS",sep=""))
write.table(ssp_variable,"ssp_variable_temp.txt",row.names=FALSE,col.names=FALSE)

ssp_shape<-arrayspecs(ssp_set[which(ssp_metadata$ssp!="bauri"),],p=m,k=2)
ssp_variable<-cbind(log(ssp_metadata[,cs_metadata_col]),ssp_metadata$ssp,
                    ssp_metadata$latitude,ssp_metadata$longitude) #%>%
  .[which(ssp_metadata$ssp!="bauri"),]

writeland.tps(ssp_shape,paste("ssp_",view,"_shape2.TPS",sep=""))
write.table(ssp_variable,"ssp_variable_temp2.txt",row.names=FALSE,col.names=FALSE)
dim(ssp_shape)
getwd()

ssp_variable2<-cbind(log(ssp_metadata[,cs_metadata_col]),ssp_metadata$ssp,
                    ssp_metadata$latitude,ssp_metadata$longitude) %>%
  array(.,dim=c(1,4,nrow(ssp_metadata)))
ssp_PLS<-two.b.pls(ssp_variable,ssp_shape) #works, but not sure it provides what I'm looking for
summary(ssp_PLS)
print(ssp_PLS)
ssp_PLS$random.r
#try Claude's function
#probably want to standardize variables first, as tpsPLS does
ssp_variable_scaled<-scale(ssp_variable,center=FALSE,scale=TRUE)
?scale
pls1<-pls(ssp_variable_scaled,PCAssp$x[,c(1:PCc)])
pls1$D/sum(pls1$D)


pls2<-pls(ssp_variable_scaled[which(ssp_metadata$ssp!="bauri"),],PCAssp$x[,c(1:PCc)][which(ssp_metadata$ssp!="bauri"),])
pls2$D/sum(pls2$D)

pls(ssp_variable_scaled[which(ssp_metadata$ssp!="carolina"),],PCAssp$x[,c(1:PCc)][which(ssp_metadata$ssp!="carolina"),])$F1
pls(ssp_variable_scaled[which(ssp_metadata$ssp!="triunguis"),],PCAssp$x[,c(1:PCc)][which(ssp_metadata$ssp!="triunguis"),])$F1
pls(ssp_variable_scaled[which(ssp_metadata$ssp!="major"),],PCAssp$x[,c(1:PCc)][which(ssp_metadata$ssp!="major"),])$F1



pls1<-pls(ssp_variable_scaled[which(ssp_metadata$ssp=="bauri"|ssp_metadata$ssp=="major"),],
          ssp_set[which(ssp_metadata$ssp=="bauri"|ssp_metadata$ssp=="major"),])

pls(ssp_variable_scaled[which(ssp_metadata$ssp=="bauri"|ssp_metadata$ssp=="carolina"),],
          ssp_set[which(ssp_metadata$ssp=="bauri"|ssp_metadata$ssp=="carolina"),])$F1

pls1<-pls(ssp_variable_scaled[which(ssp_metadata$ssp=="triunguis"|ssp_metadata$ssp=="major"),],
          ssp_set[which(ssp_metadata$ssp=="triunguis"|ssp_metadata$ssp=="major"),])
pls1<-pls(ssp_variable_scaled[which(ssp_metadata$ssp=="triunguis"|ssp_metadata$ssp=="carolina"),],
          ssp_set[which(ssp_metadata$ssp=="triunguis"|ssp_metadata$ssp=="carolina"),])

pls(ssp_variable_scaled[which(ssp_metadata$ssp=="bauri"|ssp_metadata$ssp=="triunguis"),],
    ssp_set[which(ssp_metadata$ssp=="bauri"|ssp_metadata$ssp=="triunguis"),])

colnames(iris)
####f3.4
pls<-function(M1, M2)
{p1<-dim(M1)[2]; p2<-dim(M2)[2]; n<-dim(M1)[1]
sM12<-svd(var(cbind(M1,M2))[1:p1, (p1+1):(p1+p2)])
vM12<-var(cbind(M1,M2))[1:p1, (p1+1):(p1+p2)]
vM21<-var(cbind(M1,M2))[(p1+1):(p1+p2), 1:p1]
v11<-var(M1)
v22<-var(M2)
D<-sM12$d; F1<-sM12$u; F2<-sM12$v
Rv<-sum(diag(vM12%*%vM21))/sqrt(sum(diag(v11%*%v11))*sum(diag(v22%*%v22)))
list(Rv=Rv, F1=F1, F2=F2, D=D)}


# SSP: Bayesian Clustering? -------
library(mclust)
class<-as.character(ssp_metadata$ssp)
X<-PCAssp$x[,1:PCc]
# clPairs(PCAssp$x[,1:PCc],class)
BIC<-mclustBIC(X,G=1:4)
plot(BIC)
summary(BIC)
mod1<-Mclust(X,x=BIC)
summary(mod1)
table(class, mod1$classification)
mod1dr<-MclustDR(mod1,normalized=FALSE)
summary(mod1dr)
# plot(mod1dr,what="pairs")
plot(mod1dr, what = "boundaries", ngrid = 200)

plot(ssp_metadata$longitude,ssp_metadata$latitude,
     pch=21,bg=ssp_col[mod1$classification],cex=2)

# SSP: Assignments -----
# # # Use CVAGen8 for CVA-based assignments tests with jackknife
# # # write tables for the two groups
# write.table(file=paste("cva_binary_",view,".txt",sep=""),cbind(ssp_all,ssp_all_metadata[,cs_metadata_col]),col.names=FALSE,row.names=FALSE,sep='\t')
# write.table(file=paste("cva_binary_grp_",view,".txt",sep=""),as.numeric(binary_all),col.names=FALSE,row.names=FALSE)
# write group ID's to compare to traditional ssp
# write.table(file=paste("cva_ssp_",view,"_grp.txt",sep=""),as.numeric(ssp_all_metadata$ssp),col.names=FALSE,row.names=FALSE)

classdat<-read.csv("CVA_SSP/cva_s_dor_detal.csv")
ncol(classdat)
head(classdat)
#GroupSymbol
cbind(as.character(ssp_all_metadata$ssp),classdat$OrdinalGroup)[c(1,10,16,67),]

x1<-table(classdat$GroupSymbol,classdat$AssignedGroup)

x1[1,1]/sum(x1[1,])
x1[2,2]/sum(x1[2,])
x1[3,3]/sum(x1[3,])
x1[4,4]/sum(x1[4,])

factor(classdat$Problevel3)
table(classdat$GroupSymbol[classdat$Problevel3=="p>0.05"],classdat$AssignedGroup[classdat$Problevel3=="p>0.05"])
x2<-table(classdat$GroupSymbol[classdat$Problevel1=="p>0.05"&
                             classdat$Problevel2=="p>0.05"&
                             classdat$Problevel3=="p>0.05"&
                             classdat$Problevel4=="p>0.05"],
      classdat$AssignedGroup[classdat$Problevel1=="p>0.05"&
                               classdat$Problevel2=="p>0.05"&
                               classdat$Problevel3=="p>0.05"&
                               classdat$Problevel4=="p>0.05"])

x2[1,1]/sum(x2[1,])
x2[2,2]/sum(x2[2,])
x2[3,3]/sum(x2[3,])
x2[4,4]/sum(x2[4,])
factor(classdat$Unique)
length(which(classdat$Unique=="Group=4"))

# SSP: bgPCA --------
bin_residuals_mean<-colMeans(ssp_set)#compute grand mean
bin_residuals_centered<-as.matrix(ssp_set)-rep(1,nrow(ssp_set))%*%t(bin_residuals_mean) #Subtract that grandmean from each row

#Calculate the group means
bin_means<-array(NA,dim=c(length(levels(binary)),ncol(ssp_set)))
for (i in 1:length(levels(binary))){
  bin_means[i,]<-colMeans(ssp_set[which(binary==levels(binary)[i]),])
}
B<-prcomp(bin_means)$rotation #eigenvectors
B2<-prcomp(bin_means)

bgPCAbin<-bin_residuals_centered%*%B #Get the scores for all the individuals on the eigenvectors of the PCA of the means

plot(bgPCAbin[,c(1,2)],cex=1.5,
     bg=ssp_col[binary],pch=21,xlab="",ylab="")

# SAMdata --------
write1<-lapply(ssp_metadata,as.numeric) %>% as.data.frame
write2<-cbind(write1,PCAssp$x[,c(1:PCc)])
write2$major_growth_rings[which(is.na(write2$major_growth_rings))]<-9999
colnames(write2)[which(colnames(write2)=="latitude")]<-"Latitude"
colnames(write2)[which(colnames(write2)=="longitude")]<-"Longitude"
write.csv(write2,paste("boxturtle_SAMdata_ssp_",view,".csv",sep=""),row.names=FALSE)
getwd()


