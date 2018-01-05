# Analyzes linear measures of box turtles. 
#Doesn't depend on view, so only run once. 

# size-rings-shape -----
# Assess relationship stability in terms of size, growth rings, and shape. 
# Create log-transformed length variable
age_measures$lnCL<-log(age_measures$CL) #take natural log of carapace length

# run ANOVA of carapace length (and growth rings) vs. shell ossification, 
aov.lnCL<-aov(age_measures$lnCL~age_measures$Ossification_score) %>% summary 
aov.MGR<-aov(age_measures$MGR~age_measures$Ossification_score) %>% summary
#report anova: F,p, df
aov.report<-matrix(c(aov.lnCL[[1]]$Df[1],aov.lnCL[[1]]$`F value`[1],aov.lnCL[[1]]$`Pr(>F)`[1],
                     aov.MGR[[1]]$Df[1],aov.MGR[[1]]$`F value`[1],aov.MGR[[1]]$`Pr(>F)`[1]),
  nrow=2,ncol=3,byrow=TRUE)
aov.report[,2]<-round(aov.report[,2],2)
row.names(aov.report)<-c("ln(Carapace Length)","Major Growth Rings")
colnames(aov.report)<-c("df","F","p")
write.csv(aov.report,"ossification-vs-CL-vs-MGR_ANOVA.csv")

#then t tests to look closer
#make subsets of ossification score pairs
os.pair<-combn(1:4,2)
#make table to report t statistics
t.report<-matrix(NA,nrow=ncol(os.pair)+1,ncol=6) %>% as.data.frame
colnames(t.report)<-rep(c("df","t","p <"),2)
for (i in 1:(ncol(t.report))){
  row.names(t.report)[i+1]<-paste("Score",os.pair[1,i],"vs",os.pair[2,i],sep="")
  t.subset<-which(age_measures$Ossification_score==os.pair[1,i]|
                    age_measures$Ossification_score==os.pair[2,i])
  t.result.lnCL<-t.test(age_measures$lnCL[t.subset]~age_measures$Ossification_score[t.subset])
  t.result.MGR<-t.test(age_measures$MGR[t.subset]~age_measures$Ossification_score[t.subset])
  t.report[i+1,]<-c(t.result.lnCL$parameter,t.result.lnCL$statistic,t.result.lnCL$p.value,
                  t.result.MGR$parameter,t.result.MGR$statistic,t.result.MGR$p.value)
}

t.report[c(2:nrow(t.report)),3]<-p.adjust(t.report[c(2:nrow(t.report)),3],
                                          method = "BH", n = nrow(t.report)-1)
t.report[c(2:nrow(t.report)),6]<-p.adjust(t.report[c(2:nrow(t.report)),6],
                                          method = "BH", n = nrow(t.report)-1)
t.report<-round(t.report,4) #round to manageable numbers
t.report[1,c(1,4)]<-c("lnCL","MGR")
write.csv(t.report,"ossification-vs-CL-vs-MGR_ttest.csv")

#  Plot relationship between growth rings and ln-transformed carapace length.
# Color and shape points by the degree of shell closure/fusion
plotsizelength<-ggplot(data=age_measures, aes(x=MGR,y=lnCL)) +
  geom_vline(xintercept = 8) +
  # geom_violin(aes(group=MGR)) +     #boxplot or violin plot
  geom_point(aes(fill=factor(Ossification_score),shape=factor(Ossification_score)),size=3,
             position=position_dodge(width=0.2)) + #plot points
  scale_fill_manual(name="Carapace Ossification", labels=c("open fontanelles",
                                                           "closed fontanelles",
                                                           "partial fusion",
                                                           "complete fusion"),
                    values=oss_col) +
  scale_shape_manual(name="Carapace Ossification", labels=c("open fontanelles",
                                                            "closed fontanelles",
                                                            "partial fusion",
                                                            "complete fusion"),
                     values = c(21,22,23,24)) +
  xlab("Number of Major Growth Rings") + 
  ylab("ln(Carapace Length)") +
  theme_classic() +
  theme(legend.position = c(0.7,0.2),legend.background=element_rect(fill="white"),
        text=element_text(size=12,family="ArialMT"))

# Save as eps and pdf...latter doesn't convert properly to former, former has ill-fitting artboard
# ggsave("ossification-vs-CL-vs-MGR.pdf",plot=plotsizelength, width = 5.4, height = 4, units = "in",family="ArialMT")
cairo_pdf("ossification-vs-CL-vs-MGR.pdf",width = 5.4, height = 4,family ="Arial")
plotsizelength
dev.off()

embed_fonts("ossification-vs-CL-vs-MGR.pdf") 
ggsave("ossification-vs-CL-vs-MGR.eps",plot=plotsizelength, width = 5.4, height = 4, units = "in",family="ArialMT")

# numbers to report, starting with strength of length:MRG correlation
subset<-which(age_measures$MGR>=8)
test1<-lm(age_measures$lnCL[subset]~age_measures$MGR[subset]) %>% summary(.) #%>% 
c(test1$coefficients[2,4],test1$adj.r.squared) %>% round(.,8)
  
which(age_measures$MGR==8) %>% age_measures$CL[.] %>% mean #%>% log
min(age_measures$CL)

which(age_measures$MGR<=8&
        age_measures$CL<=99&
        age_measures$Ossification_score>=2) %>% length/
  which(age_measures$MGR<=8&age_measures$CL<=99) %>% length

which(age_measures$MGR>=8&
        age_measures$CL>=99&
        age_measures$Ossification_score<=1) %>% length/
  which(age_measures$MGR>=8&age_measures$CL>=99) %>% length

which(age_measures$MGR<=7&
        age_measures$Ossification_score>=2) %>% length
which(age_measures$CL>=99&
        age_measures$Ossification_score<=1) %>% length

# SD: SDI ------
# Explore sexual size dimorphism
# boxplot(dim_metadata[,cs_metadata_col]~binary_dim*dim_metadata$sex,main="Centroid Size")
# boxplot(log(dim_metadata[,cs_metadata_col])~binary_dim*dim_metadata$sex,main="ln(CS)") #note difference in views
boxplot(log(dim_metadata$carapace_length)~binary_dim*dim_metadata$sex,main="ln(Carapace Length)") #different ways of measuring size give diff results
AV<-aov(log(dim_metadata$carapace_length)~dim_metadata$sex*binary_dim) %>% summary.aov
write.csv(AV[[1]],"SDI_ANOVA.csv")

# pdf("dim_ANOVA.pdf")
# boxplot(log(dim_metadata$carapace_length)~binary_dim*dim_metadata$sex,main="ln(Carapace Length)") #different ways of measuring size give diff results
# text(2,5.15,paste("sex F(",AV[[1]]$Df[1],",",AV[[1]]$Df[4],") = ",round(AV[[1]]$`F value`[1],2),", p = ",round(AV[[1]]$`Pr(>F)`[1],5),sep=""))
# text(2,5.1,paste("ssp F(",AV[[1]]$Df[2],",",AV[[1]]$Df[4],") = ",round(AV[[1]]$`F value`[2],2),", p = ",round(AV[[1]]$`Pr(>F)`[2],5),sep=""))
# dev.off()

# # Calculate SDI: Gibbons and Lovich 1992, average carapace length. Males are larger, based on boxplot.  
# (size of largest sex/size of smallest sex) +1 if males are larger, negative if males are larger
# 161.9/145.6-1 # calculation of SDI based on St. Clair 1998
# -141.4/125.3+1 # calculation of SDI based on Jones et al. 2016


calculate.SDI<-function(male.length,female.length){
  if(male.length<female.length){
    result<-(1)*female.length/male.length -1
  }
  if(male.length>female.length){
    result<-(-1)*male.length/female.length +1
  }
  return(result)
}


SDI_all<-calculate.SDI(male.length=mean(dim_metadata$carapace_length[which(dim_metadata$sex=="male")]),
                       female.length=mean(dim_metadata$carapace_length[which(dim_metadata$sex=="female")]))

# Use nonparametric modelling of confidence intervals to compare fossil sites.
# Resampling 10 specimens seems appropriate
# Create context by modelling confidence interval via resampling SD dataset
distribution<-NULL #create holder for null resampled distribution
sample1<-dim_metadata$carapace_length[which(dim_metadata$sex=="male")]
sample2<-dim_metadata$carapace_length[which(dim_metadata$sex=="female")]
replicates<-10000
for (i in 1:replicates) {						# start a for loop
  rm1<- sample(sample1, size=5, replace=TRUE) %>% mean	# take first bootstrap sample from pooled data (n=n1)
  rm2<- sample(sample2, size=5, replace=TRUE) %>% mean	# take second bootstrap sample from pooled data (n=n2)
  rdm<- calculate.SDI(male.length=rm1,
                      female.length=rm2)							# compute SDI for a given pair of bootstrap samples
  distribution<- rbind(distribution,rdm)					# store all differences in means
}									# end for loop

# SSP: size -----
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
# F: size histogram ------
#make the histogram. first, combine subspecies and fossil datasets
hist_metadata<-ssp_metadata
hist_metadata$site2<-hist_metadata$site
hist_metadata<-rbind(hist_metadata,fos_metadata)
hist_metadata$site2<-factor(as.character(hist_metadata$site2))

# plot carapace length distribution of fossils compared to modern sample
plotfossize<-ggplot(data=hist_metadata, aes(x=carapace_length,fill=site2)) +
  geom_histogram(data=subset(hist_metadata, !(site2 %in% "modern")),binwidth=10) +
  geom_histogram(data=subset(hist_metadata, site2 %in% "modern"),alpha=0.3, binwidth=10,fill=NA,color="black") +
  scale_fill_manual(name="Site",values=fos_col) +
  xlab("Carapace Length") +
  theme_classic() +
  theme(legend.position = c(0.8,0.55),legend.background=element_rect(fill="white"))
ggsave("CL_fos-vs-mod.pdf",plot=plotfossize, width = 5.4, height = 4, units = "in",family="ArialMT")
embed_fonts("CL_fos-vs-mod.pdf") 
ggsave("CL_fos-vs-mod.eps",plot=plotfossize, width = 5.4, height = 4, units = "in",family="ArialMT")

# F: SDI modelling ------
# # Calculate SDI for the three sites with bimodal distributions
SDI_VM<-1+(-1)*mean(fos_metadata$carapace_length[which(fos_metadata$site2=="melbourne")])/
  mean(fos_metadata$carapace_length[which(fos_metadata$site2=="vero")])
SDI_C<-1+(-1)*mean(fos_metadata$carapace_length[which(fos_metadata$site2=="camelot")])/
  mean(fos_metadata$carapace_length[which(fos_metadata$site2=="camelot2")])
# SDI_I<-1+ (-1)*mean(fos_metadata$carapace_length[which(fos_metadata$site2=="ingleside")])/
#   mean(fos_metadata$carapace_length[which(fos_metadata$site2=="ingleside2")])

# # Possible that these SDI values are just the result of sampling error?
# # Test possibility using two models built above
# # Make plot
# data<-data.frame(val=c(distribution,ssp_distribution),
#                  Dataset=c(rep("sex pre-identified",length(distribution)),rep("sex inferred",length(ssp_distribution))))
distribution<-as.data.frame(distribution)
colnames(distribution)<-"val"

cairo_pdf("SDI.pdf",width = 5.4, height = 4,family ="Arial")
ggplot(data=distribution, aes(x=val)) +
  geom_histogram(alpha=0.5, binwidth=.010,fill="black",color="white") +
  xlab("SDI") +
  geom_vline(xintercept=SDI_C,linetype=2,size=1)+
  geom_vline(xintercept=SDI_VM,linetype=3,size=1)+
  # geom_vline(xintercept=SDI_I,linetype=4,size=1)+
  geom_vline(xintercept=SDI_all,color="black",size=1)+
  geom_text(data=data.frame(x=c(SDI_C+0.17,SDI_VM+0.17),y=c(50,100),site=c("Camelot","Vero-Melbourne")),
                            aes(x=x,y=y,label=site))+
  theme_classic() #+
  # theme(legend.position = c(0.5,0.5),legend.background=element_rect(fill="white"))+
  # scale_color_manual(values=spp_col[c(1,2,1,2)])
dev.off()
embed_fonts("SDI.pdf")

# ggsave("SDI.eps",plot=plotSDI, width = 5.4, height = 4, units = "in",family="ArialMT")
all<-rbind(SDI_ssp,SDI_sspset) %>% cbind(.,c(SDI_all,SDI_sspset_all),
                                         rep(SDI_C,2),rep(SDI_VM,2))
colnames(all)[c(3,4,5)]<-c("all","Camelot","Vero-Melbourne")
rownames(all)<-c("dimorphism_dataset","ssp_dataset")
write.csv(all,"SDI.csv")
