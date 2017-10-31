# # Script to concatenate all relevant, clean data into single, saved table.
# insert NAs for specimens missing a view. Include metadata, coordinates, cs

#dependencies
library(geomorph)
library(dplyr)

# datadir <- "C:/Users/N.S/Dropbox/Documents/research/turtles/Thesis/Vitek_YR_PublishTerrapene/data"
datadir <- "D:/Dropbox/Documents/research/turtles/Thesis/Vitek_YR_PublishTerrapene/data"

setwd(datadir)

freeze<-ls()
# read in lateral tps data and file names, turn into 2D array matching other views
setwd("tps-intermediate")
namefile<-read.csv("boxturtle_lat_images.csv",header=FALSE)
raw<-readland.tps(file="boxturtle_lat_align.tps") #%>% gpagen #superimpose
# str(raw)
# IMP<-two.d.array(raw$coords) #make landmarks into 2D array
IMP<-two.d.array(raw)
colnames(IMP)<-paste("LatCoord",seq(from=1,to=ncol(IMP),by=1),sep="")
dim(IMP)
x<-namefile[,2]
m<-gregexpr("lateral\\\\.*_[0-9]",x,perl=TRUE)
shapenames<-regmatches(x,m) %>% unlist %>% gsub("lateral\\\\(.*)_[0-9]","\\1",.,perl=TRUE)
rownames(IMP)<-shapenames
lat_shape<-cbind(shapenames,IMP)
write.csv(lat_shape,"boxturtle_lat_coords.csv",row.names=FALSE)

rm(list = setdiff(ls(),c(freeze,"lat_shape")))
setwd("..")

# read in the posterior and lateral view coordinates, plus centroid sizes
pos_shape<-read.csv("boxturtle_pos_coords.csv") 
dor_shape<-read.csv("boxturtle_dor_coords.csv")

lat_cs<-read.table("boxturtle_lat_centroid.NTS",skip=3,col.names="lat_cs")
dor_cs<-read.table("boxturtle_dor_centroid.NTS",skip=3,col.names="dor_cs")
pos_cs<-read.table("boxturtle_pos_centroid.NTS",skip=3,col.names="pos_cs")


# bind centroid sizes to appropriate shape blocks
pos_shcs<-cbind(pos_shape,pos_cs)
dor_shcs<-cbind(dor_shape,dor_cs)
lat_shcs<-cbind(lat_shape,lat_cs)

# match rows, for specimens missing a view, enter in NA's
levels(lat_shcs$shapenames)<-c(levels(lat_shcs$shapenames),"UF5700")
levels(pos_shcs$id)<-c(levels(pos_shcs$id),"USNMV12000")
levels(dor_shcs$id)<-c(levels(dor_shcs$id),"TMM933-3039","USNMV11834",
                       "NCSM73537","SCSM91.166.1")

insert1<-c("UF5700",rep(NA,(ncol(lat_shcs)-1)))
lat_shcs<-rbind(lat_shcs[1:35,],insert1,lat_shcs[36:nrow(lat_shcs),])

insert2<-c("USNMV12000",rep(NA,(ncol(pos_shcs)-1)))
pos_shcs<-rbind(pos_shcs[1:42,],insert2,pos_shcs[43:nrow(pos_shcs),])
pos_shcs<-pos_shcs[c(1:47,61:71,48,72:73,49,74:75,50:51,76:80,52:60,81:nrow(pos_shcs)),]
pos_shcs<-pos_shcs[c(1:73,81:89,74,90:105,75,106:109,76:77,110:125,78:80,126:nrow(pos_shcs)),]
pos_shcs<-pos_shcs[c(1:123,126:128,124,129:141,125,142:nrow(pos_shcs)),]

#CHECK: If any specimen is printed from this for loop then there are errors to be fixed.
for (i in 1:nrow(lat_shcs)){
  if(as.character(lat_shcs$shapenames[i])!=as.character(pos_shcs$id[i])){
    print(c(i,as.character(lat_shcs$shapenames[i]),as.character(pos_shcs$id[i])))
  }
}

# lateral and posterior views should be standardized at this point. Bind together.
shape_latpos<-cbind(lat_shcs,pos_shcs[,2:ncol(pos_shcs)])

# repeat row standardization for dorsal view. 
insert3<-c("TMM933-3039",rep(NA,(ncol(dor_shcs)-1)))
dor_shcs<-rbind(dor_shcs[1:12,],insert3,dor_shcs[13:nrow(dor_shcs),])

insert4<-c("USNMV11834",rep(NA,(ncol(dor_shcs)-1)))
dor_shcs<-rbind(dor_shcs[1:41,],insert4,dor_shcs[42:nrow(dor_shcs),])
dor_shcs<-dor_shcs[c(1:45,54:59,46,60:61,47,62:63,48,64,49,65:102,50,103:nrow(dor_shcs)),]

insert5<-c("NCSM73537",rep(NA,(ncol(dor_shcs)-1)))
dor_shcs<-rbind(dor_shcs[1:108,],insert5,dor_shcs[109:nrow(dor_shcs),])

shape_latpos<-shape_latpos[-which(shape_latpos=="UF21116")[1],]
shape_latpos<-shape_latpos[-which(shape_latpos=="UF22161")[1],]

insert6<-c("SCSM91.166.1",rep(NA,(ncol(dor_shcs)-1)))
dor_shcs<-rbind(dor_shcs[1:368,],insert6,dor_shcs[369:nrow(dor_shcs),])

for (i in 1:nrow(shape_latpos)){
  if(as.character(shape_latpos$shapenames[i])!=as.character(dor_shcs$id[i])){
    print(c(i,as.character(shape_latpos$shapenames[i]),as.character(dor_shcs$id[i])))
  }
}

#create object of all shape coordinates & centroid sizes
shape_all<-cbind(dor_shcs,shape_latpos[,2:ncol(shape_latpos)])
levels(shape_all$id)<-c(levels(shape_all$id),"UF_ea_04600045")
shape_all$id[which(shape_all$id=="UF_ea_046000645")]<-"UF_ea_04600045"

subspecies<-ssp<-fos_set<-size_set<-ssp_set<-sex<-dim_set<-err_rep<-latitude<-longitude<-rep(NA,nrow(shape_all))
morph<-rep(1,nrow(shape_all))
spp<-rep("carolina",nrow(shape_all))
shape_all<-cbind(shape_all,fos_set,size_set,ssp_set,dim_set,subspecies,ssp,spp,sex,morph,err_rep,latitude,longitude)
levels(shape_all$fos_set)<-levels(shape_all$ssp_set)<-levels(shape_all$dim_set)<-levels(shape_all$size_set)<-c("yes","no")
levels(shape_all$spp)<-c("carolina","coahuila")
levels(shape_all$sex)<-c("male","female")

# read in dataset block
temp_set<-read.csv("FosSizessp_Pos_RData.csv",header=TRUE)
colnames(temp_set)
temp_set<-temp_set[,c(1,73:82)]
levels(shape_all$subspecies)<-levels(temp_set$Subspecies)
levels(shape_all$ssp)<-levels(temp_set$ssp)
# levels(shape_all$latitude)<-levels(temp_set$Latitude)
# levels(shape_all$longitude)<-levels(temp_set$Longitude)
temp_set$X<-gsub("_", "", temp_set$X, perl=TRUE)
temp_set$X<-gsub("/", "", temp_set$X, perl=TRUE)
# temp_set[which(temp_set$X=="UF21116"|temp_set$X=="UF22161"),]

# check and edit dataset affiliations
temp_set$X[which(!temp_set$X %in% as.character(shape_all$id))]
temp_set$X[which(temp_set$X=="UFFGS1645")]<-"UFFGSV1645"
temp_set$X[which(temp_set$X=="UFFGS1648")]<-"UFFGSV1648"
temp_set$X[which(temp_set$X=="UFFGS278")]<-"UFFGSV278"
temp_set<-temp_set[-which(temp_set$X=="NCSM62617"),]

for (i in 1:nrow(shape_all)){
  if(shape_all$id[i] %in% temp_set$X){
    r<-which(temp_set$X==shape_all$id[i])
    # print(c(r,i))
    shape_all$fos_set[i]<-tolower(as.character(temp_set$Fos_Set[r]))
    shape_all$size_set[i]<-tolower(as.character(temp_set$Size_Set[r]))
    shape_all$ssp_set[i]<-tolower(as.character(temp_set$ssp_Set[r]))
    shape_all$dim_set[i]<-tolower(as.character(temp_set$dim_set[r]))
    shape_all$sex[i]<-tolower(as.character(temp_set$sex[r]))
        shape_all$ssp[i]<-tolower(as.character(temp_set$ssp[r]))
    shape_all$subspecies[i]<-tolower(as.character(temp_set$Subspecies[r]))
    shape_all$latitude[i]<-as.character(temp_set$Latitude[r])
    shape_all$longitude[i]<-as.character(temp_set$Longitude[r])
    }
}

# length(which(temp_set$ssp_Set=="Yes"))
# length(which(shape_all2$ssp_set=="yes"))

length(which(tolower(as.character(temp_set$ssp_Set))=="yes"))
shape_all2<-shape_all #manual changes
shape_all2[,c(1,215:226)]
shape_all2[c(371:391),223]<-2 #coahuila morphotype
shape_all2[c(371:391),221]<-"coahuila" #coahuila morphotype
shape_all2[c(346:347,351:352,354:355,392:395),224]<-1 #error replicates
shape_all2[c(36,345,348:350,353,356:370,346:347,351:352,354:355,371:391),215]<-"yes" #forlorn fossils
shape_all2[c(36,345,348:350,353,356:370,346:347,351:352,354:355,371:391),216]<-"no" #forlorn fossils
shape_all2[c(36,345,348:350,353,356:370,346:347,351:352,354:355,371:391),c(217,218)]<-"no" #forlorn fossils
shape_all2[c(346:347,351:352,354:355,371:395),215]<-"no"

# [which(!temp_set$X[which(temp_set$ssp_Set=="Yes")] %in% shape_all2$id[which(shape_all2$ssp_set=="yes")])]



# read in specimen metadata
shape_all<-shape_all2

metadata_fos<-read.csv("boxturtlefos_metadata.csv")
metadata_mod<-read.csv("boxturtlemod_metadata.csv")

id<-paste(metadata_mod$Institution,metadata_mod$Number,sep="")
metadata_mod<-cbind(metadata_mod,id)

# want CL, state, locality/county,MGR,
carapace_length<-state<-county<-site<-major_growth_rings<-rep(NA,nrow(shape_all2))
levels(state)<-unique(c(tolower(levels(metadata_mod$State)),tolower(levels(metadata_fos$state))))
levels(county)<-unique(tolower(levels(metadata_mod$County)))
levels(site)<-unique(tolower(levels(metadata_fos$locality_data)))

shape_all3<-cbind(shape_all2,carapace_length,major_growth_rings,state,county,site)

# match metadata to shape, dataset block
for (i in 1:nrow(shape_all3)){
  if(shape_all3$id[i] %in% as.character(metadata_mod$id)){
    r<-which(as.character(metadata_mod$id)==shape_all3$id[i])
    shape_all3$state[i]<-tolower(as.character(metadata_mod$State[r]))
    shape_all3$county[i]<-tolower(as.character(metadata_mod$County[r]))
    shape_all3$major_growth_rings[i]<-metadata_mod$MGR[r]
    shape_all3$carapace_length[i]<-metadata_mod$CL[r]
    shape_all3$site[i]<-"modern"
  } else if(shape_all3$id[i] %in% as.character(metadata_fos$specimen)){
    r<-which(as.character(metadata_fos$specimen)==shape_all3$id[i])
    shape_all3$state[i]<-tolower(as.character(metadata_fos$state[r]))
    shape_all3$site[i]<-tolower(as.character(metadata_fos$locality_data[r]))
    shape_all3$carapace_length[i]<-metadata_fos$carapace_length_mm[r]
  }
}

###NEED TO CHANGE COUNTY/SITE TO SEPARATE THINGS

# Write final dataset to a csv
dim(shape_all3)
head(shape_all3[,c(1,215:226)])
shape_all3$err_rep[which(is.na(shape_all3$err_rep))]<-0 #last minute: finish filling out error classification

data_rep<-filter(shape_all3,err_rep==1) #set aside intentional replicates
data_uni<-filter(shape_all3,err_rep==0) 
data_uni<-data_uni[!duplicated(data_uni$id),]#remove unintentional replicates
length(which(duplicated(data_uni$id)==TRUE))
shape_all4<-rbind(data_uni,data_rep)
write.csv(shape_all4,"boxturtle_Rdata.csv",row.names=FALSE)
dim(shape_all4)

