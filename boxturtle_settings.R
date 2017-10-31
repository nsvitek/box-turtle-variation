
# create view-specific settings  
if(view=="dor"){
  shape_cols<-c(2:57)
  shape<-data_all[,shape_cols]
  cs_metadata_col<-2
  conx<-c(1,8,9,8,10,12,10,11,13,12,13,19:28,15,14,14,16,18,16,17,15,17,7,6,5,4,3,2,1,2,9,11) #for drawing outlines
  m<-length(shape_cols)/2
  PCc<-5 # Number of PCs with repeatability > 0.9
  thresh<-0.07 #dorsal view
  
  
}
if(view=="lat"){
  shape_cols<-c(59:146)
  shape<-data_all[,shape_cols]
  cs_metadata_col<-3
  conx<-c(1:2,13:40,3,41:44,4:12,1)
  m<-length(shape_cols)/2
  PCc<-9 # Number of PCs with repeatability > 0.9
  thresh<-0.085 #lateral view
  #etc
}
if(view=="pos"){
  shape_cols<-c(148:213)
  shape<-data_all[,shape_cols]
  cs_metadata_col<-4
  conx<-c(1:3,10:17,4,18:25,5,26:33,6:9,1)
  m<-length(shape_cols)/2
  PCc<-13 # Number of PCs with repeatability > 0.9
  thresh<-0.08 #posterior view
  
  #etc
}

# # Allometry: load --------------------------------------------------------------------
# # Pull the specimens meant for allometry analyses
# # Set filter criteria
# filtersize<-which(data_uni$size_set=="yes")
# # Make dataset and metadataset, remove specimens missing data for a given view
# size_metadata<-makemetadata(data_uni,shape,filtersize,metadata_columns)
# size_set<-shape[filtersize,] %>% na.omit() #make shape matrix

# Sexual Dimorphism: load -----------------------------------------------------------------
setwd(paste(datadir,"../data-output/",sep="/"))
#set filter criteria
filterdim<-which(data_uni$dim_set=="yes")

#make dataset and metadataset, remove specimens missing for that view
if(length(which(is.na(shape[filterdim,1])))>=1){
  dim_metadata<-data_uni[filterdim,metadata_columns] %>% 
    .[-which(is.na(shape[filterdim,1])),]
} else {
  dim_metadata<-data_uni[filterdim,metadata_columns]
}
dim_set<-shape[filterdim,] %>% na.omit()
dim_metadata<-droplevels(dim_metadata)

#make binary coding based on geographic/ssp analysis
binary_dim<-as.character(dim_metadata$ssp)
binary_dim[which(binary_dim!="bauri")]<-"not"
binary_dim<-as.factor(binary_dim)

dim.gpa<-arrayspecs(dim_set,p=ncol(dim_set)/2,k=2) %>% gpagen
dim.gdf<-geomorph.data.frame(dim.gpa,binary=binary_dim,
                             sex=dim_metadata$sex)
dim.gdf$Csize<-dim_metadata[,cs_metadata_col]

# Subspecies: load --------------------------------------------------------------
# pull specimens meant for ssp analyses
filterssp<-which(data_uni$ssp_set=="yes")

# make dataset and metadataset, remove specimens missing for that view
if(length(which(is.na(shape[filterssp,1])))>=1){
  ssp_metadata<-data_uni[filterssp,metadata_columns] %>% 
    .[-which(is.na(shape[filterssp,1])),]
} else {
  ssp_metadata<-data_uni[filterssp,metadata_columns]
}
ssp_set<-shape[filterssp,] %>% na.omit()

#add in sex data from assessments, match by specimen number
for (row_md in 1:nrow(ssp_metadata)){
  row_sd<-which(as.character(ssp_sd_id$id)==as.character(ssp_metadata$id[row_md]))
  ssp_metadata$sex[row_md]<-tolower(ssp_sd_id$assessex[row_sd])
}


ssp_all<-ssp_set #save for later
ssp_all_metadata<-ssp_metadata #save for later
binary_all<-as.character(ssp_all_metadata$ssp)
binary_all[which(binary_all!="bauri")]<-"not"
binary_all<-as.factor(binary_all)

#use only specimens with spatial data (latitude & longitude), can change back later
no.space<-which(is.na(ssp_metadata$latitude))
ssp_set<-ssp_set[-no.space,]
ssp_metadata<-ssp_metadata[-no.space,]
ssp_metadata$latitude<-ssp_metadata$latitude %>% as.character %>% as.numeric #%>% cut(.,100)
ssp_metadata$longitude<-ssp_metadata$longitude %>% as.character %>% as.numeric %>% abs(.) * -1

# try a binary coding, based on a few fragments of info
binary<-as.character(ssp_metadata$ssp)
binary[which(binary!="bauri")]<-"not"
binary<-as.factor(binary)

# make geomorph object for ssp dataset
ssp.gpa<-arrayspecs(ssp_set,p=ncol(ssp_set)/2,k=2) %>% gpagen
ssp.gdf<-geomorph.data.frame(ssp.gpa,binary=binary,ssp=ssp_metadata$ssp,
                             latitude=ssp_metadata$latitude,
                             longitude=ssp_metadata$longitude,
                             sex=ssp_metadata$sex)
ssp.gdf$Csize<-ssp_metadata[,cs_metadata_col]

# Fossils: load -------------------------------------------------------------
#pull specimens meant for fossil analyses
#set filter criteria
filterfos<-which(data_uni$fos_set=="yes"&
                   data_uni$id!="TMM933-3325"&
                   data_uni$id!="TMM933-5680"&
                   data_uni$id!="SCSM2004.01.247"&
                   data_uni$id!="UF23828"&
                   data_uni$id!="TMM30967-660")
#TMM933-3325 too deformed?
#TMM933-5680: it's not 100% clear to me that the shell is closed up. Excluded. 
#UF23828: same as TMM933-5680?
#SCSM2004.1.247: same as TMM933-5650
#TMM90367-660 looks a little odd: oddly fat. But is it biological or taphonomic?

#make dataset and metadataset, remove specimens missing for that view
if(length(which(is.na(shape[filterfos,1])))>=1){
  fos_metadata<-data_uni[filterfos,metadata_columns] %>% 
    .[-which(is.na(shape[filterfos,1])),]
} else {
  fos_metadata<-data_uni[filterfos,metadata_columns]
}
fos_set<-shape[filterfos,] %>% na.omit()
fos_metadata<-droplevels(fos_metadata)

dim(fos_set)
dim(fos_metadata)

# create some alternate codes
fos_metadata$site2<-fos_metadata$site
levels(fos_metadata$site2)<-c(levels(fos_metadata$site2),"holocene","camelot2","ingleside2","haile2")
fos_metadata$site2[which(fos_metadata$site2=="grove's orange midden"|
                           fos_metadata$site2=="fort center")]<-"holocene"
fos_metadata$site2[which(fos_metadata$site2=="camelot"&
                           fos_metadata$carapace_length<=150)]<-"camelot2"
# fos_metadata$site2[which(fos_metadata$site2=="ingleside"&
#                            fos_metadata$carapace_length<=155)]<-"ingleside2"
fos_metadata$site2<-droplevels(fos_metadata$site2)

# reorder sites by name
fos_metadata$site2<-factor(as.character(fos_metadata$site2))

# #reorder sites by carapace length, small to large
# site_order<-group_by(fos_metadata,site2) %>% summarise(., mean(carapace_length)) #%>% arrange(.,`mean(carapace_length)`)
# fos_metadata$site2<-factor(fos_metadata$site2,levels(fos_metadata$site2)[order(site_order$`mean(carapace_length)`)])
# levels(fos_metadata$site2)
