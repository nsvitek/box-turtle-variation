# Extra, unused code for box turtle project that might be handy late

# Code to create and look at "allometric stability" or "age independence" of shape, updated
# # SUSPENDED: Take a look at stability of allometric relationship given dataset composition
# # ringcor<-c(1:max(age_measures$MGR))
# # for (i in 1:length(viewoptions)){ #cycle through all the views
# #   view<-viewoptions[i]
# #   source(paste(scriptsdir,"boxturtle_settings.R",sep="/"))
# #   size_metadata<-makemetadata(data_uni,shape,filtersize,metadata_columns)
# #   size_set<-shape[filtersize,] %>% na.omit() #make shape matrix
# #   lCS<-log(size_metadata[,cs_metadata_col])
# #   cor1<-matrix(c(1:max(age_measures$MGR),rep(0,max(age_measures$MGR)*2)),ncol=3,byrow=FALSE)
# #   colnames(cor1)<-paste(view,c("mgr","rsmall","rlarge"),sep="_") #build empty matrix
# #   for (i in 1:max(age_measures$MGR)){  # build up dataset from "youngest" specimens  
# #     which_small<-which(size_metadata$major_growth_rings<i)
# #     if (length(which_small)<=10){
# #       print("Sample size too small. Skipping.")
# #       cor1[i,2]<-NA
# #       next
# #     }
# #     else {
# #       subsample1<-size_set[which_small,] #subsample shape data
# #       cor1[i,2]<-adonis(subsample1~lCS[which_small],permutations=1000,method="euclidean") %>% #MANOVA shape~size
# #         .[[("aov.tab")]] %>% .[1,5] #pull r^2
# #       print(paste("Calculated correlation for set",i))
# #     }
# #   }
# #   ringcor<-cbind(ringcor,cor1[,c(2:3)])
# # }
# # ringcorsmall<-melt(ringcor[,c(2,4,6)],id="ringcor")
# # plotcorsmall<-ggplot(data=ringcorsmall, aes(x=X1,y=value,color=X2,shape=X2,linetype=X2)) +
# #   geom_vline(xintercept = 9) +
# #   geom_point() +
# #   geom_line() +
# #     scale_color_manual(name="View", labels=c("dorsal","lateral","posterior"),
# #                       values=rep(view_col[1:3],2)) +
# #     scale_shape_manual(name="View", labels=c("dorsal","lateral","posterior"),
# #                        values=rep(c(16:18),2)) +
# #     scale_linetype_manual(name="View", labels=c("dorsal","lateral","posterior"),
# #                        values=rep(c(2:4),2)) +
# #     xlab("Number of Major Growth Rings") + 
# #     ylab("R-squared") +
# #     theme_minimal() +
# #     theme(legend.position = c(0.1,0.8),legend.background=element_rect(fill="white"))
# # grid.arrange(plotsizelength,plotcorsmall,plotcorlarge,ncol=1)

#update of old, Zelditch/MS way of plotting shapes at minimum and maximum size
# fit_val1<-model1$fitted.values
# small1<-matrix(fit_val[which(lCS==(min(lCS))),],m,k,byrow=TRUE)	
