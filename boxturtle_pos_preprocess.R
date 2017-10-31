#The script for taking care of missing landmarks in object-symmetrical datasets, such as turtle shells in certain views. Before you get here, you should already have landmarked your dataset, checked for correct landmarks, made your semi-landmarks, appended your curves to landmarks, deleted your overlappoing semi-landmarks, made your sliders file, done your semi-landmark superimposition, and double-checked that everything looked good. Then, you're going to have to reformat your data into the OSymm format from the IMP format. I'm going to try to make a loop for that here, and hope that it works.
#OSymm format: the each column is a specimen, each row is a landmark. But X and Y are each a column. 
#IMP:sp1x1 sp1y1 sp1x2
#IMP:sp2x1 sp2y1 sp2x2
#IMP:sp3x1 sp3y1 sp3x2
#IMP:sp4x1 sp4y1 sp4x2

#OS: sp1x1 sp1y1 sp2x1 sp2y1
#OS: sp1x2 sp1y2 sp2x2 sp2y2
#OS: sp1x3 sp1y3 sp2x3 sp2y3
#OS: sp1x4 sp1y4 sp2x4 sp2y4

#17 September 2012. Right now, the file we need is FosSize_Pos_PCA.txt. It's got 63 landmarks, 126 columns, and 146 specimens (rows). Read it in. It has no first column of specimen ID's, so you don't need to remove them. If you did, you'd use the following two lines:
#IMP2<-IMP[,2:NUMBER.OF.COLUMNS.PLUS.ONE]
#IMP<-IMP2
#7 October 2016: It's 4 years since I last used this script. Below are some of the first lines of code and loops I ever wrote. They suck. But not more than the time it would take to fix them. 

#load dependencies
library(geomorph)
library(dplyr)

datadir<-"C:/Users/N.S/Dropbox/Documents/research/turtles/Thesis/Vitek_YR_PublishTerrapene/data/tps-intermediate"
setwd(datadir)
namefile<-read.csv("boxturtle_pos_images.csv",header=FALSE)
raw<-readland.tps(file="boxturtle_pos_align.tps")
str(raw)
#make landmarks into 2D array
IMP<-two.d.array(raw)
#Headers should be "RawCoord#" (e.g., RawCoord1,RawCoord2,RawCoord3, etc.) for historical reasons :-P
colnames(IMP)<-paste("RawCoord",seq(from=1,to=ncol(IMP),by=1),sep="")
dim(IMP)

#Then, transpose, so that our rows become columns and vice versa.

tr.IMP<-t(IMP)

#Now I want to take every odd row and make a list of that, and then every even row and make a list of that. The odd rows are my x's and the even rows are my y's. To do that, I'll need a blank vector that is my number of rows. 

land<-1:126

exes.l<-list()
whys.l<-list()
for (i in 1:nrow(IMP)) {
	exes.l[[i]]<-tr.IMP[seq.int(1L, length(land),2),i]
	whys.l[[i]]<-tr.IMP[seq.int(2L,length(land),2),i]
}

#Now I have these two lists of rows. Next, we're going to do what Doug said, and write a loop that takes one vector from the x's and then one vector from the y's and then turns that into a list.

NewTable<-list()

j<-1
for (i in 1:nrow(IMP)){
	NewTable[[j]]<-exes.l[[i]]
	j<-j+1
	NewTable[[j]]<-whys.l[[i]]
	j<-j+1
}

#Now, we'll turn that into a matrix.

OS_Format<-do.call(cbind,NewTable)
dim(OS_Format)
head(OS_Format)

################################## The New Way ####################################
# Go in and put the NA's here in R instead of opening your exporting your new table, opening it in Excel, and replacing missing landmarks with NA's there. Hopefully you wrote down which landmarks are missing in a text document somewhere that you can easily open now. I think what you'd do is something like the following line for each missing landmark, if n is your specimen number and k is your landmark number. Say, for example, in specimen 1 landmark 6 is missing. You'd use OS_Format[6,1:2].

#OS_Format[k,2n-1:2n]<-NA

###  ENTER IN THE NUMBER OF REPLACEMENTS YOU NEED IN THE CORRECT PLACE FOR YOUR DATASET.

OS_Format[11,21:22]<-NA
OS_Format[24,21:22]<-NA
OS_Format[25,21:22]<-NA
OS_Format[26,21:22]<-NA
OS_Format[27,21:22]<-NA
OS_Format[5,45:46]<-NA
OS_Format[43,45:46]<-NA
OS_Format[44,45:46]<-NA
OS_Format[45,45:46]<-NA
OS_Format[46,45:46]<-NA
OS_Format[47,45:46]<-NA
OS_Format[7,59:60]<-NA
OS_Format[8,59:60]<-NA
OS_Format[10,75:76]<-NA
OS_Format[40,75:76]<-NA
OS_Format[41,75:76]<-NA
OS_Format[42,75:76]<-NA
OS_Format[43,75:76]<-NA
OS_Format[44,75:76]<-NA
OS_Format[45,75:76]<-NA
OS_Format[46,75:76]<-NA
OS_Format[47,75:76]<-NA
OS_Format[48,75:76]<-NA
OS_Format[49,75:76]<-NA
OS_Format[50,75:76]<-NA
#0
OS_Format[6,89:90]<-NA
OS_Format[37,89:90]<-NA
OS_Format[38,89:90]<-NA
#2
OS_Format[6,93:94]<-NA
OS_Format[7,93:94]<-NA
OS_Format[36,93:94]<-NA
OS_Format[37,93:94]<-NA
OS_Format[38,93:94]<-NA
OS_Format[39,93:94]<-NA
#3
OS_Format[6,95:96]<-NA
OS_Format[7,95:96]<-NA
OS_Format[37,95:96]<-NA
OS_Format[38,95:96]<-NA
OS_Format[39,95:96]<-NA
#4
OS_Format[7,97:98]<-NA
OS_Format[8,97:98]<-NA
OS_Format[32,97:98]<-NA
OS_Format[33,97:98]<-NA
OS_Format[34,97:98]<-NA
OS_Format[35,97:98]<-NA
#5
OS_Format[6,99:100]<-NA
OS_Format[7,99:100]<-NA
OS_Format[8,99:100]<-NA
OS_Format[38,99:100]<-NA
OS_Format[39,99:100]<-NA
#6
OS_Format[8,101:102]<-NA
OS_Format[37,101:102]<-NA
#7
OS_Format[7,103:104]<-NA
#8
OS_Format[7,105:106]<-NA
#9
OS_Format[7,107:108]<-NA
#10
OS_Format[6,109:110]<-NA
OS_Format[7,109:110]<-NA
OS_Format[35,109:110]<-NA
OS_Format[36,109:110]<-NA
OS_Format[37,109:110]<-NA
OS_Format[38,109:110]<-NA
OS_Format[39,109:110]<-NA
#11
OS_Format[6,111:112]<-NA
OS_Format[7,111:112]<-NA
OS_Format[34,111:112]<-NA
OS_Format[35,111:112]<-NA
OS_Format[36,111:112]<-NA
OS_Format[37,111:112]<-NA
OS_Format[38,111:112]<-NA
OS_Format[39,111:112]<-NA
#12
OS_Format[11,113:114]<-NA
OS_Format[12,113:114]<-NA
OS_Format[55,113:114]<-NA
OS_Format[56,113:114]<-NA
OS_Format[57,113:114]<-NA
OS_Format[58,113:114]<-NA
OS_Format[59,113:114]<-NA
OS_Format[60,113:114]<-NA
OS_Format[61,113:114]<-NA
OS_Format[62,113:114]<-NA
OS_Format[63,113:114]<-NA
#13
OS_Format[4,115:116]<-NA
OS_Format[5,115:116]<-NA
OS_Format[6,115:116]<-NA
OS_Format[7,115:116]<-NA
OS_Format[16,115:116]<-NA
OS_Format[17,115:116]<-NA
OS_Format[18,115:116]<-NA
OS_Format[19,115:116]<-NA
OS_Format[20,115:116]<-NA
OS_Format[21,115:116]<-NA
OS_Format[22,115:116]<-NA
OS_Format[23,115:116]<-NA
OS_Format[24,115:116]<-NA
OS_Format[25,115:116]<-NA
OS_Format[26,115:116]<-NA
OS_Format[27,115:116]<-NA
OS_Format[28,115:116]<-NA
OS_Format[29,115:116]<-NA
OS_Format[30,115:116]<-NA
OS_Format[31,115:116]<-NA
OS_Format[32,115:116]<-NA
OS_Format[33,115:116]<-NA
OS_Format[34,115:116]<-NA
OS_Format[35,115:116]<-NA
OS_Format[36,115:116]<-NA
OS_Format[37,115:116]<-NA
OS_Format[38,115:116]<-NA
OS_Format[39,115:116]<-NA
#14
OS_Format[8,117:118]<-NA

#ssp data +2*146=+292
OS_Format[6,327:328]<-NA
OS_Format[7,327:328]<-NA
OS_Format[6,331:332]<-NA
OS_Format[7,331:332]<-NA
OS_Format[8,339:340]<-NA
OS_Format[12,345:346]<-NA
OS_Format[6,351:352]<-NA
OS_Format[7,351:352]<-NA
OS_Format[8,351:352]<-NA
OS_Format[32,351:352]<-NA
OS_Format[33,351:352]<-NA
OS_Format[34,351:352]<-NA
OS_Format[35,351:352]<-NA
OS_Format[36,351:352]<-NA
OS_Format[37,351:352]<-NA
OS_Format[38,351:352]<-NA
OS_Format[39,351:352]<-NA
OS_Format[6,379:380]<-NA
OS_Format[7,379:380]<-NA
OS_Format[12,375:376]<-NA
OS_Format[14,375:376]<-NA
OS_Format[9,431:432]<-NA
OS_Format[12,461:462]<-NA
OS_Format[13,461:462]<-NA
OS_Format[14,461:462]<-NA
OS_Format[15,461:462]<-NA
OS_Format[14,469:470]<-NA
OS_Format[13,563:564]<-NA
OS_Format[6,579:580]<-NA
OS_Format[7,579:580]<-NA
OS_Format[8,579:580]<-NA
OS_Format[6,605:606]<-NA
OS_Format[7,605:606]<-NA
OS_Format[12,609:610]<-NA
OS_Format[6,619:620]<-NA
OS_Format[6,619:620]<-NA
OS_Format[6,619:620]<-NA
OS_Format[6,673:674]<-NA
OS_Format[7,673:674]<-NA
OS_Format[6,675:676]<-NA
OS_Format[32,677:678]<-NA
OS_Format[33,677:678]<-NA
OS_Format[34,677:678]<-NA
OS_Format[35,677:678]<-NA
OS_Format[36,677:678]<-NA
OS_Format[37,677:678]<-NA
OS_Format[38,677:678]<-NA
OS_Format[39,677:678]<-NA

#holo data: original numbers + 346  [+(n*2) + [2n-1:2n]]
hcount<-346
h5<-c(5,5,6,7,8,24:39)	 # - n5 lm 5, 6, 7, 8, 24-39
h9<-c(9,5,6,7,8,24:39)	 # - n9 lm 5, 6, 7, 8, 24-39
h10<-c(10,5,6,7,8)	 # - n10 lm 5, 6, 7, 8
h11<-c(11,5,6,7,8)	 # - n11 lm 5, 6, 7, 8
OS_Format[h5[-1],(c((h5[1])*2-1,(h5[1])*2)+hcount)]<-NA
OS_Format[h9[-1],(c((h9[1])*2-1,(h9[1])*2)+hcount)]<-NA
OS_Format[h10[-1],(c((h10[1])*2-1,(h10[1])*2)+hcount)]<-NA
OS_Format[h11[-1],(c((h11[1])*2-1,(h11[1])*2)+hcount)]<-NA

#Scoahuila data: original numbers + 357  [+(n*2) + [2n-1:2n]]
scount<-357
s2<-c(2,4,5,6,7,16:39)	 # - n2: lm 4, 5, 6, 7, 16:39
s3<-c(3,8)	 # - n3: lm 8
s7<-c(7,13:15)	 # - n7: lm 13, 14, 15
s14<-c(14,7,8)	 # - n14: lm 7, 8
s19<-c(19,14,15)	 # - n19: lm 14, 15
OS_Format[s2[-1],(c((s2[1])*2-1,(s2[1])*2)+scount)]<-NA
OS_Format[s3[-1],(c((s3[1])*2-1,(s3[1])*2)+scount)]<-NA
OS_Format[s7[-1],(c((s7[1])*2-1,(s7[1])*2)+scount)]<-NA
OS_Format[s14[-1],(c((s14[1])*2-1,(s14[1])*2)+scount)]<-NA
OS_Format[s19[-1],(c((s19[1])*2-1,(s19[1])*2)+scount)]<-NA

XX<-OS_Format

#If this were the case, you might not even have to write a new file at this point. You could just keep going. If you don't want to do it this way, go do it the old way.
################################## The Old Way #####################################
#write.table(OS_Format, file.choose(new=TRUE))

#Okay. Assuming you got that working, here's what you'll have to do next: go into the new file by hand, find the coordinates that are missing, and replace the values with "NA". Save that. Then, you'll have to take a look at the OSymm.R file, and alter it to fit your particular data set. Here, for the turtle shell in posterior view, we're going to need to rename the midline, right, and left coordinates. Midline are just going to be 1, 2, and 3. Right and left will follow a pairs file that you hopefully made (try Pos_pairs.txt). Now, read that whole function into R. I'll try to copy and paste all that below. 

#XX <- read.table(file.choose(), row.names=1, skip=1)
#####################################################################################
Fixed<-list()
j<-2
for (i in 1:nrow(IMP)){
  q<-i*2-1
  X <- as.matrix(XX[,q:j])
  midline <- c("RawCoord1","RawCoord3", "RawCoord5")
  right <- c("RawCoord7","RawCoord9","RawCoord11","RawCoord13", "RawCoord15","RawCoord17","RawCoord31","RawCoord33","RawCoord35","RawCoord37","RawCoord39","RawCoord41","RawCoord43","RawCoord45","RawCoord47","RawCoord49","RawCoord51","RawCoord53","RawCoord55","RawCoord57","RawCoord59","RawCoord61","RawCoord63","RawCoord65","RawCoord67","RawCoord69","RawCoord71","RawCoord73","RawCoord75","RawCoord77")
  left <- c("RawCoord19","RawCoord21","RawCoord23","RawCoord25","RawCoord27","RawCoord29","RawCoord79","RawCoord81","RawCoord83","RawCoord85","RawCoord87","RawCoord89","RawCoord91","RawCoord93","RawCoord95","RawCoord97","RawCoord99","RawCoord101","RawCoord103","RawCoord105","RawCoord107","RawCoord109","RawCoord111","RawCoord113","RawCoord115","RawCoord117","RawCoord119","RawCoord121","RawCoord123","RawCoord125")
  ncl <- ncol(X)

  OSymm <- function(X, midline, right, left) {
	ncl <- ncol(X)
	Xr <- cbind(X[,-ncl], -X[,ncl])
	Xow <- Xo <- rbind(X[c(midline, right, left),])
	Xrw <- Xr <- rbind(Xr[c(midline, left, right),])
	rownames(Xrw) <- rownames(Xr) <- rownames(X)
	Xo[which(is.na(Xr))] <- NA
	Xr[which(is.na(Xo))] <- NA
	mo <- matrix(apply(Xo, 2, mean, na.rm=TRUE), byrow=TRUE, nr=nrow(Xow), nc=ncol(Xow))
	mr <- matrix(apply(Xr, 2, mean, na.rm=TRUE), byrow=TRUE, nr=nrow(Xrw), nc=ncol(Xrw))
	Xrwc <- Xrw-mr
	SVD <- svd(t(na.omit(Xr-mr)) %*% na.omit(Xo-mo))
	L <- diag(SVD$d)
	S <- ifelse(L<0, -1, L)
	S <- ifelse(L>0, 1, L)
	RM <- SVD$v %*% S %*% t(SVD$u)
	Xrot <- (Xow-mo) %*% RM                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                 
	SC <- apply(array(c(Xrwc,Xrot), dim=c(nrow(Xrot),ncol(Xrot),2), dimnames=list(rownames(Xrot),colnames(Xrot))), 1:2, mean, na.rm=TRUE)
	Xrot[which(is.na(Xrot))] <- Xrwc[which(is.na(Xrot))]
	Xrwc[which(is.na(Xrwc))] <- Xrot[which(is.na(Xrwc))]
	list(rec.orig=Xrot, symmconf=SC, rec.ref=Xrwc)}
	
  Fix<-OSymm(X,midline,right,left)
  Fixed[[i]]<-Fix$symmconf
  j<-j+2
}


Fixed_Missing<-do.call(cbind,Fixed)
dim(Fixed_Missing)
which(is.na(Fixed_Missing))
#write.table(Fixed_Missing, file.choose(new=TRUE))

#Then, you will want to sort the rows back in their original order. You can't just use the sort function because it won't sort the numbers properly, so I'm going to it by hand here, because I'm too lazy to open excel.

row.names(Fixed_Missing)
sortFM1<-Fixed_Missing[1:9,]
sortFM2<-Fixed_Missing[34:39,]
sortFM3<-Fixed_Missing[10:33,]
sortFM4<-Fixed_Missing[40:63,]
sortFM<-rbind(sortFM1,sortFM2,sortFM3,sortFM4)
row.names(sortFM)

#Now, get table back into IMP format.

sea<-1:ncol(Fixed_Missing)
exes.l<-list()
whys.l<-list()
for (i in 1:63) {
	exes.l[[i]]<-sortFM[i,seq.int(1L,length(sea),by=2)]
	whys.l[[i]]<-sortFM[i,seq.int(2L,length(sea),by=2)]
}

NewTable<-list()
j<-1
for (i in 1:63){
	NewTable[[j]]<-exes.l[[i]]
	j<-j+1
	NewTable[[j]]<-whys.l[[i]]
	j<-j+1
}

Back<-data.frame(matrix(unlist(NewTable),ncol=126,))
Back2<-as.matrix(Back)
dim(Back2)

x<-namefile[,2]
m<-gregexpr("posterior\\\\.*_[0-9]",x,perl=TRUE)
shapenames<-regmatches(x,m) %>% unlist %>% gsub("posterior\\\\(.*)_[0-9]","\\1",.,perl=TRUE)
rownames(Back2)<-shapenames
write.table(Back2, "boxturtle_pos_IMP.txt",row.names=FALSE,col.names=FALSE,sep="  ")
write.csv(Back2, "boxturtle_pos_morphoj.csv",row.names=TRUE,col.names=FALSE)

# # # # # # # # # # # # # # # # # # # # # # # # # # # #
# # # # # # # # # # # # # # # # # # # # # # # # # # # #
# # # # # # # # # # # # # # # # # # # # # # # # # # # #
# # # # # # # # # # # # # # # # # # # # # # # # # # # #
# # # # # # # # # # # # # # # # # # # # # # # # # # # #
# At this point, need to go into MorphoJ.
# Procrustes fit then export the symmetric component 
# [input: _imagej.csv; output: _symm.txt]
# # # # # # # # # # # # # # # # # # # # # # # # # # # #
# # # # # # # # # # # # # # # # # # # # # # # # # # # #
# # # # # # # # # # # # # # # # # # # # # # # # # # # #
# # # # # # # # # # # # # # # # # # # # # # # # # # # #
# # # # # # # # # # # # # # # # # # # # # # # # # # # #

getwd()
symmdat<-read.table("boxturtle_pos_symm.txt",sep = "\t",header=TRUE)
dim(symmdat)
str(symmdat)

#take out the coordinates of one side. 
#for posterior: delete Symm19-30, Symm79-126; 
half<-symmdat[-1,-(c(19:30,79:126)+1)]

# superimpose landmarks
# sh<-half[,-1] %>% as.matrix %>%  arrayspecs(.,ncol(.)/2,2)  %>% gpagen #2 for 2 dimensions
# sh<-two.d.array(sh$coords)
# half[,-1]<-sh

#add column names
colnames(half)<-c("id",paste("PosCoord",seq(from=1,to=(ncol(half)-1),by=1),sep=""))
write.csv(half,"boxturtle_pos_coords.csv",row.names=FALSE)
dim(half)
head(half)
