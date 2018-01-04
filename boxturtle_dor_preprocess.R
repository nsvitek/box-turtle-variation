#The script for taking care of missing landmarks in object-symmetrical datasets, such as turtle shells in certain views. Before you get here, you should already have landmarked your dataset, checked for correct landmarks, made your semi-landmarks, appended your curves to landmarks, deleted your overlappoing semi-landmarks, made your sliders file, done your semi-landmark superimposition, and double-checked that everything looked good. Then, you're going to have to reformat your data into the OSymm format from the IMP format (try using). I'm going to try to make a loop for that here, and hope that it works.
#OSymm format: the each column is a specimen, each row is a landmark. But X and Y are each a column. 
#IMP:sp1x1 sp1y1 sp1x2
#IMP:sp2x1 sp2y1 sp2x2
#IMP:sp3x1 sp3y1 sp3x2
#IMP:sp4x1 sp4y1 sp4x2

#OS: sp1x1 sp1y1 sp2x1 sp2y1
#OS: sp1x2 sp1y2 sp2x2 sp2y2
#OS: sp1x3 sp1y3 sp2x3 sp2y3
#OS: sp1x4 sp1y4 sp2x4 sp2y4

#7 October 2016: It's 4 years since I last used this script. Below are some of the first lines of code and loops I ever wrote. They suck. But not more than the time it would take to fix them. 

#load dependencies
library(geomorph)
library(dplyr)

# datadir<-"C:/Users/N.S/Dropbox/Documents/research/turtles/Thesis/Vitek_YR_PublishTerrapene/data/tps-intermediate"
datadir <- "D:/Dropbox/Documents/research/turtles/Thesis/Vitek_YR_PublishTerrapene/data/tps-intermediate"

setwd(datadir)
namefile<-read.csv("boxturtle_dor_images.csv",header=FALSE)
raw<-readland.tps(file="boxturtle_dor_align.tps")
str(raw)
#make landmarks into 2D array
IMP<-two.d.array(raw)
#Headers should be "RawCoord#" (e.g., RawCoord1,RawCoord2,RawCoord3, etc.) for historical reasons :-P
colnames(IMP)<-paste("RawCoord",seq(from=1,to=ncol(IMP),by=1),sep="")
dim(IMP)
#Then, transpose, so that our rows become columns and vice versa.
tr.IMP<-t(IMP)

#Now I want to take every odd row and make a list of that, and then every even row and make a list of that. The odd rows are my x's and the even rows are my y's. To do that, I'll need a blank vector that is my number of columns. 
dim(tr.IMP)
land<-1:nrow(tr.IMP)

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
# Go in and put the NA's here in R instead of opening your exporting your new table, opening it in Excel, and replacing missing landmarks with NA's there. Hopefully you wrote down which landmarks are missing in a text document somewhere that you can easily open now. 
# I think what you'd do is something like the following line for each missing landmark, 
#if n is your specimen number and 
#k is your landmark number:
#OS_Format[k,2n-1:2n]<-NA
#For example, in specimen 1 landmark 6 is missing. You'd use OS_Format[6,1:2].

###  ENTER IN THE NUMBER OF REPLACEMENTS YOU NEED IN THE CORRECT PLACE FOR YOUR DATASET.

#Fossil Datest, x = landmark, y = specimen * [2n-1:2n]
OS_Format[31,21:22]<-NA
OS_Format[32,21:22]<-NA
OS_Format[33,21:22]<-NA
OS_Format[34,21:22]<-NA
OS_Format[35,21:22]<-NA
OS_Format[36,21:22]<-NA
OS_Format[22,23:24]<-NA
OS_Format[23,23:24]<-NA
OS_Format[24,23:24]<-NA
OS_Format[40,23:24]<-NA
OS_Format[41,23:24]<-NA
OS_Format[19,47:48]<-NA
OS_Format[20,47:48]<-NA
OS_Format[22,47:48]<-NA
OS_Format[23,47:48]<-NA
OS_Format[24,47:48]<-NA
OS_Format[40,47:48]<-NA
OS_Format[41,47:48]<-NA
OS_Format[42,47:48]<-NA
OS_Format[43,47:48]<-NA
OS_Format[44,47:48]<-NA
OS_Format[8,51:52]<-NA
OS_Format[9,51:52]<-NA
OS_Format[10,51:52]<-NA
OS_Format[11,51:52]<-NA
OS_Format[12,51:52]<-NA
OS_Format[13,51:52]<-NA
OS_Format[30,51:52]<-NA
OS_Format[31,51:52]<-NA
OS_Format[32,51:52]<-NA
OS_Format[14,57:58]<-NA
OS_Format[15,57:58]<-NA
OS_Format[16,57:58]<-NA
OS_Format[17,57:58]<-NA
OS_Format[18,57:58]<-NA
OS_Format[38,57:58]<-NA
OS_Format[39,57:58]<-NA
OS_Format[10,59:06]<-NA
OS_Format[11,59:60]<-NA
OS_Format[23,69:70]<-NA
OS_Format[24,69:70]<-NA
OS_Format[40,69:70]<-NA
OS_Format[25,81:82]<-NA
OS_Format[26,81:82]<-NA
OS_Format[47,81:82]<-NA
OS_Format[48,81:82]<-NA
OS_Format[49,81:82]<-NA



#Size data: original numbers + 43*2 = +86

OS_Format[8,87:88]<-NA
OS_Format[10,87:88]<-NA
OS_Format[19,89:90]<-NA
OS_Format[21,89:90]<-NA
OS_Format[23,89:90]<-NA
OS_Format[16,91:92]<-NA
OS_Format[18,91:92]<-NA
OS_Format[8,97:98]<-NA
OS_Format[40,97:98]<-NA
OS_Format[41,97:98]<-NA
OS_Format[42,97:98]<-NA
OS_Format[43,97:98]<-NA
OS_Format[44,97:98]<-NA
OS_Format[45,97:98]<-NA
OS_Format[46,97:98]<-NA
OS_Format[47,97:98]<-NA
OS_Format[48,97:98]<-NA
OS_Format[49,97:98]<-NA
OS_Format[19,99:100]<-NA
OS_Format[20,99:100]<-NA
OS_Format[21,99:100]<-NA
OS_Format[22,99:100]<-NA
OS_Format[23,99:100]<-NA
OS_Format[14,101:102]<-NA
OS_Format[15,101:102]<-NA
OS_Format[16,101:102]<-NA
OS_Format[17,101:102]<-NA
OS_Format[10,103:104]<-NA
OS_Format[11,103:104]<-NA
OS_Format[12,103:104]<-NA
OS_Format[14,103:104]<-NA
OS_Format[25,245:246]<-NA
OS_Format[40,245:246]<-NA
OS_Format[41,245:246]<-NA
OS_Format[42,245:246]<-NA
OS_Format[43,245:246]<-NA
OS_Format[44,245:246]<-NA
OS_Format[45,245:246]<-NA
OS_Format[46,245:246]<-NA
OS_Format[47,245:246]<-NA
OS_Format[48,245:246]<-NA
OS_Format[49,245:246]<-NA

#ssp data: original numbers + 144*2 = +288 [2n-1:2n]
OS_Format[18,289:290]<-NA
OS_Format[14,339:340]<-NA
OS_Format[15,339:340]<-NA
OS_Format[14,421:422]<-NA
OS_Format[16,421:422]<-NA
OS_Format[8,457:458]<-NA
OS_Format[9,457:458]<-NA
OS_Format[10,457:458]<-NA
OS_Format[11,457:458]<-NA
OS_Format[12,457:458]<-NA
OS_Format[13,457:458]<-NA
OS_Format[25,457:458]<-NA
OS_Format[27,457:458]<-NA
OS_Format[12,491:492]<-NA
OS_Format[13,491:492]<-NA
OS_Format[14,493:494]<-NA
OS_Format[15,493:494]<-NA
OS_Format[16,493:494]<-NA
OS_Format[17,493:494]<-NA
OS_Format[10,493:494]<-NA
OS_Format[12,493:494]<-NA
OS_Format[14,493:494]<-NA
OS_Format[35,499:500]<-NA
OS_Format[36,499:500]<-NA
OS_Format[37,499:500]<-NA
OS_Format[38,499:500]<-NA
OS_Format[39,499:500]<-NA
OS_Format[8,509:510]<-NA
OS_Format[9,509:510]<-NA
OS_Format[10,509:510]<-NA
OS_Format[27,521:522]<-NA
OS_Format[18,561:562]<-NA
OS_Format[19,557:558]<-NA
OS_Format[8,633:634]<-NA
OS_Format[27,637:638]<-NA

#holo data: original specimen numbers + 345  [ + [2n-1:2n]]
hcount<-345
h1<-c(1,21,22)	 # - n1 lm 21, 22
h2<-c(2,8,9)	 # - n2 lm 8, 9
h3<-c(3,21,22)	 # - n3 lm 21, 22
h5<-c(5,12,13,14,15,16,17,18,30:39)	 # - n5 lm 12, 13, 14, 15, 16, 17, 18, the left semilandmark curve
OS_Format[h1[-1],(c((h1[1])*2-1,(h1[1])*2)+hcount)]<-NA
OS_Format[h2[-1],(c((h2[1])*2-1,(h2[1])*2)+hcount)]<-NA
OS_Format[h3[-1],(c((h3[1])*2-1,(h3[1])*2)+hcount)]<-NA
OS_Format[h5[-1],(c((h5[1])*2-1,(h5[1])*2)+hcount)]<-NA

#Scoahuila data: original numbers + 356  [+(n*2) + [2n-1:2n]]
scount<-356
s2<-c(2,12,13)	 # -n2: lm 12, 13
s3<-c(3,18)	 # -n3: lm 18
s6<-c(6,23,24)	 # -n6: lm 23, 24
s7<-c(7,25,26,27,28,29)	 # -n7: lm 25, 26, 27, 28, 29
s9<-c(9,25:29,40:49)	 # -n9: lm 25, 26, 27, 27, 29, the last 10 landmarks (semilandmark curve)
s10<-c(10,8:11)	 # -n10 lm 8, 9, 10, 11
s11<-c(11,23,24)	 # -n11 lm 23, 24
s12<-c(12,14,21)	 # -n12 lm 14, 21
s34<-c(34,10,11)	 # -n34 lm 10, 11
OS_Format[s2[-1],(c((s2[1])*2-1,(s2[1])*2)+scount)]<-NA
OS_Format[s3[-1],(c((s3[1])*2-1,(s3[1])*2)+scount)]<-NA
OS_Format[s6[-1],(c((s6[1])*2-1,(s6[1])*2)+scount)]<-NA
OS_Format[s7[-1],(c((s7[1])*2-1,(s7[1])*2)+scount)]<-NA
OS_Format[s9[-1],(c((s9[1])*2-1,(s9[1])*2)+scount)]<-NA
OS_Format[s10[-1],(c((s10[1])*2-1,(s10[1])*2)+scount)]<-NA
OS_Format[s11[-1],(c((s11[1])*2-1,(s11[1])*2)+scount)]<-NA
OS_Format[s12[-1],(c((s12[1])*2-1,(s12[1])*2)+scount)]<-NA
OS_Format[s34[-1],(c((s34[1])*2-1,(s34[1])*2)+scount)]<-NA

XX<-OS_Format

Fixed<-list()
j<-2
for (i in 1:nrow(IMP)){
  q<-i*2-1
  X <- as.matrix(XX[,q:j])
  #Now, you've got to specify which coordinates are midline, right, and left. If you want to do it in your head, the name for each coordinate with landmark number n is 2n-1. Currently, setup is for the dorsal view. 
  midline <- c("RawCoord1","RawCoord3", "RawCoord5","RawCoord7","RawCoord9","RawCoord11","RawCoord13")
  right <- c("RawCoord15","RawCoord17","RawCoord19","RawCoord21","RawCoord23","RawCoord25","RawCoord27","RawCoord29","RawCoord31","RawCoord33","RawCoord35","RawCoord59","RawCoord61","RawCoord63","RawCoord65","RawCoord67","RawCoord69","RawCoord71","RawCoord73","RawCoord75","RawCoord77")
  left <- c("RawCoord37","RawCoord39","RawCoord41","RawCoord43","RawCoord45","RawCoord47","RawCoord49","RawCoord51","RawCoord53","RawCoord55","RawCoord57","RawCoord79","RawCoord81","RawCoord83","RawCoord85","RawCoord87","RawCoord89","RawCoord91","RawCoord93","RawCoord95","RawCoord97")
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

#Then, you will want to sort the rows back in their original order. 
#You can't just use the sort function because it won't sort the numbers properly, 
#so I'm going to it by hand here. 18 August 2012: set up for dorsal view.
row.names(Fixed_Missing)
sortFM1<-Fixed_Missing[1:18,]
sortFM2<-Fixed_Missing[19:28,]
sortFM3<-Fixed_Missing[29:39,]
sortFM4<-Fixed_Missing[40:49,]
sortFM<-rbind(sortFM1,sortFM3,sortFM2,sortFM4)
row.names(sortFM)

#Now, get table back into IMP format.
sea<-1:(nrow(IMP)*2)
exes.l<-list()
whys.l<-list()

###  CHANGE THE X IN 1:X TO YOUR NUMBER OF LANDMARKS
for (i in 1:nrow(Fixed_Missing)) {
	exes.l[[i]]<-sortFM[i,seq.int(1L,length(sea),by=2)]
	whys.l[[i]]<-sortFM[i,seq.int(2L,length(sea),by=2)]
}

NewTable<-list()
j<-1

###  CHANGE THE X IN 1:X TO YOUR NUMBER OF LANDMARKS
for (i in 1:nrow(Fixed_Missing)){
	NewTable[[j]]<-exes.l[[i]]
	j<-j+1
	NewTable[[j]]<-whys.l[[i]]
	j<-j+1
}

Back<-data.frame(matrix(unlist(NewTable),ncol=98))
Back2<-as.matrix(Back)
dim(Back2)
head(Back2)

x<-namefile[,2]
m<-gregexpr("dorsal\\\\.*_[0-9]",x,perl=TRUE)
shapenames<-regmatches(x,m) %>% unlist %>% gsub("dorsal\\\\(.*)_[0-9]","\\1",.,perl=TRUE)
rownames(Back2)<-shapenames
write.table(Back2, "boxturtle_dor_IMP.txt",row.names=FALSE,col.names=FALSE,sep="  ")
write.csv(Back2, "boxturtle_dor_morphoj.csv",row.names=TRUE,col.names=FALSE)

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
symmdat<-read.table("boxturtle_dor_symm.txt",sep = "\t",header=TRUE)

#take out the coordinates of one side. 
#for dorsal: delete Symm 37-58, 79-98
half<-symmdat[-1,-(c(37:58,79:98)+1)]

# superimpose landmarks
# sh<-half[,-1] %>% as.matrix %>%  arrayspecs(.,ncol(.)/2,2)  %>% gpagen #2 for 2 dimensions
# sh<-two.d.array(sh$coords)
# 
# half[,-1]<-sh

colnames(half)<-c("id",paste("DorCoord",seq(from=1,to=(ncol(half)-1),by=1),sep=""))
write.csv(half,"boxturtle_dor_coords.csv",row.names=FALSE)
dim(half)
head(half)


# tps1<-arrayspecs(Back2,ncol(Back2)/2,2)
# dim(tps1)
# plot(tps1[,,1],asp=1)
# dor_pairs<-matrix(c(8,19,
#                     9,20,
#                     10,21,
#                     11,22,
#                     12,23,
#                     13,24,
#                     14,25,
#                     15,26,
#                     16,27,
#                     17,28,
#                     18,29,
#                     30,40,
#                     31,41,
#                     32,42,
#                     33,43,
#                     34,44,
#                     35,45,
#                     36,46,
#                     37,47,
#                     38,48,
#                     39,49),ncol=2,byrow=TRUE)
# dor_pairs
# bilat<-bilat.symmetry(tps1,ind=c(1:dim(tps1)[3]),object.sym=TRUE,
#                       land.pairs=dor_pairs)

