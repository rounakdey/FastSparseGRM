calcSparseGRM<-function(prefix.in,file.score,file.train,file.seg,no_pcs=10,num_threads=0,prefix.out,degree=4,block.size=5000,max.related.block=5000)
{
options(stringsAsFactors=F)

if(file.train=="") file.train=file.score
if(max.related.block>2^16)
{
	max.related.block=65536
	print("Maximum related block size is too large, using 65536 as the maximum related block size")
}
print(paste("Sparse GRM calculation started at",Sys.time()))

score<-read.table(file.score,header=F)
if(no_pcs==0)	no_pcs=ncol(score)-2
fam<-read.table(paste0(prefix.in,".fam"),header=F)
train<-read.table(file.train,header=F)
trainIDs<-paste0(train[,1],"_",train[,2])
scoreIDs<-paste0(score[,1],"_",score[,2])
famIDs<-paste0(fam[,1],"_",fam[,2])
matchIDs<-match(famIDs,scoreIDs)
X<-as.matrix(cbind(1,score[matchIDs,c(2+1:no_pcs)]))
in.train<-which(famIDs %in% trainIDs)-1
X.in.train<-X[in.train+1,]
nullmat<-X.in.train%*%solve(t(X.in.train)%*%X.in.train)
bim<-read.table(paste0(prefix.in,".bim"),header=F)
p<-nrow(bim)
rm(bim)

kin0 <- read.table(file.seg, header = TRUE)
kin0<-removeHigherDegree(kin0,degree)
if(c("FID") %in% colnames(kin0))
{
	kin0$ID1<-paste0(kin0$FID,"_",kin0$ID1)
	kin0$ID2<-paste0(kin0$FID,"_",kin0$ID2)
} else {
	kin0$ID1<-paste0(kin0$FID1,"_",kin0$ID1)
	kin0$ID2<-paste0(kin0$FID2,"_",kin0$ID2)
}

gr<-igraph::graph_from_data_frame(kin0[,c("ID1","ID2")])
comp <- igraph::components(gr)
mem <- comp$membership
maxblock.seg<-max(comp$csize)

if(max(comp$csize)>max.related.block){
	if(!("PropIBD" %in% colnames(kin0)))
	{
		colnames(kin0)[which(colnames(kin0)=="Kinship")]<-"PropIBD"
		kin0$PropIBD=kin0$PropIBD*2
	}
	kin0kins<-kin0$PropIBD/2
	kin0<-kin0[,c("ID1","ID2")]
	print(paste0("Maximum related block size (",maxblock.seg,") is too large based on KING output. GRM will be estimated in two iterations."))
} else {
	blocks<-list()
	if(comp$no>0)
	{
		for(i in 1:comp$no)
		{
			ids<-names(mem[mem==i])
			blocks[[i]]<-t(combn(ids,2))
		}
	}
	kin0<-as.data.frame(do.call("rbind",blocks))
	colnames(kin0)<-c("ID1","ID2")
}
kinself<-data.frame(ID1=famIDs,ID2=famIDs)
kin0<-rbind(kinself,kin0)

print(paste("Estimating sparse GRM at",nrow(kin0),"coordinates"))

ptm<-proc.time()

readfile(fname=paste0(prefix.in,".bed"),length(famIDs),rep(0,p),num_threads)
setmeans(in.train,noscale=1)
beta<-postmultiply(nullmat,in.train)

#Calculate Kinship Here
Kinship<-matrix(0,nrow(kin0),2)
rep=0
repeat{
beginInd=rep*block.size
endInd=min((rep+1)*block.size,p)
print(c(beginInd,endInd))
Kinship_now=calcRel(beginInd, endInd, beta[beginInd:(endInd-1)+1,], X, match(kin0$ID1,famIDs)-1, match(kin0$ID2,famIDs)-1)
Kinship=Kinship+Kinship_now
rep=rep+1
if(endInd==p) break
}
Kinship=Kinship[,1]/Kinship[,2]

ptm1<-proc.time()
tmdiff<-round(ptm1[3]-ptm[3],0)
print(paste("Kinships calculated in",tmdiff,"seconds"))

kin0<-cbind(kin0,Kinship)
colnames(kin0)<-c("ID1","ID2","Kinship")

#Perform thresholding again
kin0noself<-kin0[-c(1:nrow(kinself)),]
kin1<-kin0noself[which(kin0noself$Kinship>=2^-(degree+1.5)),]

degree.reset=0
newthresh=2^-(degree+1.5)
repeat{
gr<-igraph::graph_from_data_frame(kin1)
comp <- igraph::components(gr)
mem <- comp$membership
maxblock.kin<-max(comp$csize)
if(maxblock.kin>max.related.block)
{
	kin1<-kin1[-which(kin1$Kinship==min(kin1$Kinship)),]
	newthresh<-min(kin1$Kinship)
} else {
	break
}
degree.reset=1
}
if(degree.reset==1) print(paste0("Maximum related block size too large, sparsity threshold modified to ",min(kin1$Kinship)))

if(maxblock.seg>max.related.block)
{
	blocks<-list()
	if(comp$no>0)
	{
		for(i in 1:comp$no)
		{
			ids<-names(mem[mem==i])
			blocks[[i]]<-t(combn(ids,2))
		}
	}
	kin1<-as.data.frame(do.call("rbind",blocks))
	colnames(kin1)<-c("ID1","ID2")

	#####difference between kin1 and kin0
	kin0joint.ids<-paste0(kin0noself$ID1,kin0noself$ID2)
	kin0joint.ids.rev<-paste0(kin0noself$ID2,kin0noself$ID1)
	kin1joint.ids<-paste0(kin1$ID1,kin1$ID2)
	inkin0<-c(which(kin1joint.ids %in% kin0joint.ids),which(kin1joint.ids %in% kin0joint.ids.rev))
	if(length(inkin0)>0)	kin1<-kin1[-inkin0,]

	print(paste("Estimating sparse GRM at",nrow(kin1),"more coordinates"))

	ptm<-proc.time()

	#Calculate Kinship Here
	Kinship<-matrix(0,nrow(kin1),2)
	rep=0
	repeat{
	beginInd=rep*block.size
	endInd=min((rep+1)*block.size,p)
	print(c(beginInd,endInd))
	Kinship_now=calcRel(beginInd, endInd, beta[beginInd:(endInd-1)+1,], X, match(kin1$ID1,famIDs)-1, match(kin1$ID2,famIDs)-1)
	Kinship=Kinship+Kinship_now
	rep=rep+1
	if(endInd==p) break
	}
	Kinship=Kinship[,1]/Kinship[,2]

	ptm1<-proc.time()
	tmdiff<-round(ptm1[3]-ptm[3],0)
	print(paste("Kinships calculated in",tmdiff,"seconds"))

	kin1<-cbind(kin1,Kinship)
	colnames(kin1)<-c("ID1","ID2","Kinship")
	kin0<-rbind(kin0,kin1)

	#Perform thresholding again
	kin0noself<-kin0[-c(1:nrow(kinself)),]
	kin1<-kin0noself[which(kin0noself$Kinship>=newthresh),]

	gr<-igraph::graph_from_data_frame(kin1)
	comp <- igraph::components(gr)
	mem <- comp$membership
}


singletons<-setdiff(famIDs,names(mem))

blocks<-list()
block.ids<-list()
if(comp$no>0)
{
	print(paste0(length(mem)," relatives in ",comp$no," clusters, largest cluster size ",max(comp$csize)))
	print(paste0(length(singletons)," singletons present"))
	for(i in 1:comp$no)
	{
		ids<-names(mem[mem==i])
		blocks[[i]]<-t(combn(ids,2))
		block.ids[[i]]<-ids
	}
}
kin1<-as.data.frame(do.call("rbind",blocks))
colnames(kin1)<-c("ID1","ID2")
kin0rev<-kin0[-c(1:nrow(kinself)),c(2,1,3)]
colnames(kin0rev)<-c("ID1","ID2","Kinship")
kin0<-rbind(kin0,kin0rev)
kin1<-rbind(kinself,kin1)

if(length(singletons)>0)        block.ids[[comp$no + 1]] <- singletons

GRMIDs<-unlist(block.ids)
kin0<-merge(kin0,kin1,all=FALSE)

sGRM <- sparseMatrix(i=match(kin0$ID1,GRMIDs), j=match(kin0$ID2,GRMIDs), x=kin0$Kinship,symmetric=TRUE)
colnames(sGRM)<-rownames(sGRM)<-GRMIDs
print(paste("Sparse GRM calculation completed at",Sys.time()))

save(sGRM,file=paste(prefix.out,".sGRM.RData",sep=""))

#Save Kinships in KING format
kin0<-kin0[-which(kin0$ID1==kin0$ID2),]
kinID1<-fam[match(kin0$ID1,famIDs),1:2]
kinID2<-fam[match(kin0$ID2,famIDs),1:2]
InfType<-rep("UN",nrow(kin0))
InfType[which(kin0$Kinship>=2^-5.5 & kin0$Kinship<2^-4.5)]<-"4th"
InfType[which(kin0$Kinship>=2^-4.5 & kin0$Kinship<2^-3.5)]<-"3rd"
InfType[which(kin0$Kinship>=2^-3.5 & kin0$Kinship<2^-2.5)]<-"2nd"
InfType[which(kin0$Kinship>=2^-2.5 & kin0$Kinship<2^-1.5)]<-"1st"
InfType[which(kin0$Kinship>=2^-1.5)]<-"Dup/MZ"

kin0<-cbind(kinID1,kinID2,kin0$Kinship,InfType)
colnames(kin0)<-c("FID1","ID1","FID2","ID2","Kinship","InfType")
write.table(kin0,paste(prefix.out,".kins",sep=""),row.names=F,col.names=T,quote=F)
deletedata()
deletetranspose()

}

