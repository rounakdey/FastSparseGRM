removeHigherDegree<-function(kin0,degree)
{
	setUN<-which(kin0$InfType=="UN")
	if(length(setUN)>0)	kin0<-kin0[-setUN,]
	if(degree==3)
	{
		set3<-which(kin0$InfType=="4th")
		if(length(set3)>0)	kin0<-kin0[-set3,]
	}
	if(degree==2)
	{
		set2<-which(kin0$InfType=="3rd")
		if(length(set2)>0)	kin0<-kin0[-set2,]
	}
	if(degree==1)
	{
		set1<-which(kin0$InfType=="2nd")
		if(length(set1)>0)	kin0<-kin0[-set1,]
	}
	return(kin0)
}

getRange<-function(nObj, nPart, index)
{
	chunkdiv=nObj %/% nPart
	chunkrem=nObj %% nPart;
	IndRange=rep(0,2)

	if(index<=(nPart-chunkrem))
	{
		IndRange[1]=(index-1)*chunkdiv+1;
		IndRange[2]=index*chunkdiv;
	} else {
		IndRange[1]=(nPart-chunkrem)*chunkdiv+(index-nPart+chunkrem-1)*(chunkdiv+1)+1;
		IndRange[2]=(nPart-chunkrem)*chunkdiv+(index-nPart+chunkrem)*(chunkdiv+1);
	}
	
	return(IndRange)
}


getDivergence<-function(file.seg,prefix.in,num_threads,divThresh=-2^5.5,degree=4,nRandomSNPs=10000,prefix.out)
{

fam <- read.table(paste0(prefix.in,".fam"),header=F)
colnames(fam)<-c("FID","IID","PID","MID","Sex","Pheno")
nSamp<-nrow(fam)

print("Calculating ancestry divergence for the samples in the seg file")
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

indexInkin0<-which(paste0(fam$FID,"_",fam$IID) %in% unique(c(kin0$ID1,kin0$ID2)))-1

bim <- read.table(paste0(prefix.in,".bim"),header=F)
nSNPs<-nrow(bim)
rm(bim)
if(nRandomSNPs>nSNPs | nRandomSNPs==0)	nRandomSNPs<-nSNPs
RandomSNPs<-sample(1:nSNPs,nRandomSNPs,replace=FALSE)
RandomSNPs<-RandomSNPs[order(RandomSNPs)]

RandomLags<-c(RandomSNPs[1],diff(RandomSNPs))-1

readfile(fname=paste0(prefix.in,".bed"),nSamp,RandomLags,num_threads)
deletedata()
divres<-calculateDivergence(indexInkin0,divThresh)
deletetranspose()
divout<-cbind(FID=fam$FID[indexInkin0+1],IID=fam$IID[indexInkin0+1],divship=divres)
write.table(divout,paste0(prefix.out,".div"),col.names=F,row.names=F,quote=F)

}

selectUnrel <- function(kin0,div,iids,id.include)
{
	unrels<-iids
	if(!("PropIBD" %in% colnames(kin0)))	colnames(kin0)[which(colnames(kin0)=="Kinship")]<-"PropIBD"
	kin0<-kin0[,c("ID1","ID2","PropIBD")]
	kin0rev<-kin0[,c("ID2","ID1","PropIBD")]
	colnames(kin0rev)<-c("ID1","ID2","PropIBD")
	kin0<-rbind(kin0,kin0rev)
	reltab<-aggregate(PropIBD ~ ID1, kin0, length)
	colnames(reltab)<-c("ID","NRel")
	IDincluded<-which(reltab$ID %in% id.include)
	if(length(IDincluded)>0)	reltab<-reltab[-IDincluded,]
	reltab<-merge(reltab,div,all.x=TRUE,all.y=FALSE,by="ID")
	if(length(which(is.na(reltab)))>0)
	{
		print("file.div does not have all the IDs present in file.seg")
		return(0)
	}
	TotalKinship<-aggregate(PropIBD ~ ID1,kin0,sum)
	colnames(TotalKinship)<-c("ID","TotalKinship")
	reltab<-merge(reltab,TotalKinship,all=FALSE,by="ID")
	reltab<-reltab[order(-reltab$NRel,reltab$Div,reltab$TotalKinship),]

	ntot<-nrow(reltab)
	lastProg<-0

	repeat
	{
		if(reltab$NRel[1]!=max(reltab$NRel))	reltab<-reltab[order(-reltab$NRel,reltab$Div,reltab$TotalKinship),]

		removeid<-reltab$ID[1]
		linked.ids<-which(reltab$ID %in% kin0$ID2[which(kin0$ID1==removeid)])
		reltab$NRel[linked.ids]<-reltab$NRel[linked.ids]-1
		reltab<-reltab[-c(1,which(reltab$NRel==0)),,drop=FALSE]
		unrels<-unrels[-which(unrels==removeid)]

		progress<-(ntot-nrow(reltab))*100/ntot
		currentProg=progress%/%5
		if(lastProg!=currentProg)	print(paste0(round(progress,0),"% completed"))
		lastProg<-currentProg

		if(nrow(reltab)==0) break
	}
	return(unrels)
}


extractUnrelated<-function(file.seg,file.div,prefix.in,degree=4,file.include,prefix.out)
{

fam <- read.table(paste0(prefix.in,".fam"),sep=" ",header=F)
colnames(fam)<-c("FID","IID","PID","MID","Sex","Pheno")

if(file.include!="")
{
	id.include<-read.table(file.include,header=F)
	if(ncol(id.include)==1)	id.include<-cbind(id.include,id.include)
} else {
	id.include=NULL
}

print("Extracting the set of unrelated samples")
kin0 <- read.table(file.seg, header = TRUE)
kin0<-removeHigherDegree(kin0,degree)
div <- read.table(file.div, header = FALSE)
colnames(div)<-c("ID","Div")

unrel_set <- selectUnrel(kin0,div,fam$IID,id.include)

unrel_set<-fam[which(fam$IID %in% unrel_set),c("FID","IID")]
write.table(unrel_set,paste0(prefix.out,".unrels"),sep=" ",col.names=F,row.names=F,quote=F)
}

