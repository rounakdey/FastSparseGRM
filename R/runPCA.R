
drpca<-function(prefix.in,unrel_ids,p,num_threads,nd=1,L=20,no_iter=10)
{
options(stringsAsFactors=F)

print(paste("PCA started at",Sys.time()))

fam<-read.table(paste0(prefix.in,".fam"),header=F)
fam_ids<-paste0(fam[,1],"_",fam[,2])
unrelIndex<-which(fam_ids%in%unrel_ids)-1
relIndex<-c(0:(length(fam_ids)-1))[-(unrelIndex+1)]	#C array format
all_ids<-c(fam_ids[unrelIndex+1],fam_ids[relIndex+1])	#C array format
nunrel<-length(unrelIndex)
nrel<-length(relIndex)
print(paste(nrel,"relateds,",nunrel,"unrelateds, and",p,"variants"))
nsamp<-nunrel+nrel

ptm<-proc.time()
#Read Data and Calculate Allele Frequencies
readfile(fname=paste0(prefix.in,".bed"),nsamp,rep(0,p),num_threads)
setmeans(unrelIndex)

ptm1<-proc.time()
tmdiff<-round(ptm1[3]-ptm[3],0)
print(paste("Reading Completed in",tmdiff,"seconds"))

H<-NULL
init<-matrix(rnorm(nunrel*L),nunrel,L)	#Initialize eigenvector
rep<-1
repeat{
	x<-postmultiply(init,unrelIndex)
	y<-premultiply(x,unrelIndex)
	H<-cbind(H,x)
	for(i in 1:ncol(y))		y[,i]<-y[,i]/norm(as.matrix(y[,i]),type='f')
	print(paste(rep,"iterations completed"))
	if(rep==no_iter) break
	init=y
	rep=rep+1
}

ptm2<-proc.time()
tmdiff<-round(ptm2[3]-ptm1[3],0)
print(paste("Iterations Completed in",tmdiff,"seconds"))

init<-svd(H,nv=0)$u
rm(H)

T<-premultiply(init,unrelIndex)
svd.T<-svd(T,nu=nd,nv=0)
rm(T)

u.res<-svd.T$u
v.res<-postmultiply(u.res,unrelIndex)
for(j in 1:nd)	v.res[,j]<-v.res[,j]/svd.T$d[j]
d.res<-svd.T$d[1:nd]

ptm3<-proc.time()
tmdiff<-round(ptm3[3]-ptm2[3],0)
print(paste("PCA on the Unrelated Samples Completed in",tmdiff,"seconds"))

if(nrel>0)
{
	u.res.rel<-premultiply(v.res,relIndex)
	for(j in 1:nd)	u.res.rel[,j]<-u.res.rel[,j]/svd.T$d[j]
	u.res<-rbind(u.res,u.res.rel)
}

deletedata()
deletetranspose()

#Match IDs with the fam file
matchID<-match(fam_ids,all_ids)
fam.file<-read.table(paste0(prefix.in,".fam"),header=F)
u.res<-cbind(fam.file[,1:2],u.res[matchID,])

print(paste("PC scores calculated at",Sys.time()))
return(list(d=d.res,u=u.res,v=v.res))
}

runPCA<-function(prefix.in,file.unrels,no_pcs=20,num_threads,prefix.out,no_iter=10)
{

if(file.unrels=="")
{
	print("No file.unrels provides. All samples are considered unrelated")
	unrels <- read.table(paste0(prefix.in,".fam"),header=F)
	unrel_ids<-paste0(unrels[,1],"_",unrels[,2])
} else {
	unrels <- read.table(file.unrels,sep=" ",header=F)
	if(ncol(unrels)==1)	unrels<-cbind(unrels,unrels)
	unrel_ids<-paste0(unrels[,1],"_",unrels[,2])
}

bim<-read.table(paste0(prefix.in,".bim"),header=F)
p=nrow(bim)
rm(bim)

out<-drpca(prefix.in,unrel_ids,p,num_threads,nd=no_pcs,L=2*no_pcs,no_iter=no_iter)

write.table(out$d^2,paste(prefix.out,".eval",sep=""),row.names=F,col.names=F,quote=F)
write.table(out$v,paste(prefix.out,".evec",sep=""),row.names=F,col.names=F,quote=F)
write.table(out$u,paste(prefix.out,".score",sep=""),row.names=F,col.names=F,quote=F)
}

