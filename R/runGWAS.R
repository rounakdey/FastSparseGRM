fitGWASNullModel<-function(prefix.fam,file.pheno,file.cov,pheno.col,cov.col,id.col="IID",file.GRM.RData,GRM.object.name="sGRM",GRM.ID.format="FID_IID",family,link,groups,method="REML",method.optim="AI",maxiter=500,tol=1e-5,taumin=1e-5,taumax=1e5,tauregion=10,verbose=TRUE,prefix.out)
{
	options(stringsAsFactors=F)
	if(link!="")
	{
		family<-eval(parse(text=family))(link=link)
	} else {
		family<-eval(parse(text=family))()
	}
	ID.format <- try(match.arg(GRM.ID.format, c("FID_IID","IID")))
	if(class(ID.format) == "try-error") stop("Error: GRM.ID.format should be one of the following: FID_IID, IID")

	fam<-read.table(paste0(prefix.fam,".fam"),header=F)
	famIDs<-paste0(fam[,1],"_",fam[,2])

	if(cov.col=="")
	{
		cov.col=NULL
	} else {
		cov.col<-trimws(unlist(strsplit(cov.col,"[,]")))
	}
	if(groups=="")
	{
		groups=NULL
	} else {
		groups<-trimws(unlist(strsplit(groups,"[,]")))
	}
	pheno.col<-trimws(unlist(strsplit(pheno.col,"[,]")))


	print("Reading Phenotypes")
	pheno<-read.table(file.pheno,header=T)
	pheno<-pheno[,which(colnames(pheno) %in% c(id.col,pheno.col,groups))]


	if(file.cov!="")
	{
		print("Reading Covariates")
		cov<-read.table(file.cov,header=T)
		cov<-cov[,which(colnames(cov) %in% c(id.col,cov.col,groups))]
		pheno<-merge(pheno,cov,by=id.col,all=FALSE)
	}

	phenoIDs<-pheno[,which(colnames(pheno) == id.col)]
	pheno<-pheno[which(phenoIDs %in% famIDs),]

	if(family$family=="binomial")
	{
		for(i in 1:length(pheno.col))
		{
			colid<-which(colnames(pheno)==pheno.col[i])
			phenocur<-pheno[,colid]
			nouniq<-length(unique(phenocur[which(!is.na(phenocur))]))
			if(nouniq!=2)
			{
				print(paste0(pheno.col," does not have two unique values, even though binomial family is selected. This will throw an error when GMMAT is called."))
			} else {
				tab<-sort(table(phenocur),decreasing=T)
				val_large<-names(tab)[1]
				val_small<-names(tab)[2]
				pheno[which(phenocur==val_large),colid]<-0
				pheno[which(phenocur==val_small),colid]<-1
			}
		}
	}

	covlist<-paste0(cov.col,collapse="+")

	print("Reading GRM")
	load(file.GRM.RData)
	GRM<-get(GRM.object.name)
	rm(list=GRM.object.name)
	if(GRM.ID.format=="IID")	colnames(GRM)<-rownames(GRM)<-paste0(rownames(GRM),"_",rownames(GRM))
	
	for(i in 1:length(pheno.col))
	{
		print(paste0("Fitting ",pheno.col[i]," started at ",Sys.time()))
		fixed<-paste0(pheno.col[i],"~",covlist)
		model.fit<-try({GMMAT::glmmkin(fixed=fixed,data=pheno,kins=GRM,id=id.col,random.slope=NULL,groups=groups,family=family,method=method,method.optim=method.optim,maxiter=maxiter,tol=tol,taumin=taumin,taumax=taumax,tauregion=tauregion,verbose=verbose)})
		if(class(model.fit) != "try-error")
		{
			save(model.fit,file=paste0(prefix.out,".",pheno.col[i],".Null.RData"))
			print(paste0("Fitting ",pheno.col[i]," completed at ",Sys.time()))
		} else {
			print(paste0("Re-fitting ",pheno.col[i]," using tau=0 at ",Sys.time()))
			model.fit<-try({GMMAT::glmmkin(fixed=fixed,data=pheno,kins=NULL,id=id.col,random.slope=NULL,groups=groups,family=family,method=method,method.optim=method.optim,maxiter=maxiter,tol=tol,taumin=taumin,taumax=taumax,tauregion=tauregion,verbose=verbose)})
			if(class(model.fit) != "try-error")
			{
				save(model.fit,file=paste0(prefix.out,".",pheno.col[i],".Null.RData"))
				print(paste0("Fitting ",pheno.col[i]," completed at ",Sys.time()))
			} else {
				write.table(model.fit[1],paste0(prefix.out,".",pheno.col[i],".err"),row.names=F,col.names=F,quote=F)
				print(paste0("Fitting ",pheno.col[i]," encountered error at ",Sys.time()))
			}
		}
	}
}
