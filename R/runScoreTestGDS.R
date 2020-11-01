runScoreTestGDS<-function(file.null.Rdata,Null.obj.name,prefix.gds,is.dosage,only.dosage,prefix.out,file.samp,miss.cutoff,minMAF,maxMAF,missing.method,nSNPsperbatch,genomem,num_threads=1)
{
	options(stringsAsFactors=F)
	load(file.null.RData)
	model.fit=get(Null.obj.name)
	MAF.range=c(minMAF,maxMAF)
	if(file.samp!="")
	{
		samp<-read.table(file.samp,header=F)
		if(ncol(samp)==1) sample.id<-samp[,1] else sample.id<-paste0(samp[,1],"_",samp[,2])
		null.model.id<-model.fit$id_include
		sample.id.order<-match(sample.id,null.model.id)
		sample.id.order[which(is.na(sample.id.order))]<-0
	} else {
		sample.id.order=NULL
	}
	if(genomem>0)	nSNPsperbatch=min(10000,max(floor(genomem*(2^20)/(8*length(model.fit$id_include))),2))

	if(num_threads>1)	Sys.setenv(MKL_NUM_THREADS = 1)

	if(is.dosage & only.dosage)
	{
		glmm.score.gds.nogeno(obj=model.fit, infile=paste0(prefix.gds,".gds"), outfile=paste0(prefix.out,".out"), select = sample.id.order, MAF.range = MAF.range, nperbatch = nSNPsperbatch, ncores = num_threads)
	} else {
		glmm.score.gds(obj=model.fit, infile=paste0(prefix.gds,".gds"), outfile=paste0(prefix.out,".out"), select = sample.id.order, MAF.range = MAF.range, miss.cutoff = miss.cutoff, missing.method = missing.method, nperbatch = nSNPsperbatch, ncores = num_threads,is.dosage=is.dosage)
	}
}


glmm.score.gds <- function(obj, infile, outfile, center = T, select = NULL, MAF.range = c(1e-7, 0.5), miss.cutoff = 1, missing.method = "impute2mean", nperbatch = 100, ncores = 1, is.dosage)
{
	if(class(obj) != "glmmkin") stop("Error: obj must be a class glmmkin object!")
	if(any(duplicated(obj$id_include))) {
		J <- Matrix(sapply(unique(obj$id_include), function(x) 1*(obj$id_include==x)), sparse = TRUE)
		res <- as.vector(as.matrix(crossprod(J, obj$scaled.residuals)))
		if(!is.null(obj$P)) obj$P <- as.matrix(crossprod(J, crossprod(obj$P, J)))
		else {
			obj$Sigma_iX <- crossprod(J, obj$Sigma_iX)
			obj$Sigma_i <- forceSymmetric(crossprod(J,crossprod(obj$Sigma_i,J)))
			obj$Sigma_i <- Matrix(obj$Sigma_i, sparse = TRUE)
		}
		rm(J)
	} else res <- obj$scaled.residuals
	miss.method <- try(match.arg(missing.method, c("impute2mean", "omit")))
	if(class(miss.method) == "try-error") stop("Error: missing.method should be one of the following: impute2mean, omit!")
	miss.method <- substr(miss.method, 1, 1)

	        if(Sys.info()["sysname"] == "Windows" && ncores > 1) {
	 		warning("The package doMC is not available on Windows... Switching to single thread...")
			ncores <- 1
		}
	        ncores <- min(c(ncores, parallel::detectCores(logical = TRUE)))
	        gds <- SeqArray::seqOpen(infile)
		sample.id <- SeqArray::seqGetData(gds, "sample.id")
		if(is.null(select)) {
			if(any(is.na(match(unique(obj$id_include), sample.id)))) warning("Check your data... Some id_include in obj are missing in sample.id of infile!")
			select <- match(sample.id, unique(obj$id_include))
			select[is.na(select)] <- 0
			if(all(select == 0)) stop("Error: id_include in obj does not match sample.id in infile!")
		}
		if(length(select) != length(sample.id)) stop("Error: number of individuals in select does not match infile!")
		select2 <- select[select > 0]
		if(any(duplicated(select2)) || max(select2) > length(unique(obj$id_include))) stop("Error: select is a vector of orders, individuals not in obj should be coded 0!")
		res <- res[select2]
		if(!is.null(obj$P)) obj$P <- obj$P[select2, select2]
		else {
			obj$Sigma_iX <- obj$Sigma_iX[select2, , drop = FALSE]
			obj$Sigma_i <- obj$Sigma_i[select2, select2]
		}
		rm(select2)
    		variant.idx.all <- SeqArray::seqGetData(gds, "variant.id")
	        SeqArray::seqClose(gds)
		p.all <- length(variant.idx.all)
		if(ncores > 1) {
			require(doMC)
			doMC::registerDoMC(cores = ncores)
			p.percore <- (p.all-1) %/% ncores + 1
			n.p.percore_1 <- p.percore * ncores - p.all
			foreach(b = 1:ncores, .inorder=FALSE, .options.multicore = list(preschedule = FALSE, set.seed = FALSE)) %dopar% {
				variant.idx <- if(b <= n.p.percore_1) variant.idx.all[((b-1)*(p.percore-1)+1):(b*(p.percore-1))] else variant.idx.all[(n.p.percore_1*(p.percore-1)+(b-n.p.percore_1-1)*p.percore+1):(n.p.percore_1*(p.percore-1)+(b-n.p.percore_1)*p.percore)]
				p <- length(variant.idx)
				gds <- SeqArray::seqOpen(infile)
				SeqArray::seqSetFilter(gds, sample.id = sample.id[select > 0], verbose = FALSE)
				rm(sample.id); rm(select)
				nbatch.flush <- (p-1) %/% 100000 + 1
				ii <- 0
				for(i in 1:nbatch.flush) {
		                        gc()
		        		tmp.variant.idx <- if(i == nbatch.flush) variant.idx[((i-1)*100000+1):p] else variant.idx[((i-1)*100000+1):(i*100000)]
					SeqArray::seqSetFilter(gds, variant.id = tmp.variant.idx, verbose = FALSE)
					MISSRATE <- SeqVarTools::missingGenotypeRate(gds, margin = "by.variant")
					AF <- 1 - SeqVarTools::alleleFrequency(gds)
					include <- (MISSRATE <= miss.cutoff & ((AF >= MAF.range[1] & AF <= MAF.range[2]) | (AF >= 1-MAF.range[2] & AF <= 1-MAF.range[1])))
					if(sum(include) == 0) next
					ii <- ii + 1
					tmp.variant.idx <- tmp.variant.idx[include]
					tmp.p <- length(tmp.variant.idx)
					SeqArray::seqSetFilter(gds, variant.id = tmp.variant.idx, verbose = FALSE)
					SNP <- SeqArray::seqGetData(gds, "annotation/id")
					SNP[SNP == ""] <- NA
					out <- data.frame(SNP = SNP, CHR = SeqArray::seqGetData(gds, "chromosome"), POS = SeqArray::seqGetData(gds, "position"))
					rm(SNP)
					alleles.list <- strsplit(SeqArray::seqGetData(gds, "allele"), ",")
					out$REF <- unlist(lapply(alleles.list, function(x) x[1]))
					out$ALT <- unlist(lapply(alleles.list, function(x) paste(x[-1], collapse=",")))
					out$MISSRATE <- MISSRATE[include]
					out$AF <- AF[include]
					rm(alleles.list, include)
					tmp.out <- lapply(1:((tmp.p-1) %/% nperbatch + 1), function(j) {
						tmp2.variant.idx <- if(j == (tmp.p-1) %/% nperbatch + 1) tmp.variant.idx[((j-1)*nperbatch+1):tmp.p] else tmp.variant.idx[((j-1)*nperbatch+1):(j*nperbatch)]
						SeqArray::seqSetFilter(gds, variant.id = tmp2.variant.idx, verbose = FALSE)
						if(is.dosage)	geno <- SeqVarTools::imputedDosage(gds, dosage.field="DS", use.names = FALSE)
						else geno <- SeqVarTools::altDosage(gds, use.names = FALSE)
        					N <- nrow(geno) - colSums(is.na(geno))
						if(center) geno <- scale(geno, scale = FALSE)
						miss.idx <- which(is.na(geno))
						if(length(miss.idx)>0) {
							geno[miss.idx] <- if(!center & missing.method == "impute2mean") colMeans(geno, na.rm = TRUE)[ceiling(miss.idx/nrow(geno))] else 0
						}
						SCORE <- as.vector(crossprod(geno, res))
						if(!is.null(obj$P)) VAR <- diag(crossprod(geno, crossprod(obj$P, geno)))
						else {
							GSigma_iX <- crossprod(geno, obj$Sigma_iX)
							VAR <- diag(crossprod(geno, crossprod(obj$Sigma_i, geno)) - tcrossprod(GSigma_iX, tcrossprod(GSigma_iX, obj$cov)))
						}
						PVAL <- ifelse(VAR>0, pchisq(SCORE^2/VAR, df=1, lower.tail=FALSE), NA)
						return(rbind(N, SCORE, VAR, PVAL))
					})
					tmp.out <- matrix(unlist(tmp.out), ncol = 4, byrow = TRUE, dimnames = list(NULL, c("N", "SCORE", "VAR", "PVAL")))
					out <- cbind(out[,c("SNP","CHR","POS","REF","ALT")], tmp.out[,"N"], out[,c("MISSRATE","AF")], tmp.out[,c("SCORE","VAR","PVAL")])
					names(out)[6] <- "N"
					rm(tmp.out)
					if(b == 1) {
				     	        write.table(out, outfile, quote=FALSE, row.names=FALSE, col.names=(ii == 1), sep="\t", append=(ii > 1), na=".")
					} else {
			     	                write.table(out, paste0(outfile, "_tmp.", b), quote=FALSE, row.names=FALSE, col.names=FALSE, sep="\t", append=(ii > 1), na=".")
					}
					rm(out)
				}
	        		SeqArray::seqClose(gds)
			}
			for(b in 2:ncores) {
			      	system(paste0("cat ", outfile, "_tmp.", b, " >> ", outfile))
				unlink(paste0(outfile, "_tmp.", b))
			}
		} else {
			variant.idx <- variant.idx.all
			rm(variant.idx.all)
			p <- length(variant.idx)
			gds <- SeqArray::seqOpen(infile)
			SeqArray::seqSetFilter(gds, sample.id = sample.id[select > 0], verbose = FALSE)
			rm(sample.id); rm(select)
			nbatch.flush <- (p-1) %/% 100000 + 1
			ii <- 0
			for(i in 1:nbatch.flush) {
		                gc()
		        	tmp.variant.idx <- if(i == nbatch.flush) variant.idx[((i-1)*100000+1):p] else variant.idx[((i-1)*100000+1):(i*100000)]
				tmp.p <- length(tmp.variant.idx)
				SeqArray::seqSetFilter(gds, variant.id = tmp.variant.idx, verbose = FALSE)
				MISSRATE <- SeqVarTools::missingGenotypeRate(gds, margin = "by.variant")
				AF <- 1 - SeqVarTools::alleleFrequency(gds)
				include <- (MISSRATE <= miss.cutoff & ((AF >= MAF.range[1] & AF <= MAF.range[2]) | (AF >= 1-MAF.range[2] & AF <= 1-MAF.range[1])))
				if(sum(include) == 0) next
				ii <- ii + 1
				tmp.variant.idx <- tmp.variant.idx[include]
				tmp.p <- length(tmp.variant.idx)
				SeqArray::seqSetFilter(gds, variant.id = tmp.variant.idx, verbose = FALSE)
				SNP <- SeqArray::seqGetData(gds, "annotation/id")
				SNP[SNP == ""] <- NA
				out <- data.frame(SNP = SNP, CHR = SeqArray::seqGetData(gds, "chromosome"), POS = SeqArray::seqGetData(gds, "position"))
				rm(SNP)
				alleles.list <- strsplit(SeqArray::seqGetData(gds, "allele"), ",")
				out$REF <- unlist(lapply(alleles.list, function(x) x[1]))
				out$ALT <- unlist(lapply(alleles.list, function(x) paste(x[-1], collapse=",")))
				out$MISSRATE <- MISSRATE[include]
				out$AF <- AF[include]
				rm(alleles.list, include)
				tmp.out <- lapply(1:((tmp.p-1) %/% nperbatch + 1), function(j) {
					tmp2.variant.idx <- if(j == (tmp.p-1) %/% nperbatch + 1) tmp.variant.idx[((j-1)*nperbatch+1):tmp.p] else tmp.variant.idx[((j-1)*nperbatch+1):(j*nperbatch)]
					SeqArray::seqSetFilter(gds, variant.id = tmp2.variant.idx, verbose = FALSE)
					if(is.dosage)	geno <- SeqVarTools::imputedDosage(gds, dosage.field="DS", use.names = FALSE)
					else geno <- SeqVarTools::altDosage(gds, use.names = FALSE)
       					N <- nrow(geno) - colSums(is.na(geno))
					if(center) geno <- scale(geno, scale = FALSE)
					miss.idx <- which(is.na(geno))
					if(length(miss.idx)>0) {
						geno[miss.idx] <- if(!center & missing.method == "impute2mean") colMeans(geno, na.rm = TRUE)[ceiling(miss.idx/nrow(geno))] else 0
					}
					SCORE <- as.vector(crossprod(geno, res))
					if(!is.null(obj$P)) VAR <- diag(crossprod(geno, crossprod(obj$P, geno)))
					else {
						GSigma_iX <- crossprod(geno, obj$Sigma_iX)
						VAR <- diag(crossprod(geno, crossprod(obj$Sigma_i, geno)) - tcrossprod(GSigma_iX, tcrossprod(GSigma_iX, obj$cov)))
					}
					PVAL <- ifelse(VAR>0, pchisq(SCORE^2/VAR, df=1, lower.tail=FALSE), NA)
					return(rbind(N, SCORE, VAR, PVAL))
				})
				tmp.out <- matrix(unlist(tmp.out), ncol = 4, byrow = TRUE, dimnames = list(NULL, c("N", "SCORE", "VAR", "PVAL")))
				out <- cbind(out[,c("SNP","CHR","POS","REF","ALT")], tmp.out[,"N"], out[,c("MISSRATE","AF")], tmp.out[,c("SCORE","VAR","PVAL")])
				names(out)[6] <- "N"
				rm(tmp.out)
			     	write.table(out, outfile, quote=FALSE, row.names=FALSE, col.names=(ii == 1), sep="\t", append=(ii > 1), na=".")
				rm(out)
			}
	        	SeqArray::seqClose(gds)
		}
		return(invisible(NULL))

}


glmm.score.gds.nogeno <- function(obj, infile, outfile, center = T, select = NULL, MAF.range = c(1e-7, 0.5), nperbatch = 100, ncores = 1, is.dosage)
{
	if(class(obj) != "glmmkin") stop("Error: obj must be a class glmmkin object!")
	if(any(duplicated(obj$id_include))) {
		J <- Matrix(sapply(unique(obj$id_include), function(x) 1*(obj$id_include==x)), sparse = TRUE)
		res <- as.vector(as.matrix(crossprod(J, obj$scaled.residuals)))
		if(!is.null(obj$P)) obj$P <- as.matrix(crossprod(J, crossprod(obj$P, J)))
		else {
			obj$Sigma_iX <- crossprod(J, obj$Sigma_iX)
			obj$Sigma_i <- forceSymmetric(crossprod(J,crossprod(obj$Sigma_i,J)))
			obj$Sigma_i <- Matrix(obj$Sigma_i, sparse = TRUE)
		}
		rm(J)
	} else res <- obj$scaled.residuals

	        if(Sys.info()["sysname"] == "Windows" && ncores > 1) {
	 		warning("The package doMC is not available on Windows... Switching to single thread...")
			ncores <- 1
		}
	        ncores <- min(c(ncores, parallel::detectCores(logical = TRUE)))
	        gds <- SeqArray::seqOpen(infile)
		sample.id <- SeqArray::seqGetData(gds, "sample.id")
		if(is.null(select)) {
			if(any(is.na(match(unique(obj$id_include), sample.id)))) warning("Check your data... Some id_include in obj are missing in sample.id of infile!")
			select <- match(sample.id, unique(obj$id_include))
			select[is.na(select)] <- 0
			if(all(select == 0)) stop("Error: id_include in obj does not match sample.id in infile!")
		}
		if(length(select) != length(sample.id)) stop("Error: number of individuals in select does not match infile!")
		select2 <- select[select > 0]
		if(any(duplicated(select2)) || max(select2) > length(unique(obj$id_include))) stop("Error: select is a vector of orders, individuals not in obj should be coded 0!")
		res <- res[select2]
		if(!is.null(obj$P)) obj$P <- obj$P[select2, select2]
		else {
			obj$Sigma_iX <- obj$Sigma_iX[select2, , drop = FALSE]
			obj$Sigma_i <- obj$Sigma_i[select2, select2]
		}
		rm(select2)
    		variant.idx.all <- SeqArray::seqGetData(gds, "variant.id")
	        SeqArray::seqClose(gds)
		p.all <- length(variant.idx.all)
		if(ncores > 1) {
			require(doMC)
			doMC::registerDoMC(cores = ncores)
			p.percore <- (p.all-1) %/% ncores + 1
			n.p.percore_1 <- p.percore * ncores - p.all
			foreach(b = 1:ncores, .inorder=FALSE, .options.multicore = list(preschedule = FALSE, set.seed = FALSE)) %dopar% {
				variant.idx <- if(b <= n.p.percore_1) variant.idx.all[((b-1)*(p.percore-1)+1):(b*(p.percore-1))] else variant.idx.all[(n.p.percore_1*(p.percore-1)+(b-n.p.percore_1-1)*p.percore+1):(n.p.percore_1*(p.percore-1)+(b-n.p.percore_1)*p.percore)]
				p <- length(variant.idx)
				gds <- SeqArray::seqOpen(infile)
				SeqArray::seqSetFilter(gds, sample.id = sample.id[select > 0], verbose = FALSE)
				rm(sample.id); rm(select)
				nbatch.flush <- (p-1) %/% nperbatch + 1
				ii <- 0
				for(i in 1:nbatch.flush) {
		          	 	gc()
		      		  	tmp.variant.idx <- if(i == nbatch.flush) variant.idx[((i-1)*nperbatch+1):p] else variant.idx[((i-1)*nperbatch+1):(i*nperbatch)]
					tmp.p <- length(tmp.variant.idx)
					SeqArray::seqSetFilter(gds, variant.id = tmp.variant.idx, verbose = FALSE)
					geno <- SeqVarTools::imputedDosage(gds, dosage.field="DS", use.names = FALSE)
					AF<-colMeans(geno,na.rm=T)/2
					include <- (((AF >= MAF.range[1] & AF <= MAF.range[2]) | (AF >= 1-MAF.range[2] & AF <= 1-MAF.range[1])))
					if(sum(include) == 0) next
					ii <- ii + 1
					geno<-geno[,include,drop=FALSE]
					tmp.variant.idx <- tmp.variant.idx[include]
					tmp.p <- length(tmp.variant.idx)
					SeqArray::seqSetFilter(gds, variant.id = tmp.variant.idx, verbose = FALSE)
					SNP <- SeqArray::seqGetData(gds, "annotation/id")
					SNP[SNP == ""] <- NA
					out <- data.frame(SNP = SNP, CHR = SeqArray::seqGetData(gds, "chromosome"), POS = SeqArray::seqGetData(gds, "position"))
					rm(SNP)
					alleles.list <- strsplit(SeqArray::seqGetData(gds, "allele"), ",")
					out$REF <- unlist(lapply(alleles.list, function(x) x[1]))
					out$ALT <- unlist(lapply(alleles.list, function(x) paste(x[-1], collapse=",")))
					out$AF <- AF[include]
					rm(alleles.list, include)

    					N <- nrow(geno) - colSums(is.na(geno))
					if(center) geno <- scale(geno, scale = FALSE)
					miss.idx <- which(is.na(geno))
					if(length(miss.idx)>0) {
						geno[miss.idx] <- if(!center & missing.method == "impute2mean") colMeans(geno, na.rm = TRUE)[ceiling(miss.idx/nrow(geno))] else 0
					}
					SCORE <- as.vector(crossprod(geno, res))
					if(!is.null(obj$P)) VAR <- diag(crossprod(geno, crossprod(obj$P, geno)))
					else {
						GSigma_iX <- crossprod(geno, obj$Sigma_iX)
						VAR <- diag(crossprod(geno, crossprod(obj$Sigma_i, geno)) - tcrossprod(GSigma_iX, tcrossprod(GSigma_iX, obj$cov)))
					}
					PVAL <- ifelse(VAR>0, pchisq(SCORE^2/VAR, df=1, lower.tail=FALSE), NA)
					out <- cbind(out[,c("SNP","CHR","POS","REF","ALT")], N, AF=out[,c("AF")], SCORE,VAR,PVAL)
					names(out) <- c("SNP","CHR","POS","REF","ALT","N","AF","SCORE","VAR","PVAL")
					if(b == 1) {
				     	        write.table(out, outfile, quote=FALSE, row.names=FALSE, col.names=(ii == 1), sep="\t", append=(ii > 1), na=".")
					} else {
			     	                write.table(out, paste0(outfile, "_tmp.", b), quote=FALSE, row.names=FALSE, col.names=FALSE, sep="\t", append=(ii > 1), na=".")
					}
					rm(out)
				}
	        		SeqArray::seqClose(gds)
			}
			for(b in 2:ncores) {
			      	system(paste0("cat ", outfile, "_tmp.", b, " >> ", outfile))
				unlink(paste0(outfile, "_tmp.", b))
			}
		} else {
			variant.idx <- variant.idx.all
			rm(variant.idx.all)
			p <- length(variant.idx)
			gds <- SeqArray::seqOpen(infile)
			SeqArray::seqSetFilter(gds, sample.id = sample.id[select > 0], verbose = FALSE)
			rm(sample.id); rm(select)
			nbatch.flush <- (p-1) %/% nperbatch + 1
			ii <- 0
			for(i in 1:nbatch.flush) {
				print(paste0("Processing batch ",i,"/",nbatch.flush))
				ptm<-proc.time()[3]
		                gc()
		        	tmp.variant.idx <- if(i == nbatch.flush) variant.idx[((i-1)*nperbatch+1):p] else variant.idx[((i-1)*nperbatch+1):(i*nperbatch)]
				tmp.p <- length(tmp.variant.idx)
				SeqArray::seqSetFilter(gds, variant.id = tmp.variant.idx, verbose = FALSE)
				geno <- SeqVarTools::imputedDosage(gds, dosage.field="DS", use.names = FALSE)
				AF<-colMeans(geno,na.rm=T)/2
				include <- (((AF >= MAF.range[1] & AF <= MAF.range[2]) | (AF >= 1-MAF.range[2] & AF <= 1-MAF.range[1])))
				print(paste0("Number of variants in MAF range: ",sum(include)))

				if(sum(include) == 0) next
				ii <- ii + 1
				geno<-geno[,include,drop=FALSE]
				tmp.variant.idx <- tmp.variant.idx[include]
				tmp.p <- length(tmp.variant.idx)
				SeqArray::seqSetFilter(gds, variant.id = tmp.variant.idx, verbose = FALSE)
				SNP <- SeqArray::seqGetData(gds, "annotation/id")
				SNP[SNP == ""] <- NA
				out <- data.frame(SNP = SNP, CHR = SeqArray::seqGetData(gds, "chromosome"), POS = SeqArray::seqGetData(gds, "position"))
				rm(SNP)
				alleles.list <- strsplit(SeqArray::seqGetData(gds, "allele"), ",")
				out$REF <- unlist(lapply(alleles.list, function(x) x[1]))
				out$ALT <- unlist(lapply(alleles.list, function(x) paste(x[-1], collapse=",")))
				out$AF <- AF[include]
				rm(alleles.list, include)

    				N <- nrow(geno) - colSums(is.na(geno))
				if(center) geno <- scale(geno, scale = FALSE)
				miss.idx <- which(is.na(geno))
				if(length(miss.idx)>0) {
					geno[miss.idx] <- if(!center & missing.method == "impute2mean") colMeans(geno, na.rm = TRUE)[ceiling(miss.idx/nrow(geno))] else 0
				}
				SCORE <- as.vector(crossprod(geno, res))
				if(!is.null(obj$P)) VAR <- diag(crossprod(geno, crossprod(obj$P, geno)))
				else {
					GSigma_iX <- crossprod(geno, obj$Sigma_iX)
					VAR <- diag(crossprod(geno, crossprod(obj$Sigma_i, geno)) - tcrossprod(GSigma_iX, tcrossprod(GSigma_iX, obj$cov)))
				}
				PVAL <- ifelse(VAR>0, pchisq(SCORE^2/VAR, df=1, lower.tail=FALSE), NA)
				tmdiff<-proc.time()[3]-ptm
				print(paste0("P value calculation completed in ", tmdiff," seconds, writing to: ",outfile))

				out <- cbind(out[,c("SNP","CHR","POS","REF","ALT")], N, AF=out[,c("AF")], SCORE,VAR,PVAL)
				names(out) <- c("SNP","CHR","POS","REF","ALT","N","AF","SCORE","VAR","PVAL")
			     	write.table(out, outfile, quote=FALSE, row.names=FALSE, col.names=(ii == 1), sep="\t", append=(ii > 1), na=".")
				rm(out)
			}
	        	SeqArray::seqClose(gds)
		}
		return(invisible(NULL))
}

