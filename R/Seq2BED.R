
Seq2SNP_wFilter<-function(gds.file,min.AVGDP,filterCat,prefix.bed)
{
	options(stringsAsFactors=F)

	gds<-seqOpen(gds.file)

	nsamp<- length(seqGetData(gds, "sample.id"))
	filter <- seqGetData(gds, "annotation/filter")
	AVGDP <- seqGetData(gds, "annotation/info/AVGDP")
	AC <- seqGetData(gds, "annotation/info/AC")
	MAF<- AC/(2*nsamp)
	MAF<-pmin(MAF,1-MAF)
	AllFilter <- filter == filterCat & AVGDP >= min.AVGDP & isSNV(gds,biallelic=TRUE) & MAF > 0.01

	variant.id <- seqGetData(gds, "variant.id")
	variant.sel<- variant.id[AllFilter]
	seqSetFilter(gds, variant.sel=NULL, sample.sel=NULL, variant.id=variant.sel,verbose=TRUE)

	SNPgds.file<-paste0(prefix.bed,".gds")
	seqGDS2SNP(gds, SNPgds.file, verbose=TRUE)
	seqClose(gds)

}

SNP2BED_wFilter<-function(SNPgds.file,min.MAF,max.miss,prefix.bed,removeSNPGDS=TRUE)
{
	options(stringsAsFactors=F)

	SNPgds<-snpgdsOpen(SNPgds.file)
	SNP.select <- snpgdsSelectSNP(SNPgds, maf=min.MAF, missing.rate=1-max.miss)
	snpgdsGDS2BED(SNPgds, prefix.bed, snp.id=SNP.select, snpfirstdim=FALSE, verbose=TRUE)
	snpgdsClose(SNPgds)
	if(removeSNPGDS==TRUE)	unlink(SNPgds.file)

	####Use SNP IDs in the format Chr:Pos_RefA/AltA
	bim<-read.table(paste0(prefix.bed,".bim"),header=F)
	bim[,2]<-paste0(bim[,1],":",bim[,4],"_",bim[,5],"/",bim[,6])
	write.table(bim,paste0(prefix.bed,".bim"),row.names=F,col.names=F,quote=F)
}


