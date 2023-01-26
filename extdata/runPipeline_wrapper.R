library('optparse')
options(stringsAsFactors=F)


option_list <- list(
  make_option("--prefix.in", type="character", default="",
    help="Prefix to the input PLINK bed file used for PCA and SparseGRM calculation (it should be accompanied by a fam and a bim file with the same prefix). Pruned and MAF filtered bed files are preferable."),
  make_option("--prefix.in.unfiltered", type="character", default="",
    help="Prefix to the input PLINK bed file used for KING (it should be accompanied by a fam and a bim file with the same prefix). To run KING, pruning and MAF filtering are not needed. If not provided, then prefix.in will be used."),
  make_option("--KING.executable", type="character", default="king",
    help="Address to the KING executable file."),
  make_option("--num_threads", type="integer", default=0,
    help="Number of threads within each core to use (if 0, then automatically calculated)"),
  make_option("--degree", type="integer", default=4,
    help="Consider samples related up to which degree (max 4, default 4)"),
  make_option("--divThresh", type="numeric", default=-0.02209709,
    help="Threshold for calculating ancestry divergence among unrelated samples (default -0.02209709)"),
  make_option("--nRandomSNPs", type="integer", default=0,
    help="Number of random SNPs to calculate ancestry divergence (default 0, which uses all available SNPs in the plink file)"),
  make_option("--file.include", type="character", default="",
    help="IDs for samples to always include in the unrelated set. The first column should be FID and second column should be IID. If only one column present, then that will be treated as both FID and IID"),
  make_option("--no_pcs", type="integer", default=10,
    help="Number of PCs to calculate and adjust (default 10)"),
  make_option("--no_iter", type="integer", default=10,
    help="Number of iterations for randomized PCA (default 10)"),
  make_option("--block.size", type="integer", default=5000,
    help="SNP block size to load in memory"),
  make_option("--max.related.block", type="integer", default=5000,
    help="Maximum allowed size of a related sample block"),
  make_option("--tempDir", type="character", default="",
    help="Directory to store intermediate files. If not provided, then default temporary directory of the system will be used."),
  make_option("--deleteTemp", type="logical", default=TRUE,
    help="Whether to delete the temporary files or not."),
  make_option("--prefix.out", type="character", default="",
    help="Prefix to the output file")
)

parser <- OptionParser(usage="%prog [options]", option_list=option_list)

args <- parse_args(parser, positional_arguments = 0)
opt <- args$options
print(opt)

prefix.in <- opt$prefix.in
prefix.in.unfiltered <- opt$prefix.in.unfiltered
KING.executable <- opt$KING.executable
num_threads <- opt$num_threads
degree <- opt$degree
divThresh <- opt$divThresh
nRandomSNPs <- opt$nRandomSNPs
file.include <- opt$file.include
no_pcs <- opt$no_pcs
no_iter <- opt$no_iter
block.size <- opt$block.size
max.related.block <- opt$max.related.block
tempDir <- opt$tempDir
deleteTemp <- opt$deleteTemp
prefix.out <- opt$prefix.out


getPCadjSparseGRM<-function(prefix.in,prefix.in.unfiltered,KING.executable,num_threads,degree,divThresh,nRandomSNPs,file.include,no_pcs,no_iter,block.size,max.related.block,tempDir,deleteTemp,prefix.out)
{
	library('FastSparseGRM')
	print(sessionInfo())

	if(prefix.in.unfiltered=="")	prefix.in.unfiltered<-prefix.in
	if(tempDir=="")	tempDir<-tempdir()
	prefix.out.basename<-basename(prefix.out)
	prefix.temp.out<-paste0(tempDir,"/",prefix.out.basename)

	print("Running KING IBD segment analysis")
	system(paste0(KING.executable," -b ",prefix.in.unfiltered,".bed --ibdseg --degree 4 --cpus ",num_threads," --prefix ",prefix.out=prefix.temp.out))

	print("Calculating ancestry divergence")
	getDivergence(file.seg=paste0(prefix.temp.out,".seg"),prefix.in,num_threads,divThresh,degree,nRandomSNPs,prefix.temp.out)

	print("Extracting unrelated samples")
	extractUnrelated(file.seg=paste0(prefix.temp.out,".seg"),file.div=paste0(prefix.temp.out,".div"),prefix.in,degree,file.include,prefix.out=prefix.temp.out)

	print("Calculating PC scores")
	runPCA(prefix.in,file.unrels=paste0(prefix.temp.out,".unrels"),no_pcs,num_threads,prefix.out,no_iter)

	print("Calculating Sparse GRM")
	calcSparseGRM(prefix.in,file.score=paste0(prefix.out,".score"),file.train=paste0(prefix.temp.out,".unrels"),file.seg=paste0(prefix.temp.out,".seg"),no_pcs,num_threads,prefix.out,degree,block.size,max.related.block)

	if(deleteTemp==TRUE)
	{
		print("Cleaning up intermediate files")
		unlink(paste0(prefix.temp.out,c(".seg",".div",".unrels",".segments.gz","allsegs.txt")))
		unlink(paste0(prefix.out,c(".evec",".eval")))
	}

}


getPCadjSparseGRM(prefix.in,prefix.in.unfiltered,KING.executable,num_threads,degree,divThresh,nRandomSNPs,file.include,no_pcs,no_iter,block.size,max.related.block,tempDir,deleteTemp,prefix.out)

