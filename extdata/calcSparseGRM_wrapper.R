library('optparse')

options(stringsAsFactors=F)


option_list <- list(
  make_option("--prefix.in", type="character", default="",
    help="Prefix to the input PLINK bed file (it should be accompanied by a fam and a bim file with the same prefix)"),
  make_option("--file.score", type="character",default="",
    help="Full path to the file containing the PC scores"), 
  make_option("--file.train", type="character",default="",
    help="Full path to the file containing the IDs of the training samples using which the individual specific allele frequencies will be calculated. The first two columns need to represent FID and IID respectively."), 
  make_option("--file.seg", type="character",default="",
    help="Full path to ibdseg file from KING"), 
  make_option("--no_pcs", type="integer", default=0,
    help="Number of PCs to use"),
  make_option("--num_threads", type="integer", default=0,
    help="Number of threads within each core to use (if 0, then automatically calculated)"),
  make_option("--block.size", type="integer", default=5000,
    help="SNP block size to load in memory"),
  make_option("--max.related.block", type="integer", default=5000,
    help="Maximum allowed size of a related sample block"),
  make_option("--degree", type="integer", default=4,
    help="Estimate kinships for samples related up to which degree (max 4, default 4)"),
  make_option("--prefix.out", type="character", default="",
    help="Prefix to the output files")
)

parser <- OptionParser(usage="%prog [options]", option_list=option_list)

args <- parse_args(parser, positional_arguments = 0)
opt <- args$options
print(opt)

file.score <- opt$file.score
file.train <- opt$file.train
file.seg<-opt$file.seg
prefix.in <- opt$prefix.in
num_threads <- opt$num_threads
block.size <- opt$block.size
no_pcs <- opt$no_pcs
prefix.out <- opt$prefix.out
degree <- opt$degree
max.related.block <- opt$max.related.block

library('FastSparseGRM')
calcSparseGRM(prefix.in,file.score,file.train,file.seg,no_pcs,num_threads,prefix.out,degree,block.size,max.related.block)





