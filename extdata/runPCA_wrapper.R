library('optparse')
options(stringsAsFactors=F)

option_list <- list(
  make_option("--prefix.in", type="character", default="",
    help="Prefix to the input PLINK bed file (it should be accompanied by a fam and a bim file with the same prefix)"),
  make_option("--file.unrels", type="character",default="",
    help="Full path to the file containing unrelated sample ids. It should have FIDs in the first column and IIDs in the second column. If only one column is provided, then that will be used as both FID and IID. If empty, all samples will be considered unrelated"), 
  make_option("--no_pcs", type="integer", default=10,
    help="Number of PCs in the output (default 10)"),
  make_option("--num_threads", type="integer", default=0,
    help="Number of threads within each core to use (if 0, then automatically calculated)"),
  make_option("--prefix.out", type="character", default="",
    help="Prefix to the output files"),
  make_option("--no_iter", type="integer", default=10,
    help="Number of iterations for randomized PCA (default 10)")
)

parser <- OptionParser(usage="%prog [options]", option_list=option_list)

args <- parse_args(parser, positional_arguments = 0)
opt <- args$options
print(opt)

file.unrels <- opt$file.unrels
prefix.in <- opt$prefix.in
no_pcs <- opt$no_pcs
no_cores <- opt$no_cores
num_threads <- opt$num_threads
prefix.out <- opt$prefix.out
no_iter <- opt$no_iter

library('FastSparseGRM')
print(sessionInfo())
runPCA(prefix.in,file.unrels,no_pcs,num_threads,prefix.out,no_iter)
