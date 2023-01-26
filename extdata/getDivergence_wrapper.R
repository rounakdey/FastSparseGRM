library('optparse')
options(stringsAsFactors=F)


option_list <- list(
  make_option("--prefix.in", type="character", default="",
    help="Prefix to the input PLINK bed file (it should be accompanied by a fam and a bim file with the same prefix)"),
  make_option("--file.seg", type="character",default="",
    help="Full path to ibdseg file from KING"),
  make_option("--num_threads", type="integer", default=0,
    help="Number of threads within each core to use (if 0, then automatically calculated)"),
  make_option("--degree", type="integer", default=4,
    help="Consider samples related up to which degree (max 4, default 4)"),
  make_option("--divThresh", type="numeric", default=-0.02209709,
    help="Threshold for calculating ancestry divergence among unrelated samples (default -0.02209709)"),
  make_option("--nRandomSNPs", type="integer", default=0,
    help="Number of random SNPs to calculate ancestry divergence (default 0, which uses all available SNPs in the plink file)"),
  make_option("--prefix.out", type="character", default="",
    help="Prefix to the output file")
)

parser <- OptionParser(usage="%prog [options]", option_list=option_list)

args <- parse_args(parser, positional_arguments = 0)
opt <- args$options
print(opt)

file.seg <- opt$file.seg
prefix.in <- opt$prefix.in
num_threads <- opt$num_threads
degree <- opt$degree
divThresh <- opt$divThresh
nRandomSNPs <- opt$nRandomSNPs
prefix.out <- opt$prefix.out

library('FastSparseGRM')
print(sessionInfo())
getDivergence(file.seg,prefix.in,num_threads,divThresh,degree,nRandomSNPs,prefix.out)

