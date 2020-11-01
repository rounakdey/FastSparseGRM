library('optparse')

options(stringsAsFactors=F)


option_list <- list(
  make_option("--prefix.in", type="character", default="",
    help="Prefix to the input PLINK fam file"),
  make_option("--file.seg", type="character",default="",
    help="Full path to the ibdseg/Kinship file from KING"),
  make_option("--file.div", type="character",default="",
    help="Full path to the ancestry divergence file"), 
  make_option("--degree", type="integer", default=4,
    help="Consider samples related up to which degree (max 4, default 4)"),
  make_option("--file.include", type="character", default="",
    help="IDs for samples to always include in the unrelated set. The first column should be FID and second column should be IID. If only one column present, then that will be treated as both FID and IID"),
  make_option("--prefix.out", type="character", default="",
    help="Prefix to the output file")
)

parser <- OptionParser(usage="%prog [options]", option_list=option_list)

args <- parse_args(parser, positional_arguments = 0)
opt <- args$options
print(opt)

file.seg <- opt$file.seg
file.div <- opt$file.div
prefix.in <- opt$prefix.in
degree <- opt$degree
prefix.out <- opt$prefix.out
file.include <- opt$file.include

library('FastSparseGRM')
extractUnrelated(file.seg,file.div,prefix.in,degree,file.include,prefix.out)
