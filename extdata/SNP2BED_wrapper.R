library('optparse')
options(stringsAsFactors=F)

option_list <- list(
  make_option("--SNPgds.file", type="character",default="",
    help="SNPRelate GDS file name"),
  make_option("--min.MAF", type="double",default="0.05",
    help="Minimum MAF"),
  make_option("--max.miss", type="double",default="0.05",
    help="Maximum missingness"),
  make_option("--removeSNPGDS", type="logical",default=FALSE,
    help="Whether to remove the SNPGDS file"),
  make_option("--prefix.bed", type="character",default="",
    help="Prefix of the intermediate plink BED file output")
)

parser <- OptionParser(usage="%prog [options]", option_list=option_list)

args <- parse_args(parser, positional_arguments = 0)
opt <- args$options
print(opt)


SNPgds.file<-opt$SNPgds.file
min.MAF<-opt$min.MAF
max.miss<-opt$max.miss
removeSNPGDS<-opt$removeSNPGDS
prefix.bed<-opt$prefix.bed


library('FastSparseGRM')
print(sessionInfo())
SNP2BED_wFilter(SNPgds.file,min.MAF,max.miss,prefix.bed,removeSNPGDS=removeSNPGDS)


