library('optparse')
options(stringsAsFactors=F)

option_list <- list(
  make_option("--gds.file", type="character",default="",
    help="SeqArray GDS or aGDS file name"),
  make_option("--min.AVGDP", type="integer",default="10",
    help="Minimum average depth"),
  make_option("--filterCat", type="character",default="PASS",
    help="Category of variants to select based on the filter field"),
  make_option("--min.MAF", type="double",default="0.05",
    help="Minimum MAF"),
  make_option("--max.miss", type="double",default="0.05",
    help="Maximum missingness"),
  make_option("--removeSNPGDS", type="logical",default=TRUE,
    help="Whether to remove the intermediate SNPGDS"),
  make_option("--prefix.bed", type="character",default="",
    help="Prefix of the intermediate plink BED file output")
)

parser <- OptionParser(usage="%prog [options]", option_list=option_list)

args <- parse_args(parser, positional_arguments = 0)
opt <- args$options
print(opt)


gds.file<-opt$gds.file
min.AVGDP<-opt$min.AVGDP
filterCat<-opt$filterCat
min.MAF<-opt$min.MAF
max.miss<-opt$max.miss
removeSNPGDS<-opt$removeSNPGDS
prefix.bed<-opt$prefix.bed


library('FastSparseGRM')
print(sessionInfo())
Seq2SNP_wFilter(gds.file,min.AVGDP,filterCat,prefix.bed)
SNPgds.file<-paste0(prefix.bed,".gds")
SNP2BED_wFilter(SNPgds.file,min.MAF,max.miss,prefix.bed,removeSNPGDS=removeSNPGDS)


