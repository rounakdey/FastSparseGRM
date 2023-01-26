# FastSparseGRM
FastSparseGRM is an R package that efficiently calculates genetic principal components (PCs) and the ancestry-adjusted sparse genetic relatedness matrix (GRM). It accounts for population heterogeneity using genetic PCs which are automatically calculated as part of the pipeline. The genetic PCs can be used as fixed effect covariates to account for the population stratification and the sparse GRM can be used to model the random effects to account for the sample relatedness in a mixed effects phenotype-genotype association testing model.

FastSparseGRM utilizes multithreading for efficient computation and can benefit greatly from using a large number of CPU cores.
   
## How to install FastSparseGRM
Install FastSparseGRM from source code:

    git clone https://github.com/rounakdey/FastSparseGRM.git
    R CMD INSTALL FastSparseGRM

FastSparseGRM also requires the softwares [KING](https://kingrelatedness.com/) and [plink](https://www.cog-genomics.org/plink/) for some of the steps in its pipeline. Executables of KING (v2.2) and plink (v1.9) for 64-bit Linux are supplied for convenience in the *extdata* folder. Use of local installations of KING (version >=v2.1.6) and plink (version >= v1.9) are also okay.

**To use the functions in FastSparseGRM, please copy all the files from the *extdata* folder to the working directory. These files provide the KING and plink executables, as well as wrapper functions to run the analysis.**

## Input format
FastSparseGRM takes genotype data in plink BED format as input. To assist the analysis of modern whole-genome sequencing (WGS) data typically saved in the [genomic data structure (GDS)](https://academic.oup.com/bioinformatics/article/33/15/2251/3072873) format, we also provides scripts to convert from [SeqArray GDS](https://bioconductor.org/packages/release/bioc/html/SeqArray.html) and [SNP-GDS](http://bioconductor.org/packages/release/bioc/html/SNPRelate.html) formats to the plink BED format.

### Convert from SeqArray GDS format (newer, more common format) to BED.

    R CMD BATCH --vanilla '--args --gds.file <gds.file> --min.AVGDP <min.AVGDP: default 10> --filterCat <filterCat: default PASS> --min.MAF <min.MAF: default 0.05> --max.miss <max.miss: default 0.05> --removeSNPGDS <removeSNPGDS: default TRUE> --prefix.bed <prefix.out>' Seq2BED_wrapper.R Seq2BED_wrapper.Rout

To convert from SeqArray GDS format, the GDS file must have the following nodes: `annotation/filter` (PASS, FAIL or any other filtering category), `annotation/info/AVGDP` (Average sequencing depth), and `annotation/info/AC` (Allele counts).

### Convert from SNP-GDS format (older format) to BED.

    R CMD BATCH --vanilla '--args --SNPgds.file <SNPgds.file> --min.MAF <min.MAF: default 0.05> --max.miss <max.miss: default 0.05> --removeSNPGDS <removeSNPGDS: default FALSE> --prefix.bed <prefix.out>' SNP2BED_wrapper.R SNP2BED_wrapper.Rout

### Details about the arguments
+ `gds.file`: Address of the SeqArray GDS input file.
+ `SNPgds.file`: Address of the SNP-GDS input file.
+ `filterCat`: Filtering category to be selected. Variants with this category (default is PASS) in the `annotation/filter` node will be included. Only one category is allowed.
+ `min.AVGDP`: Minimum average sequencing depth for a variant to be included. Default is 10.
+ `min.MAF`: Minimum minor allele frequency for a variant to be included. Default is 0.05.
+ `max.miss`: Maximum missingness for a variant to be included. Default is 0.05.
+ `removeSNPGDS`: Whether to remove the SNP-GDS intermediate/input file. Default is TRUE for Seq2BED_wrapper.R and FALSE for SNP2BED_wrapper.R.
+ `prefix.bed`: Prefix to the output BED file.

### Preparing genotype data for the analysis pipeline
**FastSparseGRM pipeline requires data from all 22 autosomes in a single BED file.** So, if your BED files are chromosome-specific, now would be the time to merge them into one single BED file. Here is an example script of how to do this on *bash* using *plink*,

    rm mergeBED.list; for i in {1..22}; do echo chr${i} >> mergeBED.list; done;
    plink --merge-list mergeBED.list --make-bed --out chrall

Here, the chromosome-specific BED files are `chr1.bed`, `chr2.bed`, ... and the merged BED output is `chrall.bed`.

**Next, some of the steps in the FastSparseGRM pipeline requires pruned genotypes.** It is okay to use unpruned genotypes there, but that will result in substantially high computation cost. Here is an example script of how to perform pruning on *bash* using *plink*,

    plink --bfile chrall --indep-pairwise 50 5 0.1 --out chrall.prunedlist;
    plink --bfile chrall --extract chrall.prunedlist.prune.in --make-bed --out chrall_pruned

Here, the unpruned genotype file is `chrall.bed`, and pruned genotype file is `chrall_pruned.bed`.

## How to run the FastSparseGRM pipeline
There are two ways of running the analysis pipeline.

### Option 1. Run the pipeline step-by-step.
This gives the user more control over memory and CPU usage in each step. **This option is suggested for analyzing very large dataset.**

Step 1. Run KING (version >=2.1.6)

    king -b <unfiltered.bedfile> --ibdseg --degree <degree: default 4> --cpus <n_cpus> --prefix <output.king>

Step 2. Get ancestry divergence estimates

    R CMD BATCH --vanilla '--args --prefix.in <prefix.bedfile> --file.seg <output.king> --num_threads <n_cpus> --degree <degree: default 4> --divThresh <divThresh: default -0.02209709> --nRandomSNPs <nRandomSNPs: default 0> --prefix.out <output.divergence>' getDivergence_wrapper.R getDivergence.Rout

Step 3. Extract unrelated samples

    R CMD BATCH --vanilla '--args --prefix.in <prefix.bedfile> --file.seg <output.king> --degree <degree: default 4> --file.div <output.divergence> --file.include <file.include: default ""> --prefix.out <output.unrelated>' extractUnrelated_wrapper.R extractUnrelated.Rout

Step 4. Run PCA

    R  CMD BATCH --vanilla '--args --prefix.in <prefix.bedfile> --file.unrels <output.unrelated> --prefix.out <output.pca> --no_pcs <no_pcs: default 20> --num_threads <n_cpus> --no_iter <no_iter: default 10>' runPCA_wrapper.R runPCA.Rout

Step 5. Calculate Sparse GRM

    R CMD BATCH --vanilla '--args --prefix.in <prefix.bedfile> --prefix.out <output.sparseGRM> --file.train <output.unrelated> --file.score <output.pca> --file.seg <output.king> --num_threads <n_cpus> --no_pcs <no_pcs: default 20> --block.size <block.size: default 5000> --max.related.block <max.related.block: default 5000> --KINGformat.out <KINGformat.out: default FALSE> --degree <degree: default 4>' calcSparseGRM_wrapper.R calcSparseGRM.Rout

### Option 2. Run the entire pipeline with one wrapper function.
The user has to specify the same memory and CPU threads for all the steps, and this option does not have the option to fine-tune every step of the computation. **This option is intended for smaller sample size, and is not suggested for very large dataset.**

    R CMD BATCH --vanilla '--args --prefix.in <prefix.bedfile> --prefix.in.unfiltered <prefix.unfiltered.bedfile> --prefix.out <output.sparseGRM> --KING.executable <king.executable.path: default king> --degree <degree: default 4> --num_threads <n_cpus> --no_pcs <no_pcs: default 20>' runPipeline_wrapper.R runPipeline.Rout

### Details about the arguments
+ `prefix.bedfile`: Prefix of the BED file used for calculating ancestry divergences, PC scores and the Sparse GRM. This file can be pruned and MAF filtered (suggested to use MAF > 0.05 and LD pruned variants). Around 200,000 variants are good enough for this input.
+ `prefix.unfiltered.bedfile` (`unfiltered.bedfile`): Prefix (full path) of the BED file used for running KING IBD segment analysis. No MAF filtering or pruning is required for this file when analyzing smaller datasets. For larger datasets (like UK Biobank or TOPMed), using MAF filtering (MAF > 0.05) and LD pruning is okay, and can save substantial computation time. No not need to specify `prefix.unfiltered.bedfile` in Option 2 if you want to use the same file for both running KING and PCA/Sparse GRM calculation.
+ `king.executable.path`: Full path to the KING executable file. Default is "king", which means KING is already installed in the system, or available at the working directory.
+ `degree`: After what degree of relatedness should the algorithm consider a pair of samples to be unrelated. Default is 4.
+ `n_cpus`: Number of CPU threads to use. Using `n_cpus` >= 20 can substantially speed up computation.
+ `no_pcs`: Number of genetic PCs to calculate/use. Default is 20.
+ `divThresh`: Threshold to use for ancestry divergence calculation. Default is -2^-(4+1.5).
+ `nRandomSNPs`: Number of random SNPs to use for ancestry divergence calculation. Default is 0, which means to use all the SNPs available. But the user can specify a smaller number for faster computation in large datasets (suggested >= 10000).
+ `file.include`: IDs for samples to always be included in the unrelated set. The first column should be FID and second column should be IID. If only one column present, then that will be treated as both FID and IID. Default is "", i.e., skipping this argument if you do not have such a list of samples to always include as unrelateds.
+ `no_iter`: Number of power method iterations to run for the randomized PCA algorithm. Default is 10.
+ `block.size`: Size of the SNP blocks to read at-a-time when calculation the sparse GRM. This is to control the memory usage. Larger block size utilizes more memory, but results in faster computation. Default is 5000.
+ `max.related.block`: Maximum size of a related family block. The sparsity threshold will be adjusted accordingly. Default is 5000.
+ `KINGformat.out`: Whether to output the entries of the sparse GRM in KING format (FID1 ID1 FID2 ID2 Kinship). Default is FALSE.
+ `output.king` (\*.seg)/`output.divergence` (\*.div)/`output.unrelated` (\*.unrels)/`output.pca` (\*.score): Intermediate outputs for each of the steps.
+ `output.sparseGRM`: Prefix of the final sparse GRM output file. The output will be in \*.RData format. If `KINGformat.out` is TRUE, then it will also output the results in KING format in a \*.kins file.
