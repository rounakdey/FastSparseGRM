# FastSparseGRM
Efficiently calculate ancestry-adjusted sparse GRM

## Example pipeline
There are two ways of running the analysis pipeline.

### Option 1. Run the pipeline step-by-step.
This gives the user more control over memory and CPU usage in each step. This option is suggested for analyzing very large dataset.

Step 1. Run KING (version >=2.1.6) - KING suggests not using pruned SNPs

    king -b <unfiltered.bedfile> --ibdseg --degree 4 --cpus <n_cpus> --prefix <output.king>

Step 2. Get ancestry divergence estimates - From this step onwards, one can use pruned SNPs to save computation costs

    R CMD BATCH --vanilla '--args --prefix.in <prefix.bedfile> --file.seg <output.king> --num_threads <n_cpus> --degree 4 --nRandomSNPs 0 --prefix.out <output.divergence>' getDivergence_wrapper.R getDivergence.Rout

Step 3. Extract unrelated samples

    R CMD BATCH --vanilla '--args --prefix.in <prefix.bedfile> --file.seg <output.king> --degree 4 --file.div <output.divergence> --prefix.out <output.unrelated>' extractUnrelated_wrapper.R extractUnrelated.Rout

Step 4. Run PCA

    R  CMD BATCH --vanilla '--args --prefix.in <prefix.bedfile> --file.unrels <output.unrelated> --prefix.out <output.pca> --no_pcs 20 --num_threads <n_cpus>' runPCA_wrapper.R runPCA.Rout

Step 5. Calculate Sparse GRM

    R CMD BATCH --vanilla '--args --prefix.in <prefix.bedfile> --prefix.out <output.sparseGRM> --file.train <output.unrelated> --file.score <output.pca.score> --file.seg <output.king> --num_threads <n_cpus> --no_pcs 20' calcSparseGRM_wrapper.R calcSparseGRM.Rout

### Option 2. Run the entire pipeline with one wrapper function.
The user has to specify the same memory and CPU threads for all the steps. This option is not suggested for very large dataset.

    R CMD BATCH --vanilla '--args --prefix.in <prefix.bedfile> --prefix.in.unfiltered <prefix.unfiltered.bedfile> --prefix.out <output.sparseGRM> --KING.executable <king.executable.path> --degree 4 --num_threads <n_cpus> --no_pcs 20' runPipeline_wrapper.R runPipeline.Rout

### Details about the arguments
+ *prefix.bedfile*: Prefix of the BED file used for calculating PC scores and Sparse GRM. This file can be pruned and MAF filtered (suggested to use MAF>0.05 and LD pruned variants).
+ *prefix.unfiltered.bedfile* (*unfiltered.bedfile*): Prefix (full path) of the BED file used for running KING IBD segment analysis. No MAF filtering or pruning is required (and suggested as well) for this file. No not need to specify *prefix.unfiltered.bedfile* in Option 2 if you want to use the same file for both running KING and PCA/Sparse GRM calculation.
+ *n_cpus*: Number of CPU threads to use.
+ *king.executable.path*: Full path to the KING executable file.
+ *output.king/output.divergence/output.unrelated/output.pca*: Intermediate outputs for each of the steps.
+ *output.sparseGRM*: Prefix of the final sparse GRM output file.
