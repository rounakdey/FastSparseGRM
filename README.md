# FastSparseGRM
Efficiently calculate ancestry-adjusted sparse GRM

## Example pipeline
Step 1. Run KING (version >=2.1.6) - KING suggests not using pruned SNPs

    king -b <bedfile> --ibdseg --degree 4 --cpus <n_cpus> --prefix <output.king>

Step 2. Get ancestry divergence estimates - From this step onwards, one can use pruned SNPs to save computation costs

    R CMD BATCH --vanilla '--args --prefix.in <prefix.bedfile> --file.seg <output.king> --num_threads <n_cpus> --degree 4 --nRandomSNPs 0 --prefix.out <output.divergence>' getDivergence_wrapper.R getDivergence.Rout

Step 3. Extract unrelated samples

    R CMD BATCH --vanilla '--args --prefix.in <prefix.bedfile> --file.seg <output.king> --degree 4 --file.div <output.divergence> --prefix.out <output.unrelated>' extractUnrelated_wrapper.R extractUnrelated.Rout

Step 4. Run PCA

    R  CMD BATCH --vanilla '--args --prefix.in <prefix.bedfile> --file.unrels <output.unrelated> --prefix.out <output.pca> --no_pcs 20 --num_threads 30' runPCA_wrapper.R runPCA.Rout

Step 5. Calculate Sparse GRM

    R CMD BATCH --vanilla '--args --prefix.in <prefix.bedfile> --prefix.out <output.sparseGRM> --file.train <output.unrelated> --file.score <output.pca.score> --file.seg <output.king> --block.size 50000 --num_threads 30 --no_pcs 20' calcSparseGRM_wrapper.R calcSparseGRM.Rout

