 # ewastools 1.6

 * `read_idats`: extract Sentrix barcode and position from .idat files, stored in `meta` data.table
 * `call_genotypes`: use pre-specified parameters for SNP mixture model by default, i.e., `learn=FALSE`
 * Update recommedation of pre-processing workflow: call `mask()` after `check_sex()`.