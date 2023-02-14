# ewastools 1.7.2
* Added support for EPICv2 chip

# ewastools 1.7.1
* Added HRS cell composition model

# ewastools 1.7
* Saliva cell composition prediction
* Test `detectionP.minfi` for compatibility
* Drop `normalize`

# ewastools 1.6
* `read_idats`: extract Sentrix barcode and position from .idat files, stored in `meta` data.table (92eb8df289a473fa0c06471be197db1953e39efd0)
* `call_genotypes`: use pre-specified parameters for SNP mixture model by default, i.e., `learn=FALSE` (507abceeb8db572d72089e423678f3cd5a7f6769)
* Update recommedation of pre-processing workflow: call `mask()` after `check_sex()` (eb54d1b0d74d2b0231db6e03a8c60710abd42e09)
* Add some unit tests
* Fix: some bisulfite conversion control probes were ignored due to inconsistent naming (24233ffcbea4337b9a7d4d3060b03b484d5ba01c)