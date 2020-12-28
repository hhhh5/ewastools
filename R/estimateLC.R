#' Cell composition
#'
#' Estimation of cell proportions for blood and saliva using the Houseman algorithm
#'
#' @rdname cell_composition
#'
#' @param meth Matrix of beta values
#' @param ref Choice of reference dataset: available options are `Reinius` [2]`, Bakulski` [3],
#' `deGoede` [4], `Gervin`[5], `Lin` [6], `Mill` [GSE103541], `Salas` [7], `Lolipop` [8] or
#' combinations of them, concatenated by `+`, e.g. `Reinius+Lin` (you might have to flip first and
#' second name)). Furthermore, the option `saliva` and `salivaEPIC` are available [9].
#' @param constrained Force that all cell proportions sum up to 1.
#'
#' @return Estimated cell proportions B-lymphocytes, CD4 T-cells, CD8 T-cells, granulocytes, monocytes,  natural killer cells (and nucleated red blood cells) using the Houseman algorithm [1]. Models were trained on various reference datasets of purified cell types.
#'
#' @references{Houseman EA, et al. DNA methylation arrays as surrogate measures of cell mixture distribution. BMC bioinformatics. 2012 Dec 1;13(1):86.}
#' @references{Reinius LE, et al. Differential DNA methylation in purified human blood cells: implications for cell lineage and studies on disease susceptibility. PloS one. 2012 Jul 25;7(7):e41361.}
#' @references{Bakulski KM, et al. DNA methylation of cord blood cell types: applications for mixed cell birth studies. Epigenetics. 2016 May 3;11(5):354-62.}
#' @references{de Goede OM, et al. Nucleated red blood cells impact DNA methylation and expression analyses of cord blood hematopoietic cells. Clinical epigenetics. 2015 Dec;7(1):95.}
#' @references{Gervin K, et al. Cell type specific DNA methylation in cord blood: A 450K-reference data set and cell count-based validation of estimated cell type composition. Epigenetics. 2016 Sep 1;11(9):690-8}
#' @references{Gervin K, et al. Systematic evaluation and validation of reference and library selection methods for deconvolution of cord blood DNA methylation data. bioRxiv. 2019 March; DOI:10.1101/570457}
#' @references{Salas LA, et al. An optimized library for reference-based deconvolution of whole-blood biospecimens assayed using the Illumina HumanMethylationEPIC BeadArray. Genome biology. 2018 Dec;19(1):64.}
#' @references{Heiss JA, et al. Training a model for estimating leukocyte composition using whole-blood DNA methylation and cell counts as reference. Epigenomics. 2017 Jan;9(1):13-20.}
#' @references{Middleton LYM, et al. Saliva cell type DNA methylation reference panel for epidemiology studies in children. 2020 Sep;}
#'
#' @export
#'

estimateLC = function(meth,ref,constrained=FALSE){
    
    J = ncol(meth)

    coefs = read.table(system.file(paste0("data/",ref,".txt"),package="ewastools"))
    coefs = as.matrix(coefs)
    n_celltypes = ncol(coefs)

    markers = match(rownames(coefs),rownames(meth))
    EST = sapply(1:J,function(j){
        tmp = meth[markers,j]
        i = !is.na(tmp)

        if(constrained == FALSE){
            return(
                quadprog::solve.QP(
                 t(coefs[i,]) %*% coefs[i,]
                ,t(coefs[i,]) %*% tmp[i]
                ,diag(n_celltypes)
                ,rep(0,n_celltypes)
            )$sol)
        }else{
            return(
                quadprog::solve.QP(
                 t(coefs[i,]) %*% coefs[i,]
                ,t(coefs[i,]) %*% tmp[i]
                ,cbind(1,diag(n_celltypes))
                ,c(1,rep(0,n_celltypes))
                ,meq=1
            )$sol)
        }
        })
    EST = t(EST)
    colnames(EST) = colnames(coefs)
    EST = data.table(EST)

    return(EST)
}
