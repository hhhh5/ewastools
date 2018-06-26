#' Estimate leukcoyte composition
#' 
#' @export
#' @param meth Matrix of beta values
#' @param ref Choice of reference dataset: available options are `gervin` [2], `reinius` [3], `goede` [4], `bakulski` [5] or specific combinations of them. Current default is `bakulski+goede`.
#' @return Estimated cell proportions B-lymphocytes, CD4 T-cells, CD8 T-cells, granulocytes, monocytes,  natural killer cells (and nucleated red blood cells) using the Houseman algorithm [1]. Models were trained on various reference datasets of purified cell types. 
#' @references{Houseman EA, Accomando WP, Koestler DC, Christensen BC, Marsit CJ, Nelson HH, Wiencke JK, Kelsey KT. DNA methylation arrays as surrogate measures of cell mixture distribution. BMC bioinformatics. 2012 Dec 1;13(1):86.}
#' @references{Gervin K, Page CM, Aass HC, Jansen MA, Fjeldstad HE, Andreassen BK, Duijts L, van Meurs JB, van Zelm MC, Jaddoe VW, Nordeng H. Cell type specific DNA methylation in cord blood: A 450K-reference data set and cell count-based validation of estimated cell type composition. Epigenetics. 2016 Sep 1;11(9):690-8}
#' @references{Reinius LE, Acevedo N, Joerink M, Pershagen G, Dahlén SE, Greco D, Söderhäll C, Scheynius A, Kere J. Differential DNA methylation in purified human blood cells: implications for cell lineage and studies on disease susceptibility. PloS one. 2012 Jul 25;7(7):e41361.}
#' @references{de Goede OM, Razzaghian HR, Price EM, Jones MJ, Kobor MS, Robinson WP, Lavoie PM. Nucleated red blood cells impact DNA methylation and expression analyses of cord blood hematopoietic cells. Clinical epigenetics. 2015 Dec;7(1):95.}
#' @references{Bakulski KM, Feinberg JI, Andrews SV, Yang J, Brown S, L. McKenney S, Witter F, Walston J, Feinberg AP, Fallin MD. DNA methylation of cord blood cell types: applications for mixed cell birth studies. Epigenetics. 2016 May 3;11(5):354-62.}
#' 
estimateLC <- function(meth,ref="bakulski+goede"){
    
    J = ncol(meth)

    coefs = read.table(system.file(paste0("data/",ref,".txt"),package="ewastools"))
    coefs = as.matrix(coefs)
    n_celltypes = ncol(coefs)

    markers = match(rownames(coefs),rownames(meth))
    EST = sapply(1:J,function(j){
        tmp = meth[markers,j]
        i = !is.na(tmp)
        quadprog::solve.QP(
             t(coefs[i,]) %*% coefs[i,]
            ,t(coefs[i,]) %*% tmp[i]
            ,diag(n_celltypes)
            ,rep(0,n_celltypes)
        )$sol
        })
    EST = t(EST)
    colnames(EST) = colnames(coefs)
    EST = data.table(EST)

    return(EST)
}
