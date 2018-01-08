#' Estimate leukcoyte composition
#' 
#' @export
#' @return Estimated cell proportions B-lymphocytes, CD4 T-cells, CD8 T-cells, Granulocytes, Monocytes and Natural killer cells using the Houseman algorithm. The model was trained on a dataset of both peripheral and cord blood, in order to apply to study populations of broad age range.
#' @references{Houseman EA, Accomando WP, Koestler DC, Christensen BC, Marsit CJ, Nelson HH, Wiencke JK, Kelsey KT. DNA methylation arrays as surrogate measures of cell mixture distribution. BMC bioinformatics. 2012 Dec 1;13(1):86.}
#' 
estimateLC <- function(meth){
    
    J = ncol(meth)

    coefs = read.table(system.file("data/blood_coefs.txt",package="ewastools"))
    coefs = as.matrix(coefs)

    markers = match(rownames(coefs),rownames(meth))
    EST = sapply(1:J,function(j){
        tmp = meth[markers,j]
        i = !is.na(tmp)
        quadprog::solve.QP(
             t(coefs[i,]) %*% coefs[i,]
            ,t(coefs[i,]) %*% tmp[i]
            ,diag(6)
            ,rep(0,6)
        )$sol
        })
    EST = t(EST)
    colnames(EST) = colnames(coefs)
    EST = data.table(EST)


    return(EST)
}
