#' @import data.table
NULL

#' @title Functions for preprocessing
#' @author Jonathan Heiss
#' @name Preprocessing
#' @rdname Preprocessing
#' @param raw Output of calling \code{\link{read_idats}}
#' @param tissue Optional. If set to "blood", predefined reference values are used for normalization. Recommened when cell proportions are to be estimated.
#'
#' @return A modified \code{raw} object with dye-bias corrected intensities using RELIC for \code{correct_dye_bias}. A matrix of beta-values, either normalized (for \code{normalize}) or not (for \code{dont_normalize}).
#' 
#' @references{Xu Z, Langie SA, De Boever P, Taylor JA, Niu L. RELIC: a novel dye-bias correction method for Illumina Methylation BeadChip. BMC Genomics. 2017 Jan 3;18(1):4.}
#' @references{Heiss JA, Brenner H. Between-array normalization for 450K data. Frontiers in Genetics. 2015;6.}
#' 
NULL

#' @rdname Preprocessing
#' @export
#' 
correct_dye_bias <- function(raw){
    
    if(!all(c('manifest','M','U','controls','ctrlG','ctrlR')%in%names(raw))) stop('Invalid argument')
    
    i1g = raw$manifest[channel=='Grn' ,]
    i2  = raw$manifest[channel=='Both',]

    ## Normalization control probes
    # They are paired (G with A, C with T)
    Ai = raw$controls[group=='NORM_A'][order(name)]$index
    Ti = raw$controls[group=='NORM_T'][order(name)]$index
    Ci = raw$controls[group=='NORM_C'][order(name)]$index
    Gi = raw$controls[group=='NORM_G'][order(name)]$index

    J = ncol(raw$M)

    for(j in 1:J){
        # Regress intensities in the red color channel on those in the green color channel
        # Relation is linear on the log scale
        x = log(raw$ctrlG[c(Gi,Ci),j])
        y = log(raw$ctrlR[c(Ai,Ti),j])

        keep = !is.na(y) & !is.na(x) & is.finite(x) & is.finite(y)
        x = x[keep]; y = y[keep]

        # Theil Sen robust linear regression (simple regression/single predictor only)
        m = mblm::mblm(y~x,repeated=FALSE)

        i = i2$index
        raw$M[i,j] = exp(stats::coef(m)[1] + log(raw$M[i,j]) * stats::coef(m)[2])

        i = i1g$index
        raw$M[i,j] = exp(stats::coef(m)[1] + log(raw$M[i,j]) * stats::coef(m)[2])
        raw$U[i,j] = exp(stats::coef(m)[1] + log(raw$U[i,j]) * stats::coef(m)[2])

        raw$oobG$U[,j] = exp(stats::coef(m)[1] + log(raw$oobG$U[,j]) * stats::coef(m)[2])
        raw$oobG$M[,j] = exp(stats::coef(m)[1] + log(raw$oobG$M[,j]) * stats::coef(m)[2])
    }

    return(raw)
}

#' @rdname Preprocessing
#' 
correct_dye_bias2 = function (raw) 
{
    # Experimental version of dye-bias correction. I found that red~green differs for 
    # ... G~A and C~T. They are however not independent. I have not understood this relation yet.
    if (!all(c("manifest", "M", "U", "controls", "ctrlG", "ctrlR") %in% names(raw)))  stop("Invalid argument")

    i1g = raw$manifest[channel == "Grn", ]
    i2 = raw$manifest[channel == "Both", ]

    Ai = raw$controls[group == "NORM_A"][order(name)]$index; Ai = raw$ctrlR[Ai,]
    Gi = raw$controls[group == "NORM_G"][order(name)]$index; Gi = raw$ctrlG[Gi,]
    Ti = raw$controls[group == "NORM_T"][order(name)]$index; Ti = raw$ctrlR[Ti,]
    Ci = raw$controls[group == "NORM_C"][order(name)]$index; Ci = raw$ctrlG[Ci,]
    
    J = ncol(raw$M)

    for(j in 1:J){
    
        x = log(Gi[,j])
        y = log(Ai[,j])
        keep = !is.na(y) & !is.na(x) & is.finite(x) & is.finite(y)
        x = x[keep]
        y = y[keep]
        m = mblm::mblm(y ~ x, repeated = FALSE)

        i = i2$index
        raw$M[i, j] = exp(stats::coef(m)[1] + log(raw$M[i, j]) * stats::coef(m)[2])
    
        x = log(Ci[,j])
        y = log(Ti[,j])
        keep = !is.na(y) & !is.na(x) & is.finite(x) & is.finite(y)
        x = x[keep]
        y = y[keep]
        m = mblm::mblm(y ~ x, repeated = FALSE)
    
        i = i1g$index
        raw$U[i, j] = exp(stats::coef(m)[1] + log(raw$U[i, j]) * stats::coef(m)[2])
        raw$M[i, j] = exp(stats::coef(m)[1] + log(raw$M[i, j]) * stats::coef(m)[2])

    }

    return(raw)
}

#' @rdname Preprocessing
#' @export
#'
dont_normalize <- function(raw){

    if(!all(c('manifest','M','U','meta')%in%names(raw))) stop('Invalid argument')

    with(raw,{

        M[M<1] = 1
        U[U<1] = 1
        
        meth = M/(M+U)
        
        rownames(meth) = manifest$probe_id
        colnames(meth) = meta$sample_id

        return(meth)
    })
}
