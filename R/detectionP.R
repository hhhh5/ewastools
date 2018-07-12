#' @title  Detection p-values
#' @description Compute detection p-values. p-values are based on the distribution of the intensities of the negative control probes or out-of-band intensities. \code{detP_threshold} generates a plot showing the number of undetected Y chromosome probes among male and female subjects for various p-value thresholds, in order to empirically choose a threshold. Finally, \code{mask} is masking all probes with detection p-values below the specified threshold. \code{detectionP.minfi} provides an implementation for \code{RGChannelSet} objects as used in the \code{minfi} package.
#' @author Jonathan A. Heiss
#' @param raw Output of calling \code{\link{read_idats}}, must include component \code{detP} for \code{mask} and \code{detP_threshold}.
#' @param threshold p-value threshold (arithmetic scale) above which oberservations are set to NA.
#' @param male/female Indices of male and female subjects
#' @param rgSet minfi rgSet object 
#' 
#' @return For \code{detectionP}, a modified \code{raw} object with a \code{detP} component, a matrix of detection p-values, added.
#' @return For \code{detectionP.minfi} a matrix of detection p-values.
#' @return For \code{mask}, a modified \code{raw} object, with undetected probes set to \code{NA}.
#'
#' @rdname detectionP
#' @export
#' 
detectionP <- function(raw){

    if(!all(c('manifest','M','U','oobG','oobR')%in%names(raw))) stop('Invalid argument')

    with(raw,{

        bkgR = oobR$U + oobR$M
        bkgG = oobG$U + oobG$M

        muR = apply(bkgR,2,median,na.rm=TRUE)
        sdR = apply(bkgR,2,mad   ,na.rm=TRUE)

        muG = apply(bkgG,2,median,na.rm=TRUE)
        sdG = apply(bkgG,2,mad   ,na.rm=TRUE)

        detP = matrix(NA_real_,nrow=nrow(U),ncol=ncol(U))

         i = manifest[channel=='Red' ,index]
            for(j in 1:ncol(M)) 
                detP[i,j] = pnorm(U[i,j]+M[i,j],mean=muR[j],sd=sdR[j],lower.tail=FALSE)

        i = manifest[channel=='Grn' ,index]
            for(j in 1:ncol(M)) 
                detP[i,j] = pnorm(U[i,j]+M[i,j],mean=muG[j],sd=sdG[j],lower.tail=FALSE)
        
        i = manifest[channel=='Both',index]
            for(j in 1:ncol(M))
                detP[i,j] = pnorm(U[i,j]+M[i,j],mean=0.5*muR[j]+0.5*muG[j],sd=sqrt(sdR[j]^2+sdG[j]^2)/2,lower.tail=FALSE)

        raw$detP = detP

        return(raw)
    })
}

#' @rdname detectionP
#' @export
#'
mask <- function(raw,threshold){

    if(!all(c('M','U','detP')%in%names(raw))) stop('Invalid argument')
        
    raw$U[raw$detP>threshold] = NA_real_
    raw$M[raw$detP>threshold] = NA_real_

    return(raw)
}

#' @rdname detectionP
#' @export
#'
eval_detP_cutoffs = function(raw,males=NULL,females=NULL){

    if(is.null(males) | is.null(females)) stop('Please specify the column indices for male and female subjects')
    if(!'detP'%in%names(raw)) stop('detP component missing')

    chrY = raw$manifest[chr=='Y',index]
    chrY = raw$detP[chrY,]

    cutoffs = c(1,0.5,0.1,0.05,0.01,0.001,0.0001)

    tmp = sapply(cutoffs,function(t){ colSums(chrY>t,na.rm=TRUE) })
    males   = apply(tmp[males  ,],2,quantile,prob=0.9)
    females = apply(tmp[females,],2,quantile,prob=0.1)

    plot  (-log10(cutoffs),females,ylim=c(0,nrow(chrY)),ylab='Chr Y # undetected ',xlab='p-value cutoff',xaxt="n")
    points(-log10(cutoffs),males,pch=3)
    axis(1,at=-log10(cutoffs),labels=cutoffs)
    legend('topleft',pch=c(3,1),legend=c('Male 90% Quantile','Female 10% Quantile'))
    invisible(NULL)
}

#' @rdname detectionP
#' @export
#'
detectionP.minfi <- function(rgSet) {
    minfi:::.isRGOrStop(rgSet)
    locusNames <- getManifestInfo(rgSet, "locusNames")
    detP <- matrix(NA_real_, ncol = ncol(rgSet), nrow = length(locusNames),
                   dimnames = list(locusNames, colnames(rgSet)))

    r <- minfi::getRed(rgSet)
    g <- minfi::getGreen(rgSet)

    oob = minfi::getOOB(rgSet)
    rMu <- matrixStats::colMedians(oob$Red)
    rSd <- matrixStats::colMads(oob$Red)

    gMu <- matrixStats::colMedians(oob$Grn)
    gSd <- matrixStats::colMads(oob$Grn)


    TypeII <- minfi::getProbeInfo(rgSet, type = "II")
    TypeI.Red <- minfi::getProbeInfo(rgSet, type = "I-Red")
    TypeI.Green <- minfi::getProbeInfo(rgSet, type = "I-Green")

    for (i in 1:ncol(rgSet)) {   
        ## Type I Red
        intensity <- r[TypeI.Red$AddressA, i] + r[TypeI.Red$AddressB, i]
        detP[TypeI.Red$Name, i] <- pnorm(intensity, mean=rMu[i]*2, sd=rSd[i]*sqrt(2),lower.tail=FALSE)
        ## Type I Green
        intensity <- g[TypeI.Green$AddressA, i] + g[TypeI.Green$AddressB, i]
        detP[TypeI.Green$Name, i] <- pnorm(intensity, mean=gMu[i]*2, sd=gSd[i]*sqrt(2),lower.tail=FALSE)
        ## Type II
        intensity <- r[TypeII$AddressA, i] + g[TypeII$AddressA, i]
        detP[TypeII$Name, i] <- pnorm(intensity, mean=rMu[i]+gMu[i], sd=sqrt(rSd[i]^2+gSd[i]^2),lower.tail=FALSE)
    }
    
    detP
}
