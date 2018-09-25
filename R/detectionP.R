#' @title  Detection p-values
#' @description Compute detection p-values. p-values are based on the distribution of the intensities of the negative control probes or the U (M) intensities observed for completely methylated (unmethylated) probes, respectively. \code{detP_threshold} generates a plot showing the number of undetected Y chromosome probes among male and female subjects for various p-value thresholds, in order to empirically choose a threshold. Finally, \code{mask} is masking all probes with detection p-values below the specified threshold. \code{detectionP.minfi} provides an implementation for \code{RGChannelSet} objects as used in the \code{minfi} package.
#' @author Jonathan A. Heiss
#' @param raw Output of calling \code{\link{read_idats}}, must include component \code{detP} for \code{mask} and \code{detP_threshold}.
#' @param threshold p-value threshold (arithmetic scale) above which oberservations are set to NA.
#' @param male/female Indices of male and female subjects
#' @param rgSet minfi rgSet object 
#' 
#' @return For \code{detectionP} a modified \code{raw} object with a \code{detP} component, a matrix of detection p-values, added.
#' @return For \code{detectionP.minfi} and \code{detectionP.neg} a matrix of detection p-values.
#' @return For \code{mask}, a modified \code{raw} object, with undetected probes set to \code{NA}.
#'

#' @rdname detectionP
#'
summits = function (beta)
{
    d <- density(beta,bw=0.01,na.rm=TRUE)

    l = which(d$x <0.4)
    u = which(d$x >0.6)

    l = l[which.max(d$y[l])]
    u = u[which.max(d$y[u])]

    d$x[c(l,u)]
}

#' @rdname detectionP
#' @export
#'
detectionP <- function(raw){

    if(!all(c('manifest','M','U')%in%names(raw))) stop('Invalid argument')

    with(raw,{

        detP = matrix(NA_real_,nrow=nrow(U),ncol=ncol(U))

        iR = manifest[channel=="Red" ,index]
        iG = manifest[channel=="Grn" ,index]
        i2 = manifest[channel=="Both",index]
        
        for(j in 1:ncol(M)){

            beta = M[,j]/(U[,j]+M[,j])
            
            sR = summits(beta[iR])
            sG = summits(beta[iG])

            bkgU = (beta[iR]-sR[2]) %>% abs %>% order %>% head(1000)
            bkgM = (beta[iR]-sR[1]) %>% abs %>% order %>% head(1000)

            bkgU = iR[bkgU]
            bkgM = iR[bkgM]

            muUR = median(U[bkgU,j],na.rm=TRUE)
            muMR = median(M[bkgM,j],na.rm=TRUE)

            sdUR = mad(U[bkgU,j],na.rm=TRUE)
            sdMR = mad(M[bkgM,j],na.rm=TRUE)

            bkgU = (beta[iG]-sG[2]) %>% abs %>% order %>% head(1000)
            bkgM = (beta[iG]-sG[1]) %>% abs %>% order %>% head(1000)

            bkgU = iG[bkgU]
            bkgM = iG[bkgM]

            muUG = median(U[bkgU,j],na.rm=TRUE)
            muMG = median(M[bkgM,j],na.rm=TRUE)

            sdUG = mad(U[bkgU,j],na.rm=TRUE)
            sdMG = mad(M[bkgM,j],na.rm=TRUE)

            detP[iR,j] = pnorm(U[iR,j]+M[iR,j],mean=muUR+muMR,sd=sqrt(sdUR^2+sdMR^2),lower.tail=FALSE)
            detP[iG,j] = pnorm(U[iG,j]+M[iG,j],mean=muUG+muMG,sd=sqrt(sdUG^2+sdMG^2),lower.tail=FALSE)
            detP[i2,j] = pnorm(U[i2,j]+M[i2,j],mean=muUR+muMG,sd=sqrt(sdUR^2+sdMG^2),lower.tail=FALSE)

        }
        
        raw$detP = detP
        return(raw)
    })
}

#' @rdname detectionP
#' @export
#'
detectionP.neg <- function(raw){
  
  with(raw,{

    bkgR = bkgG = controls[group=='NEGATIVE',index] 

    bkgR = ctrlR[bkgR,,drop=FALSE]
    bkgG = ctrlG[bkgG,,drop=FALSE]

    muG = apply(bkgG,2,median,na.rm=TRUE)
    sdG = apply(bkgG,2,mad   ,na.rm=TRUE)

    muR = apply(bkgR,2,median,na.rm=TRUE) 
    sdR = apply(bkgR,2,mad   ,na.rm=TRUE) 

    detP = matrix(NA_real_,nrow=nrow(U),ncol=ncol(U))

    i = manifest[channel=='Red' ,index] 
        for(j in 1:ncol(M))  
            detP[i,j] = pnorm(U[i,j]+M[i,j],mean=2*muR[j],sd=sqrt(2)*sdR[j],lower.tail=FALSE,log.p=TRUE) 
     
    i = manifest[channel=='Grn' ,index] 
        for(j in 1:ncol(M))  
            detP[i,j] = pnorm(U[i,j]+M[i,j],mean=2*muG[j],sd=sqrt(2)*sdG[j],lower.tail=FALSE,log.p=TRUE) 
     
    i = manifest[channel=='Both',index] 
        for(j in 1:ncol(M)) 
            detP[i,j] = pnorm(U[i,j]+M[i,j],mean=muR[j]+muG[j],sd=sqrt(sdR[j]^2+sdG[j]^2),lower.tail=FALSE,log.p=TRUE) 

    detP/log(10)
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
