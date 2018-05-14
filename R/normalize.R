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

    Ai = raw$controls[group=='NORM_A'][order(name)]$index
    Ti = raw$controls[group=='NORM_T'][order(name)]$index
    Ci = raw$controls[group=='NORM_C'][order(name)]$index
    Gi = raw$controls[group=='NORM_G'][order(name)]$index

    J = ncol(raw$M)

    for(j in 1:J){
        x = log(raw$ctrlG[c(Gi,Ci),j])
        y = log(raw$ctrlR[c(Ai,Ti),j])

        keep = !is.na(y) & !is.na(x) & is.finite(x) & is.finite(y)
        x = x[keep]; y = y[keep]

        m = mblm::mblm(y~x,repeated=FALSE)

        i = i2$index
        raw$M[i,j] = exp(coef(m)[1] + log(raw$M[i,j]) * coef(m)[2])

        i = i1g$index
        raw$M[i,j] = exp(coef(m)[1] + log(raw$M[i,j]) * coef(m)[2])
        raw$U[i,j] = exp(coef(m)[1] + log(raw$U[i,j]) * coef(m)[2])

        raw$oobG$U[,j] = exp(coef(m)[1] + log(raw$oobG$U[,j]) * coef(m)[2])
        raw$oobG$M[,j] = exp(coef(m)[1] + log(raw$oobG$M[,j]) * coef(m)[2])
    }

    return(raw)
}

#' @rdname Preprocessing
#' 
correct_dye_bias2 = function (raw) 
{
    if (!all(c("manifest", "M", "U", "controls", "ctrlG", "ctrlR") %in% names(raw)))  stop("Invalid argument")

    i1g = raw$manifest[channel == "Grn", ]
    i1r = raw$manifest[channel == "Red" & next_base=="T", ]
    i2 = raw$manifest[channel == "Both", ]

    Ai = raw$controls[group == "NORM_A"][order(name)]$index; Ai = raw$ctrlR[Ai,]
    Gi = raw$controls[group == "NORM_G"][order(name)]$index; Gi = raw$ctrlG[Gi,]
    Ti = raw$controls[group == "NORM_T"][order(name)]$index; Ti = raw$ctrlR[Ti,]
    Ci = raw$controls[group == "NORM_C"][order(name)]$index; Ci = raw$ctrlG[Ci,]
    
    J = ncol(raw$M)

    mm = sapply(1:J,function(j){
    
        x = log(Gi[,j])
        y = log(Ai[,j])
        keep = !is.na(y) & !is.na(x) & is.finite(x) & is.finite(y)
        x = x[keep]
        y = y[keep]
        m1 = mblm::mblm(y ~ x, repeated = FALSE)
    
        x = log(Ci[,j])
        y = log(Ti[,j])
        keep = !is.na(y) & !is.na(x) & is.finite(x) & is.finite(y)
        x = x[keep]
        y = y[keep]
        m2 = mblm::mblm(y ~ x, repeated = FALSE)
    
        c(coef(m1),coef(m2))

    })

    for(j in 1:J){
        Ci[,j] = exp(mm[1,j]) * Ci[,j]^mm[2,j]
        Gi[,j] = exp(mm[3,j]) * Gi[,j]^mm[4,j]

        (log(Ci[,j]) + log(Ti[,j])) %>% divide_by(2) %>% mean(na.rm=T) -> f1
        log(Ai[,j]) %>% mean(na.rm=T) -> f2
        (f2/f1) %>% exp -> f

        i = i2$index
        raw$M[i, j] = exp(mm[3,j]) * raw$M[i, j]^mm[4,j]

        i = i1g$index
        raw$M[i, j] = f * exp(mm[1,j]) * raw$M[i, j]^mm[2,j]
        raw$U[i, j] = f * exp(mm[1,j]) * raw$U[i, j]^mm[2,j]
        
        i = i1r$index
        raw$M[i, j] = f * raw$M[i, j]
        raw$U[i, j] = f * raw$U[i, j]
    }

    return(raw)
}

#' @rdname Preprocessing
#' @export
#'
normalize <- function(raw,tissue=''){

    if(!all(c('manifest','M','U','N','V','controls','ctrlG','ctrlR','meta')%in%names(raw))) stop('Invalid argument')

    with(raw,{
        J = ncol(M)
        if(J<2) stop('More than one sample required to perform "between-array" normalization')

        i1g = manifest[channel=='Grn' ,]
        i1r = manifest[channel=='Red' ,]
        i2  = manifest[channel=='Both',]

        M[M<1] = 1
        U[U<1] = 1

        ### compute the reference values
        if(tissue=='blood'){

            hk = fread(system.file('data/blood_norm.txt',package='ewastools'))
            hk = merge(hk,manifest[,list(probe_id,channel,index)])

            hk_1g = hk[channel=='Grn' ,]
            hk_1r = hk[channel=='Red' ,]
            hk_2  = hk[channel=='Both',]

            rug = hk_1g$refU
            rmg = hk_1g$refM
            rur = hk_1r$refU
            rmr = hk_1r$refM

        }else{

            hk = scan(system.file('data/hkcpgs.txt',package='ewastools'),what='',quiet=TRUE)
            hk = manifest[probe_id%in%hk]

            hk_1g = hk[channel=='Grn' ,]
            hk_1r = hk[channel=='Red' ,]
            hk_2  = hk[channel=='Both',]

            rmg = log(apply(M[hk_1g$index,],1,median,na.rm=TRUE))
            rug = log(apply(U[hk_1g$index,],1,median,na.rm=TRUE))
            rmr = log(apply(M[hk_1r$index,],1,median,na.rm=TRUE))
            rur = log(apply(U[hk_1r$index,],1,median,na.rm=TRUE))
        }

        cat('[Correcting intensity-dependent bias]\n')
        pb <- txtProgressBar(min=0,max=J,style=3)

        ### correct intensity-dependent bias
        for(j in 1:J){
            # green channel
            x = log(c(M[hk_1g$index,j],U[hk_1g$index,j]))
            y = x-c(rmg,rug)

            omit = !is.na(y)
            x = x[omit]
            y = y[omit]

            f = loess(y~x,span=.4,degree=1,family='symmetric')

            x = M[i1g$index,j]; x = x * exp(-predict(f,log(x))); M[i1g$index,j] <- ifelse(is.na(x),M[i1g$index,j],x)
            x = U[i1g$index,j]; x = x * exp(-predict(f,log(x))); U[i1g$index,j] <- ifelse(is.na(x),U[i1g$index,j],x)
            x = M[ i2$index,j]; x = x * exp(-predict(f,log(x))); M[ i2$index,j] <- ifelse(is.na(x),M[ i2$index,j],x)

            # red channel
            x = log(c(M[hk_1r$index,j],U[hk_1r$index,j]))
            y = x-c(rmr,rur)

            omit = !is.na(y)
            x = x[omit]
            y = y[omit]

            f = loess(y~x,span=.4,degree=1,family='symmetric')

            x = M[i1r$index,j]; x = x * exp(-predict(f,log(x))); M[i1r$index,j] <- ifelse(is.na(x),M[i1r$index,j],x)
            x = U[i1r$index,j]; x = x * exp(-predict(f,log(x))); U[i1r$index,j] <- ifelse(is.na(x),U[i1r$index,j],x)
            x = U[ i2$index,j]; x = x * exp(-predict(f,log(x))); U[ i2$index,j] <- ifelse(is.na(x),U[ i2$index,j],x)
            setTxtProgressBar(pb, j)
        }
        close(pb)

        meth = log(M/U)
        rm(M,U,rmg,rmr,rug,rur)

        ### correct methylation-dependent bias
        cat('[Correcting methylation-dependent bias]\n')
        pb <- txtProgressBar(min=0,max=J,style=3)

        if(tissue=='blood'){
            rb = hk_2$refU
        }else{
            rb = apply(meth[hk_2$index,],1,median,na.rm=TRUE)
        }

        for(j in 1:J){
            x = meth[hk_2$index,j]
            y = x - rb

            omit = !is.na(y)
            x = x[omit]
            y = y[omit]

            f = loess(y~x,span=.4,family='symmetric',degree=1)
            x = meth[i2$index,j]
            x = x - predict(f,x)
            meth[i2$index,j] <- ifelse(is.na(x),meth[i2$index,j],x)
            setTxtProgressBar(pb,j)
        }

        close(pb)
        rm(rb)
   
        meth = exp(meth)
        meth = meth/(meth+1)

        rownames(meth) = manifest$probe_id
        colnames(meth) = meta$sample_id

        return(meth)
        })
}

#' @rdname Preprocessing
#' @export
#'
dont_normalize <- function(raw){

    if(!all(c('manifest','M','U','N','V','meta')%in%names(raw))) stop('Invalid argument')

    with(raw,{

        M[M<1] = 1
        U[U<1] = 1
        
        meth = M/(M+U)
        
        rownames(meth) = manifest$probe_id
        colnames(meth) = meta$sample_id

        return(meth)
    })
}
