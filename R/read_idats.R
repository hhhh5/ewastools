#' @import data.table
NULL

#' Read .idat files
#' 
#' @export
#' 
#' @param idat_files Character vector of relative or absolute filepaths, but without the suffixes '_Grn.idat' and '_Red.idat'
#' @param quiet If TRUE, a progress bar is shown.
#' 
#' @return A list containing
#' \item{manifest}{A data.table describing the probes}
#' \item{M}{intensities of targeting methylated sequences}
#' \item{U}{intensities of targeting unmethylated sequences}
#' \item{N}{number of beads from which average intensities in M were derived}
#' \item{V}{number of beads from which average intensities in U were derived}
#' \item{ctrlG}{Intensities of the control probes in the green color channel}
#' \item{ctrlR}{Intensities of the control probes in the red color channel}
#' \item{meta}{A data.table containing unique sample_ids and metadata}
#' 
#' @examples
#' \dontrun{
#' read_idats('9976861004_R01C01')
#' }
#'
read_idats <- function(idat_files,quiet=FALSE){

    J = length(idat_files)

    zipped = file.exists(paste0(idat_files[1],"_Grn.idat.gz"))
    if(zipped){ suffix = ".idat.gz" } else suffix = ".idat"

    ex = file.exists( paste0(idat_files,"_",rep(c("Grn","Red"),each=J),suffix))
    if(!all(ex)) stop("Some .idat files are missing")

    # read itensities
    idat_order = illuminaio::readIDAT(paste0(idat_files[1],"_Grn",suffix))$MidBlock
    P = length(idat_order) # number of features
    print(P)

    # create annotation
    if(P==1052641)
    {
        platform="EPIC"
        manifest = data.table::copy(ewastools:::manifest_epic)
        controls = data.table::copy(ewastools:::controls_epic)
    }
    else if(P==622399)
    { 
        platform="450K"
        manifest = data.table::copy(ewastools:::manifest_450K)
        controls = data.table::copy(ewastools:::controls_450K)
    }
    else
    {
        stop("Unknown platform")
        platform="UNK"
    }

    manifest[channel!="Both",Ui:=match(addressU,idat_order)]
    manifest[channel!="Both",Mi:=match(addressM,idat_order)]
    manifest[channel=="Both",Ui:=match(addressU,idat_order)]
    manifest[channel=="Both",Mi:=Ui]
    
    controls[,i:=match(address,idat_order)]

    manifest[,index:=1L:.N]
    controls[,index:=1L:.N]

    manifest[channel=="Grn",OOBi:=1:.N]
    manifest[channel=="Red",OOBi:=1:.N]

    ### indices of probes by probe type and channel channel
    i1g = manifest[channel=="Grn" ]
    i1r = manifest[channel=="Red" ]
    i2  = manifest[channel=="Both"]

    sample_ids = strsplit(x=idat_files,split="/")
    sample_ids = sapply(sample_ids,tail,n=1L)

    M = U = matrix(NA_real_   ,nrow=nrow(manifest),ncol=J) # methylated (M) and unmethylated (U) signal intensities
    S = T = matrix(NA_integer_,nrow=nrow(manifest),ncol=J) # standard deviations
    N = V = matrix(NA_integer_,nrow=nrow(manifest),ncol=J) # number of beads underlying methylated (N) and unmethylated (V) signal intensities
    ctrlG = ctrlR = matrix(NA_real_,nrow=nrow(controls),ncol=J) # signal intensities of control probes
    oobG = list(M=matrix(NA_real_,nrow=nrow(i1r),ncol=J),U=matrix(NA_real_,nrow=nrow(i1r),ncol=J))
    oobR = list(M=matrix(NA_real_,nrow=nrow(i1g),ncol=J),U=matrix(NA_real_,nrow=nrow(i1g),ncol=J))


    if(!quiet) pb <- txtProgressBar(min=0,max=2*J,style=3)

    barcodes = rep(NA_character_,J)
    dates    = rep(NA_character_,J)

    ### read the intensities of the green channel
    for(j in 1:J){
        tmp = illuminaio::readIDAT(paste0(idat_files[j],"_Grn",suffix))
        if(!identical(tmp$MidBlock,idat_order)) stop("Different versions of .idat files")
        
        barcodes[j] = tmp$Barcode

        # This information is sometimes not recorded
        if(nrow(tmp$RunInfo)>1) dates[j] = tmp$RunInfo[2,1]

        tmp = tmp$Quants
        means  = tmp[,"Mean"]
        sds    = tmp[,"SD"]
        nbeads = tmp[,"NBeads"]
        
        M[i1g$index,j] = means [i1g$Mi]
        S[i1g$index,j] = sds   [i1g$Mi]
        N[i1g$index,j] = nbeads[i1g$Mi]
        
        U[i1g$index,j] = means [i1g$Ui]
        T[i1g$index,j] = sds   [i1g$Ui]
        V[i1g$index,j] = nbeads[i1g$Ui]

        M[i2$index ,j] = means [ i2$Mi]
        S[i2$index ,j] = sds   [ i2$Mi]
        N[i2$index ,j] = nbeads[ i2$Mi]

        ctrlG[,j]    = means[controls$i]

        oobG$M[i1r$OOBi,j] = means[i1r$Mi]
        oobG$U[i1r$OOBi,j] = means[i1r$Ui]

        if(!quiet) setTxtProgressBar(pb, j)
    }

    ### the red channel
    for(j in 1:J){
        tmp = illuminaio::readIDAT(paste0(idat_files[j],"_Red",suffix))
        if(!identical(tmp$MidBlock,idat_order)) stop("Different versions of .idat files")
        
        tmp = tmp$Quants
        means  = tmp[,"Mean"]
        sds    = tmp[,"SD"]
        nbeads = tmp[,"NBeads"]

        M[i1r$index,j] = means [i1r$Mi]
        S[i1r$index,j] = sds   [i1r$Mi]
        N[i1r$index,j] = nbeads[i1r$Mi]

        U[i1r$index,j] = means [i1r$Ui]
        T[i1r$index,j] = sds   [i1r$Ui]
        V[i1r$index,j] = nbeads[i1r$Ui]

        U[i2$index ,j] = means [ i2$Ui]
        T[i2$index ,j] = sds   [ i2$Ui]
        V[i2$index ,j] = nbeads[ i2$Ui]

        ctrlR[,j]    = means[controls$i]

        oobR$M[i1g$OOBi,j] = means[i1g$Mi]
        oobR$U[i1g$OOBi,j] = means[i1g$Ui]

        if(!quiet) setTxtProgressBar(pb, j+J)
    }

    if(!quiet) close(pb)

    M[N==0] = NA
    U[V==0] = NA

    meta = data.table(
         sample_id = sample_ids
        ,date = as.IDate(dates,"%m/%d/%Y %r")
        ,time = as.ITime(dates,"%m/%d/%Y %r")
        )

    raw = list(
         platform=platform
        ,manifest=manifest
        ,M=M,S=S,U=U
        ,N=N,T=T,V=V
        ,controls=controls
        ,ctrlG=ctrlG,ctrlR=ctrlR
        ,oobG=oobG,oobR=oobR
        ,meta=meta
    )

    return(raw)
}

#' Drop samples from raw data.
#' @author Jonathan A. Heiss
#'
#' @export
#'
#' @param raw Output of calling \code{\link{read_idats}}
#' @param j Indices of the samples to drop
#'
#' @return A modified \code{raw} object
#'
drop_samples <- function(raw,j=NULL){

    raw$M = raw$M[,-j,drop=FALSE]; raw$N = raw$N[,-j,drop=FALSE]
    raw$U = raw$U[,-j,drop=FALSE]; raw$V = raw$V[,-j,drop=FALSE]

    raw$ctrlG = raw$ctrlG[,-j,drop=FALSE]
    raw$ctrlR = raw$ctrlR[,-j,drop=FALSE]

    raw$oobG$M = raw$oobG$M[,-j,drop=FALSE]
    raw$oobG$U = raw$oobG$U[,-j,drop=FALSE]
    raw$oobR$M = raw$oobR$M[,-j,drop=FALSE]
    raw$oobR$U = raw$oobR$U[,-j,drop=FALSE]

    raw$meta = raw$meta[-j,,drop=FALSE]

    if("detP"%in%names(raw)){
        raw$detP = raw$detP[,-j,drop=FALSE]
    }

    return(raw)
}
