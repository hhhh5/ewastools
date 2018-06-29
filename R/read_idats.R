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

    zipped = !file.exists(paste0(idat_files,"_Grn.idat"))
    suffix = rep(".idat",times=J)
    suffix[zipped] = ".idat.gz"

    ex = file.exists(paste0(idat_files,"_Grn",suffix)) & file.exists(paste0(idat_files,"_Red",suffix))
    if(!all(ex)) stop("Some .idat files are missing")

    # read itensities
    P = illuminaio::readIDAT(paste0(idat_files[1],"_Grn",suffix[1]))$nSNPsRead
    print(P)

    # create annotation
    if(P %in% c(1051815,1051943,1052641))
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

    setkeyv(manifest,c("chr","mapinfo"))
    manifest[,index:=1L:.N]
    controls[,index:=1L:.N]

    manifest[channel=="Grn",OOBi:=1:.N]
    manifest[channel=="Red",OOBi:=1:.N]

    M = U = matrix(NA_real_   ,nrow=nrow(manifest),ncol=J) # methylated (M) and unmethylated (U) signal intensities
    S = T = matrix(NA_integer_,nrow=nrow(manifest),ncol=J) # standard deviations
    N = V = matrix(NA_integer_,nrow=nrow(manifest),ncol=J) # number of beads underlying methylated (N) and unmethylated (V) signal intensities
    ctrlG = ctrlR = matrix(NA_real_,nrow=nrow(controls),ncol=J) # signal intensities of control probes
    ctrlN = matrix(NA_integer_,nrow=nrow(controls),ncol=J)
    oobG = list(M=matrix(NA_real_,nrow=manifest[channel=="Red",.N],ncol=J),U=matrix(NA_real_,nrow=manifest[channel=="Red",.N],ncol=J))
    oobR = list(M=matrix(NA_real_,nrow=manifest[channel=="Grn",.N],ncol=J),U=matrix(NA_real_,nrow=manifest[channel=="Grn",.N],ncol=J))


    if(!quiet) pb <- txtProgressBar(min=0,max=J,style=3)

    barcodes = rep(NA_character_,J)
    dates    = rep(NA_character_,J)

    for(j in 1:J){

        red = illuminaio::readIDAT(paste0(idat_files[j],"_Red",suffix[j]))
        grn = illuminaio::readIDAT(paste0(idat_files[j],"_Grn",suffix[j]))

        idat_order = red$MidBlock
        if(!identical(idat_order,grn$MidBlock)) stop("Red and green .idat files do not agree!")

        barcodes[j] = red$Barcode

        # This information is sometimes not recorded
        if(nrow(red$RunInfo)>1) dates[j] = red$RunInfo[2,1]

        manifest[,Ui:=match(addressU,idat_order)]
        manifest[,Mi:=match(addressM,idat_order)]
        controls[, i:=match(address ,idat_order)]

        setindexv(manifest,"channel")
        manifest["Both",Mi:=Ui,on="channel"]

        ### Type I Red probes
        i = manifest["Red",on="channel"]

        U[ i$index,j ] = red$Quants[ i$Ui,1 ] # Mean
        T[ i$index,j ] = red$Quants[ i$Ui,2 ] # SD
        V[ i$index,j ] = red$Quants[ i$Ui,3 ] # NBeads
    
        M[ i$index,j ] = red$Quants[ i$Mi,1 ]
        S[ i$index,j ] = red$Quants[ i$Mi,2 ]
        N[ i$index,j ] = red$Quants[ i$Mi,3 ]

        oobG$U[ i$OOBi,j ] = grn$Quants[ i$Ui,1 ]
        oobG$M[ i$OOBi,j ] = grn$Quants[ i$Mi,1 ]

        ### Type I Green probes
        i = manifest["Grn",on="channel"]

        U[ i$index,j ] = grn$Quants[ i$Ui,1 ]
        T[ i$index,j ] = grn$Quants[ i$Ui,2 ]
        V[ i$index,j ] = grn$Quants[ i$Ui,3 ]

        M[ i$index,j ] = grn$Quants[ i$Mi,1 ]
        S[ i$index,j ] = grn$Quants[ i$Mi,2 ]
        N[ i$index,j ] = grn$Quants[ i$Mi,3 ]

        oobR$U[ i$OOBi,j ] = red$Quants[ i$Ui,1 ]
        oobR$M[ i$OOBi,j ] = red$Quants[ i$Mi,1 ]
        
        ### Type II probes
        i = manifest["Both",on="channel"]

        U[ i$index,j ] = red$Quants[ i$Ui,1 ]
        T[ i$index,j ] = red$Quants[ i$Ui,2 ]
        V[ i$index,j ] = red$Quants[ i$Ui,3 ]

        M[ i$index,j ] = grn$Quants[ i$Mi,1 ]
        S[ i$index,j ] = grn$Quants[ i$Mi,2 ]
        N[ i$index,j ] = grn$Quants[ i$Mi,3 ]

        ctrlR[ controls$index,j ] = red$Quants[ controls$i,1 ]
        ctrlG[ controls$index,j ] = grn$Quants[ controls$i,1 ]
        ctrlN[ controls$index,j ] = red$Quants[ controls$i,3 ]

        if(!quiet) setTxtProgressBar(pb, j)

    }

    if(!quiet) close(pb)

    M[N==0] = NA
    U[V==0] = NA

    ctrlG[ctrlN==0] = NA
    ctrlR[ctrlN==0] = NA

    sample_ids = strsplit(x=idat_files,split="/")
    sample_ids = sapply(sample_ids,tail,n=1L)

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
        ,ctrlG=ctrlG,ctrlR=ctrlR,ctrlN=ctrlN
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
