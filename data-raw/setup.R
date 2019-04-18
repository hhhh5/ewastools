library(data.table)
library(stringi)
library(magrittr)
library(parallel)
library(purrr)
library(forcats)
library(ewastools)
library(strict)
library(minfi)


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
    
    raw$detP = detP/log(10)
    return(raw)
  })
}


# -------------------------------- Purified blood cells 

# 450K datasets

#--------------------------------------------------------------
# de Goede (GSE68456) https://doi.org/10.1186/s13148-015-0129-6
#--------------------------------------------------------------

dir.create("GSE68456")

### Select datasets by GSE accession 
goede = "https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE68456&targ=gsm&form=text&view=brief"
goede %<>% map(readLines) %>% unlist
goede = split( goede, cumsum(goede %like% "^\\^SAMPLE = GSM") )

names(goede) =  map(goede,1) %>% stri_match_first(regex="GSM\\d+")

goede %<>% imap(function(s,acc){
  s = strsplit(s,split=" = ",fixed=TRUE)  
  data.table(gsm=acc,variable=map_chr(s,1),value=map_chr(s,2))
})

goede %<>% rbindlist

goede = goede[variable %chin% c("!Sample_characteristics_ch1","!Sample_supplementary_file","!Sample_title")]

i = goede[variable == "!Sample_characteristics_ch1",which=TRUE]
ch = goede$value[i] %>% stri_split(fixed=": ")
goede$variable[i] = map_chr(ch,1)
goede$value   [i] = map_chr(ch,2)
rm(ch)

goede[variable == "!Sample_title",variable:="title"    ]
goede[variable == "Sex"          ,variable:="sex"      ]
goede[variable == "cell type"    ,variable:="cell_type"]
goede = goede[!variable %in% c("facs strategy","tissue")]

# Find the URLs pointing to the two .idat files
goede[variable == "!Sample_supplementary_file" & value %like% "_Red\\.idat",variable:="red"]
goede[variable == "!Sample_supplementary_file" & value %like% "_Grn\\.idat",variable:="grn"]

goede = dcast(goede,gsm ~ variable,fill=NA)
goede[,sex:=factor(sex,levels=c("M","F"),labels=c("m","f"))]
goede[,file:=paste0("GSE68456","/",gsm)]
goede[,cell_type:=fct_recode(cell_type,B="B cell",CD4="CD4 T cell",CD8="CD8 T cell",GR="Granulocyte",MO="Monocyte",NK="NK cell",nRBC="nRBC")]
goede[,donor:=stri_extract_last(title,regex="ind\\d+")]

map2(goede$red, goede$file %s+% "_Red.idat.gz", ~ download.file(.x,.y) ) %>% invisible
map2(goede$grn, goede$file %s+% "_Grn.idat.gz", ~ download.file(.x,.y) ) %>% invisible

goede = goede[,list(study="goede",cell_type,donor,sex,file)]

#-----------------------------------------------------
# Gervin https://doi.org/10.1080/15592294.2016.1214782
#-----------------------------------------------------

gervin = read.csv("gervin/sampleSheet_cordblood.csv",skip=7)
setDT(gervin)
 
gervin[,cell_type:=fct_recode(cellType,B="cd19",CD4="cd4",CD8="cd8",GR="gran",MO="cd14",WB="wb")]
gervin[,donor:=stri_extract_first(Sample_Name,regex="^\\d+")]
gervin[,file:=paste0("gervin","/",Sentrix_ID,"_",Sentrix_Position)]

gervin = gervin[,list(study="gervin",cell_type,donor,file)]

#----------------------------------------------------
# Reinius https://doi.org/10.1371/journal.pone.004136
#----------------------------------------------------

reinius = fread("reinius/sample_sheet_IDAT.csv")

reinius[,donor:=stri_sub(Sample,-3,-1)]
reinius[,`Chip Row Pos`:=stri_replace_all_fixed(`Chip Row Pos`,"O","0")] 
reinius[,`Chip CO position`:=stri_replace_all_fixed(`Chip CO position`,"O","0")] 
reinius[,file:=paste0("reinius/",`Chip#ID`,"_",`Chip Row Pos`,`Chip CO position`)] 
reinius$cell_type = fct_recode(reinius$Type,MO="CD14+ Monocytes",NK="CD56+ NK-cells",CD8="CD8+ T-cells",CD4="CD4+ T-cells",GR="Granulocytes",B="CD19+ B-cells",WB="Whole blood",EO="Eosinophils",NE="Neutrophils")
reinius$sex = "m"
reinius = reinius[,list(study="reinius",cell_type,donor,sex,file)]

#----------------------------------------------------

purified = rbind(goede,gervin,reinius,use.names=TRUE,fill=TRUE)
purified %<>% droplevels
purified[,j:=1:.N]

beta = purified$file %>% read_idats(.) %>% correct_dye_bias(.) %>% ewastools::detectionP(.) %>% ewastools::mask(.,0.05) %>% dont_normalize

#----------------------------------------------------
# Bakulski
#----------------------------------------------------

load("bakulski/FlowSorted.CordBlood.450k.rda") 
 
pheno = pData(FlowSorted.CordBlood.450k)
pheno %<>% as.data.frame %>% setDT
J = nrow(pheno)

pheno$cell_type = fct_recode(pheno$CellType,MO="Mono",NK="NK",CD8="CD8T",CD4="CD4T",GR="Gran",B="Bcell",WB="WholeBlood",nRBC="nRBC")
pheno[,sex:=factor(Sex,levels=c("M","F"),labels=c("m","f"))]
pheno = pheno[,list(study="bakulski",cell_type,donor=Individual.ID,sex,file=NA,j=(nrow(purified)+1):(nrow(purified)+J))]

purified = rbind(purified,pheno)


r = getRed  (FlowSorted.CordBlood.450k)
g = getGreen(FlowSorted.CordBlood.450k)
idat_order = rownames(r)

meth = list()
meth$platform = "450K"
manifest = copy(ewastools:::manifest_450K)
controls = copy(ewastools:::controls_450K)

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

M = U = matrix(NA_real_   ,nrow=nrow(manifest),ncol=J) # methylated (M) and unmethylated (U) signal intensities
ctrlG = ctrlR = matrix(NA_real_,nrow=nrow(controls),ncol=J) # signal intensities of control probes
oobG = list(M=matrix(NA_real_,nrow=nrow(i1r),ncol=J),U=matrix(NA_real_,nrow=nrow(i1r),ncol=J))
oobR = list(M=matrix(NA_real_,nrow=nrow(i1g),ncol=J),U=matrix(NA_real_,nrow=nrow(i1g),ncol=J))

### the red channel
M[i1r$index,]     = r[i1r$Mi,]
U[i1r$index,]     = r[i1r$Ui,]
U[i2$index ,]     = r[ i2$Ui,]
ctrlR             = r[controls$i,]
oobR$M[i1g$OOBi,] = r[i1g$Mi,]
oobR$U[i1g$OOBi,] = r[i1g$Ui,]

### green channel
M[i1g$index,]     = g[i1g$Mi,]
U[i1g$index,]     = g[i1g$Ui,]
M[i2$index ,]     = g[ i2$Mi,]
ctrlG             = g[controls$i,]
oobG$M[i1r$OOBi,] = g[i1r$Mi,]
oobG$U[i1r$OOBi,] = g[i1r$Ui,]

M[M==0] = NA
U[U==0] = NA

ctrlG[ctrlG==0] = NA
ctrlR[ctrlR==0] = NA

meth$manifest = manifest
meth$M = M
meth$U = U
meth$controls = controls
meth$ctrlR = ctrlR
meth$ctrlG = ctrlG
meth$oobR = oobR
meth$oobG = oobG
meth$meta = data.table(sample_id="BB" %s+% 1:104)
rm(M,U,oobG,oobR,controls,ctrlR,ctrlG,manifest,pheno)

beta2 = meth %>% correct_dye_bias(.) %>% ewastools::detectionP(.) %>% ewastools::mask(.,0.05) %>% dont_normalize(.)

beta = cbind(beta,beta2)

# --------------------------------

train_model = function(studies,cell_types,output){
    
    train = purified[study %in% studies & cell_type %in% cell_types]
    train_beta = na.omit(beta[,train$j])

    markers = list() 

    for(ct in cell_types){ 
        cat(ct,'\n') 
     
        j = train$cell_type == ct
     
        tmp = apply(train_beta,1,function(x){ 
            if(!any(is.na(x[j]))) 
            {  
                tmp = t.test(x[j],x[!j],var.equal=T) 
                return(c(tmp$p.value,tmp$estimate[1]-tmp$estimate[2])) 
            }else{ 
                return(c(NA,NA)) 
            } 
        }) 
     
        i = which(p.adjust(tmp[1,])<0.05) 
        o = order(tmp[2,i],na.last=NA)
        markers[[ct]] = c(i[head(o,50)],i[tail(o,50)]) 
    } 

    markers %<>% unlist %>% unique 

    coefs = sapply(cell_types,function(ct){ 
        j = train$cell_type == ct
        rowMeans(train_beta[markers,j],na.rm=TRUE) 
    }) 
    rownames(coefs) = rownames(train_beta)[markers] 
    colnames(coefs) = cell_types 

    write.table(coefs,file=output,row.names=TRUE) 
}

### First model using only the Gervin dataset
train_model("gervin" ,c("GR","MO","B","CD4","CD8","NK"),"../data/gervin.txt")
train_model("reinius",c("GR","MO","B","CD4","CD8","NK"),"../data/reinius.txt")
train_model(c("gervin","reinius"),c("GR","MO","B","CD4","CD8","NK"),"../data/gervin+reinius.txt")
train_model("goede",c("GR","MO","B","CD4","CD8","NK","nRBC"),"../data/goede.txt")
train_model("bakulski",c("GR","MO","B","CD4","CD8","NK","nRBC"),"../data/bakulski.txt")
train_model(c("bakulski","goede"),c("GR","MO","B","CD4","CD8","NK","nRBC"),"../data/bakulski+goede.txt")
train_model(c("bakulski","gervin"),c("GR","MO","B","CD4","CD8","NK"),"../data/bakulski+gervin.txt")
train_model(c("bakulski","reinius"),c("GR","MO","B","CD4","CD8","NK"),"../data/bakulski+reinius.txt")
train_model(c("goede","reinius"),c("GR","MO","B","CD4","CD8","NK"),"../data/goede+reinius.txt")

# --------------------------------
# EPIC reference dataset https://doi.org/10.1186/s13059-018-1448-7

dir.create("GSE110554")

### Select datasets by GSE accession 
salas = "https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE110554&targ=gsm&form=text&view=brief"
salas %<>% map(readLines) %>% unlist
salas = split( salas, cumsum(salas %like% "^\\^SAMPLE = GSM") )

names(salas) =  map(salas,1) %>% stri_match_first(regex="GSM\\d+")

salas %<>% imap(function(s,acc){
  s = strsplit(s,split=" = ",fixed=TRUE)  
  data.table(gsm=acc,variable=map_chr(s,1),value=map_chr(s,2))
})

salas %<>% rbindlist

salas = salas[variable %chin% c("!Sample_characteristics_ch1","!Sample_supplementary_file")]

i = salas[variable == "!Sample_characteristics_ch1",which=TRUE]
ch = salas$value[i] %>% stri_split(fixed=": ")
salas$variable[i] = map_chr(ch,1)
salas$value   [i] = map_chr(ch,2)
rm(ch)

# Find the URLs pointing to the two .idat files
salas[variable == "!Sample_supplementary_file" & value %like% "_Red\\.idat",variable:="red"]
salas[variable == "!Sample_supplementary_file" & value %like% "_Grn\\.idat",variable:="grn"]

salas[,variable:=fct_recode(variable,B="bcell",CD4="cd4t",CD8="cd8t",NE="neu",MO="mono",NK="nk")]
salas = dcast(salas,gsm ~ variable,fill=NA)
salas = salas[,list(gsm,B,CD4,CD8,MO,NE,NK,purity,cell_type=`cell type`,red,grn)]

salas[,cell_type:=fct_recode(cell_type,B="Bcell",CD4="CD4T",CD8="CD8T",NE="Neu",MO="Mono",NK="NK")]
salas[,file:=paste0("GSE110554","/",gsm)]

map2(salas$red, salas$file %s+% "_Red.idat.gz", ~ download.file(.x,.y) ) %>% invisible
map2(salas$grn, salas$file %s+% "_Grn.idat.gz", ~ download.file(.x,.y) ) %>% invisible

beta = salas %$% file %>% read_idats %>% correct_dye_bias %>% ewastools::detectionP(.) %>% ewastools::mask(.,0.05) %>% dont_normalize

purified = salas
purified$study = "salas"
purified[,j:=1:.N]

# --------------------------------

train_model("salas",c("B","CD4","CD8","MO","NE","NK"),"../data/salas.txt")



