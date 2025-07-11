library(data.table)
library(stringi)
library(magrittr)
library(parallel)
library(forcats)
library(ewastools)
library(purrr)

correct_dye_bias = function(raw){
    
    if(!all(c("manifest","M","U","controls","ctrlG","ctrlR")%in%names(raw))) stop("Invalid argument")
    
    i1g = raw$manifest[channel=="Grn" ,]
    i2  = raw$manifest[channel=="Both",]

    Ai = raw$controls[group=="NORM_A"][order(name)]$index
    Ti = raw$controls[group=="NORM_T"][order(name)]$index
    Ci = raw$controls[group=="NORM_C"][order(name)]$index
    Gi = raw$controls[group=="NORM_G"][order(name)]$index

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
    }

    return(raw)
}

# -------------------------------- Purified blood cells 

load("FlowSorted.CordBloodCombined.450k.rda") 

pheno1 = minfi::pData(FlowSorted.CordBloodCombined.450k)
pheno1 %<>% as.data.table
J = nrow(pheno1)

pheno1[,cell_type:=fct_recode(CellType,`B`="Bcell",`CD4`="CD4T",`CD8`="CD8T",`GR`="Gran",`MO`="Mono",`NK`="NK",`WBC`="WBC")]
pheno1[,sex:=fct_recode(Sex,`m`="M",`f`="F")]
pheno1[,sample_id:=paste0(Slide,"_",Array)]

pheno1 = pheno1[,.(
     sample_id
    ,donor=SampleID
    ,sex
    ,cell_type
    ,tissue="Blood"
    ,study=Study
    )]

r = minfi::getRed  (FlowSorted.CordBloodCombined.450k)
g = minfi::getGreen(FlowSorted.CordBloodCombined.450k)
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
meth$meta = copy(pheno1[,.(sample_id)])

rm(M,U,oobG,oobR,controls,ctrlR,ctrlG,manifest,r,g,i2,i1g,i1r,idat_order,FlowSorted.CordBloodCombined.450k)

beta1 =
    meth %>%
    ewastools::detectionP.neg(.) %>%
    ewastools::mask(-2) %>%
    correct_dye_bias %>%
    dont_normalize

rm(meth)

#----------------------------------------------------
# Reinius https://doi.org/10.1371/journal.pone.004136
#----------------------------------------------------

reinius = fread("reinius/sample_sheet_IDAT.csv")

reinius[,donor:=stri_sub(Sample,-3,-1)]
reinius[,`Chip Row Pos`:=stri_replace_all_fixed(`Chip Row Pos`,"O","0")] 
reinius[,`Chip CO position`:=stri_replace_all_fixed(`Chip CO position`,"O","0")] 
reinius[,sample_id:=paste0(`Chip#ID`,"_",`Chip Row Pos`,`Chip CO position`)] 
reinius$cell_type = fct_recode(reinius$Type,MO="CD14+ Monocytes",NK="CD56+ NK-cells",CD8="CD8+ T-cells",CD4="CD4+ T-cells",GR="Granulocytes",B="CD19+ B-cells",WB="Whole blood",EO="Eosinophils",NE="Neutrophils")
reinius$sex = "m"
reinius = reinius[,.(sample_id,donor,sex,cell_type,tissue="Blood",study="Reinius")]
reinius = reinius[cell_type %in% c("MO","NK","CD8","GR","CD4","B")]
reinius = reinius[donor != "105"]
reinius = reinius[! (donor=="160" & cell_type %in% c("CD8","CD4"))]
reinius = reinius[! (donor=="261" & cell_type %in% c("NK","B"))]
reinius = reinius[! (donor=="218" & cell_type=="NK")]

reinius[,file:=paste0("reinius/",sample_id)]

beta2 =
    reinius$file %>%
    ewastools::read_idats(.) %>%
    ewastools::detectionP.neg(.) %>%
    ewastools::mask(-2) %>%
    correct_dye_bias %>%
    dont_normalize

#----------------------------------------------------

### EPIC datasets

# https://doi.org/10.1186/s13059-018-1448-7
if(!dir.exists("GSE110554")) dir.create("GSE110554")

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

salas = dcast(salas,gsm ~ variable,fill=NA)
salas = salas[,.(
     sample_id=gsm
    ,sex=factor(Sex,levels=c("M","F"),labels=c("m","f"))
    ,cell_type = fct_recode(`cell type`,B="Bcell",CD4="CD4T",CD8="CD8T",GR="Neu",MO="Mono",NK="NK")
    ,tissue = "Blood"
    ,study = "Salas"
    ,red,grn
    )]

# map2(salas$red, salas$file %s+% "_Red.idat.gz", ~ download.file(.x,.y) ) %>% invisible
# map2(salas$grn, salas$file %s+% "_Grn.idat.gz", ~ download.file(.x,.y) ) %>% invisible

salas$red = NULL
salas$grn = NULL

salas[,file:=paste0("GSE110554","/",sample_id)]
salas = salas[cell_type!="MIX"]

# --------------------------------
if(!dir.exists("GSE103541")) dir.create("GSE103541")

### Select datasets by GSE accession 
mill = "https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE103541&targ=gsm&form=text&view=brief"
mill %<>% map(readLines) %>% unlist
mill = split( mill, cumsum(mill %like% "^\\^SAMPLE = GSM") )

names(mill) =  map(mill,1) %>% stri_match_first(regex="GSM\\d+")

mill %<>% imap(function(s,acc){
  s = strsplit(s,split=" = ",fixed=TRUE)  
  data.table(gsm=acc,variable=map_chr(s,1),value=map_chr(s,2))
})

mill %<>% rbindlist

mill = mill[variable %chin% c("!Sample_characteristics_ch1","!Sample_supplementary_file")]

i = mill[variable == "!Sample_characteristics_ch1",which=TRUE]
ch = mill$value[i] %>% stri_split(fixed=": ")
mill$variable[i] = map_chr(ch,1)
mill$value   [i] = map_chr(ch,2)
rm(ch)

# Find the URLs pointing to the two .idat files
mill[variable == "!Sample_supplementary_file" & value %like% "_Red\\.idat",variable:="red"]
mill[variable == "!Sample_supplementary_file" & value %like% "_Grn\\.idat",variable:="grn"]

mill = mill[value!="blood"] # variable `cell type` appears twice
mill = dcast(mill,gsm ~ variable,fill=NA)

mill = mill[,.(
     sample_id = gsm
    ,donor = `patient number`
    ,cell_type = fct_recode(`cell type`,B="B-cells",CD4="CD4 T-cells",CD8="CD8 T-cells",GR="Granulocytes",MO="Monocytes")
    ,tissue = "Blood"
    ,study = "Mill"
    ,red,grn
    )]

# map2(mill$red, mill$file %s+% "_Red.idat.gz", ~ download.file(.x,.y) ) %>% invisible
# map2(mill$grn, mill$file %s+% "_Grn.idat.gz", ~ download.file(.x,.y) ) %>% invisible

mill$red = NULL
mill$grn = NULL

# drop problematic samples
mill = mill[!sample_id %in% c("GSM2773348","GSM2773349","GSM2773350")]
mill[,file:=paste0("GSE103541/",sample_id)]

#--------------------------------------
# Saliva (Kelly Bakulski)

bakulski = fread("saliva/clean_saliva_samples.csv")
bakulski = bakulski[cell_type%in%c("Leukocytes","Epithelial cells")]

pheno3 = rbind(salas,mill,bakulski,use.names=TRUE,fill=TRUE); rm(salas,mill,bakulski)

beta3 =
    pheno3$file %>%
    ewastools::read_idats(.) %>%
    ewastools::detectionP.neg(.) %>%
    ewastools::mask(-2) %>%
    correct_dye_bias %>%
    dont_normalize

#--------------------------------------
### MERGER

common = paste0("beta",1:3) %>% map(get) %>% map(rownames) %>% reduce(intersect)
common = intersect(common,ewastools:::manifest_450K[!chr%in%c("X","Y")]$probe_id)

beta = cbind(
     beta1[ match(common,rownames(beta1)) ,]
    ,beta2[ match(common,rownames(beta2)) ,]
    ,beta3[ match(common,rownames(beta3)) ,]
    )

rm(beta1,beta2,beta3)

pheno = rbind(pheno1,reinius,pheno3,use.names=TRUE,fill=TRUE)
rm(common,i,J,pheno1,reinius,pheno3)


j = pheno[cell_type!="WBC",which=TRUE]
beta = beta[,j]
pheno = pheno[j]
pheno[,j:=1:.N]
pheno %<>% droplevels


table(pheno[,.(study,cell_type)])
#           cell_type
# study       B CD4 CD8 GR MO NK nRBC Leukocytes Epithelial cells
#   Bakulski 11   9   2 11 14 14    4         18               22
#   deGoede   7   7   6  7 12  6    7          0                0
#   Gervin   11  11  11 11  8 11    0          0                0
#   Lin      13  14  14 14 14 14    0          0                0
#   Mill     28  28  28 29 29  0    0          0                0
#   Reinius   4   4   4  5  5  3    0          0                0
#   Salas     6   7   6  6  6  6    0          0                0


train_model = function(train,output){
    
    train_beta = na.omit(beta[,train$j])
    cell_types = unique(train$cell_type)

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

train = copy(pheno[study=="Bakulski" & tissue=="Saliva"])
train_model(train,output="saliva.txt")
