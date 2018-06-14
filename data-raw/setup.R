library(data.table)
library(stringi)
library(magrittr)
library(parallel)
library(purrr)
library(forcats)
library(ewastools)

# -------------------------------- EPIC chip manifest

### CSV contains both 'normal' and control probes. Create two separate tables for them (split at line 865927)

manifest_epic = fread("MethylationEPIC_v-1-0_B3.csv"
    ,skip="IlmnID",header=TRUE,nrows=865918,integer64="character",sep=",",sep2=";")

manifest_epic = manifest_epic[,list(
     probe_id=IlmnID
    ,addressU=as.integer(AddressA_ID)
    ,addressM=as.integer(AddressB_ID)
    ,channel=Color_Channel
    ,next_base=Next_Base
    ,chr=CHR
    ,mapinfo=MAPINFO
)]

manifest_epic[                          ,probe_type:="cg"]
manifest_epic[substr(probe_id,1,2)=="ch",probe_type:="ch"]
manifest_epic[substr(probe_id,1,2)=="rs",probe_type:="rs"]

manifest_epic[channel=="",channel:="Both"]

controls_epic = fread("MethylationEPIC_v-1-0_B3.csv",skip=865927,header=FALSE)
controls_epic = controls_epic[,1:4]
names(controls_epic) = c("address","group","channel","name")

# -------------------------------- 450K chip manifest

### CSV contains both 'normal' and control probes. Create two separate tables for them (split at line 865927)

manifest_450K = fread("HumanMethylation450_15017482_v1-2.csv"
    ,skip="IlmnID",header=TRUE,nrows=485577,integer64="character",sep=",",sep2=";")

manifest_450K = manifest_450K[,list(
     probe_id=IlmnID
    ,addressU=AddressA_ID
    ,addressM=AddressB_ID
    ,channel=Color_Channel
    ,next_base=Next_Base
    ,chr=CHR
    ,mapinfo=MAPINFO
)]

manifest_450K[                          ,probe_type:="cg"]
manifest_450K[substr(probe_id,1,2)=="ch",probe_type:="ch"]
manifest_450K[substr(probe_id,1,2)=="rs",probe_type:="rs"]

manifest_450K[channel=="",channel:="Both"]

controls_450K = fread("HumanMethylation450_15017482_v1-2.csv",skip=485586,header=FALSE)
controls_450K = controls_450K[,1:4]
names(controls_450K) = c("address","group","channel","name")


save(manifest_epic,controls_epic,manifest_450K,controls_450K,file="../R/sysdata.rda",compress="xz")

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

# -------------------------------- 

purified = rbind(goede,gervin,reinius,use.names=TRUE,fill=TRUE)
purified %<>% droplevels
purified[,j:=1:.N]

meth = read_idats(purified$file) %>% correct_dye_bias %>% detectionP %>% mask(0.05)
beta = meth %>% dont_normalize

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



for(ct in cell_types){ 
    cat(ct,'\n') 
 
    j = purified$cell_type == ct 
 
    tmp = apply(beta,1,function(x){ 
        if(!any(is.na(x[j]))) 
        {  
            tmp = t.test(x[j],x[!j],var.equal=T) 
            return(c(tmp$p.value,tmp$estimate[1]-tmp$estimate[2])) 
        }else{ 
            return(c(NA,NA)) 
        } 
    }) 
 
    i = which(p.adjust(tmp[1,])<0.05) 
    o = order(tmp[2,i]) 
    markers[[ct]] = c(i[head(o,50)],i[tail(o,50)]) 
} 

markers %<>% unlist %>% unique 

coefs = sapply(cell_types,function(ct){ 
    j = purified$cell_type == ct 
    rowMeans(beta[markers,j],na.rm=TRUE) 
}) 
rownames(coefs) = rownames(beta)[markers] 
colnames(coefs) = cell_types 

write.table(coefs,file="../data/blood_coefs.txt",row.names=TRUE) 

# --------------------------------
