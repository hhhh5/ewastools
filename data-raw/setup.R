library(data.table)
library(stringi)
library(magrittr)
library(parallel)
library(purrr)
library(forcats)

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

dir.create("idats")

# Download meta data
meta = "https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE68456&targ=gsm&form=text&view=brief"
meta = readLines(meta)

# Split into individual samples
meta = split(meta,cumsum(meta %like% "^\\^SAMPLE = GSM"))

# Extract GSM accessions
names(meta) = map(meta,1) %>% stri_match_first(regex="GSM\\d+")

# Parse meta data
imap(meta,function(s,acc){
    s = strsplit(s,split=" = ",fixed=TRUE)  
    data.table(gsm=acc,variable=map_chr(s,1),value=map_chr(s,2))
}) -> meta

meta = rbindlist(meta)

# Keep only information on sample characteristics and supplementary files
meta = meta[variable %chin% c("!Sample_characteristics_ch1","!Sample_supplementary_file")]
i = meta[variable == "!Sample_characteristics_ch1",which=TRUE]
ch = meta$value[i] %>% stri_split(fixed=": ")
meta$variable[i] = map_chr(ch,1)
meta$value   [i] = map_chr(ch,2)
rm(ch)

# Find the URLs pointing to the two .idat files
meta[variable == "!Sample_supplementary_file" & value %like% "_Red\\.idat",variable:="red"]
meta[variable == "!Sample_supplementary_file" & value %like% "_Grn\\.idat",variable:="grn"]

# Reshape data.table from long to wide format
meta = dcast(meta, gsm ~ variable)

# Download .idat files
# map2(meta$red, "idats/" %s+% meta$gsm %s+% "_Red.idat.gz", ~ download.file(.x,.y) )
# map2(meta$grn, "idats/" %s+% meta$gsm %s+% "_Grn.idat.gz", ~ download.file(.x,.y) )
meta$red = NULL; meta$grn = NULL

meta[,sex:=factor(Sex,levels=c("M","F"),labels=c("m","f"))]
meta$cell_type = fct_recode(meta$`cell type`,MO="Monocyte",NK="NK cell",CD8="CD8 T cell",CD4="CD4 T cell",GR="Granulocyte",B="B cell",nRBC="RBC nRBC")
meta[,file:= "idats/" %s+% gsm]
meta$tissue = "Cord Blood"
meta = meta[,list(cell_type,tissue,sex,file)]

cord_blood = meta
 
reinius = fread("reinius/sample_sheet_IDAT.csv") 
reinius[,`Chip Row Pos`:=stri_replace_all_fixed(`Chip Row Pos`,"O","0")] 
reinius[,`Chip CO position`:=stri_replace_all_fixed(`Chip CO position`,"O","0")] 
reinius[,file:=paste0("reinius/",`Chip#ID`,"_",`Chip Row Pos`,`Chip CO position`)] 
reinius$cell_type = fct_recode(reinius$Type,MO="CD14+ Monocytes",NK="CD56+ NK-cells",CD8="CD8+ T-cells",CD4="CD4+ T-cells",GR="Granulocytes",B="CD19+ B-cells")
reinius$tissue = "Blood (adult)"
reinius$sex = "m"
reinius = reinius[,list(cell_type,tissue,sex,file)]

purified = rbind(cord_blood,reinius); rm(cord_blood,reinius) 
purified = purified[cell_type%in%c('MO','NK','CD8','CD4','GR','B')]

# -------------------------------- 

meth = read_idats(purified$file)
beta = meth %>% correct_dye_bias %>% dont_normalize

cell_types = unique(purified$cell_type) 

markers = list() 

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
