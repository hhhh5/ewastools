library(data.table)

# -------------------------------- EPIC V2 chip manifest

### NEW EPIC array added by Costanza L. Vallerga
### CSV contains both 'normal' and control probes. Create two separate tables for them (split at line 937055)

manifest_epic_v2 = fread("EPIC-8v2-0_A1.csv",skip="IlmnID",header=TRUE,nrows=937055,integer64="character",sep=",",sep2=";")

manifest_epic_v2 = manifest_epic_v2[,list(
     probe_id=IlmnID
    ,addressU=as.integer(AddressA_ID)
    ,addressM=as.integer(AddressB_ID)
    ,channel=Color_Channel
    ,next_base=Next_Base
    ,chr=CHR
    ,mapinfo=MAPINFO
    ,strand=factor(Strand_FR)
)]

manifest_epic_v2[                          ,probe_type:="cg"]
manifest_epic_v2[substr(probe_id,1,2)=="ch",probe_type:="ch"]
manifest_epic_v2[substr(probe_id,1,2)=="nv",probe_type:="nv"]
manifest_epic_v2[substr(probe_id,1,2)=="rs",probe_type:="rs"]

manifest_epic_v2[channel=="",channel:="Both"]

controls_epic_v2 = fread("EPIC-8v2-0_A1.csv",skip=937056,header=FALSE)
controls_epic_v2 = controls_epic_v2[,1:4]
names(controls_epic_v2) = c("address","group","channel","name")

# -------------------------------- EPIC chip manifest b5

### Added the latest version of the EPIC v1 (b5)
### CSV contains both 'normal' and control probes. Create two separate tables for them (split at line 865918) --> 635 controls

manifest_epic = fread("infinium-methylationepic-v-1-0-b5-manifest-file.csv",skip="IlmnID",header=TRUE,nrows=865918,integer64="character",sep=",",sep2=";")

manifest_epic = manifest_epic[,list(
     probe_id=IlmnID
    ,addressU=as.integer(AddressA_ID)
    ,addressM=as.integer(AddressB_ID)
    ,channel=Color_Channel
    ,next_base=Next_Base
    ,chr=CHR
    ,mapinfo=MAPINFO
    ,strand=factor(Strand)
)]

manifest_epic[                          ,probe_type:="cg"]
manifest_epic[substr(probe_id,1,2)=="ch",probe_type:="ch"]
manifest_epic[substr(probe_id,1,2)=="rs",probe_type:="rs"]

manifest_epic[channel=="",channel:="Both"]

controls_epic = fread("infinium-methylationepic-v-1-0-b5-manifest-file.csv",skip=865927,header=FALSE,fill=TRUE)
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
    ,strand=factor(Strand)
)]

manifest_450K[                          ,probe_type:="cg"]
manifest_450K[substr(probe_id,1,2)=="ch",probe_type:="ch"]
manifest_450K[substr(probe_id,1,2)=="rs",probe_type:="rs"]

manifest_450K[channel=="",channel:="Both"]

controls_450K = fread("HumanMethylation450_15017482_v1-2.csv",skip=485586,header=FALSE)
controls_450K = controls_450K[,1:4]
names(controls_450K) = c("address","group","channel","name")

save(manifest_epic,controls_epic,manifest_450K,controls_450K,file="../R/sysdata.rda",compress="xz")
