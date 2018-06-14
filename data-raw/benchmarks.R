library(data.table)
library(magrittr)
library(purrr)
library(stringi)
library(car)
library(ewastools)


samples = "https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE100197&targ=gsm&form=text&view=brief"
samples = readLines(samples)
samples = split(samples,cumsum(samples %like% "^\\^SAMPLE = GSM") )

# Extract GSM accessions
names(samples) = map(samples,1) %>% stri_match_first(regex="GSM\\d+")

# Parse meta data
imap(samples,function(s,acc){
	s = strsplit(s,split=" = ",fixed=TRUE)	
	data.table(gsm=acc,variable=map_chr(s,1),value=map_chr(s,2))
}) -> samples

samples = rbindlist(samples)

# Keep only information on sample characteristics and supplementary files
samples = samples[variable %chin% c("!Sample_characteristics_ch1","!Sample_supplementary_file","!Sample_title")]
i = samples[variable == "!Sample_characteristics_ch1",which=TRUE]
ch = samples$value[i] %>% stri_split(fixed=": ")
samples$variable[i] = map_chr(ch,1)
samples$value   [i] = map_chr(ch,2)
rm(ch)

# Find the URLs pointing to the two .idat files
samples[variable=="!Sample_supplementary_file" & value %like% "_Red\\.idat",variable:="red"]
samples[variable=="!Sample_supplementary_file" & value %like% "_Grn\\.idat",variable:="grn"]

samples[variable=="!Sample_title",variable:="name"]

samples[,variable:=chartr(old=" ",new="_",variable)]

# Reshape data.table from long to wide format
samples = dcast(samples, gsm ~ variable)

setnames(samples,"450k_plate","plate")
setnames(samples,"450k_sentrix_id","chip")
setnames(samples,"450k_sentrix_position","well")

samples[subject_id %like% "^P[ML]\\d+(r\\d?)$",subject_id:=stri_extract(subject_id,regex="P[ML]\\d+")]

samples = samples[gsm %in% c("GSM2674469","GSM2674488","GSM2674479","GSM2674490","GSM2674440","GSM2674492","GSM2674494","GSM2674495")]

map2(samples$red, samples$gsm %s+% "_Red.idat.gz", ~ download.file(.x,.y) ) %>% invisible
map2(samples$grn, samples$gsm %s+% "_Grn.idat.gz", ~ download.file(.x,.y) ) %>% invisible
samples$red = NULL; samples$grn = NULL

meth =
	read_idats(samples$gsm,quiet=TRUE) %>%
	correct_dye_bias %>%
	detectionP





# test set (GSE77797)  https://doi.org/10.1186/s12859-016-0943-7

dir.create("GSE77797")

lc_test = "https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE77797&targ=gsm&form=text&view=brief"
lc_test = readLines(lc_test)
lc_test = split(lc_test,cumsum(lc_test %like% "^\\^SAMPLE = GSM") )

# Extract GSM accessions
names(lc_test) = map(lc_test,1) %>% stri_match_first(regex="GSM\\d+")

# Parse meta data
imap(lc_test,function(s,acc){
	s = strsplit(s,split=" = ",fixed=TRUE)	
	data.table(gsm=acc,variable=map_chr(s,1),value=map_chr(s,2))
}) -> lc_test

lc_test = rbindlist(lc_test)

# Keep only information on sample characteristics and supplementary files
lc_test = lc_test[variable %chin% c("!Sample_characteristics_ch1","!Sample_supplementary_file","!Sample_title")]

i = lc_test[variable == "!Sample_characteristics_ch1",which=TRUE]
ch = lc_test$value[i] %>% stri_split(fixed=": ")
lc_test$variable[i] = map_chr(ch,1)
lc_test$value   [i] = map_chr(ch,2)
rm(ch)

lc_test[variable == "!Sample_title",variable:="title"    ]
lc_test[,variable:= recode(variable,"'b cell (%)'='B';'cd4+ t cell (%)'='CD4';'cd8+ t cell (%)'='CD8';'granulocyte (%)'='GR';'monocyte (%)'='MO';'natural killer cell (%)'='NK'")]

# Find the URLs pointing to the two .idat files
lc_test[variable == "!Sample_supplementary_file" & value %like% "_Red\\.idat",variable:="red"]
lc_test[variable == "!Sample_supplementary_file" & value %like% "_Grn\\.idat",variable:="grn"]

lc_test = dcast(lc_test,gsm ~ variable,fill=NA)

lc_test[,file:=paste0("GSE77797","/",gsm)]

# map2(lc_test$red, lc_test$file %s+% "_Red.idat.gz", ~ download.file(.x,.y) ) %>% invisible
# map2(lc_test$grn, lc_test$file %s+% "_Grn.idat.gz", ~ download.file(.x,.y) ) %>% invisible
lc_test$grn = NULL; lc_test$red = NULL

lc_test$NK  %<>% as.numeric
lc_test$GR  %<>% as.numeric
lc_test$MO  %<>% as.numeric
lc_test$CD8 %<>% as.numeric
lc_test$CD4 %<>% as.numeric
lc_test$B   %<>% as.numeric

meth = lc_test %$% file %>% read_idats %>% correct_dye_bias %>% detectionP %>% mask(0.05)
beta = dont_normalize(meth)

LC1 = estimateLC(beta,ref="reinius")
LC2 = estimateLC(beta,ref="goede")
LC3 = estimateLC(beta,ref="gervin")
LC4 = estimateLC(beta,ref="gervin+reinius")

cell_types = c("B","CD4","CD8","GR","MO","NK")

sapply(cell_types,function(ct) cor(LC1[,ct],lc_test[[ct]]))
sapply(cell_types,function(ct) cor(LC2[,ct],lc_test[[ct]]))
sapply(cell_types,function(ct) cor(LC3[,ct],lc_test[[ct]]))
sapply(cell_types,function(ct) cor(LC4[,ct],lc_test[[ct]]))
