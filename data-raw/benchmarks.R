library(data.table)
library(magrittr)
library(purrr)
library(stringi)
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

