lvls_chr37 = GenomeInfoDb::getChromInfoFromUCSC("hg19") $ chrom
lvls_chr38 = GenomeInfoDb::getChromInfoFromUCSC("hg38") $ chrom

# Add msytery chromosome "chr0"
lvls_chr37 = c(lvls_chr37, "chr0")
lvls_chr38 = c(lvls_chr38, "chr0")

lvls_base = c("A", "G", "C", "T")

MANIFESTS = list()
CONTROLS  = list()

# -------------------------------- EPIC V2 chip manifest
# Added by Costanza L. Vallerga
# Use version A2 of manifest
# CSV contains both 'normal' and control probes. Create two separate tables for them.

MANIFESTS$EPICv2 =
    "EPIC-8v2-0_A2.csv" |>
    readr::read_csv(
        skip = 7, # Column names are at line 8
        col_names = TRUE,
        col_types = "ccicicfffffffffficfcccicccccccccccccccccccccclcicfi",
        col_select = c(
            # IlmnID                             # cg25324105_BC11
            probe_id = Name,                     # cg25324105
            addressU = AddressA_ID,              # 1754126
            # AlleleA_ProbeSeq                   # ATTTATAAAC...
            addressM = AddressB_ID,              # 99753217
            # AlleleB_ProbeSeq                   # GTTTATAAA...
            next_base = Next_Base,               # A
            channel = Color_Channel,             # Red
            # col                                # R
            probe_type = Probe_Type,             # cg
            # strand_FR = Strand_FR,             # F
            # strand_TB = Strand_TB,             # B
            # strand_CO = Strand_CO,             # C
            # Infinium_Design                    # 1
            probe_design = Infinium_Design_Type, # I
            chr38 = CHR,                         # chr19
            mapinfo38 = MAPINFO,                 # 37692358
            # Species                            # Human
            # Genome_Build                       # GRCh38
            # Source_Seq                         # GTTTGTGGG...
            # Forward_Sequence                   # CGGTTCCGC...GGC[CG]ACGTGCT...
            # Top_Sequence                       # AGCAG...ACGT[CG]GCCGC...
            probe_rep = Rep_Num,                 # 1
            # UCSC_RefGene_Group                 # TSS200;TSS200
            # UCSC_RefGene_Name                  # ZNF781;ZNF781
            # UCSC_RefGene_Accession             # NR_173332.1;NR_173331.1
            # UCSC_CpG_Islands_Name              # chr19:37691892-37692426
            # Relation_to_UCSC_CpG_Island        # Island
            # GencodeV41_Group                   # TSS200;TSS200
            # GencodeV41_Name                    # ENSG00000120784.17;ENSG00000120784.17
            # GencodeV41_Accession               # ENST00000587199.5;ENST00000589676.5
            # Phantom5_Enhancers                 # NA
            # HMM_Island                         # NA
            # Regulatory_Feature_Name            # 19:38182603-38183622
            # Regulatory_Feature_Group           # Unclassified_Cell_type_specific
            # 450k_Enhancer                      # TRUE
            # DMR                                # DMR
            # DNase_Hypersensitivity_NAME        # EH38D5131145
            # Encode_CisReg_Site                 # EH38E3303729|DNase-only;EH38E3303729|Low-DNase
            # Encode_CisReg_Site_Evid            # 1;18;46;22
            # OpenChromatin_NAME                 # Het;Quies;TssA
            # OpenChromatin_Evidence_Count       # 261;34;540;408;120;186;50;42;12;3;7
            # Methyl450_Loci                     # cg25324105
            # Methyl27_Loci                      # NA
            # EPICv1_Loci                        # cg25324105
            # Manifest_probe_match               # TRUE
            # SNP_ID                             # NA
            # SNP_DISTANCE                       # NA
            # SNP_MinorAlleleFrequency           # NA
            chr37 = CHR_37,                      #
            mapinfo37 = MAPINFO_37               #
        ),
        n_max = 937055) # This number is documented in line 6 of the manifest

MANIFESTS$EPICv2 =
    MANIFESTS$EPICv2 |>
    dplyr::mutate(
        chr37 = stringr::str_c("chr", chr37), # Don't use paste0() as it doesn't handle NA
        chr38 = stringr::str_c("chr", chr38), # Don't use paste0() as it doesn't handle NA
        # ensure proper order of levels
        chr37 = forcats::fct_relevel(chr37, !!!lvls_chr37), 
        chr38 = forcats::fct_relevel(chr38, !!!lvls_chr38),
        probe_design = forcats::fct_relevel(probe_design, "I", "II"),
        channel = dplyr::if_else(probe_design == "II", "Both", channel),
        channel = forcats::fct_relevel(channel, "Red", "Green", "Both"),
        probe_type = forcats::fct_relevel(probe_type, "cg", "ch", "rs", "nv"),
        next_base = forcats::fct_relevel(next_base, !!!lvls_base))|>
    pointblank::col_vals_in_set(chr37, c(lvls_chr37, NA)) |>
    pointblank::col_vals_in_set(chr38, c(lvls_chr38, NA)) |>
    pointblank::row_count_match(937055)

CONTROLS$EPICv2 =
    "EPIC-8v2-0_A1.csv" |>
    readr::read_csv(
        skip = 7 + 937055 + 2, #
        col_names = c("address", "group", "channel", "name"),
        col_types = "iccc",
        col_select = 1:4) |>
    pointblank::row_count_match(635)

# -------------------------------- EPICv1
# Use version b5 of manifest
# CSV contains both 'normal' and control probes. Create two separate tables for them.

MANIFESTS$EPICv1 =
    "infinium-methylationepic-v-1-0-b5-manifest-file.csv" |>
    readr::read_csv(
        skip = 7, # Column names are at line 8
        col_names = TRUE,
        col_types = "ccicicfffcfficfccccccccccccccccccccccccllficid?lfiif",
        col_select = c(
            # IlmnID                              # cg07881041
            probe_id = Name,                      # cg07881041
            addressU = AddressA_ID,               # 0085713262
            # AlleleA_ProbeSeq                    # CTACAAATA...
            addressM = AddressB_ID,               # NA
            # AlleleB_ProbeSeq                    # NA
            probe_design = Infinium_Design_Type,  # II
            next_base = Next_Base,                # NA
            channel = Color_Channel,              # NA
            # Forward_Sequence                    # CTGCACGC...TAA[CG]CAT...AGGTG
            # Genome_Build                        # 37
            chr37 = CHR,                          # 19
            mapinfo37 = MAPINFO,                  # 5236016
            # SourceSeq                           # TGCAGGTG...
            # strand37 = Strand,                  # R
            # UCSC_RefGene_Name                   # PTPRS;PTPRS;PTPRS;PTPRS
            # UCSC_RefGene_Accession              # NM_130855;NM_002850;NM_130854;NM_130853
            # UCSC_RefGene_Group                  # Body;Body;Body;Body
            # UCSC_CpG_Islands_Name               # chr19:5237294-5237669
            # Relation_to_UCSC_CpG_Island         # N_Shore
            # Phantom4_Enhancers                  # NA
            # Phantom5_Enhancers                  # NA
            # DMR                                 # NA
            # 450k_Enhancer                       # NA
            # HMM_Island                          # NA
            # Regulatory_Feature_Name             # NA
            # Regulatory_Feature_Group            # NA
            # GencodeBasicV12_NAME                # NA
            # GencodeBasicV12_Accession           # NA
            # GencodeBasicV12_Group               # NA
            # GencodeCompV12_NAME                 # NA
            # GencodeCompV12_Accession            # NA
            # GencodeCompV12_Group                # NA
            # DNase_Hypersensitivity_NAME         # NA
            # DNase_Hypersensitivity_Evidence_Count NA
            # OpenChromatin_NAME                  # NA
            # OpenChromatin_Evidence_Count        # NA
            # TFBS_NAME                           # NA
            # TFBS_Evidence_Count                 # NA
            # Methyl27_Loci                       # NA
            # Methyl450_Loci                      # TRUE
            # Chromosome_36                       # 19
            # Coordinate_36                       # 5187016
            # SNP_ID                              # rs187313142
            # SNP_DISTANCE                        # 18
            # SNP_MinorAlleleFrequency            # 0.000200
            # Random_Loci                         # NA
            # MFG_Change_Flagged                  # FALSE
            chr38     = CHR_hg38,                 # chr19
            mapinfo38 = Start_hg38,               # 5236004
            #End_hg38                             # 5236006
            strand38  = Strand_hg38               # +
        ),
        n_max = 865918) # This number is documented in line 6 of the manifest

MANIFESTS$EPICv1 =
    MANIFESTS$EPICv1 |>
    dplyr::mutate(
        probe_rep = 1L, # There are no replicate probes for EPIC v1
        probe_type = substr(probe_id, 1L, 2L), # Column not in original manifest
        chr37 = stringr::str_c("chr", chr37), # Don't use paste0() as it doesn't handle NA
        # chr38 has already the "chr" prefix
        # ensure proper order of levels
        chr37 = forcats::fct_relevel(chr37, !!!lvls_chr37), 
        chr38 = forcats::fct_relevel(chr38, !!!lvls_chr38),
        probe_design = forcats::fct_relevel(probe_design, "I", "II"),
        channel = dplyr::if_else(probe_design == "II", "Both", channel),
        channel = forcats::fct_relevel(channel, "Red", "Green", "Both"),
        probe_type = forcats::fct_relevel(probe_type, "cg", "ch", "rs", "nv"),
        next_base = forcats::fct_relevel(next_base, !!!lvls_base))|>
    pointblank::col_vals_in_set(chr37, c(lvls_chr37, NA)) |>
    pointblank::col_vals_in_set(chr38, c(lvls_chr38, NA)) |>
    pointblank::row_count_match(865918)

CONTROLS$EPICv1 =
    "infinium-methylationepic-v-1-0-b5-manifest-file.csv" |>
    readr::read_csv(
        skip = 7 + 865918 + 2, # initial lines + loci + 2 headers
        col_names = c("address", "group", "channel", "name"),
        col_types = "iccc",
        col_select = 1:4) |>
    pointblank::row_count_match(count = 635)

# -------------------------------- 450K chip manifest
# Use version 1.2 of manifest
# CSV contains both 'normal' and control probes. Create two separate tables for them.

MANIFESTS$`450K` =
    "humanmethylation450_15017482_v1-2.csv" |>
    readr::read_csv(
        skip = 7, # Column names are at line 8
        col_names = TRUE,
        col_types = "ccicicfffcfficfifcccccccccccccccc",
        col_select = c(                                                      
            # IlmnID                             # cg00035864                                                                                                                  
            probe_id = Name,                     # cg00035864                                                                                                                  
            addressU = AddressA_ID,              # 31729416                                                                                                                    
            # AlleleA_ProbeSeq                   # AAAACACTAACAATC...                                                                          
            addressM = AddressB_ID,              # NA                                                                                                                            
            # AlleleB_ProbeSeq                   # NA                                                                                                                            
            probe_design = Infinium_Design_Type, # II                                                                                                                          
            next_base = Next_Base,               # NA                                                                                                                            
            channel = Color_Channel,             # NA                                                                                                                            
            # Forward_Sequence                   # AATCCAA...AAC[CG]AA...
            # Genome_Build                       # 37                                                                                                                          
            chr37 = CHR,                         # Y                                                                                                                           
            mapinfo37 = MAPINFO                  # 8553009                                                                                                                     
            # SourceSeq                          # AGACATTCG...                                                                          
            # Chromosome_36                      # Y                                                                                                                           
            # Coordinate_36                      # 8613009                                                                                                                    
            # Strand                             # F                                                                                                                           
            # Probe_SNPs                         # NA                                                                                                                            
            # Probe_SNPs_10                      # NA                                                                                                                            
            # Random_Loci                        # NA                                                                                                                            
            # Methyl27_Loci                      # NA                                                                                                                            
            # UCSC_RefGene_Name                  # TTTY18                                                                                                                      
            # UCSC_RefGene_Accession             # NR_001550                                                                                                                 
            # UCSC_RefGene_Group                 # TSS1500                                                                                                                     
            # UCSC_CpG_Islands_Name              # NA                                                                                                                            
            # Relation_to_UCSC_CpG_Island        # NA                                                                                                                            
            # Phantom                            # NA                                                                                                                            
            # DMR                                # NA                                                                                                                            
            # Enhancer                           # NA                                                                                                                            
            # HMM_Island                         # NA                                                                                                                            
            # Regulatory_Feature_Name            # NA                                                                                                                            
            # Regulatory_Feature_Group           # NA                                                                                                                            
            # DHS                                # NA  
        ),
        n_max = 485553 + 24) # The number mentioned in lines 6 of the manifest deviates
            # from the actual count by 24

MANIFESTS$`450K` =
    MANIFESTS$`450K` |>
    dplyr::mutate(
        probe_rep = 1L, # There are no replicate probes for EPIC v1
        probe_type = substr(probe_id, 1L, 2L), # Column not in original manifest
        chr37 = stringr::str_c("chr", chr37), # Don't use paste0() as it doesn't handle NA
        # ensure proper order of levels
        chr37 = forcats::fct_relevel(chr37, !!!lvls_chr37), 
        probe_design = forcats::fct_relevel(probe_design, "I", "II"),
        channel = dplyr::if_else(probe_design == "II", "Both", channel),
        channel = forcats::fct_relevel(channel, "Red", "Green", "Both"),
        probe_type = forcats::fct_relevel(probe_type, "cg", "ch", "rs", "nv"),
        next_base = forcats::fct_relevel(next_base, !!!lvls_base))|>
    pointblank::col_vals_in_set(chr37, c(lvls_chr37, NA)) |>
    pointblank::row_count_match(485553 + 24)

CONTROLS$`450K` =
    "humanmethylation450_15017482_v1-2.csv" |>
    readr::read_csv(
        skip = 7 + 485553 + 24 + 2, #
        col_names = c("address", "group", "channel", "name"),
        col_types = "iccc",
        col_select = 1:4) |>
    pointblank::row_count_match(850)

save(MANIFESTS, CONTROLS, file="../R/sysdata.rda",compress="xz")
