lvls_chr36 = GenomeInfoDb::getChromInfoFromUCSC("hg18") $ chrom
lvls_chr37 = GenomeInfoDb::getChromInfoFromUCSC("hg19") $ chrom
lvls_chr38 = GenomeInfoDb::getChromInfoFromUCSC("hg38") $ chrom
lvls_mm10  = GenomeInfoDb::getChromInfoFromUCSC("mm10") $ chrom

# Add mystery chromosome "chr0"
lvls_chr37 = c(lvls_chr37, "chr0")
lvls_chr38 = c(lvls_chr38, "chr0")
lvls_mm10  = c(lvls_mm10,  "chr0")

lvls_pdsign  = c("I", "II")
lvls_base = c("A", "G", "C", "T")
lvls_ptypes = c("cg", "ch", "rs", "nv")
lvls_channel = c("Red", "Grn", "Both")

MANIFESTS = list()
CONTROLS  = list()

fp_manifest = list(
    `27K`  = "humanmethylation27_270596_v1-2.csv",
    `450K` = "humanmethylation450_15017482_v1-2.csv",
    EPICv1 = "infinium-methylationepic-v-1-0-b5-manifest-file.csv",
    EPICv2 = "EPIC-8v2-0_A2.csv",
    MSA    = "MSA-48v1-0_20102838_B1.csv",
    MOUSE  = "MouseMethylation-12v1-0_A2.csv")

# -------------------------------- 27K chip manifest
# Use version 1.2 of manifest
# CSV contains both 'normal' and control probes. Create two separate tables for them.
# SNP probes (rs) are listed among the control probes with the two beads on separate rows

MANIFESTS$`27K` =
    fp_manifest$`27K` |>
    readr::read_csv(
        skip = 7, # Column names are at line 8
        col_names = TRUE,
        col_types = "ccficicffifffffccffifcccccccclcc",
        col_select = c(
             ilmn_id = IlmnID            # cg00000292
            ,probe_id = Name             # cg00000292
            # IlmnStrand                 # TOP
            ,addressU = AddressA_ID      # 990370
            # AlleleA_ProbeSeq
            ,addressM = AddressB_ID      # 6660678
            # AlleleB_ProbeSeq
            # GenomeBuild                # 36
            ,chr36 = Chr
            ,mapinfo36 = MapInfo
            # Ploidy                     # diploid
            # Species                    # Homo sapiens
            # Source                     # NCBI:RefSeq
            # SourceVersion              # 36.1     
            # SourceStrand               # TOP
            # SourceSeq
            # TopGenomicSeq
            ,next_base = Next_Base       # T
            ,channel = Color_Channel     # Red
            # TSS_Coordinate
            # Gene_Strand
            # Gene_ID
            # Symbol
            # Synonym
            # Accession
            # GID
            # Annotation
            # Product
            # Distance_to_TSS
            # CPG_ISLAND
            # CPG_ISLAND_LOCATIONS
            # MIR_CPG_ISLAND,MIR_NAMES
        ),
        n_max = 27578)

MANIFESTS$`27K` =
    MANIFESTS$`27K` |>
    dplyr::mutate(
        probe_rep = 1L, # There are no replicate probes for 27K
        probe_type = "cg", # 27K chip used only type "cg"
        chr36 = stringr::str_c("chr", chr36), # Don't use paste0() as it doesn't handle NA
        # ensure proper order of levels
        chr36 = forcats::fct_relevel(chr36, !!!lvls_chr36),
        probe_design = forcats::fct("I", levels = c("I", "II")),
        channel = dplyr::if_else(probe_design == "II", "Both", channel),
        channel = forcats::fct_relevel(channel, !!!lvls_channel),
        probe_type = forcats::fct_relevel(probe_type, !!!lvls_ptypes),
        next_base = forcats::fct_relevel(next_base, !!!lvls_base))

SNPS27K = fp_manifest$`27K` |>
    readr::read_csv(
        skip = 27699, # Column names are at line 8
        col_names = c("address", "group", "name", "probe_id", NA),
        col_types = "iccc_",
        col_select = 1:4) |>
    pointblank::col_vals_equal(group, "Genotyping") |>
    dplyr::select(address, probe_id) |>
    tidyr::separate_wider_delim(probe_id, delim = "_", names = c("probe_id", "bead")) |>
    dplyr::mutate(bead = dplyr::recode_values(bead, "A" ~ "M", "B" ~ "U", unmatched = "error")) |>
    tidyr::pivot_wider(id_cols = probe_id,
        names_from = bead, names_prefix = "address", names_sep = "", values_from = address) |>
    dplyr::mutate(probe_type = forcats::fct("rs", levels = !!lvls_ptypes))

df_missing_channel = tibble::tibble(
    probe_id = c("rs1019916", "rs10457834", "rs1416770", "rs1941955", "rs2125573", "rs2235751",
      "rs2521373", "rs264581", "rs2804694", "rs2959823", "rs5931272", "rs6546473", 
      "rs739259",  "rs798149",  "rs845016",  "rs866884"),
    channel = c("Red", "Red", "Red", "Red", "Red", "Grn", "Grn", "Red", "Red", "Red", "Red",
        "Grn", "Red", "Grn", "Grn", "Red"))

SNPS27K = dplyr::inner_join(SNPS27K, df_missing_channel)

MANIFESTS$`27K` = dplyr::bind_rows(MANIFESTS$`27K`, SNPS27K)

MANIFESTS$`27K` = MANIFESTS$`27K` |>
    pointblank::col_vals_in_set(chr36, c(lvls_chr36, NA)) |>
    pointblank::row_count_match(27578 + 16) |>
    invisible()

CONTROLS$`27K` =
    fp_manifest$`27K` |>
    readr::read_csv(
        skip = 7 + 27578 + 2, # initial lines + loci + deviation + 2 headers
        col_names = c("address", "group", "channel", "name", NA),
        col_types = "iccc_",
        col_select = 1:4) |>
    dplyr::filter(!group %in% c("pACYC174", "phiX174", "pUC19", "Genotyping")) |>
    dplyr::mutate(
        group = stringr::str_to_upper(group),
        group = dplyr::replace_values(group,
            "SPECIFICITY" ~ "SPECIFICITY I",
            "BISULFITE CONVERSION" ~ "BISULFITE CONVERSION I",
            "NORMALIZATION-RED" ~ "NORM_A",
            "NORMALIZATION-GREEN" ~ "NORM_G")) |>
    pointblank::row_count_match(94)

setDT(CONTROLS$`27K`)

CONTROLS$`27K`[group == "STAINING", name := stringr::str_replace(name, "Bgnd", "Bkg")]
CONTROLS$`27K`[group == "TARGET REMOVAL", name := "Target Removal 1"]
CONTROLS$`27K`[group == "BISULFITE CONVERSION I", name := stringr::str_replace(name, "BS conversion ", "BS Conversion I-")]
CONTROLS$`27K`[group == "SPECIFICITY I", name := stringr::str_replace(name, "GT mismatch (\\d)(.+$)", "GT Mismatch \\1 \\2")]
CONTROLS$`27K`[group == "NORM_A", `:=`(channel = "Red",  name = stringr::str_replace(name, " Red ", "_A"))]
CONTROLS$`27K`[group == "NORM_G", `:=`(channel = "Blue", name = stringr::str_replace(name, " Green ", "_G"))]

rm(SNPS27K)

# -------------------------------- 450K chip manifest
# Use version 1.2 of manifest
# CSV contains both 'normal' and control probes. Create two separate tables for them.

MANIFESTS$`450K` =
    fp_manifest$`450K` |>
    readr::read_csv(
        skip = 7, # Column names are at line 8
        col_names = TRUE,
        col_types = "ccicicfffcfficfifcccccccccccccccc",
        col_select = c(
            ilmn_id = IlmnID,                    # cg00035864
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
        probe_rep    = 1L, # There are no replicate probes for EPIC v1
        probe_type   = substr(probe_id, 1L, 2L), # Column not in original manifest
        chr37        = stringr::str_c("chr", chr37), # Don't use paste0() as it doesn't handle NA
        chr37        = forcats::fct_relevel(chr37, !!!lvls_chr37),
        probe_design = forcats::fct_relevel(probe_design, !!!lvls_pdsign),
        channel      = dplyr  ::if_else(probe_design == "II", "Both", channel),
        channel      = forcats::fct_relevel(channel,    !!!lvls_channel),
        probe_type   = forcats::fct_relevel(probe_type, !!!lvls_ptypes),
        next_base    = forcats::fct_relevel(next_base,  !!!lvls_base))|>
    pointblank::col_vals_in_set(chr37, c(lvls_chr37, NA)) |>
    pointblank::row_count_match(485553 + 24)

CONTROLS$`450K` =
    fp_manifest$`450K` |>
    readr::read_csv(
        skip = 7 + 485553 + 24 + 2, # initial lines + loci + deviation + 2 headers
        col_names = c("address", "group", "channel", "name"),
        col_types = "iccc",
        col_select = 1:4) |>
    pointblank::row_count_match(850)

# -------------------------------- EPICv1
# Use version b5 of manifest
# CSV contains both 'normal' and control probes. Create two separate tables for them.

MANIFESTS$EPICv1 =
    fp_manifest$EPICv1 |>
    readr::read_csv(
        skip = 7, # Column names are at line 8
        col_names = TRUE,
        col_types = "ccicicfffcfficfccccccccccccccccccccccccllficid?lfiif",
        col_select = c(
            ilmn_id = IlmnID,                     # cg07881041
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
            # End_hg38                            # 5236006
            # Strand_hg38                         # +
        ),
        n_max = 865918) # This number is documented in line 6 of the manifest

MANIFESTS$EPICv1 =
    MANIFESTS$EPICv1 |>
    dplyr::mutate(
        probe_rep    = 1L, # There are no replicate probes for EPIC v1
        probe_type   = substr(probe_id, 1L, 2L), # Column not in original manifest
        chr37        = stringr::str_c("chr", chr37), # Don't use paste0() as it doesn't handle NA
        chr37        = forcats::fct_relevel(chr37, !!!lvls_chr37),
        chr38        = forcats::fct_relevel(chr38, !!!lvls_chr38), # chr38 already prefixed
        probe_design = forcats::fct_relevel(probe_design, !!!lvls_pdsign),
        channel      = dplyr  ::if_else(probe_design == "II", "Both", channel),
        channel      = forcats::fct_relevel(channel,    !!!lvls_channel),
        probe_type   = forcats::fct_relevel(probe_type, !!!lvls_ptypes),
        next_base    = forcats::fct_relevel(next_base,  !!!lvls_base))|>
    pointblank::col_vals_in_set(chr37, c(lvls_chr37, NA)) |>
    pointblank::col_vals_in_set(chr38, c(lvls_chr38, NA)) |>
    pointblank::row_count_match(865918)

CONTROLS$EPICv1 =
    fp_manifest$EPICv1 |>
    readr::read_csv(
        skip = 7 + 865918 + 2, # initial lines + loci + 2 headers
        col_names = c("address", "group", "channel", "name"),
        col_types = "iccc",
        col_select = 1:4) |>
    pointblank::row_count_match(count = 635)

# -------------------------------- EPIC V2 chip manifest
# Added by Costanza L. Vallerga
# Use version A2 of manifest
# CSV contains both 'normal' and control probes. Create two separate tables for them.

MANIFESTS$EPICv2 =
    fp_manifest$EPICv2 |>
    readr::read_csv(
        skip = 7, # Column names are at line 8
        col_names = TRUE,
        col_types = "ccicicfffffffffficfcccicccccccccccccccccccccclcicfi",
        col_select = c(
            ilmn_id = IlmnID,                    # cg25324105_BC11
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
        chr37        = forcats::fct_relevel(chr37, !!!lvls_chr37),
        chr38        = forcats::fct_relevel(chr38, !!!lvls_chr38),
        probe_design = forcats::fct_relevel(probe_design, "I", "II"),
        channel      = dplyr  ::if_else(probe_design == "II", "Both", channel),
        channel      = forcats::fct_relevel(channel, !!!lvls_channel),
        probe_type   = forcats::fct_relevel(probe_type, !!!lvls_ptypes),
        next_base    = forcats::fct_relevel(next_base, !!!lvls_base))|>
    pointblank::col_vals_not_null(columns =
        c(ilmn_id, probe_id, addressU, channel, probe_type, probe_design, probe_rep)) |>
    pointblank::rows_distinct(columns = ilmn_id) |>
    pointblank::rows_distinct(columns = addressU) |>
    pointblank::rows_distinct(columns = addressM,
        preconditions = \(df) dplyr::filter(df, !is.na(addressM))) |>
    pointblank::col_vals_in_set(probe_design, c("I", "II")) |>
    pointblank::col_vals_in_set(chr37, c(lvls_chr37, NA)) |>
    pointblank::col_vals_in_set(chr38, c(lvls_chr38, NA)) |>
    pointblank::row_count_match(937055) |>
    invisible()

CONTROLS$EPICv2 =
    fp_manifest$EPICv2 |>
    readr::read_csv(
        skip = 7 + 937055 + 2, #
        col_names = c("address", "group", "channel", "name"),
        col_types = "iccc",
        col_select = 1:4) |>
    pointblank::row_count_match(635)

# -------------------------------- Methylation Screening Array (MSA)
# Use version B1 of manifest
# CSV contains both 'normal' and control probes. Create two separate tables for them.

MANIFESTS$MSA =
    fp_manifest$MSA |>
    readr::read_csv(
        skip = 7, # Column names are at line 8
        col_names = TRUE,
        col_types = "ccicicfffffffffficfcccicccccccccccccccccccccclcicfi",
        col_select = c(
            ilmn_id = IlmnID,                    # cg25324105_BC11
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
            probe_rep = Rep_Num,                 # 1
            chr38 = CHR,                         # chr19
            mapinfo38 = MAPINFO,                 # 37692358
            # Species                            # Human
            # Genome_Build                       # GRCh38
            # Source_Seq                         # GTTTGTGGG...
            # Forward_Sequence                   # CGGTTCCGC...GGC[CG]ACGTGCT...
            # Top_Sequence                       # AGCAG...ACGT[CG]GCCGC...
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
            # DNase_Hypersensitivity_NAME        # EH38D5131145
            # Encode_CisReg_Site                 # EH38E3303729|DNase-only;EH38E3303729|Low-DNase
            # Encode_CisReg_Site_Evid            # 1;18;46;22
            # OpenChromatin_NAME                 # Het;Quies;TssA
            # OpenChromatin_Evidence_Count       # 261;34;540;408;120;186;50;42;12;3;7
            # EPICv2_Locus_Match
            # EPICv1_Locus_Match
            # Methyl450_Locus_Match
            # EPICv2_ProbeSeq_Match
            # EPICv1_ProbeSeq_Match
            # Methyl450_ProbeSeq_Match
            # SNP_ID                             # NA
            # SNP_DISTANCE                       # NA
            # SNP_MinorAlleleFrequency           # NA
            chr37 = CHR_GRCh37,                  #
            mapinfo37 = MAPINFO_GRCh37           #
        ),
        n_max = 281798) # This number is documented in line 6 of the manifest

MANIFESTS$MSA =
    MANIFESTS$MSA |>
    dplyr::mutate(
        chr37        = stringr::str_c("chr", chr37), # Don't use paste0() as it doesn't handle NA
        chr38        = stringr::str_c("chr", chr38), # Don't use paste0() as it doesn't handle NA
        chr37        = forcats::fct_relevel(chr37, !!!lvls_chr37),
        chr38        = forcats::fct_relevel(chr38, !!!lvls_chr38),
        probe_design = forcats::fct_relevel(probe_design, !!!lvls_pdsign),
        channel      = dplyr  ::if_else(probe_design == "II", "Both", channel),
        channel      = forcats::fct_relevel(channel, !!!lvls_channel),
        probe_type   = forcats::fct_relevel(probe_type, !!!lvls_ptypes),
        next_base    = forcats::fct_relevel(next_base, !!!lvls_base))|>
    pointblank::col_vals_not_null(columns =
        c(ilmn_id, probe_id, addressU, channel, probe_type, probe_design, probe_rep)) |>
    pointblank::rows_distinct(columns = ilmn_id) |>
    pointblank::rows_distinct(columns = addressU) |>
    pointblank::rows_distinct(columns = addressM,
        preconditions = \(df) dplyr::filter(df, !is.na(addressM))) |>
    pointblank::col_vals_in_set(probe_design, c("I", "II")) |>
    pointblank::col_vals_in_set(chr37, c(lvls_chr37, NA)) |>
    pointblank::col_vals_in_set(chr38, c(lvls_chr38, NA)) |>
    pointblank::row_count_match(281798) |>
    invisible()

CONTROLS$MSA =
    fp_manifest$MSA |>
    readr::read_csv(
        skip = 7 + 281798 + 2, #
        col_names = c("address", "group", "channel", "name"),
        col_types = "iccc",
        col_select = 1:4) |>
    pointblank::row_count_match(2511)

# -------------------------------- Mouse chip manifest
# Manifest from May 25, 2021
# CSV contains both 'normal' and control probes. Create two separate tables for them.

MANIFESTS$MOUSE =
    fp_manifest$MOUSE |>
    readr::read_csv(
        skip = 7, # Column names are at line 8
        col_names = TRUE,
        col_types = "ccicicfffffffffficfcccicccccccccccccccccccccclcicfi",
        col_select = c(
            ilmn_id = IlmnID,                    # cg25324105_BC11
            probe_id = Name,                     # cg25324105
            addressU = AddressA_ID,              # 1754126
            # AlleleA_ProbeSeq                   # ATTTATAAAC...
            addressM = AddressB_ID,              # 99753217
            # AlleleB_ProbeSeq                   # GTTTATAAA...
            next_base = Next_Base,               # A
            channel = Color_Channel,             # Red
            # Col                                # F/R
            probe_type = Probe_Type,             # cg (ctl, cg, ch, my, rp, rs)
            # strand_FR = Strand,                # F
            # strand_TB = Strand_TB,             # B
            # strand_CO = Strand_CO,             # C
            probe_design = Infinium_Design_Type, # 1/2
            probe_rep = Rep_Num,                 # 1
            chr = CHR,                           # 
            mapinfo = MAPINFO,                   # 75515719
            # Species                            # Mus musculus
            # Genome_Build                       # mm10
            # SourceSeq                          # GTTTGTGGG...
            # Forward_Sequence                   # CGGTTCCGC...GGC[CG]ACGTGCT...
            # Top_Sequence                       # AGCAG...ACGT[CG]GCCGC...
            # Genome_Build_NCBI                  #
            # N_Shelf                            #
            # N_Shore                            #
            # CpG_Island                         #
            # CpG_Island_chrom                   #
            # CpG_Island_chromStart              #
            # CpG_Island_chromEnd                #
            # CpG_Island_length                  #
            # CpG_Island_cpgNum                  #
            # CpG_Island_gcNum                   #
            # CpG_Island_perCpg                  #
            # CpG_Island_perGc                   #
            # CpG_Island_obsExp                  #
            # S_Shore                            #
            # S_Shelf                            #
            # MFG_Change_Flagged                 #
        ),
        n_max = 287050) # This number is documented in line 6 of the manifest

MANIFESTS$MOUSE =
    MANIFESTS$MOUSE |>
    dplyr::mutate(
        chr = stringr::str_c("chr", chr), # Don't use paste0() as it doesn't handle NA
        chr = forcats::fct_recode(chr, "chrM" = "chrMT"),
        # ensure proper order of levels
        chr = forcats::fct_relevel(chr, !!!lvls_mm10),
        probe_design = forcats::fct_recode(probe_design, "I" = "1", "II" = "2"),
        channel = dplyr::if_else(probe_design == "II", "Both", channel),
        channel = forcats::fct_relevel(channel, !!!lvls_channel),
        probe_type = forcats::fct_relevel(probe_type, "cg", "ch", "rs", "rp", "mu"), # nv
        next_base = forcats::fct_relevel(next_base, !!!lvls_base, "Y"))|>
    pointblank::col_vals_in_set(chr, c(lvls_mm10, NA)) |>
    pointblank::row_count_match(287050)

CONTROLS$MOUSE =
    fp_manifest$MOUSE |>
    readr::read_csv(
        skip = 7 + 287050 + 2, #
        col_names = c("address", "group", "channel", "name"),
        col_types = "iccc",
        col_select = 1:4) |>
    pointblank::row_count_match(642)

setDT(CONTROLS$MOUSE)


# Move species tag to dedicated column
CONTROLS$MOUSE[, species := stringr::str_extract(name, pattern = "_(HSA|MUS)$", group = 1)]
CONTROLS$MOUSE[, name := stringr::str_remove(name, pattern = "_(HSA|MUS)$")]

# Harmonize names of control probes
CONTROLS$MOUSE[group == "RESTORATION", name := "Restore"]
CONTROLS$MOUSE[group == "HYBRIDIZATION", name := stringr::str_replace(name, "HYB-(.+)$", "Hyb (\\1)")]
CONTROLS$MOUSE[group == "EXTENSION", name := stringr::str_replace(name, "EXT-(.)", "Extension (\\1)")]
CONTROLS$MOUSE[group == "STAINING", name := stringr::str_replace(name, "STN-(DNP|Biotin)-(.+)", "\\1 (\\2)")]
CONTROLS$MOUSE[group == "TARGET REMOVAL", name := stringr::str_replace(name, "^TRM-(\\d)$", "Target Removal \\1")]
CONTROLS$MOUSE[group == "BISULFITE CONVERSION I", name := stringr::str_replace(name, "^BS1-(.+)$", "BS Conversion I-\\1")]
CONTROLS$MOUSE[group == "BISULFITE CONVERSION II", name := stringr::str_replace(name, "^BS2-(\\d+)$", "BS Conversion II-\\1")]
CONTROLS$MOUSE[group == "SPECIFICITY I", name := stringr::str_replace(name, "^SP1-(PM|MM)(\\d)$", "GT Mismatch \\2 (\\1)")]
CONTROLS$MOUSE[group == "SPECIFICITY II", name := stringr::str_replace(name, "^SP2-(\\d)", "Specificity \\1")]
CONTROLS$MOUSE[group == "NON-POLYMORPHIC", name := stringr::str_replace(name, "^NPM-([ACGT])", "NP (\\1)")]
CONTROLS$MOUSE[group == "NON-POLYMORPHIC", name := stringr::str_replace(name, "^NPM-(\\d)", "NP (G) \\1")]
CONTROLS$MOUSE[group == "NEGATIVE", name := stringr::str_replace(name, "^NEG-(\\d+)$", "Negative \\1")]
CONTROLS$MOUSE[group %like% "^NORM_[ACGT]$", name := stringr::str_replace(name, "^NRM-([ACGT])-(\\d+)$", "Norm_\\1\\2")]

# --------------------------------

save(MANIFESTS, CONTROLS, file="../R/sysdata.rda",compress="xz")
