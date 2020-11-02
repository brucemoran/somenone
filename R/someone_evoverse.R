#' Parse per-sample VCF and CN inputs, creating a list for evoverse input
#' evoverse is github.com/caravagna/mobster, cnaqc, viber, ctree etc.
#' Format:
#' chr, from, to, ref, alt, is_driver, driver_label, FILTER, DP, NV, VAF
#' ref/alt are the bases, DP total depth and NV No. reads at Variant, VAF = NV/DP
#' mut_rdata is master_consensus_shared from somatic_n-of-1 workflow
#' For somatic_n-of-1 workflow, CN is based on FACETS inputs
#' Format:
#' chr, from, to, Major, minor (not sic)
#' Purity also required as float
#'
#' @param mut_rdata is RData file with shared mutations in 2+ samples
#'  [1] "seqnames"             "start"                "end"
#'  [4] "width"                "strand"               "rowname"
#'  [7] "Consequence"          "IMPACT"               "SYMBOL"
#' [10] "HGVSc"                "HGVSp"                "HGVSp1"
#' [13] "CLIN_SIG"             "samples_n"            "sampleIDs"
#' [16] <sampleIDn.AD>       "sampleIDn.AD.1"     "sampleIDn.AF"
#' @param cn_pattern is a pattern that matches FACETS input filenames
#' @param pp_pattern is a pattern that matches polidy/purity table from FACETS
#' @param which_genome is GRCh37 or GRCh38
#' @return none, writes a tsv for input to Chimaera
#' @export

mut_cn_to_evoverse <- function(mut_rdata, cn_pattern, pp_pattern, which_genome, tag){

  ##genome
  if(which_genome %in% c("GRCh37", "hg19")){
    which_genome <- "hg19"
  } else {
    which_genome <- "GRCh38"
  }

  ##make tibble of Mutation format
  load(mut_rdata)
  mut_df <- as.data.frame(gr_master_consensus_all[[2]])
  mut_df$ref <- unlist(apply(mut_df, 1, function(f){
    strsplit(f["rowname"], "[_/]")[[1]][2]
  }))
  mut_df$alt <- unlist(apply(mut_df, 1, function(f){
    strsplit(f["rowname"], "[_/]")[[1]][3]
  }))

  ##add driver status for alll HIGH, MODERATE impacts
  mut_df$is_driver <- FALSE
  mut_df$is_driver[mut_df$IMPACT %in% c("HIGH", "MODERATE")] <- TRUE
  mut_df$driver_label <- FALSE
  mut_df$driver_label[mut_df$is_driver] <- paste0(mut_df$SYMBOL[mut_df$is_driver], ":", mut_df$HGVSp1[mut_df$is_driver])
  colnames(mut_df)[c(1,2,3)] <- c("chr", "from", "to")

  ##make df list of CN with Mutations
  sampleIDs <- sort(unique(unlist(strsplit(mut_df$sampleIDs, "\\,"))))

  ##parse out each sample and make correct format
  mut_df_list <- parse_sample_mut_list(sampleIDs, mut_df)

  ##parse out CN files and makke correct format
  cn_df_list <- parse_sample_cn_list(cn_pattern)

  purity_list <- parse_pur_list(pp_pattern)

  lapply(sampleIDs, function(s){
    readr::write_tsv(mut_df_list[[s]], file = paste0(s, ".mut.evoverse.tsv"))
    readr::write_tsv(cn_df_list[[s]], file = paste0(s, ".cna.evoverse.tsv"))
    data <- list(mut_df_list[[s]], cn_df_list[[s]], as.numeric(purity_list[s]), which_genome)
    names(data) <- c("mut", "cna", "purity", "reference")
    assigned_name <- paste0(s, "_evoverse_data")
    assign(assigned_name, value = data)
    save_file <- paste0(s, "_evoverse_data.RData")
    save(list = assigned_name, file = save_file)
  })
}

run_multivar_mobster_viber_ctree <- function(evo_pattern){
  evov_files <- dir(pattern = evo_pattern)
  viber_list <- lapply(evov_files, function(f){
      load(f)
      ff <- gsub(".RData", "", f)
      print(paste0("Working on: ", ff))
      ds <- get(ff)
      dataset_muts <- ds[[1]]
      run_mobster(dataset_muts)
    })
}
#' Parse mut_df and output as list per sample
#'
#' @param sampleIDS is vector of sample IDs to parse out depth data
#' @param mut_df is data frame of data with samples <sampleID[1].AD, .AD.1, .AF
#' @param vaf_filter is level of VAF to filter, default 5%
#' @return list of correct format per sample
#' @export

parse_sample_mut_list <- function(sampleIDs, mut_df, vaf_filter = 0.05){

  ##put muts in a list
  mut_list <- lapply(sampleIDs, function(s){

    ##create FILTER, DP, NV, VAF for sampleIDs[s]
    mut_sampleID <- dplyr::select(.data = mut_df, dplyr::starts_with(!!s))
    VAF <- dplyr::select(.data = mut_df, dplyr::starts_with(paste0(!!s, ".AF")))
    ADS <- dplyr::select(.data = mut_df, dplyr::starts_with(paste0(!!s, ".AD")))
    ADS[is.na(ADS)] <- 0
    ADS$DP <- as.numeric(ADS[,1]) + as.numeric(ADS[,2])
    mut_df$DP <- ADS$DP
    mut_df$NV <- ADS[,2]
    mut_df$VAF <- as.numeric(mut_df$NV) / as.numeric(mut_df$DP)
    mut_df$VAF[is.na(mut_df$VAF)] <- 0
    mut_df$FILTER <- mut_df$VAF > vaf_filter
    wanted <- c("chr", "from", "to", "ref", "alt", "is_driver", "driver_label", "FILTER", "DP", "NV", "VAF")
    muts <- dplyr::select(.data = mut_df, !!wanted)
    return(muts)
  })
  names(mut_list) <- sampleIDs
  return(mut_list)
}

#' Parse CN files and output as list
#'
#' @param cn_pattern is a pattern that matches CN input filenames
#' @return list of correct format per sample CNs
#'         chr, from, to, Major, minor (not sic)
#' @export

parse_sample_cn_list <- function(cn_pattern){

  ##files are in current working dir
  cn_files <- dir(pattern = cn_pattern)

  ##get sampleIDs as first element of split on dot
  cn_names <- unlist(lapply(cn_files, function(f){
    strsplit(f, "\\.")[[1]][1]
  }))

  ##put cns in a list
  cn_list <- lapply(cn_files, function(f){

    ##read CN data, make GRanges
    r1 <- readr::read_tsv(f)
    sample_id <- strsplit(f, "\\.")[[1]][1]

    print(paste0("Working on: ", sample_id))

    wanted <- c("chrom", "start", "end", "tcn.em", "lcn.em")
    cns <- r1[,wanted]
    colnames(cns) <- c("chr", "from", "to", "Major", "minor")

    return(cns)
  })
  names(cn_list) <- cn_names
  return(cn_list)
}

#' Parse purity files and output as list
#'
#' @param pp_pattern is a pattern that matches polidy/purity table from FACETS
#' @return list of correct format per sample
#' @export

parse_pur_list <- function(pp_patern){

  ##files are in current working dir
  cn_files <- dir(pattern = cn_pattern)

  ##get sampleIDs as first element of split on dot
  cn_names <- unlist(lapply(cn_files, function(f){
    strsplit(f, "\\.")[[1]][1]
  }))

  ##files are in current working dir
  pp_files <- dir(pattern = pp_pattern)

  ##get sampleIDs as first element of split on dot
  pp_names <- unlist(lapply(pp_files, function(f){
    strsplit(f, "\\.")[[1]][1]
  }))
  names(pp_files) <- pp_names

  pur_list <- lapply(pp_names, function(s){
    purity <- readr::read_tsv(pp_files[s])[[2]]
  })
  names(pur_list) <- pp_names
  return(pur_list)
}

#' Run CNAqc

run_cnaqc <- function(dataset){

  cnaqc_obj <- CNAqc::init(
                snvs = dataset[[1]],
                cna = dataset[[2]],
                purity = dataset[[3]],
                ref = dataset[[4]])

  # Set7_mapped_calls <- CNAqc::map_calls(
  #               CNA_calls = dataset,
  #               mutation_calls = Set7_mutations,
  #               purities = Set7_purity,
  #               samples = Set7_samples)
}

run_mobster <- function(dataset_muts){
  dataset_muts_f <- dplyr::filter(.data = dataset_muts, FILTER == TRUE)
  m_fit <- mobster::mobster_fit(dataset_muts_f)
  cm_best <- Clusters(m_fit$best)
  b_fit <- mobster::mobster_fit(cm_best)
  return(b_fit$best$data)
}
