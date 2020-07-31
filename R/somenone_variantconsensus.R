#' somtic_variant_consensus functions
#'
#' Primary function to call others and produce plots, table
#' @param germline_id ID for germline_id sample
#' @param vep_vcf_pattern pattern to match VEP annotated VCFs
#' @param raw_vcf_pattern pattern to match raw unannotated, unfiltered VCFs unfiltered
#' @param included_order oredering of samples for plotting
#' @param name_callers named variant callers to use primary to random selection
#' @param tag is a string to tag output files
#' @return vector of single-letter HGVS protein IDs
#' @export

variant_consensus <- function(germline_id, vep_vcf_pattern, raw_vcf_pattern = "raw.vcf", included_order = NULL, name_callers = NULL, tag = NULL) {

  options(stringAsFactors = FALSE)

  ##parse VCFs
  ##all should exist in current dir, all output to same dir
  vcf_vec <- dir(pattern = vep_vcf_pattern)
  vcf_list <- vcf_vec[grep(paste0(vep_vcf_pattern, "$"), vcf_vec)]
  raw_vec <- dir(pattern = raw_vcf_pattern)
  raw_list <- raw_vec[grep(paste0(raw_vcf_pattern, "$"), raw_vec)]
  input_list <- list(vcf_list, raw_list)
  vcf_ext <- gsub("-","_",
                 base::paste(stringr::str_split(
                                gsub("\\.vcf","",
                                vcf_list[[1]]),
                                "\\.")[[1]][-c(1,2)],
                      collapse="."))

  ##create order from vcf_list
  if(!is.null(included_order)){
      included_order <- stringr::str_split(included_order, ",")[[1]]
  } else {
      included_order <- unique(unlist(lapply(vcf_list, function(f){
                          stringr::str_split(f, "\\.")[[1]][1]
                        })))
  }

  ##operate over vep, raw VCFs
  for(x in 1:2) {
    in_vec <- input_list[[x]]
    callers <- unique(unlist(lapply(in_vec, function(f){
                stringr::str_split(f, "\\.")[[1]][2]
              })))
    out_ext <- gsub("-","_",
                   base::paste(stringr::str_split(
                                  gsub("\\.vcf","",
                                  in_vec[[1]]),
                                  "\\.")[[1]][-c(1,2)],
                        collapse="."))
    print(paste0("Working on: ", out_ext))

    ##test if RData already exists, then make if not
    ##uses both VEP annotated and raw VCFs
    dir_test <- dir(pattern = paste0(out_ext, ".RData"))
    if(length(dir_test) == 0){
      rdata_gr_list(in_vec, germline_id, callers, out_ext, raw_vcf_pattern)
    } else {
      print(paste0("Found: ", out_ext, ".RData, this is used"))
    }
  }

  ##load set of RData GRanges
  print("Loading GRanges...")
  rdata_in <- dir(pattern = ".RData")
  loaded_gr <- lapply(rdata_in, function(x){
    load(x, envir = .GlobalEnv)
  })
  names(loaded_gr) <- unlist(loaded_gr)

  ##set GRnages into lists
  raw_list <- get(loaded_gr$raw)
  var_list <- get(loaded_gr[[vcf_ext]])

  ##callers used
  callers <- names(var_list)
  two_callers <- c(callers[1], callers[2])

  ##ensure some variants in each sample to check on
  any_vars <- unlist(lapply(var_list, function(calr){
                      lapply(calr, length)
                    }))

  if(all(any_vars > 0)){

    ##get GRanges superset for HIGH, MODERATE IMPACTS from VEP
    grsuper_plot_list <- gr_super_alt_plot(var_list,
                                           raw_list,
                                           two_callers,
                                           impacts = c("HIGH", "MODERATE"),
                                           tag,
                                           included_order)

    ##run to get all impacts, print but not plot
    ##get GRanges superset per mutype
    grsuper_plot_list <- gr_super_alt_plot(var_list,
                                           raw_list,
                                           two_callers,
                                           impacts = c("HIGH", "MODERATE", "MODIFIER", "LOW"),
                                           tag,
                                           included_order)
  } else {
    print("No varaints found in one or more callers, please check and exclude")
  }
}

#' run function to make list of GRanges per caller
#'
#' @param in_vec vector of VCFs to parse
#' @param germline_id ID for germline_id sample
#' @param callers the variant callers for which VCFs are present in dir
#' @param out_ext the extension after caller in filename of VCF
#' @param raw_vcf_pattern pattern to match raw unannotated, unfiltered VCFs unfiltered
#' @return none, creates an RData object
#' @export

rdata_gr_list <- function(in_vec, germline_id, callers, out_ext, raw_vcf_pattern) {

  gr_list <- as.list(unique(callers))
  names(gr_list) <- unique(callers)

  for (callern in 1:length(unique(callers))){
    caller <- unique(callers)[callern]
    print(paste0("Caller: ", caller))

    ##parse VCFs, raw or VEP annotated
    raw_pattern <- gsub("\\.vcf", "", raw_vcf_pattern)
    if(raw_pattern %in% out_ext){
      ##parse VCFs from use_list based on caller, into gr_list
      gr_list[[caller]] <- lapply(in_vec[grep(caller, in_vec)], function(f){
        print(paste0("Parsing: ", f))
        suppressWarnings(vcf_parse_gr(f, germline_id))
      })
    } else {
      gr_list[[caller]] <- lapply(in_vec[grep(caller, in_vec)], function(f){
        print(paste0("Parsing: ", f))
        suppressWarnings(vcf_vep_ann_parse_soma_gr(f, germline_id))
      })
    }
    samps <- unlist(lapply(in_vec, function(f) {
      stringr::str_split(f, "\\.")[[1]][1]
    }))
    names(gr_list[[caller]]) <- samps[grep(caller, in_vec)]
  }

  ##assign output, save to dir
  assigned_name <- paste0(out_ext)
  assign(assigned_name, value = gr_list)
  save_file <- paste0(out_ext, ".RData")
  save(list = assigned_name, file = save_file)
}

#' Parses for VCFs into GRanges object
#'
#' @param vcf_in VCF path input
#' @param germline_id ID for germline_id sample
#' @return vector of single-letter HGVS protein IDs
#' @export

vcf_parse_gr <- function(vcf_in, germline_id){

  vcf <- VariantAnnotation::readVcf(vcf_in)
  if(dim(vcf)[1] != 0){
    gr <- suppressWarnings(customProDB::InputVcf(vcf_in))

    ##parse info
    infor <- VariantAnnotation::info(VariantAnnotation::header(vcf))

    ##somatic
    som_name <- names(gr)[names(gr) != germline_id]
    som <- gr[[som_name]]
    GenomeInfoDb::seqinfo(som) <- GenomeInfoDb::seqinfo(vcf)[GenomeInfoDb::seqlevels(som)]

    ##ensure an AF is there, pisces has VF instead (thanks pisces dev=D)
    if(! "AF" %in% names(S4Vectors::mcols(som))) {
      ad <- as.numeric(unlist(S4Vectors::mcols(som)["AD"]))
      ad1 <- as.numeric(unlist(S4Vectors::mcols(som)["AD.1"]))
      tot <- ad + ad1
      S4Vectors::mcols(som)$AF <- ad1 / tot
    }
    GenomeInfoDb::seqinfo(som) <- GenomeInfoDb::seqinfo(vcf)[GenomeInfoDb::seqlevels(som)]
    return(som)
  } else {
    print("No variants found")
    return(GenomicRanges::GRanges())
  }
}

#' Parses for VCFs annotated by VEP into GRanges object
#'    takes CANONICAL annotation or first if no CANONICAL
#'
#' @param vcf_in VCF path input
#' @param germline_id ID for germline_id sample
#' @return vector of single-letter HGVS protein IDs
#' @export

vcf_vep_ann_parse_soma_gr  <- function(vcf_in, germline_id){

  ##for single sample within a single VCF
  vcf <- suppressWarnings(VariantAnnotation::readVcf(vcf_in))

  if(dim(vcf)[1] != 0){
    gr <- suppressWarnings(customProDB::InputVcf(vcf_in))

    ##parse info
    infor <- VariantAnnotation::info(VariantAnnotation::header(vcf))

    ##VEP annotation naming
    ann_names <- unlist(stringr::str_split(infor[rownames(infor)=="ANN",]$Description, "\\|"))

    ##somatic
    som_name <- names(gr)[names(gr) != germline_id]
    som <- gr[[som_name]]
    GenomeInfoDb::seqinfo(som) <- GenomeInfoDb::seqinfo(vcf)[GenomeInfoDb::seqlevels(som)]

    ##annotation by CANONICAL, and add to mcols
    som_ann_df <- t(as.data.frame(lapply(som$ANN, function(ff){
      ffu <- unique(unlist(ff))
      ffuret <- unlist(lapply(stringr::str_split(ffu,"\\|"), function(fff){
        if(fff[ann_names == "CANONICAL"] == "YES"){
          if(length(fff)!=length(ann_names)){
            lengExtra <- length(ann_names) - length(fff)
            fff <-c(fff, rep("", lengExtra))
          }
          return(fff)
        }
      }))[1:length(ann_names)]
      if(length(ffuret)>0){
        return(ffuret[1:length(ann_names)])
      }
      else{
        return(rep("", length(ann_names)))
      }
    })))
    colnames(som_ann_df) <- ann_names

    if(sum(dim(som_ann_df)) != 0){
      S4Vectors::values(som) <- cbind(as.data.frame(S4Vectors::mcols(som)), som_ann_df)
      som$ANN <- NULL
    }
    som <- unique(som)
    return(som)
  }
  else{
    print("No variants found")
    return(GenomicRanges::GRanges())
  }
}

#' Create GRanges 'superset' from GRanges per-sample, per-caller variants
#'    contains all variants found per sample in any caller
#' @param var_list is a nested list of [[caller]][[samples1..n]]
#' @param name_callers is a set of two callers used for initial screening
#' @param impacts VEP impacts (one or combination of "HIGH", "MODERATE", "MODIFIER", "LOW")
#' @return GRanges 'superset' of all callers, and samples therein
#' @export

gr_super_set <- function(var_list, name_callers, impacts){

  print("Super-setting GRanges for IMPACTs:")
  print(impacts)
  ##set wanted mcols
  mcols_want <- c("AD", "AD.1", "AF", "Consequence", "IMPACT", "SYMBOL", "HGVSc", "HGVSp", "CLIN_SIG")

  if(length(name_callers) != 2){
    print("Require only 2 callers, no more and no less!")
  } else {

  ##set up output
  gr_super <- as.list(names(var_list[[1]]))
  caller_names <- names(var_list)
  samps <-  names(var_list[[1]])

  ##read first entry
  call_1 <- var_list[[name_callers[1]]]
  call_2 <- var_list[[name_callers[2]]]

  #exclude MT, GL
  seqwant <- c(seq(from=1, to=22, by=1), "X")

  ##iterate over samples in call_1, call_2 (calls from the two named callers)
  ##in the same sample
  if(length(call_1) > 1){
    for (x in seq_along(call_1)){
      print(paste0("Working on: ", names(call_1)[x]))

      calls_1 <- call_1[[x]]
      calls_2 <- call_2[[x]]
      GenomeInfoDb::seqlevels(calls_1, pruning.mode = "coarse") <- seqwant
      GenomeInfoDb::seqlevels(calls_2, pruning.mode = "coarse") <- seqwant

      ##test all wanted mcols exist, rename if "VF" not "AF" (Pisces)
      for(y in 1:2){
        if(length(mcols_want[mcols_want %in% names(S4Vectors::mcols(call_1[[y]]))]) != length(mcols_want)){
          gsub("VF","AF", names(S4Vectors::mcols(call_1[[y]])))
        }
        if(length(mcols_want[mcols_want %in% names(S4Vectors::mcols(call_1[[y]]))]) != length(mcols_want)){
          gsub("VF","AF", names(S4Vectors::mcols(call_1[[y]])))
        }
      }
      calls_1$HGVSp1 <- sub_hgvsp(calls_1$HGVSp)
      calls_2$HGVSp1 <- sub_hgvsp(calls_2$HGVSp)

      ##sets of calls_1, and those different between calls_1 and calls_2
      ##therefore represent all possible variants per sample from any caller
      gr_11 <- calls_1[calls_1$IMPACT %in% impacts, names(S4Vectors::mcols(calls_1)) %in% mcols_want]
      gr_22 <- calls_2[calls_2$IMPACT %in% impacts, names(S4Vectors::mcols(calls_2)) %in% mcols_want]
      gr_12 <- suppressWarnings(IRanges::subsetByOverlaps(gr_22, gr_11))
      if(length(gr_11) != 0 & length(gr_12) != 0){
        S4Vectors::mcols(gr_11) <- S4Vectors::mcols(gr_11)[, mcols_want]
        S4Vectors::mcols(gr_12) <- S4Vectors::mcols(gr_12)[, mcols_want]
        ##add 1 and difference (the superset of variants)
        gr_super[[x]] <- suppressWarnings(c(gr_11, gr_12))
        if(length(gr_11) != 0 & length(gr_12) == 0){
          gr_super[[x]] <- suppressWarnings(c(gr_11, gr_12))
        }
      }
    }
  } else {
      ##only one sample, return this
      calls_1 <- call_1[[1]]
      calls_2 <- call_2[[1]]

      ##test all wanted mcols exist, rename if "VF" not "AF" (Pisces)
      if(length(mcols_want[mcols_want %in% names(S4Vectors::mcols(call_1[[1]]))]) != length(mcols_want)){
          gsub("VF", "AF", names(S4Vectors::mcols(call_1[[1]])))
      }
      if(length(mcols_want[mcols_want %in% names(S4Vectors::mcols(call_2[[1]]))]) != length(mcols_want)){
          gsub("VF", "AF", names(S4Vectors::mcols(call_2[[1]])))
      }
      calls_1$HGVSp1 <- sub_hgvsp(calls_1$HGVSp)
      calls_2$HGVSp1 <- sub_hgvsp(calls_2$HGVSp)

      ##sets of call_1, 2 and the difference
      gr_11 <- calls_1[calls_1$IMPACT %in% impacts, names(S4Vectors::mcols(calls_1)) %in% mcols_want]
      gr_22 <- calls_2[calls_2$IMPACT %in% impacts, names(S4Vectors::mcols(calls_2)) %in% mcols_want]
      gr_12 <- suppressWarnings(dplyr::setdiff(gr_22, gr_11))

      if(length(gr_11)!=0 & length(gr_12) != 0){
        S4Vectors::mcols(gr_11) <- S4Vectors::mcols(gr_11)[, mcols_want]
        S4Vectors::mcols(gr_12) <- S4Vectors::mcols(gr_12)[, mcols_want]
      }

      ##add 1 and difference (the superset of variants)
      gr_super[[1]] <- suppressWarnings(c(gr_11, gr_12))
    }
    names(gr_super) <- samps
    return(gr_super)
  }
}

#' Find consensus of at least two callers using GRanges 'superset'
#'
#' @param var_list is a nested list of [[caller]][[samples1..n]]
#' @param gr_super is a GRanges superset from gr_super_set()
#' @param tag is a string used to tag output files
#' @return GRanges object of consensus per sample
#' @export

at_least_two <- function(var_list, gr_super, tag) {

  ##set TMB, no longer used as PCGR determines TMB
  # if(is.null(tmb)){
  #   tmb <- "null"
  # }

  ##set vars
  callers <- names(var_list)
  samps <- names(var_list[[1]])

  ##iterate over list of callers
  if(length(samps) > 1){
    gr_plots <- lapply(seq_along(samps), function(x) {
      samp <- samps[x]
      print(paste0("Working on: ", samp))
      combn_at_least_two(var_list, gr_super, callers, samp, tag)
    })
  } else {
    samp <- samps[1]
    print(paste0("Only one sample: ", samp))
    gr_plots <- combn_at_least_two(var_list, callers, samp, tag)
  }
  names(gr_plots) <- samps
  return(gr_plots)
}

#' Find consensus of at least two callers using GRanges superset
#'
#' @importFrom rlang :=
#' @param var_list is a nested list of [[caller]][[samples1..n]]
#' @param gr_super is a GRanges superset from gr_super_set()
#' @param callers are the named variant callers
#' @param samp is the sample being operated on
#' @param tag is a string to tag output files
#' @return GRanges object of consensus per sample
#' @export

combn_at_least_two <- function(var_list, gr_super, callers, samp, tag) {

  r <- POS <- AD <- AD1 <- ADSUM <- ID <- REF <- ALT <- QUAL <- FILTER <- INFO <- FORMAT <- sampleID <- NULL
  ##all possible combinations of intersects of callers
  ##output to new GRanges object
  up_l1 <- apply(t(utils::combn(length(callers), m = 2)), 1, function(xx){
    print(paste(callers[xx[1]], " vs. ", callers[xx[2]]))
    gr_1 <- var_list[[names(var_list)[xx[1]]]][[samp]]
    gr_2 <- var_list[[names(var_list)[xx[2]]]][[samp]]
    ##found in both
    suppressWarnings(dplyr::intersect(gr_2, gr_1))
  })

  ##with more than 2 variant callers you get output of more than one GRanges
  ##with only two callers, you get one comparison
  if(length(up_l1) > 1){
    for(xx in 2:length(up_l1)){
      up_l2 <- suppressWarnings(c(up_l1[[1]], up_l1[[xx]]))
    }
    gr_plot <- suppressWarnings(gr_super[[samp]][gr_super[[samp]] %in% unique(up_l2)])
  } else {
    gr_plot <- dplyr::intersect(gr_super[[samp]], unique(up_l1[[1]]))
  }

  #output files
  file_out <- paste0(samp, ".", tag, ".consensus.tsv")
  vcf_out <- paste0(samp, ".", tag, ".pcgr.all.tsv.vcf")

  # if(tmb == "snv"){
  #     tmb_out <- exome_tumour_mutation_burden(gr_plot)
  #     file_out <- paste0(samp, ".", tag, ".TMB_", tmb_out, "_SNV-Mb.consensus.tab")
  # }
  readr::write_tsv(as.data.frame(gr_plot), path = file_out)

  ##VCF output
  ###CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT
  ##1	33238563	.	G	A	.	PASS	TVAF=0.145833333333;TDP=49;	GT:DPC:DPT:ADC:ADT	0/1:.:DP:.,.:AD,AD.1
  vcf_grp <- tibble::as_tibble(as.data.frame(gr_plot), rownames="r")
  vcf_grps <- tidyr::separate(vcf_grp,
                              r,
                              c("#CHROM", "POS", "REF", "ALT"),
                              sep = "[:_\\/]")
  vcf_grpr <- dplyr::rename(vcf_grps, AD1='AD.1')
  vcf_grpm <- dplyr::mutate(vcf_grpr,
                            POS = as.integer(POS),
                            ID = ".",
                            QUAL = ".",
                            INFO = ".",
                            FILTER = "PASS",
                            FORMAT = "GT:DPC:DPT:ADC:ADT",
                            ADSUM = as.numeric(AD) + as.numeric(AD1),
                            sampleID = paste0("0/1:.:",
                                              ADSUM,
                                              ":.,.:",
                                              AD,
                                              ",",
                                              AD1))
    vcf_grpms <- dplyr::select(vcf_grpm, '#CHROM', POS, ID, REF, ALT, QUAL, FILTER, INFO,  FORMAT, sampleID)
    vcf_gr_plot <- dplyr::rename(vcf_grpms, !!samp:="sampleID")

  readr::write_tsv(as.data.frame(vcf_gr_plot), path = vcf_out)
  return(gr_plot)
}

#' Generate raw allele frequencies from list of GRanges
#'
#' @param raw_list is a nested list of raw calls [[caller]][[samples1..n]]
#' @param comb_gr is GRanges object combined with all samples
#' @return data.frame object of raw variant call allele frequencies
#' @export

raw_afs <- function(raw_list, comb_gr){
  ##ff contains all GR of each sample for the caller
  as.data.frame(lapply(raw_list, function(raw_samp_list){
    afs <- rep(as.numeric(0), length(comb_gr))
    ##per sample, test which variants are found in entire set
    lapply(seq_along(raw_samp_list), function(raw_samp){
      fff <- raw_samp_list[[raw_samp]]
      GenomeInfoDb::seqlevels(fff,  pruning.mode="coarse") <- GenomeInfoDb::seqlevels(comb_gr)
      ffs <- sort(fff)
      ffsi <- ffs[ffs %in% comb_gr]
      afs[names(comb_gr) %in% names(ffs)] <- as.numeric(S4Vectors::mcols(ffsi)$AF)
      return(afs)
    })
  }))
}

#' Create plot of shared variants among samples
#'
#' @param plot_list is a nested list of plot data [[caller]][[samples1..n]]
#' @param raw_list is a nested list of raw calls [[caller]][[samples1..n]]
#' @param tag is a string to tag output files
#' @param included_order oredering of samples for plotting
#' @return none, plots PDF and writes out tsv files
#' @export

plot_consensus_list <- function(plot_list, raw_list, tag, included_order){

  row_sum <- sum_01 <- row_sum_01 <- NULL

  print("Determining shared variants for plotting...")
  ##remove hyphens from names
  samps <- gsub("-","_", names(plot_list))
  included_order <- gsub("-", "_", included_order)

  ##combined set of all samples variants
  comb_gr <- suppressWarnings(unique(do.call("c", unname(plot_list))))
  GenomeInfoDb::seqlevels(comb_gr) <- sort(GenomeInfoDb::seqlevels(comb_gr))
  comb_gr <- sort(comb_gr)
  comb_df <- as.data.frame(comb_gr)

  ##labels for plot
  hgvsp <- unlist(lapply(comb_gr$HGVSp1,function(f){
    strsplit(f,"\\.")[[1]][3]
  }))
  uniq_labels <- paste(names(comb_gr),
                     comb_df$SYMBOL,
                     comb_df$Consequence,
                     hgvsp, sep=":")

  ##take comb_gr positions, then query raw calls for them
  ##allows 'rescue' of those falling out because of filters
  ##use raw allele frequencies
  plot_df_raw_afs <- raw_afs(raw_list, comb_gr)

  ##rename based on callers.samp
  colnames(plot_df_raw_afs) <- unlist(lapply(names(raw_list), function(f){
    paste(f, samps, sep=".")
  }))
  ##take maximum AF from the two callers
  plot_df_raw_max <- do.call(cbind, lapply(samps, function(ss){
    apply(plot_df_raw_afs[, grep(ss, colnames(plot_df_raw_afs))], 1, max)
  }))

  plot_df_raw_max <- as.data.frame(plot_df_raw_max)
  colnames(plot_df_raw_max) <- samps
  rownames(plot_df_raw_max) <- uniq_labels

  ##remove those with all 0 frequency
  plot_df_raw_order <- plot_df_raw_max[rowSums(plot_df_raw_max)!=0,]

  ##rowSums for arranging
  ##gtz function to create '!= 0' -> 1, else 0
  ##sum those to determine how many are shared (sum_01 > 1)
  ##arrange by sum_01, then row_sum
  gtz <- function(x) {
    ifelse(x > 0, 1, 0)
  }
  pdf_rawr <- tibble::as_tibble(plot_df_raw_order, rownames = "label")
  pdf_rawm <- dplyr::mutate(pdf_rawr, row_sum = rowSums(pdf_rawr[-1]))
  pdf_rawi <- dplyr::mutate_if(pdf_rawm, is.numeric, dplyr::funs("01" = gtz))
  pdf_raws <- dplyr::select(pdf_rawi, -row_sum_01)
  pdf_rawm <- dplyr::mutate(pdf_raws, sum_01 = rowSums(dplyr::select(pdf_raws, tidyselect::ends_with("_01"))))
  pdf_rawd <- as.data.frame(pdf_rawm)
  pdf_rawr <- dplyr::arrange(pdf_rawd, sum_01, row_sum)
  print(colnames(pdf_rawr))
  plot_df_raw_order_ar <- dplyr::select(pdf_rawr, 1, samps, sum_01)

  ##those shared between samples
  plot_df_sh <- plot_df_raw_order_ar[, c("label", included_order, "sum_01")]
  plot_df_sha <- dplyr::filter(plot_df_sh, sum_01 > 1)
  plot_df_shar <- dplyr::select(plot_df_sha, -sum_01)
  plot_df_shared <- tibble::column_to_rownames(plot_df_shar, "label")

  ##write out those in plotDF
  write_consensus_all(plot_list = plot_list,
                      plot_df = plot_df_shared,
                      tag = tag,
                      included_order = included_order,
                      cons = "shared")
  write_consensus_all(plot_list = plot_list,
                      plot_df = plot_df_raw_order,
                      tag = tag,
                      included_order = included_order,
                      cons="all")

  plot_vec <- plot_tag <- c()

  if(dim(plot_df_shared)[1] == 0){
    print("No shared variants found...")
  } else {
    ##shared variants
    plot_vec <- plot_df_shared
    plot_tag <- "shared"
    print(paste0("Plotting ", dim(plot_df_shared)[1], " shared variants..."))
  }

  if(dim(plot_vec)[1] != 0){
    plot_labels <- rep("", dim(plot_vec)[1])
    row_fontsize <- 1
    colz <- grDevices::colorRampPalette(c("lightgrey", "dodgerblue", "blue"))
    ##plotting and whether to use labels, size of labels
    if(dim(plot_vec)[1] < 120){
      row_fonttype = "bold"
      if(dim(plot_vec)[1] < 20){row_fontsize = 12}
      if(dim(plot_vec)[1] < 20){row_fontsize = 8}
      if(dim(plot_vec)[1] < 50){row_fontsize = 6}
      if(dim(plot_vec)[1] > 50 & dim(plot_vec)[1] < 100){row_fontsize = 4}
      if(dim(plot_vec)[1] > 100){row_fontsize = 2}
      plot_labels <- rownames(plot_vec)
    }

    grDevices::pdf(paste0(tag, ".", plot_tag, ".pdf"), onefile = F)
    pheatmap::pheatmap(plot_vec[, c(1:length(included_order))],
       breaks = seq(from = 0, to = 0.5, length.out = 101),
       color = colz(100),
       cluster_rows = FALSE,
       cluster_cols = FALSE,
       clustering_distance_rows = NA,
       cellwidth = 12,
       legend = TRUE,
       fontsize_row = row_fontsize,
       labels_row = plot_labels,
       border_color = "lightgrey",
       gaps_col = c(1:length(included_order))
    )
    grDevices::dev.off()
  }
}

#' Create two plots: all consensus, those in 2+ samples
#'
#' @importFrom rlang .data
#' @param plot_list is a nested list of plot data [[caller]][[samples1..n]]
#' @param raw_list is a nested list of raw calls [[caller]][[samples1..n]]
#' @param tag is a string to tag output files
#' @param included_order ordering of samples for plotting
#' @return data.frame object of raw variant call allele frequencies
#' @export

plot_consensus_single <- function(plot_list, raw_list, tag, included_order){

  ##remove hyphens
  samp <- names(plot_list)[1]

  ##combined set of all samples
  comb_gr <- suppressWarnings(unique(plot_list))
  GenomeInfoDb::seqlevels(comb_gr) <- sort(GenomeInfoDb::seqlevels(comb_gr))
  comb_gr <- sort(comb_gr)
  comb_df <- as.data.frame(comb_gr)

  ##labels for plot
  hgvsp <- unlist(lapply(comb_gr$HGVSp1,function(f){
    strsplit(f,"\\.")[[1]][3]
  }))
  uniq_labels <- paste(names(comb_gr),
                       comb_df$SYMBOL,
                       comb_df$Consequence,
                       hgvsp, sep=":")

  ##take those positions, then query raw calls
  ##allows 'rescue' of those falling out from arbitrary filters
  ##enough support previously to allow re-entry
  plot_df_raw_afs <- raw_afs(raw_list, comb_gr)

  ##rename based on callers.samp
  colnames(plot_df_raw_afs) <- unlist(lapply(names(raw_list), function(f){
    paste(f, samp, sep=".")
  }))

  ##take maximum AF from the two callers
  plot_df_raw_max <- do.call(cbind, lapply(samp, function(ss){
    apply(plot_df_raw_afs[, grep(ss, colnames(plot_df_raw_afs))], 1, max)
  }))
  plot_df_raw_max <- as.data.frame(plot_df_raw_max)
  colnames(plot_df_raw_max) <- samp
  rownames(plot_df_raw_max) <- uniq_labels

  ##ordering
  plot_df_raw_ordered <- tibble::as_tibble(plot_df_raw_max, rownames = "label")
  plot_df_raw_ordered <- dplyr::arrange(.data = plot_df_raw_ordered, .data[[2]])
  plot_df_raw_ordered <- base::as.data.frame(plot_df_raw_ordered)

  if(dim(plot_df_raw_ordered)[1] != 0){
    plot_vec <- data.frame(row.names = plot_df_raw_ordered[,1],
                           tag = plot_df_raw_ordered[,2])
    plot_tag <- "variants"
    plot_labels <- rep("",times=dim(plot_vec)[1])
    row_fontsize <- 1
    colz <- grDevices::colorRampPalette(c("lightgrey", "dodgerblue", "blue"))
    if(dim(plot_vec)[1] < 120){
      if(dim(plot_vec)[1] < 20){row_fontsize = 8}
      if(dim(plot_vec)[1] < 50){row_fontsize = 6}
      if(dim(plot_vec)[1] > 50 & dim(plot_vec)[1] < 100){row_fontsize = 4}
      if(dim(plot_vec)[1] > 100){row_fontsize = 2}
      plot_labels <- rownames(plot_vec)
    }

    grDevices::pdf(paste0(tag, ".", plot_tag, ".consensus.onlyOverlap.pdf"), onefile=F)
    pheatmap::pheatmap(plot_vec,
       breaks = seq(from = 0, to = 0.5, length.out = 101),
       color = colz(100),
       cluster_rows = FALSE,
       cluster_cols = FALSE,
       clustering_distance_rows = NA,
       cellwidth = 12,
       legend = TRUE,
       fontsize_row = row_fontsize,
       labels_row = plot_labels,
       labels_col = tag,
       border_color = "lightgrey"
    )
    grDevices::dev.off()
  }
}

#' Return overlapping variants for all samples
#'
#' @param plot_list list produced by at_least_two()
#' @param plot_df data.frame of plotting information
#' @param tag is a string used to tag output files
#' @param included_order ordering of samples for plotting
#' @param cons string to define consensus (shared or all)
#' @return GRanges object of all  of single-letter HGVS protein IDs
#' @export

write_consensus_all <- function(plot_list, plot_df, tag, included_order, cons){

  ##input names
  samps <- gsub("-","_",names(plot_list))

  ##split on colon, paste into label for plot_list
  strSplitFun <- function(input, sepn){
    lapply(input, function(f){
      strsplit(f, sepn)[[1]]
    })
  }
  labels_out <- unlist(lapply(strSplitFun(rownames(plot_df), ":"), function(f){
    paste0(f[1],":",f[2])
  }))

  ##parse required from plot_list
  plot_out <- lapply(plot_list, function(f){
    f[names(f) %in% labels_out]
  })

  ##test any output
  out_gr <- suppressWarnings(unique(do.call("c", unname(plot_out))))
  GenomeInfoDb::seqlevels(out_gr) <- sort(GenomeInfoDb::seqlevels(out_gr))
  out_gr <- sort(out_gr)
  out_df <- as.data.frame(out_gr)

  if(dim(out_df)[1] != 0){
    ##write output
    file_out <- paste0(tag, ".", cons, ".tsv")
    readr::write_tsv(out_df, path = file_out)
  }
}

#' Create single-letter HGVS protein annotation (VEP outputs 3-letter)
#'    make vector, gsub out aa3 for aa1
#'
#' @param in_vec vector input
#' @return vector of single-letter HGVS protein IDs
#' @export

sub_hgvsp <- function(in_vec){

  aa1 <- c("A", "R", "N", "D", "C", "Q", "E", "G", "H", "I", "L", "K", "M", "F", "P", "S", "T", "W", "Y", "V", "X", "D", "R", "C", "C", "C", "C", "C", "C", "C", "C", "H", "G", "H", "H", "H", "H", "H", "H", "D", "K", "K", "M", "K", "M", "C", "F", "Y", "S", "T")
  ##amino acid 3 letter to gsub HGVSp
  aa3 <- c("Ala","Arg", "Asn", "Asp", "Cys", "Gln", "Glu", "Gly", "His", "Ile", "Leu", "Lys", "Met", "Phe", "Pro", "Ser", "Thr", "Trp", "Tyr", "Val", "Aba", "Ash", "Cir", "Cme", "Cmt", "Csd", "Cso", "Csw", "Csx", "Cym", "Cyx", "Dde", "Glh", "Hid", "Hie", "Hip", "Hsd", "Hse", "Hsp", "Ias", "Kcx", "Lyn", "Mho", "Mly", "Mse", "Ocs", "Pff", "Ptr", "Sep", "Tpo")

  ##include * for Ter
  aa1 <-c(aa1,"*")
  aa3 <- c(aa3, "Ter")

  unlist(lapply(in_vec,function(f){
    #check matches (should be none or two)
    a3 <- aa3[!is.na(unlist(stringi::stri_match_all(f,regex=aa3)))]
    a1 <- aa1[!is.na(unlist(stringi::stri_match_all(f,regex=aa3)))]
    ##beauty:
    #https://stackoverflow.com/questions/19424709/r-gsub-pattern-vector-and-replacement-vector
    if(length(a3)>0){
      names(a1) <- a3
      stringr::str_replace_all(f,a1)
    }
    else{
      return("")
    }
  }))
}

#' Wrapper of main functions
#'
#' @param var_list is a nested list of VEP annotated [[caller]][[samples1..n]]
#' @param raw_list is a nested list of unfiltered calls [[caller]][[samples1..n]]
#' @param impacts VEP impacts (one or combination of "HIGH", "MODERATE", "MODIFIER", "LOW")
#' @param name_callers two of the variant callers
#' @param tag is a string used to tag output
#' @param included_order oredering of samples for plotting
#' @return GRanges object of all  of single-letter HGVS protein IDs
#' @export

gr_super_alt_plot <- function(var_list, raw_list, name_callers, impacts, tag, included_order) {

  ##GRanges superset
  gr_super <- gr_super_set(var_list, name_callers, impacts)

  ##get list to plot from with at least two callers supporting
  plot_list <- at_least_two(var_list, gr_super, tag)

  ##if single sample plot_list is actually a GRanges object(!)
  if(!is.list(plot_list)){
    plot_consensus_single(plot_list, raw_list, tag)
  } else {
    if(length(plot_list[[1]]) != 0 & length(plot_list[[2]]) != 0){
      plot_consensus_list(plot_list, raw_list, tag, included_order)
    }
    if(length(plot_list[[1]]) == 0 & length(plot_list[[2]]) == 0){
      print(paste0("No variants for IMPACTS: ", impacts, ", support across callers lacking"))
    }
  }
  return(list(gr_super, plot_list))
}
