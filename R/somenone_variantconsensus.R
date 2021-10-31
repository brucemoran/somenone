#' somatic_variant_consensus functions
#'
#' Primary function to call others and produce plots, table
#' @param germline_id ID for germline_id sample
#' @param vep_vcf_pattern pattern to match VEP annotated VCFs
#' @param raw_vcf_pattern pattern to match raw unannotated, unfiltered VCFs unfiltered
#' @param which_genome hg19 or hg38
#' @param included_order oredering of samples for plotting
#' @param name_callers named variant callers to use primary to random selection
#' @param tag is a string to tag output files
#' @return vector of single-letter HGVS protein IDs
#' @export

variant_consensus <- function(germline_id, vep_vcf_pattern, raw_vcf_pattern = "raw.vcf", tag = "somatic_n_of_1", which_genome, included_order = NULL, name_callers = NULL, impacts = NULL) {

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
  if(length(grep(",", included_order)) > 0){
      included_order <- stringr::str_split(included_order, ",")[[1]]
  } else {
      included_order <- unique(unlist(lapply(vcf_list, function(f){
                          stringr::str_split(f, "\\.")[[1]][1]
                        })))
  }

  ##which_genome
  if(! which_genome %in% c("hg19", "hg38")){
    if(length(grep("37", which_genome))==1){
      which_genome <- "hg19"
    } else {
      if(length(grep("19", which_genome))==1){
        which_genome <- "hg19"
      } else {
        which_genome <- "hg38"
      }
    }
  }

  ##operate over vep, raw VCFs
  callers <- unique(unlist(lapply(input_list[[1]], function(f){
              stringr::str_split(f, "\\.")[[1]][2]
            })))
  out_ext <- gsub("-","_",
                 base::paste(stringr::str_split(
                                gsub("\\.vcf","",
                                input_list[[1]]),
                                "\\.")[[1]][-c(1,2)],
                      collapse="."))

  ##lists from VEP annotated
  var_list <- somenone::rdata_gr_list(input_list[[1]],
                                      germline_id,
                                      callers,
                                      out_ext = "snv_indel.pass.vep",
                                      raw_vcf_pattern)

  ##callers used
  callers <- names(var_list)
  two_callers <- c(callers[1], callers[2])

  ##ensure some variants in each sample to check on
  any_vars <- unlist(lapply(seq_along(var_list), function(x){
                     lapply(seq_along(var_list[[x]]), function(y){
                       length(var_list[[x]][[y]])
                     })}))
  ##run twice as for some reason results differ from 1st to 2nd, but not beyond(?)
  any_vars <- unlist(lapply(seq_along(var_list), function(x){
                    lapply(seq_along(var_list[[x]]), function(y){
                      length(var_list[[x]][[y]])
                    })}))
  if(all(any_vars > 0)){

    ##define impacts from comma separated impacts input
    if(! is.null(impacts)){
      impact <- strsplit(impacts, ",")[[1]]
    } else {
      impact <- c("HIGH", "MODERATE", "MODIFIER", "LOW")
    }
    ##string defining impacts
    impact_str <- paste(unlist(lapply(impact, function(f){
          strsplit(f,"")[[1]][1]
        })), collapse = "")

    gr_super_plot_out <- somenone::gr_super_alt_plot(var_list = var_list,
                                                     name_callers = two_callers,
                                                     impacts = impact,
                                                     taga = paste0(tag, ".",  impact_str, "_impacts"),
                                                     included_order,
                                                     which_genome)

   ##get GRanges superset for HIGH, MODERATE IMPACTS from VEP
   gr_super_plot_out_list <- list(gr_super_plot_out = gr_super_plot_out,
                                  var_list = var_list,
                                  name_callers = two_callers,
                                  impacts = impact,
                                  taga = paste0(tag, ".", impact_str, "_impacts"),
                                  included_order = included_order,
                                  which_genome = which_genome, call = "gr_super_alt_plot(gr_super_alt_plot_list$var_list, gr_super_alt_plot_list$name_callers, gr_super_alt_plot_list$impacts, gr_super_alt_plot_list$taga, gr_super_alt_plot_list$included_order, gr_super_alt_plot_list$which_genome)")

   save(gr_super_plot_out_list,
        file = paste0(paste0(tag, ".", impact_str, "_impacts"), ".gr_super_alt_plot.RData"))

  } else {
    print("No variants found in one or more callers, please check and exclude")
    vcf_out <- paste0(names(var_list[[1]]), ".no_vars.impacts.pcgr.tsv.vcf")
    readr::write_tsv(data.frame(), file = vcf_out)
    file_out <- paste0(names(var_list[[1]]), ".consensus.tsv")
    readr::write_tsv(data.frame(), file = file_out)
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
        suppressWarnings(somenone::vcf_parse_gr(f, germline_id))
      })
    } else {
      gr_list[[caller]] <- lapply(in_vec[grep(caller, in_vec)], function(f){
        print(paste0("Parsing: ", f))
        suppressWarnings(somenone::vcf_vep_ann_parse_soma_gr(f, germline_id))
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
  return(gr_list)
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
  mcols_want <- c("AD", "AD.1", "AF", "Consequence", "IMPACT", "SYMBOL", "HGVSc", "HGVSp", "HGVSp1", "CLIN_SIG")

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
      gr_12 <- granges_sdin(gr_22, gr_11, "setdiff")

      S4Vectors::mcols(gr_11) <- S4Vectors::mcols(gr_11)[, mcols_want]
      S4Vectors::mcols(gr_12) <- S4Vectors::mcols(gr_12)[, mcols_want]

      ##add 1 and difference (the superset of variants)
      gr_super[[x]] <- c(gr_11, gr_12)
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
      calls_1$HGVSp1 <- somenone::sub_hgvsp(calls_1$HGVSp)
      calls_2$HGVSp1 <- somenone::sub_hgvsp(calls_2$HGVSp)

      ##sets of call_1, 2 and the difference
      gr_11 <- calls_1[calls_1$IMPACT %in% impacts, names(S4Vectors::mcols(calls_1)) %in% mcols_want]
      gr_22 <- calls_2[calls_2$IMPACT %in% impacts, names(S4Vectors::mcols(calls_2)) %in% mcols_want]
      gr_12 <- somenone::granges_sdin(gr_22, gr_11, "setdiff")

      S4Vectors::mcols(gr_11) <- S4Vectors::mcols(gr_11)[, mcols_want]
      S4Vectors::mcols(gr_12) <- S4Vectors::mcols(gr_12)[, mcols_want]

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

at_least_two <- function (var_list, gr_super, tag){
  print("Finding at-least-two callers supporting variants...")
  callers <- names(var_list)
  samps <- names(var_list[[1]])
  print("Samples available:")
  print(samps)
  gr_plot <- NULL
  gr_plots <- lapply(seq_along(samps), function(x) {
      print(paste0("Working on: ", samps[x]))
      samp <- samps[x]
      gr_samp <- gr_super[[samp]]
      up_l1 <- apply(t(utils::combn(length(callers), m = 2)),
          1, function(xx) {
              print(paste(callers[xx[1]], " vs. ", callers[xx[2]]))
              gr_1 <- var_list[[names(var_list)[xx[1]]]][[samp]]
              gr_2 <- var_list[[names(var_list)[xx[2]]]][[samp]]
              gr_in <- granges_sdin(gr_2, gr_1, "intersect")
          })
      if(length(up_l1) > 1) {
          up_l2 <- up_l1[[1]]
          for (xx in 2:length(up_l1)) {
              up_l2 <- suppressWarnings(c(up_l2, up_l1[[xx]]))
          }
          if (!is.null(names(gr_super[[samp]]))) {
              gr_plot <- gr_samp[names(gr_samp) %in%
                unique(names(up_l2))]
          }
          else {
              gr_plot <- GRanges()
          }
      } else {
          gr_plot <- granges_sdin(gr_super[[samp]], unique(up_l1[[1]]), "intersect")
      }
      if (!length(names(gr_plot)) == 0) {
        file_out <- paste0(samp, ".", tag, ".consensus.tsv")
        vcf_out <- paste0(samp, ".", tag, ".pcgr.tsv.vcf")
        readr::write_tsv(as.data.frame(gr_plot), file = file_out)
        vcf_grp <- tibble::as_tibble(as.data.frame(gr_plot))
        vcf_grp <- dplyr::mutate(vcf_grp, r = names(gr_plot))
        vcf_grps <- tidyr::separate(vcf_grp, r, c("#CHROM",
            "POS", "REF", "ALT"), sep = "[:_\\/]")
        vcf_grpr <- dplyr::rename(vcf_grps, AD1 = "AD.1")
        vcf_grpm <- dplyr::mutate(vcf_grpr, POS = as.integer(POS),
            ID = ".", QUAL = ".", INFO = ".", FILTER = "PASS",
            FORMAT = "GT:DPC:DPT:ADC:ADT", ADSUM = as.numeric(AD) +
              as.numeric(AD1), sampleID = paste0("0/1:.:",
              ADSUM, ":.,.:", AD, ",", AD1))
        vcf_grpms <- dplyr::select(vcf_grpm, "#CHROM", POS,
            ID, REF, ALT, QUAL, FILTER, INFO, FORMAT, sampleID)
        colnames(vcf_grpms)[colnames(vcf_grpms) == "sampleID"] <- samp
        readr::write_tsv(as.data.frame(vcf_grpms), file = vcf_out)
      } else {
        file_out <- paste0(samp, ".", tag, ".consensus.tsv")
        readr::write_tsv(as.data.frame(gr_plot), file = file_out)
        vcf_out <- paste0(names(var_list[[1]]), ".no_vars.", tag, ".pcgr.tsv.vcf")
        readr::write_tsv(data.frame(), file = vcf_out)
      }
      return(gr_plot)
  })
  names(gr_plots) <- samps
  return(gr_plots)
}
#' Create plot of shared variants among samples
#'
#' @param granges is a GRanges object
#' @param tag is a string to tag output files
#' @param sampleID name of sample
#' @param sample_map map included_order to new names, must be name vector where
#'        names equate to sampleID, vector element is name to change to
#' @param colours to use for colouring/shading, made into a rampPalette
#' @param plot_label_pattern match this to print label of variant
#'        (NB too vague and you will have a huge amount of labels which looks shit;
#'        currently set to show 'patho'genic)
#' @return none, plots PDF and writes out tsv files
#' @export

plot_single <- function(granges, tag, sampleID, sample_map = NULL, colours = NULL, plot_label_pattern = "patho"){

  row_sum <- sum_01 <- row_sum_01 <- NULL

  ##colouring
  if(is.null(colours)){
    colours <- c("lightgrey", "dodgerblue", "blue")
  }

  if(length(granges) == 0){
    print("No variants...")
  } else {

    ##save everything into a RData to allow rerunning with other options
    plot_single_list <- list(granges, tag, sampleID, sample_map, colours, plot_label_pattern, "plot_single(granges, tag, sampleID, sample_map, colours, plot_label_pattern)")
    names(plot_single_list) <- c("master_gr", "tag", "included_order", "sample_map", "colours", "plot_label_pattern", "plot_single_call")
    save(plot_single_list,
         file = paste0(tag, ".plot_single.RData"))

    ##remove hyphens from names
    samps <- gsub("-","_", sampleID)

    ##combined set of all samples variants
    comb_df <- as.data.frame(granges)

    ##labels for plot
    hgvsp <- unlist(lapply(granges$HGVSp1,function(f){
      strsplit(f,"\\.")[[1]][3]
    }))
    uniq_labels <- make.unique(gsub(" : NA", "",
                              gsub(" : $", "", paste(comb_df$SYMBOL,
                                                      comb_df$Consequence,
                                                      hgvsp,
                                                      comb_df$CLIN_SIG,
                                                      sep=" : "))))

    ##set up plotting
    plot_af <- as.data.frame(comb_df[,"AF"])
    if(!is.null(sample_map)){
      sampleID <- sample_map[1]
    }
    colnames(plot_af) <- paste0(sampleID, ".AF")
    rownames(plot_af) <- uniq_labels

    ##arrange rows
    plot_af <- dplyr::arrange(.data = plot_af, dplyr::across(colnames(plot_af)))
    plot_adf <- dplyr::mutate(.data = plot_af, dplyr::across(where(is.character), as.numeric))
    rownames(plot_adf) <- rownames(plot_af)
    plot_af <- plot_adf

    ##change names
    print(paste0("Plotting ", dim(plot_af)[1], " variants..."))

    ##plot_setup
    row_fontsize <- 1
    colz <- grDevices::colorRampPalette(colours)
    ##plotting and whether to use labels, size of labels
    if(dim(plot_af)[1] < 120){
      row_fonttype = "bold"
      if(dim(plot_af)[1] < 20){row_fontsize = 12}
      if(dim(plot_af)[1] < 20){row_fontsize = 8}
      if(dim(plot_af)[1] < 50){row_fontsize = 6}
      if(dim(plot_af)[1] > 50 & dim(plot_af)[1] < 100){row_fontsize = 4}
      if(dim(plot_af)[1] > 100){row_fontsize = 2}
      plot_labels <- rownames(plot_af)
    } else {
      ###only include rownames that are pathogenic
      row_fonttype = "bold"
      plot_labels <- grep(plot_label_pattern, rownames(plot_af), value = TRUE)
      if(length(plot_labels)[1] < 20){row_fontsize = 12}
      if(length(plot_labels)[1] < 20){row_fontsize = 8}
      if(length(plot_labels)[1] < 50){row_fontsize = 6}
      if(length(plot_labels)[1] > 50 & length(plot_labels)[1] < 100){row_fontsize = 4}
      if(length(plot_labels)[1] > 100){row_fontsize = 2}
    }

    grDevices::pdf(paste0(tag, ".pdf"), onefile = F)
    pheatmap::pheatmap(plot_af,
       breaks = seq(from = 0, to = 0.5, length.out = 101),
       color = colz(100),
       cluster_rows = FALSE,
       cluster_cols = FALSE,
       clustering_distance_rows = NA,
       cellwidth = 12,
       legend = TRUE,
       fontsize_row = row_fontsize,
       labels_row = plot_labels,
       border_color = "lightgrey")
    grDevices::dev.off()
  }
}

#' Create plot of shared variants among samples
#'
#' @param master_gr is a named list of GRanges object [[samples1..n]]
#' @param tag is a string to tag output files
#' @param included_order ordering of samples for plotting
#' @param sample_map map included_order to new names, must be name vector where
#'        names equate to included_order elements
#' @param colours to use for colouring/shading, made into a rampPalette, order according to low to high allele frequency
#' @param plot_label_pattern match this to print label of variant
#'        (NB too vague and you will have a huge amount of labels which looks shit;
#'        currently set to show 'patho'genic)
#' @return none, plots PDF and writes out tsv files
#' @export

plot_consensus <- function(master_gr, tag, included_order, sample_map = NULL, colours = NULL, plot_label_pattern = "patho"){

  row_sum <- sum_01 <- row_sum_01 <- NULL

  ##colouring
  if(is.null(colours)){
    colours <- c("lightgrey", "dodgerblue", "blue")
  }

  if(length(master_gr) == 0){
    print("No variants...")
  } else {

    ##save everything into a RData to allow rerunning with other options
    plot_consensus_list <- list(master_gr, tag, included_order, sample_map, colours, plot_label_pattern, "somenone::plot_consensus(master_gr, tag, included_order, sample_map, colours, plot_label_pattern)")
    names(plot_consensus_list) <- c("master_gr", "tag", "included_order", "sample_map", "colours", "plot_label_pattern", "plot_consensus_call")
    save(plot_consensus_list,
         file = paste0(tag, ".plot_consensus.RData"))

    ##remove hyphens from names
    samps <- gsub("-","_", included_order)
    included_order <- gsub("-", "_", included_order)

    ##combined set of all samples variants
    comb_gr <- master_gr
    comb_df <- as.data.frame(comb_gr)

    ##labels for plot
    hgvsp <- unlist(lapply(comb_gr$HGVSp1,function(f){
      strsplit(f,"\\.")[[1]][3]
    }))
    uniq_labels <- make.unique(gsub(" : NA", "",
                              gsub(" : $", "", paste(comb_df$SYMBOL,
                                                      comb_df$Consequence,
                                                      hgvsp,
                                                      comb_df$CLIN_SIG,
                                                      sep=" : "))))
    ##set up plotting
    if(!is.null(sample_map)){
      ##change names
      plot_af <- comb_df[, c("samples_n", paste0(included_order,".AF"))]
      colnames(plot_af) <- c("samples_n", unlist(lapply(seq_along(sample_map), function(x){
        gsub(names(sample_map)[x], sample_map[x], colnames(plot_af)[match(paste0(names(sample_map)[x], ".AF"), colnames(plot_af))])
      })))
    } else {
      plot_af <- comb_df[,  c("samples_n", paste0(included_order,".AF"))]
    }
    colnames(plot_af) <- gsub(".AF", "", colnames(plot_af))
    print(paste0("Plotting ", dim(plot_af)[1], " variants..."))

    ##remove NAs (set to 0)
    plot_af[is.na(plot_af)] <- 0
    plot_af <- sapply(plot_af, as.numeric)
    if(is.null(dim(plot_af))){
      tt <- t(as.data.frame(plot_af))
      plot_af <- tt
      rownames(plot_af) <- uniq_labels
    } else {
      rownames(plot_af) <- uniq_labels
    }

    ##order plot based on input set (shared or all)
    if(length(grep("shared", tag)) == 1){

      plot_af <- plot_af[order(rowSums(plot_af), plot_af[,2]),-1]

    } else {

      ##arrange on each colname
      plot_adf <- as.data.frame(plot_af[,-1])
      plot_af <- dplyr::arrange(.data = plot_adf, dplyr::across(colnames(plot_adf)))

    }

    ##plot_setup
    row_fontsize <- 1
    colz <- grDevices::colorRampPalette(colours)
    ##plotting and whether to use labels, size of labels
    if(is.null(dim(plot_af))){
      row_fonttype = "bold"
      row_fontsize = 12
      plot_labels <- rownames(plot_af)
    } else {
      if(dim(plot_af)[1] < 120){
        row_fonttype = "bold"
        if(dim(plot_af)[1] < 20){row_fontsize = 12}
        if(dim(plot_af)[1] < 20){row_fontsize = 8}
        if(dim(plot_af)[1] < 50){row_fontsize = 6}
        if(dim(plot_af)[1] > 50 & dim(plot_af)[1] < 100){row_fontsize = 4}
        if(dim(plot_af)[1] > 100){row_fontsize = 2}
        plot_labels <- rownames(plot_af)
      } else {
        ###only include rownames that are pathogenic
        row_fonttype = "bold"
        plot_labels <- grep(plot_label_pattern, rownames(plot_af), value = TRUE)
        if(length(plot_labels)[1] < 20){row_fontsize = 12}
        if(length(plot_labels)[1] < 20){row_fontsize = 8}
        if(length(plot_labels)[1] < 50){row_fontsize = 6}
        if(length(plot_labels)[1] > 50 & length(plot_labels)[1] < 100){row_fontsize = 4}
        if(length(plot_labels)[1] > 100){row_fontsize = 2}
      }
    }

    if(length(plot_af)>2){
      gaps_col <- c(1:length(included_order))
    } else {
      gaps_col <- NULL
    }

    grDevices::pdf(paste0(tag, ".pdf"), onefile = F)
    pheatmap::pheatmap(plot_af,
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
       gaps_col = gaps_col)
    grDevices::dev.off()
  }
}

#' Find intersect from list of SNV GRanges
#' @param gr_list list of named GRanges objects
#' @param ps_vec mcols columns to keep p(er) s(ample; appended with list element's name)
#' @param dp_vec mcols columns to d(edu)p(licate); appended as-is if all same; if varaitions, usual 'dot number' formatting applies)
#' @param tag to apply to output file (master)
#' @param which_genome hg19 or hg38
#' @return list of GRanges object of shared intersecting, and all mutations with ps_vec, dp_vec columns, and original rownames (should be unique therefore)
#' @export

master_intersect_snv_grlist <- function(gr_list, ps_vec, dp_vec, tag, which_genome){

  options(scipen=999)
  
  ##check gr is A GRanges object
  if(!as.vector(class(gr_list)) %in% c("GRangesList", "list")){
    stop("Input \'gr_list\' is not a GRangesList nor list object, retry")
  }

  ##use bedr::bedr.join.multiple.region to make a complete set and samples
  ##lapply over gr_list which creates character vector list in "chrX:1-100" format
  ##thinks it's 0-based so add 1...
  chr_list <- lapply(gr_list, function(f){
      bs <- apply(as.data.frame(unique(f)), 1, function(ff){
          ff <- unlist(ff)
          paste0("chr", ff[1], ":", as.numeric(gsub(" ", "", ff[2])), "-", as.numeric(gsub(" ", "", ff[3]))+1)
        })
        bedr::bedr.sort.region(bs)
    })
  join_chr_all_tb <- tibble::as_tibble(bedr::bedr.join.multiple.region(chr_list, build = which_genome))

  ##if two SNV are adjacent, above will join them
  ##so unjoin them
  ##remove MNVs
  mnv_index <- unlist(lapply(join_chr_all_tb$index, function(f){
    splt <- strsplit(unlist(f), ":|-")[[1]]
    if(as.numeric(splt[3])-as.numeric(splt[2]) != 1){
      return(f)
    }
  }))

  not_line <- join_chr_all_tb[!join_chr_all_tb$index %in% mnv_index,]

  if(dim(not_line)[1]>0){
    for(x in 1:dim(not_line)[1]){
      nlx <- not_line[x,]
      splt <- strsplit(unlist(nlx), ":|-")[[1]]
      seqrange <- seq.int(from = as.numeric(splt[2]), to = as.numeric(splt[3]), by = 1)
      rangeo <- c()
      for(xx in seq_along(seqrange)){
        rangeo <- c(rangeo, paste0(splt[1], ":", seqrange[xx], "-", seqrange[xx+1]))
      }
      rangeo <- rangeo[-length(rangeo)]
      join_chr_all_tb <- tibble::add_row(.data = join_chr_all_tb,
                                         index = rangeo,
                                         nlx[2],
                                         nlx[3],
                                         nlx[4],
                                         nlx[5])
    }
  }

  print("Joining into single GRanges...")
  join_chr_all_gr_tb <- tidyr::separate(data = join_chr_all_tb,
                                        col =  index,
                                        into = c("seqnames", "start", "end"),
                                        sep = "[:-]")

  join_chr_all_gr_tb <- dplyr::select(.data = join_chr_all_gr_tb,
                                      seqnames,
                                      "ranges" = start,
                                      "samples_n" = n.overlaps,
                                      "sampleIDs" = names)

  ##make seqinfo
  ##"'seqinfo' must be NULL, or a Seqinfo object, or a character vector of
  ## seqlevels, or a named numeric vector of sequence lengths"
  seqinf <- tryCatch(GenomeInfoDb::fetchExtendedChromInfoFromUCSC(which_genome),
              error = function(f){
                GenomeInfoDb::getChromInfoFromUCSC(which_genome)
              }
            )

  seqinf[,1] <- gsub("chr","",seqinf[,1])
  seqinf <- seqinf[grep("_", seqinf[,1], invert = TRUE),c(1,2)]
  seqinf <- GenomeInfoDb::Seqinfo(seqnames = seqinf[,1],
                                  seqlengths = seqinf[,2],
                                  genome = which_genome)

  join_chr_all_gr <- GenomicRanges::GRanges(seqnames = factor(gsub("chr", "", unlist(join_chr_all_gr_tb[,1]))),
                         ranges = IRanges::IRanges(start = as.numeric(unlist(join_chr_all_gr_tb[,2])),
                                                   end = as.numeric(unlist(join_chr_all_gr_tb[,2]))),
                         strand = NULL,
                         mcols = join_chr_all_gr_tb[,c(3,4)],
                         seqinfo = seqinf)

  colnames(S4Vectors::mcols(join_chr_all_gr)) <- gsub("mcols.", "", colnames(S4Vectors::mcols(join_chr_all_gr)))

  join_chr_all_gr <- unique(join_chr_all_gr)

  ##sample names
  samples <- unique(unlist(lapply(unique(join_chr_all_gr$sampleIDs), function(s){
      stringr::str_split(s, ",")[[1]]
    })))

  ##function to make master table mcols
  master_mcols <- function(gr_list, gr_master, ps_vec, dp_vec, seqinf){
    lapply(seq_along(gr_list), function(ff){
      ##first GRanges object
      gr_ff <- unique(gr_list[[ff]])
      GenomeInfoDb::seqinfo(gr_ff) <- seqinf

      ##mcols of those hits
      gr_ff_df <- S4Vectors::mcols(gr_ff[, c(ps_vec, dp_vec)])

      ##name based on list names and return
      names(gr_ff_df) <- c(paste0(names(gr_list)[ff], ".",
                                  gsub("mcols.", "", names(gr_ff_df[,ps_vec]))),
                           gsub("mcols.", "", names(gr_ff_df[,dp_vec])))

      ##make names into a col also
      if(!is.null(names(gr_ff))){
        gr_ff_df$rowname <- names(gr_ff)
      } else {
        gr_df <- as.data.frame(gr_ff)
        gr_ff_df$rowname <- names(gr_ff) <- apply(gr_df, 1, function(p){
          gsub(" ","",paste0(p[1], ":", p[2], "-", p[3]))
        })
      }
      ##create master output, same size as gr_master and with mcols from queryHits
      ##in the place where subjectHits matched
      df_master <- as.data.frame(matrix(nrow = length(gr_master),
                                        ncol = length(names(gr_ff_df))))

      colnames(df_master) <- names(gr_ff_df)

      ##where gr_ff intersects with the supplied ranges
      hits <- base::as.data.frame(GenomicRanges::findOverlaps(gr_ff, gr_master, ignore.strand = TRUE, minoverlap = 0))

      ##place 'hits' annotation as per subjectHits
      qh <- hits$queryHits

      df_master[hits[qh, 2],] <- unlist(gr_ff_df[hits[qh, 1],])

      return(df_master)
    })
  }

  ##using the above as a master GRanges object, walk through per sample
  ##create per sample df to be added to mcol
  print("Creating master GRanges...")
  S4Vectors::mcols(join_chr_all_gr) <- c(S4Vectors::mcols(join_chr_all_gr), do.call(cbind, master_mcols(gr_list, join_chr_all_gr, ps_vec, dp_vec, seqinf)))

  ##need to collapse the duplicated columns into one
  ##include a 'rowname' for unique naming

  ##is mcols prefix in mcols? remove unless also in dp_vec
  col_dp <- colnames(S4Vectors::mcols(join_chr_all_gr)) %in% dp_vec
  dp_vecr <- c(dp_vec, "rowname")
  if(length(table(col_dp))==1){
    dp_vecr <- gsub("mcols.", "", dp_vecr)
    col_dp <- colnames(S4Vectors::mcols(join_chr_all_gr)) %in% dp_vecr
  } else {
    col_dp <- colnames(S4Vectors::mcols(join_chr_all_gr)) %in% dp_vecr
  }

  ##take inverse to keep
  col_kp <- !col_dp

  ##get unique set of values for dp
  mcols_dpd <- as.data.frame(matrix(nrow = length(join_chr_all_gr),
                                    ncol = length(dp_vecr)))
  colnames(mcols_dpd) <- colnames(as.data.frame(S4Vectors::mcols(join_chr_all_gr)))[col_dp][1:(length(dp_vecr))]

  ##set up tibble to allow checking of NAs, to condense dp_vec
  dp_tb <- tibble::as_tibble(S4Vectors::mcols(join_chr_all_gr))[, col_dp]

  ##unique colnames
  col_uniq <- grep("\\.", colnames(dp_tb), invert = TRUE, value = TRUE)

  ##check these for NA
  dp_chk <- unlist(lapply(1:length(samples), function(s){
    if(s == 1){
      return(col_uniq[1])
    } else {
      return(paste0(col_uniq[1], ".", s - 1))
    }
  }))

  ##which of dp_chk are not NA (use first match for specifying mcols)
  not_isna <- !is.na(dp_tb[, dp_chk])
  dp_cd_tb <- dplyr::bind_rows(lapply(1:dim(dp_tb)[1], function(f){
     print(f)

    ff <- dp_tb[f, ]

    ##have come across two NAs which give both FALSE
    mtch_i <- match(dp_chk[match(TRUE, not_isna[f,])], colnames(dp_tb))
    mtch_tb <- dp_tb[f, mtch_i:(mtch_i+length(dp_vecr)-1)]
    colnames(mtch_tb) <- col_uniq
    return(mtch_tb)
  }))

  ##set as mcols and rename
  dp_cd_tb <- dplyr::select(.data = dp_cd_tb, rowname, dplyr::everything())

  S4Vectors::mcols(join_chr_all_gr) <- c(as.data.frame(dp_cd_tb), as.data.frame(S4Vectors::mcols(join_chr_all_gr)[col_kp]))
  names(join_chr_all_gr) <- join_chr_all_gr$rowname

  ##write output
  join_chr_kp_gr <- sort(join_chr_all_gr[join_chr_all_gr$samples_n > 1,])
  names(join_chr_kp_gr) <- join_chr_kp_gr$rowname
  adr <- as.data.frame(S4Vectors::mcols(join_chr_kp_gr))
  readr::write_tsv(adr, file = paste0(tag, ".master_consensus.tsv"))
  readr::write_tsv(as.data.frame(S4Vectors::mcols(join_chr_all_gr)),
                  file = paste0(tag, ".master_all.tsv"))
  gr_master_consensus_all <- list(join_chr_kp_gr, join_chr_all_gr)
  save(gr_master_consensus_all, file = paste0(tag, ".master_consensus_all.RData"))
  return(gr_master_consensus_all)
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
    a3 <- aa3[!is.na(unlist(stringi::stri_match_all(f, regex = aa3)))]
    a1 <- aa1[!is.na(unlist(stringi::stri_match_all(f, regex = aa3)))]
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
#' @param which_genome hg19 or hg38
#' @return GRanges object of all  of single-letter HGVS protein IDs
#' @export

gr_super_alt_plot <- function(var_list, name_callers, impacts, taga, included_order, which_genome) {

  ##GRanges superset
  gr_super <- somenone::gr_super_set(var_list, name_callers, impacts)

  ##get list to plot from with at least two callers supporting
  plot_list <- somenone::at_least_two(var_list, gr_super, taga)

  ##if single sample, make plot_list a GRanges object
  if(length(plot_list) == 1){

    if(length(plot_list[[1]]) == 0){

      print(paste0("No variants for IMPACTS: ", impacts))
    } else {
      print("Variants found, plotting...")
        somenone::plot_single(granges = plot_list[[1]], sampleID = names(plot_list), tag = paste0(names(plot_list), ".", taga, ".solo"))
    }
  } else {

    ##test for empty and rename to exclude those empty
    nz_plot_list <- lapply(seq_along(plot_list), function(f){
      if(length(plot_list[[f]]) != 0){
        return(sort(plot_list[[f]]))
      }
    })
    nm_vec <- unlist(lapply(seq_along(plot_list), function(f){
      if(length(plot_list[[f]])!=0){
        return(names(plot_list)[f])
    }}))

    ##remove NULL (where nothing was returned above)
    nz_plot_list[sapply(nz_plot_list, is.null)] <- NULL
    names(nz_plot_list) <- nm_vec

    ##create a GRanges of shared elements from list
    if(length(nz_plot_list)==0){
      print(paste0("No shared variants for IMPACTS: ", impacts, ", support across callers lacking"))

    } else {
      col_vec <- names(S4Vectors::mcols(nz_plot_list[[1]]))
      ps_vec <- col_vec[1:3]
      dp_vec <- col_vec[4:length(col_vec)]
      master_gr_list <- somenone::master_intersect_snv_grlist(gr_list = nz_plot_list,
                                                ps_vec = ps_vec,
                                                dp_vec = dp_vec,
                                                tag = taga,
                                                which_genome)

      ##based on elements in nz_plot_list, plot or do not
      ##shared
      if(length(master_gr_list[[1]]) == 0){
        print(paste0("No shared variants for IMPACTS: ", impacts, ", support across callers lacking"))
      } else {
        print("Shared variants found, plotting...")
          somenone::plot_consensus(master_gr = master_gr_list[[1]], tag = paste0(taga, ".shared"), included_order)
      }
      ##all
      if(length(master_gr_list[[2]]) == 0){
        print(paste0("No shared variants for IMPACTS: ", impacts, ", support across callers lacking"))
      } else {
        print("Shared variants found, plotting...")
          plot_consensus(master_gr = master_gr_list[[2]], tag = paste0(taga, ".all"),  included_order)
      }
    }
  }
  return(list(gr_super, plot_list))
}

#' Helper deals with setdiff, intersect to return full GRanges not just ranges
#'
#' @param granges_1 is a GRanges object
#' @param granges_2 is a GRanges object
#' @param method is 'setdiff' or 'intersect' to find difference of granges_1 vs. 2, or the intersect of the two
#' @return GRanges object based on method
#' @export

granges_sdin <- function(granges_1, granges_2, method){
  if(method == "setdiff"){
    gr_sd <- suppressWarnings(GenomicRanges::setdiff(BiocGenerics::sort(granges_1), BiocGenerics::sort(granges_2)))
    nms_12 <- paste0(GenomeInfoDb::seqnames(gr_sd), ":", BiocGenerics::start(IRanges::ranges(gr_sd)), "_")
    grp_12 <- unlist(lapply(nms_12, function(f){
                grep(f, names(granges_1), value=TRUE)
              }))
    return(granges_1[grp_12])
  } else {
    gr_sd <- suppressWarnings(GenomicRanges::intersect(BiocGenerics::sort(granges_1), BiocGenerics::sort(granges_2)))
    nms_12 <- paste0(GenomeInfoDb::seqnames(gr_sd), ":", BiocGenerics::start(IRanges::ranges(gr_sd)), "_")
    grp_12 <- unlist(lapply(nms_12, function(f){
                grep(f, names(granges_1), value=TRUE)
              }))
    return(granges_1[grp_12])
  }
}
