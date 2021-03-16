#' Write output for pairtree based on mutations in master TSV input
#' Format:
#' id,	name,	var_reads,	total_reads,	var_read_prob
#' id is sequential s0..sn
#' name is ID for the variant (chr:pos_ref>alt_HGVS)
#' the below are comma-sep
#' var_reads is reads supporting variant at variant locus
#' total_reads is all reads at variant locus
#' var_read_prob is based on CNA at locus; mutation at diploid locus is 0.5;
#' in a duplicated (n=3) locus, 1/3 etc.
#'
#' @param rdata_input is RData from process vcfGra including master_all, and
#' gr_master_consensus_all
#' @param cn_master master CNA from fctscon process in somatic_n-of-1 Nextflow pipeline)
#' @param cn_pattern pattern matching FACETS input (fit_cncf_jointsegs.tsv)
#' @param pp_pattern pattern matching FACETS polidy/purity (fit_ploidy_purity.tsv)
#' @param which_genome is GRCh37 or GRCh38
#' @param tag is string to tag output
#' @return none, writes a tsv for input to pairtree
#' @export

make_pairtree_input <- function(rdata_input, cn_master, cn_pattern, pp_pattern, which_genome, tag){

  ##genome
  if(which_genome %in% c("GRCh37", "hg19")){
    which_genome <- "hg19"
  } else {
    which_genome <- "hg38"
  }

  ##make tibble of Mutation format
  load(rdata_input)
  mut_df <- as.data.frame(gr_master_consensus_all[[2]])
  sampleIDs <- sort(unique(unlist(strsplit(mut_df$sampleIDs, "\\,"))))

  ##filter on all samples
  mut_df_filt <- mut_df[mut_df$samples_n == length(sampleIDs),]

  ##required columns
  mut_df_filt$id <- paste0("s", 0:(dim(mut_df_filt)[1]-1))
  mut_df_filt$name <- paste0("S_", 0:(dim(mut_df_filt)[1]-1))
  mut_df_filt$mut_name <- rownames(mut_df_filt)
  mut_df_filt$var_reads <- apply(dplyr::select(.data = mut_df_filt,
                                         tidyselect::ends_with("AD.1")),
                                 1,
                                 function(f){
                                    paste(f, collapse = ",")
                                  })
  mut_df_filt$total_reads <- apply(dplyr::select(.data = mut_df_filt,
                                         tidyselect::ends_with("AD")),
                                   1,
                                   function(f){
                                      paste(f, collapse = ",")
                                    })
  for(x in 1:length(mut_df_filt$total_reads)){
    trds <- as.numeric(strsplit(mut_df_filt$total_reads[x], ",")[[1]])
    vrds <- as.numeric(strsplit(mut_df_filt$var_reads[x], ",")[[1]])
    mut_df_filt$total_reads[x] <- paste(trds + vrds, collapse = ",")
  }

  ##parse out CN files and makke correct format
  cn_df <- cna_master_muts(cn_master, mut_df_filt)
  mut_df_filt$var_read_prob <- cn_df

  ##write outputs
  psm <- data.frame(id = mut_df_filt$id,
                    name = mut_df_filt$name,
                    mut_name = mut_df_filt$mut_name,
                    var_reads = mut_df_filt$var_reads,
                    total_reads = mut_df_filt$total_reads,
                    var_read_prob = mut_df_filt$var_read_prob)
  readr::write_tsv(psm, file = paste0(tag, ".pairtree.psm"))

  mut_gr <- gr_master_consensus_all[[2]][names(gr_master_consensus_all[[2]]) %in% mut_df_filt$mut_name]

  mut_cn_to_pyclone(mut_gr = mut_gr,
                    cn_pattern = cn_pattern,
                    pp_pattern = pp_pattern,
                    tag)

  ##blank JSON output for clustervars
  json_list <- list(samples = sampleIDs,
                    clusters = list(),
                    garbage = list())
  pt_json <- rjson::toJSON(json_list)
  write(pt_json, file = paste0(tag, ".in_params_fn.json"))
}

#' Parse CN file and output as list
#'
#' @param cn_master is a pattern that matches CN input filenames
#' e.g. *.ENS.facets.CNA.master.tsv from fctcons
#' @param mut_df_filt object with mut_df_filt from above
#' @return cn values overlapping mut_df_filt loci
#' @export

cna_master_muts <- function(cn_master, mut_df_filt){

  ##files are in current working dir
  cna_master <- readr::read_tsv(cn_master)
  sampleIDs <- sort(unique(unlist(strsplit(mut_df_filt$sampleIDs, "\\,"))))
  tcns <- dplyr::select(.data = cna_master, 1, 2, 3, "samples_n",
                        tidyselect::ends_with("Total_Copy_Number"))
  tcns <- tcns[tcns$samples_n == length(sampleIDs),]

  ##iterate over muts, if they fall in a cna change prob
  mut_df_filt_cn <- apply(mut_df_filt, 1, function(f){
      var_read_prob <- paste(rep("0.5", length(sampleIDs)), collapse = ",")
      is_tcn <- dplyr::filter(.data = tcns,
                              seqnames %in% as.character(f["seqnames"]),
                              start < as.numeric(f["start"]),
                              end > as.numeric(f["end"]))
      if(dim(is_tcn)[1] > 0){
        var_read_prob <- paste(1/is_tcn[c(5,6,7,8)], collapse = ",")
      }
      return(var_read_prob)
  })

  return(mut_df_filt_cn)
}

#' Parse mutation and CN input, creating a master table of PyClone(VI) input
#' Format:
#' mutation_id, sample_id, ref_counts, alt_counts, normal_cn, major_cn, minor_cn, tumour_content
#' mutation_id is a unique ID (position-based)
#' sample_id is sample
#' ref/alt_counts are read counts at Mutation, VAF the variant allele frequency
#' mut_rdata is master_consensus_shared from somatic_n-of-1 workflow
#' For somatic_n-of-1 workflow, CN is based on FACETS inputs
#'
#' @param mut_gr gr_master_consensus_all
#' @param cn_pattern pattern matching FACETS input (fit_cncf_jointsegs.tsv)
#' @param pp_pattern pattern matching FACETS polidy/purity (fit_ploidy_purity.tsv)
#' @return none, writes a tsv for input to pyclone(VI)
#' @export

mut_cn_to_pyclone <- function(mut_gr, cn_pattern, pp_pattern, tag){

  ##make tibble of Mutation format
  mcols_df <- as.data.frame(S4Vectors::mcols(mut_gr))
  mut_df <- dplyr::select(.data = mcols_df, dplyr::contains(".AD"))
  ren_df <- gsub("AD", "ref_counts", gsub("AD.1", "alt_counts", colnames(mut_df)))
  colnames(mut_df) <- ren_df
  mut_gr1 <- GenomicRanges::GRanges(mut_gr,
                                   mcols = tibble::as_tibble(mut_df, rownames = "mutation_id"))
  colnames(S4Vectors::mcols(mut_gr)) <- gsub("mcols.", "", colnames(S4Vectors::mcols(mut_gr)))

  ##make df list of CN with Mutations
  mut_cn_df_list <- parse_cn_list(mut_gr = mut_gr1,
                                  cn_pattern = cn_pattern,
                                  pp_pattern = pp_pattern)

  ##reduce into a single df
  mut_cn_df <- mut_cn_df_list[[1]]
  for (x in 2:length(mut_cn_df_list)){
    mut_cn_df <- rbind(mut_cn_df, mut_cn_df_list[[x]])
  }
  readr::write_tsv(mut_cn_df, file = paste0(tag, ".mutation_CN.pyclone.tsv"))
}

#' Parse CN files and output as list
#'
#' @param mut_gr is the GRanges of Mutation data
#' @param cn_pattern is a pattern that matches CN input filenames
#' @param pp_pattern is a pattern that matches polidy/purity table from FACETS
#' @return list of correct format per sample
#' @export

parse_cn_list <- function(mut_gr, cn_pattern, pp_pattern){

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

  ##put cns in a list
  cn_list <- lapply(cn_files, function(f){

    ##read CN data, make GRanges
    r1 <- readr::read_tsv(f)
    gr_cn <- GenomicRanges::GRanges(seqnames = r1$chrom,
                           ranges = IRanges::IRanges(start = r1$start, end = r1$end),
                           mcols = r1[,c("tcn.em", "lcn.em")])
    sampleID <- strsplit(f, "\\.")[[1]][1]
    print(paste0("Working on: ", sampleID))

    colnames(S4Vectors::mcols(gr_cn)) <- c("major_cn", "minor_cn")

    ##compare mutation and cn
    mu_cn_overlap <- suppressWarnings(as.data.frame(GenomicRanges::findOverlaps(mut_gr, gr_cn)))
    if("rowname" %in% colnames(S4Vectors::mcols(mut_gr))){
      colnames(S4Vectors::mcols(mut_gr))[colnames(S4Vectors::mcols(mut_gr)) == "rowname"] <- "mutation_id"
    }

    order_cb <- c("mutation_id", "sample_id", "ref_counts", "alt_counts", "normal_cn", "major_cn", "minor_cn", "tumour_content")

    ##read pp data
    purity <- readr::read_tsv(pp_files[sampleID])[[2]]

    mu_cn_df <- do.call(rbind, lapply(1:length(mut_gr), function(m){

      ##is mutation found in mu_cn overlap?
      matcm <- match(m, mu_cn_overlap[,1])
      if(!is.na(matcm)){

        ##if it's found, add it to the output
        cn <- mu_cn_overlap[mu_cn_overlap[,1] %in% m, 2]
        mcolm <- S4Vectors::mcols(mut_gr[m])
        colnames(mcolm) <- gsub("mcols.", "", colnames(mcolm))

        cb <- cbind(as.data.frame(mcolm[c("mutation_id", grep(sampleID, colnames(mcolm), value = TRUE))]),
                    as.data.frame(gr_cn[cn])[c(6, 7)])
      } else {

        ##mut not overlapped by cn, return 2, 1 for CN major, minor

        fn <- rev(colnames(as.data.frame(gr_cn)))[c(2,1)]
        mcolm <- S4Vectors::mcols(mut_gr[m])
        colnames(mcolm) <- gsub("mcols.", "", colnames(mcolm))
        cb <- cbind(as.data.frame(mcolm[c("mutation_id", grep(sampleID, colnames(mcolm), value = TRUE))]),
                    data.frame(major_cn = 2, minor_cn = 1))
      }

      colnames(cb)[c(2, 3)] <- c("ref_counts", "alt_counts")
      cb$sample_id <- sampleID
      cb$normal_cn <- 2
      cb$tumour_content <- purity
      if(is.na(purity)){
        cb$tumour_content <- 1
      }

      return(cb[order_cb])
    }))
    mu_cn_df[is.na(mu_cn_df)] <- 0

    return(mu_cn_df)
  })
  names(cn_list) <- cn_names
  return(cn_list)
}


#' Write JSON for pairtree
#'
#' @param pyclone_res results from pyclone-vi write-results-file
#' @param pairtree_psm pairtree psm from make_pairtree_input(), with mut_name col
#' @param tag is string to tag output
#' @return none, writes JSON for pairtree
#' @export

make_pairtree_json <- function(pyclone_res, pairtree_psm, tag){

  ##read pyclone results and SSM for pairtree
  pc_res <- readr::read_tsv(pyclone_res)
  psm <- readr::read_tsv(pairtree_psm, col_types = "cccccc")

  ##get useful infos
  samples <- sort(unique(pc_res$sample_id))
  clusters_1 <- dplyr::filter(.data = pc_res,
                             sample_id %in% samples[1])

  ##keep only diploid a la https://github.com/morrislab/pairtree/issues/3
  garbage <- psm$id[psm$var_read_prob != "0.5,0.5,0.5,0.5"]

  ##make a list of mutation ids from SSM from pyclone res
  cluster_list <- lapply(unique(clusters_1$cluster_id), function(f){
      dfc <- dplyr::filter(.data = clusters_1,
                           cluster_id %in% f)
      psm_dfc <- dplyr::inner_join(psm, dfc, by = c("mut_name" = "mutation_id"))
      return(unlist(psm_dfc$id[!psm_dfc$id %in% garbage]))
    })


  ##JSON output
  json_list <- list(samples = samples,
                    clusters = cluster_list,
                    garbage = garbage)
  pt_json <- rjson::toJSON(json_list)
  write(pt_json, file = paste0(tag, ".pairtree.json"))
  in_json_list <- list(samples = samples,
                      clusters = list(),
                      garbage = list())
  pt_json <- rjson::toJSON(in_json_list)
  write(pt_json, file = "in_params_fn.json")

  ssm <- data.frame(id = psm$id,
                    name = psm$name,
                    var_reads = psm$var_reads,
                    total_reads = psm$total_reads,
                    var_read_prob = psm$var_read_prob)
  readr::write_tsv(ssm, file = paste0(tag, ".pairtree.ssm"))

}
