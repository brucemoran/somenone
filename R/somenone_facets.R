#' Facets functions
#'
#' Standard Facets CNA call
#'
#' @import pctGCdata
#' @param input_csv path to CSV format file as per requirements of Facets
#' @return none, saves plots and writes output in current dir
#' @export

facets_cna_call <- function(input_csv){

  attachNamespace("pctGCdata")
  tumour <- stringr::str_split(input_csv, "\\.")[[1]][1]
  set.seed(1234)
  snpmat <- facets::readSnpMatrix(input_csv)
  pps <- facets::preProcSample(snpmat)
  oo <- facets::procSample(pps)
  diplogr <- oo$dipLogR
  fit <- facets::emcncf(oo)
  ploidpurdf <- data.frame(PLOIDY = round(fit$ploidy, digits=3),
                           PURITY = round(fit$purity, digits=3))
  readr::write_tsv(ploidpurdf, path = paste0(tumour,".fit_ploidy_purity.tsv"))
  readr::write_tsv(round(fit$cncf, 3), path = paste0(tumour,".fit_cncf_jointsegs.tsv"))
  pcgr_out <- tibble::tibble(Chromosome = fit$cncf$chrom,
                             Start = as.integer(fit$cncf$start),
                             End = as.integer(fit$cncf$end),
                             Segment_Mean = fit$cncf$cnlr.median)
  readr::write_tsv(pcgr_out, path = paste0(tumour,".cncf_jointsegs.pcgr.tsv"))
  grDevices::pdf(paste0(tumour, ".facets_CNA.pdf"))
    facets::plotSample(x = oo, emfit = fit, plot.type = "both", sname = tumour)
  grDevices::dev.off()
}

#' Consensus CNA estimation and plotting with Facets input
#'
#' @import pctGCdata
#' @param pattern a pattern within file names to use for parsing input
#' @param dict_file a dictionary file for the fasta used to align data
#' @param tag a string used to tag output
#' @param cgc_bed bedfile from the Cancer Gene Census from COSMIC
#' @return none, saves plots and data in current dir
#' @export

facets_cna_consensus <- function(pattern, dict_file, tag, cgc_bed = NULL) {

  ##housekeeping
  options(scipen = 999)
  options(stringAsFactors = FALSE)

  ##set genome assembly version based on dictfile
  which_genome <- "hg19"
  if(length(grep("38", dict_file))==1) {
    which_genome <- "hg38"
  }

  ##load Cancer Gene Census bed file if supplied and run process
  cgc_gr <- NULL
  if(! is.null(cgc_bed)){
    print("Working on CGC")
    cgc <- utils::read.table(cgc_bed)
    cgc_gr <- GenomicRanges::GRanges(seqnames = cgc[,1],
                                    ranges = IRanges::IRanges(cgc[,2], cgc[,3]),
                                    CGC_gene = unlist(lapply(as.vector(cgc[,4]), function(f){
                                                        strsplit(f, ";")[[1]][1]
                                                      })))
    in_list <- as.list(dir(pattern = pattern))
    out_list <- lapply(in_list, function(f){
      somenone::process_in_list(f, which_genome, cgc_gr)
    })
    output_out_list(out_list, in_list, dict_file, which_genome, tag = paste0(tag, ".CGC"))
  }

  ##using ENS annotation otherwise/also
  ##set input list based on file pattern
  print("Working on ENS")
  in_list <- as.list(dir(pattern = pattern))
  out_list <- lapply(in_list, function(f){
    process_in_list(f, which_genome, cgc_gr = NULL)
  })

  output_out_list(out_list, in_list, dict_file, which_genome, tag = paste0(tag, ".ENS"))
}

#' Processing list of Facets input
#' @param in_list of data from Facets
#' @param which_genome genoem assembly being used "hg19" or "hg38"
#' @param cgc_gr GRanges of cancer gene census from COSMIC
#' @return purity-ploidy dataframe, also write rds used
#' @export

process_in_list <- function(in_list, which_genome, cgc_gr){

  ##allow CGC, else full Ensembl annotations
  anno <- "ENS"
  if(!is.null(cgc_gr)){
    anno <- "CGC"
  }

  ##iterate over samples
  for(samp in 1:length(in_list)){
    jointseg <- in_list[[samp]]
    sampleID <- out_name <- stringr::str_split(jointseg, "\\.")[[1]][1]
    out_ext <- gsub(".tsv", "", jointseg)

    ##also read and store ploidy:
    ploidypurity <- dir(pattern = paste0(sampleID, ".fit_ploidy_purity.tsv"))
    ploidy <- as.vector(utils::read.table(ploidypurity)[2,1])
    purity <- as.vector(utils::read.table(ploidypurity)[2,2])
    pp_df <- data.frame(PLOIDY = ploidy, PURITY = purity)

    print(paste0("Working on: ", sampleID))

    ##run function to make GRangesList
    ##(jointsegs_in, which_genome, cgc_gr = NULL, anno = NULL, bsgenome = NULL)
    grl <- somenone::facets_jointsegs_parse_to_gr(jointseg, sampleID, which_genome, anno = anno, cgc_gr = cgc_gr)
  }
  return(list(grl, pp_df))
}

#' Parse jointsegs into GRanges object
#' @param jointseg filename (output from Facets) with columns:
#'    "seg", "num.mark", "nhet", "cnlr.median", "mafR", "segclust", "cnlr.median.clust", "mafR.clust", "cf.em", "tcn.em", "lcn.em"
#' @param sampleID the sample identifier, taken as the string before first dot in jointseg
#' @param which_genome the genome assembly used, "hg19" or "hg38"
#' @param cgc_gr GRanges object of cancer gene census from COSMIC
#' @param anno annotation strategy, "ENS" or "CGC"
#' @param bsgenome which version of bsgenome to use
#' @return a GRanges object with annotated jointsegs from input
#' @export

facets_jointsegs_parse_to_gr <- function(jointseg, sampleID, which_genome, anno = NULL, cgc_gr = NULL, bsgenome = NULL){

  ##bed file from: https://cancer.sanger.ac.uk/census
  ##NB requires login hence not provided (but parameters exist to supply
  ##credentials for download-references.nf in somatic_n-of-1 workflow)
  if(!is.null(cgc_gr)){
    CGCGR <- cgc_gr
    anno <- "CGC"
  }

  ##annotation source
  if(is.null(anno)){
    anno <- "ENS"
  }

  ##default genome version is hg38
  if(which_genome == "hg19"){
    bsgenome <- BSgenome.Hsapiens.UCSC.hg19::BSgenome.Hsapiens.UCSC.hg19
  } else {
    bsgenome <- BSgenome.Hsapiens.UCSC.hg38::BSgenome.Hsapiens.UCSC.hg38
  }

  ##read input
  js <- utils::read.table(jointseg, header=T)

  ##convert chr23 -> X
  if("23" %in% js$chrom){
    js$chrom[js$chrom == 23] <- "X"
  }
  gr <- GenomicRanges::GRanges(seqnames = js$chrom,
                ranges = IRanges::IRanges(start = js$start, end = js$end),
                mcols = js[, c("seg", "num.mark", "nhet", "cnlr.median", "mafR", "segclust", "cnlr.median.clust", "mafR.clust", "cf.em", "tcn.em", "lcn.em")])
  GenomeInfoDb::seqlevels(gr) <- GenomeInfoDb::mapSeqlevels(GenomeInfoDb::seqlevels(gr), "UCSC")
  GenomeInfoDb::seqinfo(gr) <- GenomeInfoDb::seqinfo(bsgenome)[GenomeInfoDb::seqlevels(gr)]
  GenomeInfoDb::seqlevelsStyle(gr) <- "NCBI"

  ##overlaps
  if(anno == "ENS"){

    ##annotation
    if(which_genome == "hg19"){
      genes <- ensembldb::genes(EnsDb.Hsapiens.v75::EnsDb.Hsapiens.v75)
      gene_symbols <- c("GRCh37_v75_SYMBOL", "count_GRCh37_v75_SYMBOLs")
    } else {
      genes <- ensembldb::genes(EnsDb.Hsapiens.v86::EnsDb.Hsapiens.v86)
      gene_symbols <- c("GRCh38_v86_SYMBOL", "count_GRCh38_v86_SYMBOLs")
    }

    GenomeInfoDb::seqlevelsStyle(genes) <- "NCBI"
    GenomeInfoDb::genome(genes) <- which_genome
    GenomeInfoDb::seqlevels(genes, pruning.mode="coarse") <- GenomeInfoDb::seqlevels(gr)
    names(gr) <- gr$mcols.seg

    hits <- as.data.frame(GenomicRanges::findOverlaps(gr, genes, ignore.strand=TRUE))
    hits$SYMBOL <- biomaRt::select(org.Hs.eg.db::org.Hs.eg.db,
                                   as.character(genes[hits$subjectHits]$entrezid),
                                   "SYMBOL")$SYMBOL
    gr$SYMBOL <- "-"
    gr$SYMBOLs <- "-"
    ##loop to collapse symbols per region
    for(x in 1:max(hits$queryHits)){
      hitsx <- sort(unique(hits$SYMBOL[hits$queryHits==x]))
      gr$SYMBOL[x] <- paste(hitsx[2:length(hitsx)], collapse=";")
      gr$SYMBOLs[x] <- length(hitsx)-1;
    }

    ##rename S4Vectors::mcols
    names(S4Vectors::mcols(gr)) <- c(names(S4Vectors::mcols(gr))[1:9], "Total_Copy_Number", "Minor_Copy_Number", gene_symbols)
    readr::write_tsv(as.data.frame(gr), path = paste0(sampleID, ".facets.CNA.ENS.tsv"))
    saveRDS(gr, file = paste0(sampleID, ".facets.CNA.ENS.RData"))
  }

  if(anno == "CGC"){
    hits <- as.data.frame(GenomicRanges::findOverlaps(gr, CGCGR, ignore.strand=TRUE))
    hits$SYMBOL <- CGCGR[hits$subjectHits]$CGC_gene
    gr$CGC_SYMBOL <- "-"
    gr$CGC_SYMBOLs <- "-"
    ##loop to collapse symbols per region
    for(x in 1:max(hits$queryHits)){
      hitsx <- as.vector(sort(unique(hits$SYMBOL[hits$queryHits==x])))
      hitsx <- hitsx[!is.na(hitsx)]
      if(length(hitsx)==0){gr$CGC_SYMBOL[x] <- NA; gr$CGC_SYMBOLs[x] <- 0}
      else{
        gr$CGC_SYMBOL[x] <- paste(hitsx[2:length(hitsx)], collapse=";")
        gr$CGC_SYMBOLs[x] <- length(hitsx)-1;
      }
    }

    ##rename S4Vectors::mcols
    names(S4Vectors::mcols(gr)) <- c(names(S4Vectors::mcols(gr))[1:9], "Total_Copy_Number", "Minor_Copy_Number", "CGC_SYMBOL", "n_CGC_SYMBOLs")
    readr::write_tsv(as.data.frame(gr), path = paste0(sampleID, ".facets.CNA.CGC.tsv"))
    saveRDS(gr, file = paste0(sampleID, ".facets.CNA.CGC.RData"))
  }

  return(gr)
}


#' Plotting and Output from 'out_list' object
#' @param out_list resulting from process_in_list
#' @param in_list input to process_in_list
#' @param dict_file dictionary file of fasta used to align BAMs
#' @param which_genome the genome assembly used, "hg19" or "hg38"
#' @param tag runID and annotation strategy, "ENS" or "CGC"
#' @return none writes output files
#' @export

output_out_list <- function(out_list, in_list, dict_file, which_genome, tag){

  ##separate elements in list object into list objects
  facets_list <- lapply(out_list, function(f){
    ff <- f[[1]]
    names(ff) <- ff$mcols.seg
    return(ff)
  })
  pp_list <- lapply(out_list, function(f){
    return(f[[2]])
  })

  samples <- names(pp_list) <- names(facets_list) <- unlist(lapply(in_list, function(f){
    stringr::str_split(f,"\\.")[[1]][1]
  }))

  ##list of only those CNA (i.e. not TCN == 2)
  cna_list <- lapply(facets_list, function(f){
    fo <- f[S4Vectors::mcols(f)$Total_Copy_Number != 2]
    sort(unique(c(fo)))
  })

  ##make bed, bedr:: and GRanges from that
  ps_vec <- c("mcols.seg",
              "mcols.cnlr.median",
              "Total_Copy_Number",
              "Minor_Copy_Number",
              rev(colnames(S4Vectors::mcols(cna_list[[1]])))[2])

  cna_master_gr <- master_intersect_cna_grlist(cna_list, ps_vec)
  cna_master_df <- as.data.frame(cna_master_gr)
  cna_master_df$width <- cna_master_df$end - cna_master_df$start
  readr::write_tsv(cna_master_df, path = paste0(tag, ".facets.CNA.master.tsv"))

  ##make CNA df
  cna_df <- c()
  for(x in seq_along(samples)){
    if(length(cna_list[[x]])>0){
      cna_dfb <- as.data.frame(cna_list[[x]], row.names = seq(1:length(cna_list[[1]])))
      cna_dfb$purity <- unlist(rep(pp_list[[x]][2], length(cna_dfb[,1])))
      cna_dfb$ploidy <- unlist(rep(pp_list[[x]][1], length(cna_dfb[,1])))
      cna_dfb$sampleID <- samples[x]
      cna_df <- rbind(cna_df, cna_dfb)
      readr::write_tsv(cna_df, path = paste0(samples[x],".facets.CNA.jointsegs.tsv"))
    } else {
      readr::write_tsv(as.data.frame(cna_list[[x]]), path = paste0(samples[x], ".facets.CNA.jointsegs.tsv"))
    }
  }

  ##remove character chromosomes NB no MT, GL...
  sqn <- gsub("chr", "", as.vector(cna_df[,1]))
  sqn[sqn=="X"] <- 23
  sqn[sqn=="Y"] <- 24
  cna_df$seqnames <- sqn

  seqlengths <- seqlengths_df(unique(cna_df[,1]), dict_file, which_genome)
  cna_df <- cna_df[cna_df$seqnames %in% seqlengths$seqnames,]
  cum_sum_add <- unlist(apply(cna_df, 1, function(x) {
              cum_sum_off <- seqlengths$cum_sum_0[rownames(seqlengths) %in% x[1]]
              return(cum_sum_off)
  }))

  ##plotting
  plot_start <- plot_end <- cna_call <- minor_call <- purity <- ploidy <- diploid <- colour <- sample <- seqlength <- cum_sum_0 <- cum_sum_centro <- NULL

  ##set maximum plotted CNA
  cna_maxd <- cna_df$Total_Copy_Number
  if(max(cna_maxd) > 8){
    cna_maxd[cna_maxd > 8] <- 8
  }

  ##colours for plotting
  colz <- c("blue", "darkblue", "black", "darkred", "firebrick3", "red2", rep("red",199))
  names(colz) <- c(seq(from=0,to=198,by=1))
  cnames <- sort(unique(cna_maxd))
  colz <- colz[is.element(names(colz), cnames)]

  ##plot data.frame
  plot_df <- data.frame(row.names=seq(from = 1,to = dim(cna_df)[1], by = 1),
                        seqnames = cna_df[,1],
                        plot_start = cna_df[,2] + (cum_sum_add - 1),
                        plot_end = (cna_df[,3] - 1) + cum_sum_add,
                        cna_call = as.numeric(cna_maxd),
                        minor_call = as.numeric(unlist(cna_df$Minor_Copy_Number)),
                        purity = round(as.numeric(cna_df$purity)),
                        ploidy = round(as.numeric(cna_df$ploidy)),
                        diploid = 2,
                        colour = plyr::mapvalues(cna_maxd, cnames, colz),
                        sample = unlist(cna_df$sampleID))

  ##plot
  ggp <-  ggplot2::ggplot() +
          ggplot2::geom_hline(data = plot_df,
                ggplot2::aes(yintercept = ploidy),
                    linetype = 1,
                    color = "purple",
                    size = 0.5) +
          ggplot2::geom_hline(data = plot_df,
                ggplot2::aes(yintercept = diploid),
                    linetype = 1,
                    color = "orange",
                    size = 0.5) +
          ggplot2::geom_vline(data = seqlengths,
                     ggplot2::aes(xintercept = cum_sum_0),
                     linetype = 1,
                     color = "dodgerblue",
                     size = 0.2) +
          ggplot2::geom_vline(data = seqlengths,
                     ggplot2::aes(xintercept = cum_sum_centro),
                     linetype = 4,
                     color = "grey",
                     size = 0.2) +
          ggplot2::geom_segment(data = plot_df,
                    ggplot2::aes(x = plot_start,
                     xend = plot_end,
                     y = cna_call,
                     yend = cna_call,
                     colour = colour,
                     size = 4.5)) +
          ggplot2::scale_colour_identity() +
          ggplot2::geom_segment(data = plot_df,
                    ggplot2::aes(x = plot_start,
                     xend = plot_end,
                     y = minor_call,
                     yend = minor_call,
                     colour = "forestgreen",
                     size = 2)) +
          ggplot2::scale_x_continuous(name = "Chromosome",
                     labels = as.vector(seqlengths$seqnames),
                     breaks = seqlength$cum_sum_centro) +
          ggplot2::scale_y_continuous(name = "Total CNA (facets tcn.em)",
                     labels = seq(from = 0, to = max(cna_maxd), by = 1),
                     breaks = seq(from = 0, to = max(cna_maxd), by = 1),
                     limits = c(0, max(cna_maxd))) +
          ggplot2::theme(axis.text.x = ggplot2::element_text(size = 5),
                     panel.grid.major = ggplot2::element_blank(),
                     panel.grid.minor = ggplot2::element_blank(),
                     legend.position="none") +
          ggplot2::facet_grid(sample ~ .)

  ggplot2::ggsave(filename = paste0(tag, ".facets_consensus.plot.pdf"), plot = ggp)
}

#' Create lengths of chromosomes for plotting
#' @param in_seqs a jointsegs filename (output from Facets) with columns:
#'    "seg", "num.mark", "nhet", "cnlr.median", "mafR", "segclust", "cnlr.median.clust", "mafR.clust", "cf.em", "tcn.em", "lcn.em"
#' @param dict_file dictionary file for fasta used in alignment
#' @param which_genome the genome assembly used, "hg19" or "hg38"
#' @return a GRanges object with annotated jointsegs from input
#' @export

seqlengths_df <- function(in_seqs, dict_file, which_genome){

  V1 <- V3 <- V4 <- V5 <- seqnames <- NULL

  ##create centromere data
  url19 <- paste0("http://hgdownload.cse.ucsc.edu/goldenpath/hg19/database/cytoBand.txt.gz")
  url38 <- paste0("http://hgdownload.cse.ucsc.edu/goldenpath/hg38/database/cytoBand.txt.gz")
  temp19 <- tempfile()
  utils::download.file(url19, temp19)
  centroms19 <- tibble::as_tibble(utils::read.table(temp19))
  temp38 <- tempfile()
  utils::download.file(url38, temp38)
  centroms38 <- tibble::as_tibble(utils::read.table(temp38, fill=TRUE))

  centromeres <- function(chr, which_genome){
    if(which_genome == "hg19"){
      centrom <- centroms19
    } else {
      centrom <- centroms38
    }
    centrout <- dplyr::mutate(centrom, seqnames = gsub("chr", "", V1))
    centrout <- dplyr::select(centrout, seqnames, -V1, dplyr::everything())
    centrout <- dplyr::filter(centrout, V5 %in% "acen")
    centrout <- dplyr::filter(centrout, V4 %in% grep("p11", V4, value=T))
    centrout <- dplyr::filter(centrout, seqnames %in% chr)
    centrout <- dplyr::select(centrout, seqnames, centromere = V3)
    return(centrout)
  }

  ##parse dict
  dict_seqs <- utils::read.table(dict_file, skip=1)
  seqlens <- data.frame(seqnames = gsub("SN:","", dict_seqs[,2]),
                        start = 1,
                        end = as.numeric(gsub("LN:","", dict_seqs[,3])),
                        stringsAsFactors = FALSE)

  ##find those chromosomes in inSeqs
  seqlens <- seqlens[is.element(seqlens[,1],
                     unique(as.vector(in_seqs))),]

  ## make centromere data from function
  ##load centromere
  centromere <- centromeres(seqlens[,1], which_genome)
  seqlengths <- dplyr::left_join(seqlens, centromere, by = "seqnames")

  ##non-numeric chr IDs are numeric as rownames!
  ##need to output a table to convert between inSeqs and newSeqs
  seqlengths$cum_sum_0 <- c(1,cumsum(as.numeric(seqlengths$end))[1:length(seqlengths[,1])-1])
  seqlengths$cum_sum_1 <- c(cumsum(as.numeric(seqlengths$end)))
  seqlengths$cum_sum_centro <- seqlengths$centromere + seqlengths$cum_sum_0

  return(seqlengths)
}

#' Processing list of Facets input
#' @param cna_df data frame created in facets_cna_consensus()
#' @param cum_sum_add data_frame created in seqlengths_df()
#' @param seqlength a data.frame made by seqlengths_df
#' @param cna_max max CNA call to plot (eveything else is squashed to this
#' @return purity-ploidy dataframe, also write rds outputs
#' @export

plot_cna_df <- function(cna_df, cum_sum_add, seqlength, cna_max = 8){

  plot_start <- plot_end <- cna_call <- minor_call <- purity <- ploidy <- diploid <- colour <- sample <- seqlength <- cum_sum_0 <- cum_sum_centro <- NULL

  ##set maximum plotted CNA
  cna_maxd <- cna_df$Total_Copy_Number
  if(max(cna_maxd) > cna_max){
    cna_maxd[cna_maxd > cna_max] <- cna_max
  }

  ##colours for plotting
  colz <- c("blue", "darkblue", "black", "darkred", "firebrick3", "red2", rep("red",199))
  names(colz) <- c(seq(from=0,to=198,by=1))
  cnames <- sort(unique(cna_maxd))
  colz <- colz[is.element(names(colz), cnames)]

  ##plot data.frame
  plot_df <- data.frame(row.names=seq(from = 1,to = dim(cna_df)[1], by = 1),
                        seqnames = cna_df[,1],
                        plot_start = cna_df[,2] + (cum_sum_add - 1),
                        plot_end = (cna_df[,3] - 1) + cum_sum_add,
                        cna_call = as.numeric(cna_maxd),
                        minor_call = as.numeric(unlist(cna_df$Minor_Copy_Number)),
                        purity = round(as.numeric(cna_df$purity)),
                        ploidy = round(as.numeric(cna_df$ploidy)),
                        diploid = 2,
                        colour = plyr::mapvalues(cna_maxd, cnames, colz),
                        sample = unlist(cna_df$sampleID))

    ##plot
    ggp <-  ggplot2::ggplot() +
            ggplot2::geom_hline(data = plot_df,
                  ggplot2::aes(yintercept = ploidy),
                      linetype = 1,
                      color = "purple",
                      size = 0.5) +
            ggplot2::geom_hline(data = plot_df,
                  ggplot2::aes(yintercept = diploid),
                      linetype = 1,
                      color = "orange",
                      size = 0.5) +
            ggplot2::geom_vline(data = seqlength,
                       ggplot2::aes(xintercept = cum_sum_0),
                       linetype = 1,
                       color = "dodgerblue",
                       size = 0.2) +
            ggplot2::geom_vline(data = seqlength,
                       ggplot2::aes(xintercept = cum_sum_centro),
                       linetype = 4,
                       color = "grey",
                       size = 0.2) +
            ggplot2::geom_segment(data = plot_df,
                      ggplot2::aes(x = plot_start,
                       xend = plot_end,
                       y = cna_call,
                       yend = cna_call,
                       colour = colour,
                       size = 4.5)) +
            ggplot2::scale_colour_identity() +
            ggplot2::geom_segment(data = plot_df,
                      ggplot2::aes(x = plot_start,
                       xend = plot_end,
                       y = minor_call,
                       yend = minor_call,
                       colour = "forestgreen",
                       size = 2)) +
            ggplot2::scale_x_continuous(name = "Chromosome",
                       labels = as.vector(seqlength$seqnames),
                       breaks = seqlength$cum_sum_centro) +
            ggplot2::scale_y_continuous(name = "Total CNA (facets tcn.em)",
                       labels = seq(from = 0, to = max(cna_maxd), by = 1),
                       breaks = seq(from = 0, to = max(cna_maxd), by = 1),
                       limits = c(0, max(cna_maxd))) +
            ggplot2::theme(axis.text.x = ggplot2::element_text(size = 5),
                       panel.grid.major = ggplot2::element_blank(),
                       panel.grid.minor = ggplot2::element_blank(),
                       legend.position="none") +
            ggplot2::facet_grid(sample ~ .)
    return(ggp)
}

#' Find intersect from list of CNA GRanges
#' @param cna_list list of named GRanges objects
#' @param ps_vec mcols columns to keep p(er) s(ample; appended with list element's name)
#' @return GRanges object of all CNA, intersecting or not
#'         includes the ranges and each mcol from ps_vec, annotated `ps_vec[x].sample`
#' @export

master_intersect_cna_grlist <- function(gr_list, ps_vec){

  ##check gr is A GRanges object
  if(!as.vector(class(gr_list)) %in% c("GRangesList", "list")){
    print("Input \'gr_list\' is not a GRangesList nor list object, retry")
    break
  }

  ##use bedr::bedr.join.multiple.region to make a complete set and samples
  ##lapply over gr_list which creates character vector list in "chrX:1-100" format
  ##thinks it's 0-based so add 1...
  chr_list <- lapply(gr_list, function(f){
      fb <- apply(as.data.frame(f), 1, function(ff){
          ff <- unlist(ff)
          paste0("chr", ff[1], ":", gsub(" ", "", ff[2]), "-", gsub(" ", "", ff[3]))
        })
      bedr::bedr.sort.region(fb)
    })
  join_chr_all_tb <- tibble::as_tibble(bedr::bedr.join.multiple.region(chr_list))

  print("Joining into single GRanges...")
  join_chr_all_gr_tb <- tidyr::separate(data = join_chr_all_tb,
                                        col =  index,
                                        into = c("seqnames", "start", "end"),
                                        sep = "[:-]")
  join_chr_all_gr_tb <- dplyr::select(.data = join_chr_all_gr_tb,
                                      seqnames,
                                      start, end,
                                      "samples_n" = n.overlaps,
                                      "sampleIDs" = names)
  join_chr_all_gr <- GenomicRanges::GRanges(seqnames = gsub("chr", "", unlist(join_chr_all_gr_tb[,1])),
                     ranges = IRanges::IRanges(start = as.numeric(unlist(join_chr_all_gr_tb[,2])),
                                               end = as.numeric(unlist(join_chr_all_gr_tb[,3]))),
                     strand = NULL,
                     mcols = join_chr_all_gr_tb[,c("samples_n", "sampleIDs")],
                     seqinfo = GenomicRanges::seqinfo(gr_list[[1]]))

  adr <- as.data.frame(join_chr_all_gr)
  GenomicRanges::width(join_chr_all_gr) <- adr$end - adr$start

  colnames(S4Vectors::mcols(join_chr_all_gr)) <- gsub("mcols.", "", colnames(S4Vectors::mcols(join_chr_all_gr)))

  ##sample names
  samples <- unique(unlist(lapply(unique(join_chr_all_gr$sampleIDs), function(s){
      stringr::str_split(s, ",")[[1]]
    })))

  ##function to make master table mcols
  master_mcols <- function(gr_list, gr_master, ps_vec){
    lapply(seq_along(gr_list), function(ff){
      ##first GRanges object
      gr_ff <- gr_list[[ff]]

      ##where gr_ff intersects with the supplied ranges
      hits <- base::as.data.frame(GenomicRanges::findOverlaps(gr_ff, gr_master, ignore.strand = TRUE, minoverlap = 2))

      ##mcols of those hits
      mchits <- S4Vectors::mcols(gr_ff[, ps_vec][hits$queryHits])

      ##name based on list names and return
      names(mchits) <- paste0(names(gr_list)[ff], ".", gsub("mcols.", "", names(mchits)))

      ##create master output, same size as gr_master and with mcols from queryHits
      ##in the place where subjectHits matched
      master_df <- as.data.frame(matrix(nrow = length(gr_master),
                                        ncol = length(names(mchits))))
      colnames(master_df) <- names(mchits)

      ##place mchits annotation as per subjectHits
      for(x in 1:length(gr_master)){
        sh <- hits$subjectHits
        if(x %in% sh){
          master_df[hits[sh %in% x, 2],] <- unlist(mchits[hits[sh %in% x, 1],])
        } else {
          next
        }
      }
      return(master_df)
    })
  }

  ##using the above as a master GRanges object, walk through per sample
  ##create per sample df to be added to mcol
  S4Vectors::mcols(join_chr_all_gr) <- c(S4Vectors::mcols(join_chr_all_gr), do.call(cbind, master_mcols(gr_list, join_chr_all_gr, ps_vec)))

  ##return GRanges with mcols per sample
  return(join_chr_all_gr)
}
