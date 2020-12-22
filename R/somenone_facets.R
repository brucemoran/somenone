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
  readr::write_tsv(ploidpurdf, file = paste0(tumour,".fit_ploidy_purity.tsv"))
  readr::write_tsv(round(fit$cncf, 3), file = paste0(tumour,".fit_cncf_jointsegs.tsv"))
  pcgr_out <- tibble::tibble(Chromosome = fit$cncf$chrom,
                             Start = as.integer(fit$cncf$start),
                             End = as.integer(fit$cncf$end),
                             Segment_Mean = fit$cncf$cnlr.median)
  readr::write_tsv(pcgr_out, file = paste0(tumour,".cncf_jointsegs.pcgr.tsv"))
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
      process_in_list(f, which_genome, cgc_gr)
    })
    output_out_list(out_list, in_list, dict_file, which_genome, tag = paste0(tag, ".CGC"), cgc_gr = cgc_gr)
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
    grl <- facets_jointsegs_parse_to_gr(jointseg, sampleID, which_genome, anno = anno, cgc_gr = cgc_gr)
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
    anno <- "CGC"
  }

  ##annotation source
  if(is.null(anno)){
    anno <- "ENS"
  }

  ##default genome version is hg38
  if(which_genome == "hg19"){
    bsgenome <- BSgenome.Hsapiens.UCSC.hg19::BSgenome.Hsapiens.UCSC.hg19
    gene_symbols <- c("GRCh37_v75_SYMBOL", "count_GRCh37_v75_SYMBOLs")
  } else {
    bsgenome <- BSgenome.Hsapiens.UCSC.hg38::BSgenome.Hsapiens.UCSC.hg38
    gene_symbols <- c("GRCh38_v86_SYMBOL", "count_GRCh38_v86_SYMBOLs")
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
    print("Annotating Ensembl genes")
    gr_anno <- anno_ens_cna(gr, which_genome)
    ##rename S4Vectors::mcols
    names(S4Vectors::mcols(gr_anno)) <- c(names(S4Vectors::mcols(gr_anno))[1:9], "Total_Copy_Number", "Minor_Copy_Number", gene_symbols)
    readr::write_tsv(as.data.frame(gr_anno), file = paste0(sampleID, ".facets.CNA.ENS.tsv"))
    saveRDS(gr_anno, file = paste0(sampleID, ".facets.CNA.ENS.RData"))
  }

  if(anno == "CGC"){
    print("Annotating Cancer Gene Census genes")
    gr_anno <- anno_cgc_cna(gr, cgc_gr, which_genome)

    ##rename S4Vectors::mcols
    names(S4Vectors::mcols(gr_anno)) <- c(names(S4Vectors::mcols(gr_anno))[1:9], "Total_Copy_Number", "Minor_Copy_Number", "CGC_SYMBOL", "n_CGC_SYMBOLs")
    readr::write_tsv(as.data.frame(gr_anno), file = paste0(sampleID, ".facets.CNA.CGC.tsv"))
    saveRDS(gr_anno, file = paste0(sampleID, ".facets.CNA.CGC.RData"))
  }

  return(gr_anno)
}

#' Annotate CNA with ENSembl IDs
#' @param gr GRanges object
#' @param which_genome the genome assembly used, "hg19" or "hg38"
#' @return GRanges object with annotation
#' @export

anno_ens_cna <- function(gr, which_genome){
    ##annotation
    if(which_genome == "hg19"){
      genes <- ensembldb::genes(EnsDb.Hsapiens.v75::EnsDb.Hsapiens.v75)
      gene_symbols <- c("GRCh37_v75_SYMBOL", "count_GRCh37_v75_SYMBOLs")
    } else {
      genes <- ensembldb::genes(EnsDb.Hsapiens.v86::EnsDb.Hsapiens.v86)
      gene_symbols <- c("GRCh38_v86_SYMBOL", "count_GRCh38_v86_SYMBOLs")
    }

    ##used
    genes <- genes[GenomeInfoDb::seqnames(GenomeInfoDb::seqinfo(genes)) %in% GenomeInfoDb::seqnames(GenomeInfoDb::seqinfo(gr)),]
    GenomeInfoDb::seqlevels(genes, pruning.mode="coarse") <- GenomeInfoDb::seqlevels(gr)
    GenomeInfoDb::seqinfo(genes) <- GenomeInfoDb::seqinfo(gr)
    hits <- as.data.frame(GenomicRanges::findOverlaps(gr, genes, ignore.strand=TRUE))

    hits$SYMBOL <- biomaRt::select(org.Hs.eg.db::org.Hs.eg.db,
                                   as.character(genes[hits$subjectHits]$entrezid),
                                   "SYMBOL")$SYMBOL
    gr$SYMBOL <- "-"
    gr$SYMBOLs <- "-"

    ##loop to collapse symbols per region
    print("collapsing gene symbol annotations per region")
    for(x in 1:max(hits$queryHits)){
      hitsx <- sort(unique(hits$SYMBOL[hits$queryHits==x]))
      gr$SYMBOL[x] <- paste(hitsx[2:length(hitsx)], collapse=";")
      gr$SYMBOLs[x] <- length(hitsx)-1;
    }

    return(gr)
}

#' Annotate CNA with cgc_gr IDs
#' @param gr GRanges object
#' @param cgc_gr cancer gene census GRanges
#' @param which_genome the genome assembly used, "hg19" or "hg38"
#' @return GRanges object with annotation
#' @export

anno_cgc_cna <- function(gr, cgc_gr, which_genome){
  hits <- as.data.frame(GenomicRanges::findOverlaps(gr, cgc_gr, ignore.strand=TRUE))
  hits$SYMBOL <- cgc_gr[hits$subjectHits]$CGC_gene
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

output_out_list <- function(out_list, in_list, dict_file, which_genome, tag, cgc_gr = NULL){

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
              "Minor_Copy_Number")

  cna_master_gr <- master_intersect_cna_grlist(cna_list, ps_vec, which_genome)

  ##annotate
  anno <- "ENS"
  if(!is.null(cgc_gr)){
    anno <- "CGC"
    cna_master_anno_gr <-  anno_cgc_cna(cna_master_gr, cgc_gr, which_genome)
  } else {
    cna_master_anno_gr <-  anno_ens_cna(cna_master_gr, which_genome)
  }

  GenomeInfoDb::seqinfo(cna_master_anno_gr) <- GenomeInfoDb::seqinfo(cna_list[[1]])
  cna_master_anno_df <- as.data.frame(cna_master_anno_gr)
  cna_master_anno_df$width <- cna_master_anno_df$end - cna_master_anno_df$start
  readr::write_tsv(cna_master_anno_df, file = paste0(tag, ".facets.CNA.master.tsv"))

  ##write output of CNA analysis to XLSX
  cna_df_list <- lapply(cna_list, as.data.frame)
  cna_df_list$master_anno_df <- cna_master_anno_df

  ##replace NAs
  cna_df_list_na <- lapply(cna_df_list, function(f){
    f[is.na(f)] <- "-"
    return(tibble::as_tibble(f))
  })

  ##summarise master gr
  summ_tb <- summarise_master(cna_master_anno_gr)
  cna_df_list_na$summary <- as.data.frame(summ_tb)

  openxlsx::write.xlsx(cna_df_list_na, file = paste0(tag, ".facets.CNA.full.xlsx"))

  ##plot
  plot_out_list(cna_list, pp_list, dict_file, which_genome, tag, samples, write_out = TRUE, max_cna_maxd = 8, sample_map = NULL)
}

#' Plotting from 'out_list' object
#' @param cna_list GRanges in list of facets results
#' @param pp_list ploidy/purity list
#' @param dict_file dictionary file of fasta used to align BAMs
#' @param which_genome the genome assembly used, "hg19" or "hg38"
#' @param tag runID and annotation strategy, "ENS" or "CGC"
#' @param samples vector of sampleIDs (use sample_map to change naming in plots)
#' @param write_out runID and annotation strategy, "ENS" or "CGC"
#' @param max_cna_maxd runID and annotation strategy, "ENS" or "CGC"
#' @param sample_map named vector, names are all in samples, elements are new names
#' @return none writes output files
#' @export

plot_out_list <- function(cna_list, pp_list, dict_file, which_genome, tag, samples, write_out = TRUE, max_cna_maxd = 8, sample_map = NULL){

  if(length(dir(pattern = paste0(tag, ".facets_consensus.plot.pdf"))) == 1){
    print(paste0("tag: ", tag, " has already been used, please use another"))
  } else {

    ##make CNA df
    cna_df_list <- lapply(seq_along(samples), function(x){
      if(length(cna_list[[x]])>0){
        print(samples[x])
        cna_dfb <- as.data.frame(cna_list[[x]], row.names = seq(1:length(cna_list[[x]])))
        cna_dfb$purity <- unlist(rep(pp_list[[x]][2], length(cna_dfb[,1])))
        cna_dfb$ploidy <- unlist(rep(pp_list[[x]][1], length(cna_dfb[,1])))
        if(is.null(sample_map)){
          cna_dfb$sampleID <- samples[x]
        } else {
          cna_dfb$sampleID <- as.vector(sample_map[samples[x]])
        }
        if(write_out == TRUE){
          readr::write_tsv(cna_dfb, file = paste0(samples[x],".facets.CNA.jointsegs.tsv"))
        }
        return(cna_dfb)
      } else {
        if(write_out == TRUE){
          readr::write_tsv(as.data.frame(cna_list[[x]]), file = paste0(samples[x], ".facets.CNA.jointsegs.tsv"))
        }
      }
    })

    cna_df <- do.call(rbind, cna_df_list)

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
    if(max(cna_maxd) > max_cna_maxd){
      cna_maxd[cna_maxd > max_cna_maxd] <- max_cna_maxd
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
    ggp <- ggplot2::ggplot() +
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
                       breaks = seqlengths$cum_sum_centro) +
           ggplot2::scale_y_continuous(name = "Total CNA (facets tcn.em)",
                       labels = seq(from = 0, to = max_cna_maxd, by = 1),
                       breaks = seq(from = 0, to = max_cna_maxd, by = 1),
                       limits = c(0, max_cna_maxd)) +
           ggplot2::theme(axis.text.x = ggplot2::element_text(size = 5),
                       panel.grid.major = ggplot2::element_blank(),
                       panel.grid.minor = ggplot2::element_blank(),
                       legend.position="none") +
           ggplot2::facet_grid(sample ~ .)

    ggplot2::ggsave(filename = paste0(tag, ".facets_consensus.plot.pdf"), plot = ggp)
  }
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

#' Find intersect from list of CNA GRanges
#' @param cna_list list of named GRanges objects
#' @param ps_vec mcols columns to keep p(er) s(ample; appended with list element's name)
#' @param which_genome the genome assembly used, "hg19" or "hg38"
#' @return GRanges object of all CNA, intersecting or not
#'         includes the ranges and each mcol from ps_vec, annotated `ps_vec[x].sample`
#' @export

master_intersect_cna_grlist <- function(gr_list, ps_vec, which_genome){

  ##check gr is A GRanges object
  if(!as.vector(class(gr_list)) %in% c("GRangesList", "list")){
    stop("Input \'gr_list\' is not a GRangesList nor list object, retry")
  }

  ##use bedr::bedr.join.multiple.region to make a complete set and samples
  ##lapply over gr_list which creates character vector list in "chrX:1-100" format
  ##thinks it's 0-based so add 1...
  if(length(gr_list) > 1){

    chr_list <- lapply(gr_list, function(f){
        fb <- apply(as.data.frame(f), 1, function(ff){
            ff <- unlist(ff)
            paste0("chr", ff[1], ":", gsub(" ", "", ff[2]), "-", gsub(" ", "", ff[3]))
          })
        bedr::bedr.sort.region(fb)
      })

    join_chr_all_tb <- tibble::as_tibble(bedr::bedr.join.multiple.region(chr_list, build = which_genome))

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
    join_chr_all_gr <- sortSeqlevels(join_chr_all_gr)
    join_chr_all_gr <- sort(join_chr_all_gr)
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
        ##NB this will almost certainly be a subset of gr_master ranges
        hits <- base::as.data.frame(GenomicRanges::findOverlaps(gr_ff, gr_master, ignore.strand = TRUE, minoverlap = 2))

        ##due to being a subset, we insert NA rows for those missing in subjectHits (master_gr)
        hits_list <- lapply(seq_along(gr_master), function(x){
          if(x %in% hits$subjectHits){
            ho <- hits[hits$subjectHits == x,]
            rownames(ho) <- x
            return(ho)
          } else {
            ho <- data.frame(queryHits = NA, subjectHits = x)
            rownames(ho) <- x
            return(ho)
          }
        })

        hits_na <- do.call(rbind, hits_list)

        ##mcols of those hits
        # mchits <- S4Vectors::mcols(gr_ff[, ps_vec][hits_na$queryHits])
        mchits <- S4Vectors::mcols(gr_ff[, ps_vec])

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
          master_df[hits[sh %in% x, 2],] <- unlist(mchits[hits[sh %in% x, 1],])
        }
        return(master_df)
      })
    }

    ##using the above as a master GRanges object, walk through per sample
    ##create per sample df to be added to mcol
    S4Vectors::mcols(join_chr_all_gr) <- c(S4Vectors::mcols(join_chr_all_gr), do.call(cbind, master_mcols(gr_list, join_chr_all_gr, ps_vec)))

    ##return GRanges with mcols per sample

  } else{
    join_chr_all_tb <- tibble::as_tibble(gr_list[[1]])
    join_chr_all_tb$sampleIDs <- rep(names(gr_list), times = dim(join_chr_all_tb)[1])
    join_chr_all_gr <- GenomicRanges::GRanges(join_chr_all_tb)

    ##return single GRanges

  }
  return(join_chr_all_gr)
}

#' Allow tunable reformatting of plots
#' @param pattern string e.g. "facets.CNA.CGC.tsv" as output by main functions
#' @param dict_file dictionary file of fasta used to align BAMs
#' @param which_genome the genome assembly used, "hg19" or "hg38"
#' @param tag write output with this tag prefix
#' @param write_out write tsv output?
#' @param max_cna_maxd maximum CNA to show on plot
#' @param sample_map named vector, names are filenames before first dot, elements are new names
#' @return none, prints a plot and writes output if flag used
#' @export

plot_from_tsv <- function(pattern, dict_file, which_genome, tag, write_out = FALSE, max_cna_maxd = 8, sample_map = NULL){
  files_in <- dir(pattern = pattern)
  out_list <- lapply(files_in, function(f){
    ff <- readr::read_tsv(f)
    colnames(ff) <- gsub("mcols.", "", colnames(ff))
    gr <- GenomicRanges::GRanges(seqnames = ff$seqnames,
                  ranges = IRanges::IRanges(start = ff$start, end = ff$end),
                  mcols = ff[, c("seg", "num.mark", "nhet", "cnlr.median", "mafR", "segclust", "cnlr.median.clust", "mafR.clust", "cf.em", "Total_Copy_Number", "Minor_Copy_Number", rev(colnames(ff))[5], rev(colnames(ff))[4])])
    return(list(gr, data.frame(PLOIDY = unique(ff$ploidy), PURITY = unique(ff$purity))))
    })
  names(out_list) <- samples <- gsub("\\.", "", gsub(pattern, "", files_in))

  facets_list <- lapply(out_list, function(f){
    ff <- f[[1]]
    names(ff) <- ff$mcols.seg
    return(ff)
  })
  pp_list <- lapply(out_list, function(f){
    return(f[[2]])
  })

  ##list of only those CNA (i.e. not TCN == 2)
  cna_list <- lapply(facets_list, function(f){
    fo <- f[S4Vectors::mcols(f)$Total_Copy_Number != 2]
    sort(unique(c(fo)))
  })

  plot_out_list(cna_list, pp_list, dict_file, which_genome, tag, samples, write_out = TRUE, max_cna_maxd, sample_map)
}

#' Summarise master table
#' @param cna_master_anno_gr annotated GRanges from master_intersect_cna_grlist + anno_cgc_cna
#' @return list of various delights per sampleIDs combinations extant
#' @export

summarise_master <- function(cna_master_anno_gr){

  cna_master_anno_tb <- tibble::as_tibble(cna_master_anno_gr)
  genome_size <- sum(GenomeInfoDb::seqlengths(GenomeInfoDb::seqinfo(cna_master_anno_gr)))

  ##find widths for each set of sampleIDs found overlapped
  width_unique_list <- lapply(unique(cna_master_anno_tb$sampleIDs), function(f){
    ff <- dplyr::filter(.data = cna_master_anno_tb, sampleIDs %in% f)
    dplyr::summarise(.data = ff, width)
  })

  ##sum those to get summary
  width_sum_unique_list <- lapply(width_unique_list, function(f){
    return(list(sum = sum(f), prop = sum(f)/genome_size))
  })

  ##genes affected by those
  symbol <- grep("_SYMBOL", colnames(cna_master_anno_tb), value = TRUE)[1]

  genes_affected_list <- lapply(unique(cna_master_anno_tb$sampleIDs), function(f){
    ff <- dplyr::filter(.data = cna_master_anno_tb, sampleIDs %in% f)
    fg <- dplyr::summarise(.data = ff, !!as.symbol(symbol))
    fgo <- unlist(na.omit(fg))
    fgo <- gsub("NA;", "", fgo)
    fgo[fgo!="-"]
  })

  #losses and gains
  summarise_tcn <- function(rowin){
    apply(rowin, 1, function(r){
      rowina <- na.omit(r)
      if(all(rowina>2)){
        return("GAIN")
      } else if(all(rowina < 2)){
        return("LOSS")
      } else {
        return("BOTH")
      }
    })
  }

  tcn_summary_list <- lapply(unique(cna_master_anno_tb$sampleIDs), function(f){
    ff <- dplyr::filter(.data = cna_master_anno_tb, sampleIDs %in% f)
    ft <- dplyr::select(.data = ff, tidyselect::ends_with("Total_Copy_Number"))
    table(summarise_tcn(ft))
  })

  ##name all the same
  names(tcn_summary_list) <- names(genes_affected_list) <- names(width_unique_list) <- names(width_sum_unique_list) <- unique(cna_master_anno_tb$sampleIDs)

  ##create table output
  genesaf <- unlist(lapply(genes_affected_list, function(f){paste(sort(f), collapse = ",")}))

  summ_tb <- tibble::tibble(sampleIDs = unique(cna_master_anno_tb$sampleIDs),
                            tcn_summary = unlist(lapply(tcn_summary_list, function(f){paste(names(f), collapse = ";")})),
                            tcn_ns = unlist(lapply(tcn_summary_list, function(f){paste(f, collapse = ";")})),
                            width_sum_total = unlist(lapply(width_sum_unique_list, function(f){f$sum})),
                            width_sum_prop = unlist(lapply(width_sum_unique_list, function(f){f$prop})),
                            n_genes_affected = unlist(lapply(genesaf, function(f){length(strsplit(f, ",")[[1]])})),
                            genes_affected = genesaf)

  return(dplyr::arrange(.data = summ_tb, desc(width_sum_prop)))
}

#' Take input of a GRanges object, and use that seqinfo to bin by bin_size
#'
#' @param gr GRanges object with seqinfo used to make bins
#' @param genome use this to add seqinfo if none extant (e.g. "hg19", "hg38")
#' @param bin_size base pairs by which to bin
#' @param mcol_in string naming mcols to use for binning (one only)
#' @param mcol_in_lim numeric on which to screen mcol_in column; if length>1, take all outside the range suplied
#' @param one_to_x logical indicating whether to use chr1..chrX for bin, as opposed to what is in input GRanges
#' @return GRanges binned on mcol_in
#' @export

bin_granges <- function(gr, genome = NULL, bin_size=NULL, mcol_in, mcol_in_lim = NULL, one_to_x = TRUE){

  ##check gr is A GRanges object
  if(!as.vector(class(gr)) %in% c("GRanges", "GRangesList", "list")){
    stop("Input \'gr\' is not a GRanges, GRangesList nor list object, retry")
  }

  ##if grList, set gr to a single entity therein for binning
  if(as.vector(class(gr)) %in% c("GRangesList", "list")){
    gri <- gr[[1]]
  } else {
    if(!as.vector(class(gr)) == "GRanges"){
      stop("Input \'gr\' list does not contain a GRanges object, retry")
    } else {
      gri <- gr
    }
  }

  ##test seqinfo
  print("Testing seqinfo on GRanges input...")
  gri <- test_seqinfo(gr = gri, genome = genome)

  ##set default
  if(is.null(bin_size)){
    print("Bin size is not specified, setting to 10MB")
    bin_size <- 10000000
  } else {
    bin_size <- as.numeric(bin_size)
  }

  ##make bin gr
  print("Making bins...")
  grb <- bin_maker(grr = gri, bin_size = bin_size, genome = genome, one_to_x = TRUE)

  ##screen based on mcols_in_lim
  if(!is.null(mcol_in_lim)){
    if(length(mcol_in_lim) == 1) {
      grf <- gri[unlist(S4Vectors::mcols(gri[,mcol_in]))!=mcol_in_lim]
    } else {
      grf <- gri[unlist(S4Vectors::mcols(gri[,mcol_in]))<mcol_in_lim[1] & unlist(S4Vectors::mcols(gri[,mcol_in]))>mcol_in_lim[2] ]
    }
  } else {
    grf <- gri
  }

  ##findoverlaps, mcols on mcol_in
  print("Finding overlaps...")
  hits <- as.data.frame(GenomicRanges::findOverlaps(grf, grb, ignore.strand=TRUE, type = "any"))
  bin_mcol_qh <- S4Vectors::mcols(grf)[hits$queryHits, mcol_in]
  S4Vectors::mcols(grb)[mcol_in] <- NA
  S4Vectors::mcols(grb)[hits$subjectHits, mcol_in] <- bin_mcol_qh

  return(grb)

  ##old
  ##N.B. that hits can include multiple 'subjectHits', i.e. bins can be hit by more that one CNA
  ##here we test if: those are the same value (and include just that value)
  ##different value with one == 2, include other
  ##fail, reduce bin size below
  # grb$x <- S4Vectors::mcols(grf[hits$queryHits, mcol_in])
  # gr$CGC_SYMBOLs <- "-"
  # ##loop to collapse symbols per region
  # for(x in 1:max(hits$queryHits)){
  #   hitsx <- as.vector(sort(unique(hits$SYMBOL[hits$queryHits==x])))
  #   hitsx <- hitsx[!is.na(hitsx)]
  #   if(length(hitsx)==0){gr$CGC_SYMBOL[x] <- NA; gr$CGC_SYMBOLs[x] <- 0}
  #   else{
  #     gr$CGC_SYMBOL[x] <- base::paste(hitsx[2:length(hitsx)], collapse=";")
  #     gr$CGC_SYMBOLs[x] <- length(hitsx)-1;
  #   }
  # }
}

#' Take input of a GRanges object, and use that seqinfo to bin by bin_size
#'
#' @param grr GRanges object with seqinfo used to make bins
#' @param bin_size base pairs by which to bin
#' @param genome use this to add seqinfo if none extant (e.g. "hg19", "hg38")
#' @param one_to_x logical indicating whether to use chr1..chrX for bin, as opposed to what is in input GRanges
#' @return GRanges bin
#' @export

bin_maker <- function(grr, bin_size = 30000, genome = NULL, one_to_x = TRUE){

  if(one_to_x == FALSE){
    ##new GRanges across that bin based on that seqinfo
    grri <- test_seqinfo(grr, genome = genome)

    sn <- GenomeInfoDb::seqnames(GenomeInfoDb::seqinfo(grri))
    sl <- GenomeInfoDb::seqlengths(GenomeInfoDb::seqinfo(grri))
    seqinf <- GenomeInfoDb::seqinfo(grri)
    genome <- as.vector(GenomeInfoDb::genome(gri)[1])

  } else {
    sn <- c(1:22, "X")
    seqinf <- GenomeInfoDb::fetchExtendedChromInfoFromUCSC(genome)
    seqinf <- seqinf[seqinf$NCBI_seqlevel %in% sn,]
    sl <- seqinf$UCSC_seqlength
    genome <- genome
  }
  ##create IRanges for seqnames and seqlengths
  ##combine into bin GRange
  print(sn)
  print(sl)
  print(seqinf)

  grbList <- lapply(seq_along(sn), function(f){
      st <- seq(from = 0, to = sl[f], by = bin_size)
      nd <- c(seq(st[2], sl[f], by = bin_size), as.vector(sl[f])+1)-1
      grsn <- GenomicRanges::GRanges(seqnames = sn[f],
                                     IRanges::IRanges(start = st, end = nd),
                                     strand = "*")
      GenomeInfoDb::genome(grsn) <- genome
      GenomeInfoDb::seqlengths(grsn) <- sl[f]
      return(grsn)
  })
  grb <- unlist(as(grbList, "GRangesList"))
  S4Vectors::mcols(grb) <- bin_size
  names(S4Vectors::mcols(grb)) <- "bin_size"
  return(grb)
}

#' Test for seqinfo

#' @param gr GRanges to test for seqinfo
#' @param genome string indicating genome to add seqinfo from if missing
#' @return GRanges gr with seqinfo from genome if required
#' @export

test_seqinfo <- function(gr, genome = NULL){
  adf <- as.data.frame(GenomeInfoDb::seqinfo(gr))[1,1]
  if(is.na(adf)){
    if(is.null(genome)){
      stop("Input has no seqinfo defined, and no genome specified, please retry with seqinfo or genome")
    } else {
      print(paste0("Using ", genome, " to get seqinfo"))
      seqinf <- GenomeInfoDb::fetchExtendedChromInfoFromUCSC(genome)
      seqinf <- seqinf[seqinf$NCBI_seqlevel %in% GenomeInfoDb::seqlevels(gr),]
      GenomeInfoDb::seqlengths(gr) <- seqinf$UCSC_seqlength
      GenomeInfoDb::genome(gr) <- genome
      return(gr)
    }
  } else {
    print("Input has seqinfo already")
    return(gr)
  }
}
