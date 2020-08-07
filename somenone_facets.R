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
  fit <- facets::emcncf(oo)
  ploidpurdf <- data.frame(PLOIDY = round(fit$ploidy, digits=3),
                           PURITY = round(fit$purity, digits=3))
  readr::write_tsv(ploidpurdf, path = paste0(tumour,".fit_ploidy-purity.tsv"))
  readr::write_tsv(round(fit$cncf, 3), path = paste0(tumour,".fit_cncf-jointsegs.tsv"))
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
  which_genome <- "19"
  if(length(grep("38", dict_file))==1) {
    which_genome <- "38"
  }

  ##load Cancer Gene Census bed file if supplied
  cgc <- c()
  cgc_gr <- GenomicRanges::GRanges()
  if(! is.null(cgc_bed)){
    cgc <- utils::read.table(cgc_bed)
    cgc_gr <-GenomicRanges::GRanges(seqnames = cgc[,1],
                                    ranges = IRanges::IRanges(cgc[,2], cgc[,3]),
                                    CGC_gene = cgc[,4])
  }

  ##set input list based on file pattern
  in_list <- as.list(dir(pattern = pattern))
  pp_list <- process_in_list(in_list, which_genome)

  ##load set of RData GRanges created above into GRangesList
  gr_in <- dir(pattern = ".RData$")
  loaded_gr <- unlist(lapply(gr_in, function(x){
    load(x, envir = .GlobalEnv)
  }))
  facets_list <- unlist(lapply(loaded_gr, function(f){get(f)}))
  samples <- names(facets_list) <- stringr::str_split(loaded_gr,"\\.")[[1]][1]

  cna_list <- lapply(facets_list, function(f){
    fo <- f[S4Vectors::mcols(f)$Total_Copy_Number != 2 & S4Vectors::mcols(f)$Total_Copy_Number != 1]
    fo <- unique(c(fo))
  })

  ##make holder df
  cna_df <- c()
  for(x in 1:length(cna_list)){
    if(length(cna_list[[x]]) > 0){
      cna_df <- as.data.frame(cna_list[[x]], row.names = seq(1:length(cna_list[[x]])))
      cna_df$sampleID <- cna_df$ploidy <- cna_df$purity <- "-"
      cna_df <- cna_df[-1,]
      break;
    }
  }

  for(x in 1:length(samples)){
    if(length(cna_list[[x]])>0){
      cna_dfb <- as.data.frame(cna_list[[x]], row.names = seq(1:length(cna_list[[1]])))
      cna_dfb$purity <- unlist(rep(pp_list[[x]][2], length(cna_dfb[,1])))
      cna_dfb$ploidy <- unlist(rep(pp_list[[x]][1], length(cna_dfb[,1])))
      cna_dfb$sampleID <- samples[x]
      cna_df <- rbind(cna_df, cna_dfb)
      utils::write.table(cna_df, file=paste0(samples[x],".facets.CNA.jointsegs.tab"), quote=FALSE,sep="\t",col=TRUE)
    } else {
      utils::write.table(as.data.frame(cna_list[[x]]), file=paste0(samples[x],".facets.CNA.jointsegs.tab"), quote=FALSE,sep="\t",col=TRUE)
    }
  }

  ##remove character chromosomes NB no MT, GL...
  sqn <- gsub("chr", "", as.vector(cna_df[,1]))
  sqn[sqn=="X"] <- 23
  sqn[sqn=="Y"] <- 24
  cna_df$seqnames <- sqn

  seqlength <- seqlengths_df(unique(cna_df[,1]), dict_file, which_genome)
  cna_df <- cna_df[cna_df$seqnames %in% seqlength$seqnames,]
  cum_sum_add <- unlist(apply(cna_df, 1, function(x) {
              cum_sum_off <- seqlength$cum_sum_0[rownames(seqlength) %in% x[1]]
              return(cum_sum_off)
  }))

  ##findOverlaps of all samples in list, write overlaps out
  if(length(facets_list) > 1){
    cna_df_sample_intersect <- Reduce(intersect, facets_list)
  } else {
    cna_df_sample_intersect <- facets_list
  }
  readr::write_tsv(cna_df_sample_intersect,
            file = paste0(tag,".facets_consensus.CNA.tsv"))

  ##plotting
  ggp <- plot_cna_df(cna_df, cum_sum_add, cna_max = 8)
  ggplot2::ggsave(filename = paste0(tag, ".facets_consensus.call.pdf"), plot = ggp)
}

#' Processing list of Facets input
#' @param in_list of data from Facets
#' @param which_genome genoem assembly being used "hg19" or "hg38"
#' @param cgc_bed bedfile of cancer gene census from COSMIC
#' @return purity-ploidy dataframe, also write rds used
#' @export

process_in_list <- function(in_list, which_genome, cgc_bed = NULL){

  ##iterate over samples
  for(samp in 1:length(in_list)){
    sample <- in_list[[samp]]
    in_name <- out_name <- stringr::str_split(sample,"\\.")[[1]][1]
    out_ext <- gsub(".tsv", "", sample)

    ##also read and store ploidy:
    ploidypurity <- dir(pattern = paste0(in_name, ".fit_ploidy-purity.tab"))
    ploidy <- as.vector(utils::read.table(ploidypurity)[2,1])
    purity <- as.vector(utils::read.table(ploidypurity)[2,2])
    pp_df <- data.frame(PLOIDY = ploidy, PURITY = purity)

    print(paste0("Working on: ", in_name))

    ##run function to make GRangesList
    ##(jointsegs_in, which_genome, cgc_gr = NULL, anno = NULL, bsgenome = NULL)
    grl <- facets_jointsegs_parse_to_gr(sample, which_genome, anno= "ENS")
    names(grl) <- in_name

    ##assign output
    if(!is.null(cgc_bed)){
      assigned_name <- paste0(out_ext,".CGC")
      assign(assigned_name, value = grl)
      save_file <- paste0(out_ext,".CGC.RData")
    } else {
      assigned_name <- paste0(out_ext,".gr")
      assign(assigned_name, value = grl)
      save_file <- paste0(out_ext,".gr.RData")
    }

    ##save to current dir
    save(list = assigned_name, file = save_file)
  }

  return(pp_df)
}

#' Parse jointsegs into GRanges object
#' @param jointsegs_in filename (output from Facets) with columns:
#'    "seg", "num.mark", "nhet", "cnlr.median", "mafR", "segclust", "cnlr.median.clust", "mafR.clust", "cf.em", "tcn.em", "lcn.em"
#' @param which_genome the genome assembly used, "hg19" or "hg38"
#' @param cgc_gr GRanges object of cancer gene census from COSMIC
#' @param anno annotation strategy, "ENS" or "CGC"
#' @param bsgenome which version of bsgenome to use
#' @return a GRanges object with annotated jointsegs from input
#' @export

facets_jointsegs_parse_to_gr <- function(jointsegs_in, which_genome, cgc_gr = NULL, anno = NULL, bsgenome = NULL){

  ##bed file from: https://cancer.sanger.ac.uk/census
  ##NB requires login hence not provided
  if(!is.null(cgc_gr)){
    CGCGR <- cgc_gr
  }

  ##annotation source
  if(is.null(anno)){
    anno <- "NONE"
  }

  ##default genome version is hg38
  if(which_genome == "hg19"){
    bsgenome <- BSgenome.Hsapiens.UCSC.hg19::BSgenome.Hsapiens.UCSC.hg19
  } else {
    bsgenome <- BSgenome.Hsapiens.UCSC.hg38::BSgenome.Hsapiens.UCSC.hg38
  }

  ##read input
  js <- utils::read.table(jointsegs_in, header=T)

  ##convert chr23 -> X
  js$chrom[js$chrom == 23] <- "X"
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
    GenomeInfoDb::genome(genes) <- "hg19"
    GenomeInfoDb::seqlevels(genes, pruning.mode="coarse") <- GenomeInfoDb::seqlevels(gr)

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
  }

  if(anno == "NONE"){
    ##rename S4Vectors::mcols
    names(S4Vectors::mcols(gr)) <- c(names(S4Vectors::mcols(gr))[1:9], "Total_Copy_Number", "Minor_Copy_Number")
  }
  return(gr)
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
  seqlengths <- data.frame(seqnames = gsub("SN:","", dict_seqs[,2]),
                           start = 1,
                           end = as.numeric(gsub("LN:","", dict_seqs[,3])))

  ##find those chromosomes in inSeqs
  seqlengths <- seqlengths[is.element(seqlengths[,1],
                           unique(as.vector(in_seqs))),]

  ## make centromere data from function
  ##load centromere
  centromere <- centromeres(seqlengths[,1], which_genome)
  seqlengths <- dplyr::left_join(seqlengths, centromere, by="seqnames")

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
#' @param cna_max max CNA call to plot (eveything else is squashed to this
#' @return purity-ploidy dataframe, also write rds outputs
#' @export

plot_cna_df <- function(cna_df, cum_sum_add, cna_max = 8){

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
