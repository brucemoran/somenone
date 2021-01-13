#' GRIDSS functions
#'
#' Wrapper to plot multi-sample VCF list (recurrent and not)
#'
#' @param vcf, the VCF file path
#' @param dict_file dictionary file for fasta used in alignment
#' @param which_genome, genome assembly used ("hg19", "hg38")
#' @param output_path, file path to output plots
#' @return none, plots recurrent and private SVs and writes output to TSV format file
#' @export

gridss_parse_plot <- function(vcf, dict_file, germline_id, which_genome = NULL, output_path = NULL){

  ##set genome as hg19 unless spec'd
  if (is.null(which_genome)){
    if(unlist(lapply(c("19", "37"), function(f){grep(f, dict_file)})) > 0){
    which_genome <- "hg19"
    }
    else{
      which_genome <- "hg38"
    }
  }

  ##set output_path based on VCF if NULL
  if (is.null(output_path)){
    output_path <- paste0(stringr::str_split(vcf, "\\.vcf")[[1]][1])
  }

  gridss_parsed_multi_list <- gridss_parse_multi_vcf(vcf = vcf,
                                                     germline_id = germline_id,
                                                     which_genome = which_genome)
  rec_df <- as.data.frame(gridss_parsed_multi_list[[1]])
  prv_df <- as.data.frame(gridss_parsed_multi_list[[2]])

  ##plotting
  if(dim(rec_df)[1] > 0){
  prep_plot_circos_sv(input_df = rec_df,
                      dict_file = dict_file,
                      which_genome = which_genome,
                      output_path = paste0(output_path, ".recurrent"))
  }

  prep_plot_circos_sv(input_df = prv_df,
                      dict_file = dict_file,
                      which_genome = which_genome,
                      output_path = paste0(output_path, ".private"))
}

#' Parse multi-sample VCF (recommended output of GRIDSS)
#'
#' @param vcf, the VCF file path
#' @param which_genome, genome assembly used ("hg19", "hg38")
#' @return a list with recurrent SVs and private SVs as elements
#' @export

gridss_parse_multi_vcf <- function(vcf, germline_id, which_genome = NULL){

  options(stringsAsFactors = FALSE)
  sampleID <- num_rows <- qualscore <- vars <- bedpe_list_filter <- n <- NULL

  ##set genome as hg19 unless spec'd
  if (is.null(which_genome)){
    which_genome <- "hg19"
  }

  ##set samplenames
  sample_names <- grep(germline_id,     VariantAnnotation::samples(VariantAnnotation::scanVcfHeader(vcf)), value = TRUE, invert = TRUE)

  ##readVcf and get per-sample GRanges
  rvcf <- VariantAnnotation::readVcf(vcf, genome = which_genome)
  bpgr <- StructuralVariantAnnotation::breakpointRanges(rvcf)
  #pairs <- StructuralVariantAnnotation::breakpointgr2pairs(gr)

  ##iterate over sample_names
  ##https://github.com/PapenfussLab/StructuralVariantAnnotation/issues/31
  bpgr_list <- lapply(sample_names, function(samp) {

    ##get sample-specific bpgr object
#    bpgr[VariantAnnotation::geno(rvcf[bpgr$sourceId])$QUAL[, samp] >= 2]
    ##older version has vcfId vs. sourceId
    bpgr[VariantAnnotation::geno(rvcf[names(bpgr)])$QUAL[, samp] >= 2]

  })

  names(bpgr_list) <- sample_names

  ##make single bedpe for all VCFs
  gridss_bedpe_list <- somenone::gridss_bedpe_list(bpgr_list, which_genome)
  bedpe_rb <- do.call(rbind, gridss_bedpe_list)

  ##find reccuring events (same at 1,2,3,4,5,6,7,8)
  bedpe_rb_rec <- dplyr::mutate_if(bedpe_rb, is.factor,as.numeric)
  bedpe_rb_rec <- dplyr::group_by_at(bedpe_rb_rec, dplyr::vars(-qualscore, -sampleID))
  bedpe_rb_rec <- dplyr::mutate(bedpe_rb_rec, num_rows = base::sum(dplyr::n()))
  bedpe_rb_rec <- dplyr::filter(bedpe_rb_rec, num_rows > 1)
  bedpe_rb_rec <- dplyr::ungroup(bedpe_rb_rec)
  bedpe_rb_rec <- dplyr::select(bedpe_rb_rec, -qualscore, -sampleID, -num_rows)
  bedpe_rb_rec <- dplyr::distinct(bedpe_rb_rec)

  ##lapply to paste sampleIDs, colour by combinations therein
  if(dim(bedpe_rb_rec)[1] > 0){
    sample_rec_list <- lapply(1:dim(bedpe_rb_rec)[1], function(f){
      lj <- dplyr::left_join(bedpe_rb_rec[f,],
                bedpe_rb,
                by=c("start1" = "start1",
                     "end1" = "end1",
                     "start2" = "start2",
                     "end2" = "end2"))
      ljs <- dplyr::select(lj, sampleID)
      lju <- unlist(ljs)
      ljs <- sort(lju)
      paste(ljs, collapse=",")
    })

    ##add to bedpe_rb_rec
    bedpe_rb_rec$sampleID <- unlist(sample_rec_list)

    ##colour in based on samples with SV
    colz_rec <- grDevices::rainbow(n = length(table(unlist(sample_rec_list))))
    names(colz_rec) <- names(table(unlist(sample_rec_list)))
    bedpe_rb_rec$colour <- colz_rec[match(bedpe_rb_rec$sampleID, names(colz_rec))]
  }

  ##non-recurrent
  bedpe_rb_pri <- dplyr::mutate_if(bedpe_rb, is.factor,as.numeric)
  bedpe_rb_pri <- dplyr::group_by_at(bedpe_rb_pri, dplyr::vars(-qualscore, -sampleID))
  bedpe_rb_pri <- dplyr::mutate(bedpe_rb_pri, num_rows = sum(dplyr::n()))
  bedpe_rb_pri <- dplyr::filter(bedpe_rb_pri, num_rows == 1)
  bedpe_rb_pri <- dplyr::ungroup(bedpe_rb_pri)
  bedpe_rb_pri <- dplyr::select(bedpe_rb_pri, -num_rows)
  colz_pri <- grDevices::rainbow(n = length(table(unlist(bedpe_rb_pri$sampleID))))
  names(colz_pri) <- names(table(unlist(bedpe_rb_pri$sampleID)))
  bedpe_rb_pri$colour <- colz_pri[match(bedpe_rb_pri$sampleID, names(colz_pri))]

  return(list(bedpe_rb_rec, bedpe_rb_pri))
}

#' Produce a list of bedpe, from a list of breakpoint GRanges,
#'    with score > qual_filt and both breakpoints within annotated genes
#' @param bpgr_list list of per-sample breakpont GRanges objects
#' @param which_genome, genome assembly used ("hg19", "hg38")
#' @param qual_filt quality score filter above which is returned
#' @return filtered bedpe list
#' @export

gridss_bedpe_list <- function(bpgr_list, which_genome, qual_filt = 2) {

  start <- end <- vcf_list <- symbol1 <- symbol2 <- NULL

  bedpe_list_filter <- lapply(seq_along(bpgr_list), function(f){

    print(paste0("Working on: ", names(bpgr_list)[f]))
    gr <- bpgr_list[[f]]

    ##remove non-found partners!
    pairs <- StructuralVariantAnnotation::breakpointgr2pairs(gr)
    grp <- gr[names(gr) %in% S4Vectors::mcols(pairs)[,"name"]]
    gr_partner <- gr[names(gr) %in% grp$partner,]

    ##force into only 1:22, X
    GenomeInfoDb::seqlevels(grp, pruning.mode = "coarse") <- GenomeInfoDb::seqlevels(grp)[c(1:23)]

    ##annotate
    grp_anno <- gridss_annotate_gr(grp, which_genome)
    gr_partner_anno <- gridss_annotate_gr(gr_partner, which_genome)
    grp_anno_use <- grp_anno[names(grp_anno) %in% gr_partner_anno$partner]
    gr_partner_anno_use <- gr_partner_anno[names(gr_partner_anno) %in% grp_anno_use$partner]

    #bedpe for input, simplest format
    bedpe1k <- data.frame(chrom1 = as.data.frame(grp_anno_use)[,1],
                          start1 = c(as.data.frame(grp_anno_use)[,2]),
                          end1 = c(as.data.frame(grp_anno_use)[,3]),
                          symbol1 = grp_anno_use$SYMBOL,
                          chrom2 = as.data.frame(gr_partner_anno_use)[,1],
                          start2 = c(as.data.frame(gr_partner_anno_use)[,2]),
                          end2 = c(as.data.frame(gr_partner_anno_use)[,3]),
                          symbol2 = gr_partner_anno_use$SYMBOL,
                          qualscore = grp_anno_use$QUAL,
                          sampleID = names(bpgr_list)[f])
    bedpe1k <- dplyr::filter(bedpe1k, ! symbol1 %in% "",! symbol2 %in% "")
    return(bedpe1k[bedpe1k$qualscore > qual_filt,])
  })
  return(bedpe_list_filter)
}

#' Annotate a GRanges object from GRIDSS
#' @param gr GRanges object
#' @param which_genome, genome assembly used ("hg19", "hg38")
#' @return annotated GRanges object
#' @export

gridss_annotate_gr <- function(gr, which_genome) {

  # annotate breakends with gene names and gene orientation
  if(which_genome == "hg19"){
    gns <- suppressMessages(GenomicFeatures::genes(TxDb.Hsapiens.UCSC.hg19.knownGene::TxDb.Hsapiens.UCSC.hg19.knownGene))
  } else {
    gns <- suppressMessages(GenomicFeatures::genes(TxDb.Hsapiens.UCSC.hg38.knownGene::TxDb.Hsapiens.UCSC.hg38.knownGene))
  }

  ##remove chr if chr in seqnames
  if(length(grep("chr", levels(GenomeInfoDb::seqnames(gr)))) == 0){
    gns <- GenomeInfoDb::renameSeqlevels(gns, sub("chr", "", GenomeInfoDb::seqlevels(gns)))
  }

  hits <- as.data.frame(GenomicRanges::findOverlaps(gr,
                                                    gns,
                                                    ignore.strand=TRUE))
  if(dim(hits)[1] > 0){
    hits$SYMBOL <- suppressMessages(biomaRt::select(org.Hs.eg.db::org.Hs.eg.db,
                                gns[hits$subjectHits]$gene_id, "SYMBOL")$SYMBOL)
    hits$gene_strand <- as.character(BiocGenerics::strand(gns[hits$subjectHits]))

    ##add annotation to grp
    gr$SYMBOL <- gr$geneStrand <- ""
    gr$SYMBOL[hits$queryHits] <- hits$SYMBOL
    gr$geneStrand[hits$queryHits] <- hits$gene_strand

    return(gr)
  }
}

#' OmicCircos plotting from data.frame from GRIDSS bedpe
#' @param input_df data.frame with chrom1, start1, end1, symbol1,
#'   chrom2, start2, end2, symbol2, sampleID, colour columns
#' @param which_genome, genome assembly used ("hg19", "hg38")
#' @param dict_file dictionary file for fasta used in alignment
#' @param output_path path to where PDF file with output plot is written
#'    N.B. two pages, one with Circos plot and one with legend
#' @return none, prints plot to file, and a file too
#' @export

prep_plot_circos_sv <- function(input_df, which_genome, dict_file, output_path){

  if(!is.data.frame(input_df)){
    print("Require data.frame input with chromsomes as factors")
  } else {

    ##write output table
    readr::write_tsv(input_df, file = paste0(output_path, ".tsv"))

    ##impose condition that over 200 SVs in input_df
    ##plot only those that are on different chromosomes
    input_df_tr_2 <- NULL

    if(dim(input_df)[1] > 50){
      print(paste0("Total of ", dim(input_df)[1], " SV found, also plotting those on different chromosomes only"))
      input_df_tr_2 <- input_df[as.vector(input_df$chrom1) != as.vector(input_df$chrom2),]
      input_df_tr_1 <- input_df[order(input_df$qualscore, decreasing = TRUE),]
    } else {
      input_df_tr_1 <- input_df[order(input_df$qualscore, decreasing = TRUE),]
    }

    ##create seqlengths df for ideogram
    url19 <- paste0("http://hgdownload.cse.ucsc.edu/goldenpath/hg19/database/cytoBand.txt.gz")
    url38 <- paste0("http://hgdownload.cse.ucsc.edu/goldenpath/hg38/database/cytoBand.txt.gz")
    tempf <- tempfile()
    if(which_genome == "hg19"){
      utils::download.file(url19, tempf)
      cytoband <- utils::read.table(tempf)
    } else {
      utils::download.file(url19, tempf)
      cytoband <- utils::read.table(tempf)
    }

    plot_circos_sv(input_df_tr_1, output_path, cytoband)
    if(!is.null(input_df_tr_2)){
      plot_circos_sv(input_df_tr_2, paste0(output_path, ".diff_chrs"), cytoband)
    }
  }
}

#' OmicCircos plotting from data.frame from GRIDSS bedpe
#' @param input_df data.frame with chrom1, start1, end1, symbol1,
#'   chrom2, start2, end2, symbol2, sampleID, colour columns
#' @param output_path path to where PDF file with output plot is written
#'    N.B. two pages, one with Circos plot and one with legend
#' @return none, prints plot to file
#' @export

plot_circos_sv <- function(input_df, output_path, cytoband){

  ##circlize requires 'chr' label on chroms
  plot_df_list <- parse_input_df(input_df)

  labels_o20 <- rbind(plot_df_list[["labels_o1"]][1:20,], plot_df_list[["labels_o2"]][1:20,])

  ##cChimerKBv4
  chim_df <- findin_chimerkb4(input_df)

  ##circize plot
  ##colours per-sample
  if(!is.null(chim_df)){
    output_path <- paste0(output_path, ".ChimerKB_V4")
  }
  grDevices::pdf(paste0(output_path, ".pdf"), width = 9, height = 9)

  ##initialise blank ideogram
  circlize::circos.initializeWithIdeogram(plotType = NULL)

  ##labels on outer track if enough space...
  if(!is.null(chim_df)){
    chim_df_list <- parse_input_df(chim_df)
    circlize::circos.genomicLabels(chim_df_list[["labels_o"]],
                         labels.column = "symbol",
                         side = "outside",
                         col = chim_df_list[["labels_o"]]$colour,
                         line_col = chim_df_list[["labels_o"]]$colour)

    readr::write_tsv(chim_df, file = paste0(output_path, ".ChimerKB_v4.tsv"))

  } else {
    if(dim(region_1)[1] < 20) {
      circlize::circos.genomicLabels(plot_df_list[["labels_o"]],
                           labels.column = "symbol",
                           side = "outside",
                           col = plot_df_list[["labels_o"]]$colour,
                           line_col = plot_df_list[["labels_o"]]$colour)
    } else {
      ##top 20 best by qualscore
      circlize::circos.genomicLabels(labels_o20,
                           labels.column = "symbol",
                           side = "outside",
                           col = labels_o20$colour,
                           line_col = labels_o20$colour)
    }
  }

  ##cytoband data
  circlize::circos.genomicIdeogram(cytoband)

  ##links
  circlize::circos.genomicLink(region1 = plot_df_list[["region_1"]],
                               region2 = plot_df_list[["region_2"]],
                               col = plot_df_list[["region_c"]],
                               lwd = plot_df_list[["log10quals"]]/1.8)

  if(!is.null(chim_df)){
    circlize::circos.genomicLink(region1 = chim_df_list[["region_1"]],
                                region2 = chim_df_list[["region_2"]],
                                col = "#00000080",
                                lwd = 10)
  }

  ##legend
  nms <- gsub(",", ", ", unique(names(region_c)))
  colz <- unique(region_c)
  lgnd <- ComplexHeatmap::Legend(at = nms,
                                 type = "lines",
                                 legend_gp = grid::gpar(col = colz, lwd = 4),
                                 title_position = "topleft",
                                 title = "Legend")

  ##plot legend on a new page to allow manipulation into final image
  grid::grid.newpage()
  grid::grid.draw(lgnd)

  grDevices::dev.off()
}

#' Download and parse ChimerKB v4 for annotating known fusions
#' @return chimerkb tibble of parsed ChimerKB v4
#' @export

download_chimerkb4 <- function(){
  urlckb <- paste0("https://www.kobic.re.kr/chimerdb_mirror/downloads?name=ChimerKB4.xlsx")
  tempf <- tempfile()
  utils::download.file(urlckb, tempf)
  chimerkb4 <- readxl::read_xlsx(tempf, sheet = 1)
  return(chimerkb4)
}

#' find ChimerKB data in user data
#' @param input_df data.frame with chrom1, start1, end1, symbol1,
#'   chrom2, start2, end2, symbol2, sampleID, colour columns
#' @return output_df rows of input_df where both fusion partners
#'   were found in ChimerKB4
#' @export

findin_chimerkb4 <- function(input_df){

  ##DL ChimerKB4
  chimerkb4 <- download_chimerkb4()

  ##unique sets of genes
  chimerkb4_uni_genes <- unique(c(unlist(chimerkb4$H_gene),
                                  unlist(chimerkb4$T_gene)))
  input_df_uni_genes <- unique(c(input_df$symbol1,
                                  input_df$symbol2))
  chimerkb4_idf_genes <- chimerkb4_uni_genes[chimerkb4_uni_genes %in% input_df_uni_genes]

  ##create test set of chimerkb4
  chimerkb4_test <- chimerkb4[chimerkb4$H_gene %in% chimerkb4_idf_genes & chimerkb4$H_gene %in% chimerkb4_idf_genes,]

  ##find lines in input_df supporting pairs in chimerkb4_test
  nn <- lapply(1:dim(input_df)[1], function(f){
    f <- input_df[f,]
    if(!is.na(f["symbol1"]) & !is.na(f["symbol2"])){
      if(f["symbol1"] %in% chimerkb4_test$H_gene){
        ct <- dplyr::filter(.data = chimerkb4_test, H_gene %in% f["symbol1"])
        if(f["symbol2"] == ct$T_gene){
          return(f)
        }
      }
      if(f["symbol2"] %in% chimerkb4_test$H_gene){
        ct <- dplyr::filter(.data = chimerkb4_test, H_gene %in% f["symbol2"])
        if(f["symbol1"] == ct$T_gene){
          return(f)
        }
      }
    }
  })

  output_df <- do.call(rbind, nn[lengths(nn)!=0])

  return(output_df)
}

#' find ChimerKB data in user data
#' @param input_df data.frame with chrom1, start1, end1, symbol1,
#'   chrom2, start2, end2, symbol2, qualscore, sampleID, colour columns
#' @return plot_input_df_list list with elements input_df (as per input), region1, region2, region_c, log10quals, labels_o1, labels_o2, labels_o
#' @export

parse_input_df <- function(input_df){

  if(length(grep("chr", input_df$chrom1[1])) == 0){
    input_df[,1] <- paste0("chr", input_df[,1])
    input_df[,5] <- paste0("chr", input_df[,5])
  }

  ##test all chrs are in each set of regions
  region_1 <- data.frame(chr = input_df[,"chrom1"],
                      start = input_df[,"start1"],
                      end = input_df[,"end1"])
  region_2 <- data.frame(region = input_df[,"chrom2"],
                      start = input_df[,"start2"],
                      end = input_df[,"end2"])
  region_1_in <- region_1[,1] %in% cytoband[,1]
  region_2_in <- region_2[,1] %in% cytoband[,1]
  region_1 <- region_1[region_1_in & region_2_in,]
  region_2 <- region_2[region_1_in & region_2_in,]
  region_c <- input_df[region_1_in & region_2_in, "colour"]
  log10quals <- log10(input_df[region_1_in & region_2_in, "qualscore"])
  labels_o1 <- input_df[region_1_in & region_2_in, c("chrom1", "start1", "end1", "symbol1", "colour")]
  labels_o2 <- input_df[region_1_in & region_2_in, c("chrom2", "start2", "end2", "symbol2", "colour")]
  colnames(labels_o1) <- colnames(labels_o2) <- c("chr1", "start", "end", "symbol", "colour")
  labels_o <- rbind(labels_o1, labels_o2)

  return(list(input_df = input_df,
              region_1 = region_1,
              region_2 = region_2,
              region_c = region_c,
              log10quals = log10quals,
              labels_o1 = labels_o1,
              labels_o2 = labels_o2,
              labels_o = labels_o))
}
