#' GRIDSS functions
#'
#' Parse multi-sample VCF (recommended output of GRIDSS)
#'
#' @param vcf, the VCF file path
#' @param which_genome, genome assembly used ("hg19", "hg38")
#' @return a list with recurrent SVs and private SVs as elements
#' @export

gridss_parse_multi_vcf <- function(vcf, which_genome = NULL){

  options(stringsAsFactors = FALSE)

  ##set genome as hg19 unless spec'd
  if (is.null(which_genome)){
    which_genome <- "hg19"
  }

  ##set samplenames
  sample_names <- VariantAnnotation::samples(VariantAnnotation::scanVcfHeader(vcf))

  ##readVcf and get per-sample GRanges
  rvcf <- VariantAnnotation::readVcf(vcf, genome = which_genome)
  bpgr <- StructuralVariantAnnotation::breakpointRanges(rvcf)
  #pairs <- StructuralVariantAnnotation::breakpointgr2pairs(gr)

  ##iterate over sample_names
  ##https://github.com/PapenfussLab/StructuralVariantAnnotation/issues/31
  bpgr_list <- lapply(sample_names, function(samp) {

    ##get sample-specific bpgr object
    bpgr[VariantAnnotation::geno(rvcf[bpgr$sourceId])$QUAL[, samp] >= 2]
  })

  names(bpgr_list) <- sample_names

  ##make single bedpe for all VCFs
  bedpe_rb <- do.call(rbind, bedpe_list_filter)

  ##find reccuring events (same at 1,2,3,4,5,6,7,8)
  bedpe_rb_rec <- bedpe_rb %>%
                  dplyr::mutate_if(is.factor,as.numeric) %>%
                  group_by_at(vars(-qualscore, -sampleID)) %>%
                  mutate(num_rows = sum(n())) %>%
                  filter(num_rows > 1) %>%
                  ungroup() %>%
                  dplyr::select(-qualscore, -sampleID, -num_rows) %>%
                  distinct()

  ##lapply to paste sampleIDs, colour by combinations therein
  sample_rec_list <- lapply(1:dim(bedpe_rb_rec)[1], function(f){
    left_join(bedpe_rb_rec[f,],
              bedpe_rb,
              by=c("start1" = "start1",
                   "end1" = "end1",
                   "start2" = "start2",
                   "end2" = "end2")) %>%
    dplyr::select(sampleID) %>%
    unlist() %>%
    sort() %>%
    paste(., collapse=",")
  })

  ##add to bedpe_rb_rec
  bedpe_rb_rec$sampleID <- unlist(sample_rec_list)

  ##colour in based on samples with SV
  colz_rec <- rainbow(n = length(table(unlist(sample_rec_list))))
  names(colz_rec) <- names(table(unlist(sample_rec_list)))
  bedpe_rb_rec$colour <- colz_rec[match(bedpe_rb_rec$sampleID, names(colz_rec))]

  ##non-recurrent
  bedpe_rb_pri <- bedpe_rb %>%
                  dplyr::mutate_if(is.factor,as.numeric) %>%
                  group_by_at(vars(-qualscore, -sampleID)) %>%
                  mutate(num_rows = sum(n())) %>%
                  filter(num_rows == 1) %>%
                  ungroup() %>%
                  dplyr::select(-num_rows)
  colz_pri <- rainbow(n = length(table(unlist(bedpe_rb_pri$sampleID))))
  names(colz_pri) <- names(table(unlist(bedpe_rb_pri$sampleID)))
  bedpe_rb_pri$colour <- colz_pri[match(bedpe_rb_pri$sampleID, names(colz_pri))]

  return(list(bedpe_rb_rec, bedpe_rb_pri))
}

#' Produce a list of bedpe, from a list of breakpoint GRanges,
#'    with score > qual_filt and both breakpoints within annotated genes
#' @param bpgr_list list of per-sample breakpont GRanges objects
#' @param qual_filt qulaity score filter above which is returned
#' @return filtered bedpe list
#' @export

gridss_bedpe_list <- function(bpgr_list, qual_filt) {

  bedpe_list_filter <- lapply(seq_along(bpgr_list), function(f){

    print(paste0("Working on: ", names(bpgr_list)[f]))
    gr <- bpgr_list[[f]]

    ##remove non-found partners!
    pairs <- StructuralVariantAnnotation::breakpointgr2pairs(gr)
    grp <- gr[names(gr) %in% mcols(pairs)[,"name"]]
    gr_partner <- gr[names(gr) %in% grp$partner,]

    ##force into only 1:22, X
    GenomeInfoDb::seqlevels(grp, pruning.mode="coarse") <- GenomeInfoDb::seqlevels(grp)[c(1:23)]

    ##annotate
    grp_anno <- gridss_annotate_gr(grp, which_genome)
    gr_partner_anno <- gridss_annotate_gr(gr_partner, which_genome)
    grp_anno_use <- grp_anno[names(grp_anno) %in% gr_partner_anno$partner]
    gr_partner_anno_use <- gr_partner_anno[names(gr_partner_anno) %in% grp_anno_use$partner]

    #bedpe for input, simplest format
    bedpe1k <- data.frame(chrom1 = GenomeInfoDb::seqnames(grp_anno_use),
                          start1 = start(grp_anno_use),
                          end1 = end(grp_anno_use),
                          symbol1 = grp_anno_use$SYMBOL,
                          chrom2 = GenomeInfoDb::seqnames(gr_partner_anno_use),
                          start2 = start(gr_partner_anno_use),
                          end2 = end(gr_partner_anno_use),
                          symbol2 = gr_partner_anno_use$SYMBOL,
                          qualscore = grp_anno_use$QUAL,
                          sampleID = names(vcf_list)[f]) %>%
                dplyr::filter(! symbol1 %in% "",! symbol2 %in% "")

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

  ##remove chr if no chr
  if(! "chr1" %in% levels(GenomeInfoDb::seqnames(gr))){
    gns <- GenomeInfoDb::renameSeqlevels(gns, sub("chr", "", GenomeInfoDb::seqlevels(gns)))
  }

  hits <- as.data.frame(GenomicRanges::findOverlaps(gr,
                                                    gns,
                                                    ignore.strand=TRUE))
  hits$SYMBOL <- suppressMessages(biomaRt::select(org.Hs.eg.db::org.Hs.eg.db,
                              gns[hits$subjectHits]$gene_id, "SYMBOL")$SYMBOL)
  hits$gene_strand <- as.character(strand(gns[hits$subjectHits]))

  ##add annotation to grp
  gr$SYMBOL <- gr$geneStrand <- ""
  gr$SYMBOL[hits$queryHits] <- hits$SYMBOL
  gr$geneStrand[hits$queryHits] <- hits$gene_strand

  # require the breakpoint to be between different genes
  # grp <- grp[grp$SYMBOL != partner(grp)$SYMBOL & !is.na(grp$SYMBOL) & !is.na(partner(grp)$SYMBOL),]
  # grp <- grp[grp$SYMBOL != "" & partner(grp)$SYMBOL != "" & !is.na(grp$SYMBOL) & !is.na(partner(grp)$SYMBOL),]

  # require the breakpoint to possibly generate a fusion transcript
  # grp$couldBeThreePrimeStart <- stringr::str_detect(grp$geneStrand,
  #                                   stringr::fixed(as.character(strand(grp))))
  # grp$couldBeFivePrimeEnd <- stringr::str_detect(grp$geneStrand,
  #            stringr::fixed(ifelse(as.character(strand(grp))=="+", "-", "+")))
  # grp <- grp[(grp$couldBeThreePrimeStart & partner(grp)$couldBeFivePrimeEnd) | (grp$couldBeFivePrimeEnd & partner(grp)$couldBeThreePrimeStart),]
  return(gr)
}

#' OmicCircos plotting from data.frame from GRIDSS bedpe
#' @param input_df data.frame with chrom1, start1, end1, symbol1,
#'   chrom2, start2, end2, symbol2, sampleID, colour columns
#' @param which_genome, genome assembly used ("hg19", "hg38")
#' @param pdf_path path to where PDF file with output plot is written
#'    N.B. two pages, one with Circos plot and one with legend
#' @return none, prints plot to file
#' @export

plot_circos_sv <- function(input_df, which_genome, pdf_path){

  if(!is.data.frame(input_df)){
    print("Require data.frame input with chromsomes as factors")
  } else {

    ##impose condition that over 200 SVs in input_df
    ##plot only those that are on different chromosomes
    if(dim(input_df)[1]>100){
      print(paste0("Total of ", dim(input_df)[1], " SV found, reducing to those on different chromosomes"))
      input_df <- input_df[as.vector(input_df$chrom1) != as.vector(input_df$chrom2),]
      print(paste0("Now using with ", dim(input_df)[1], " SV"))
    }

    ##create labelling coordinates, colouring
    label_df <- data.frame(chrom = c(unlist(input_df[,"chrom1"]), unlist(input_df[,"chrom2"])),
                          start = c(unlist(input_df[,"start1"]), unlist(input_df[,"start2"])),
                          gene = c(unlist(input_df[,"symbol1"]), unlist(input_df[,"symbol2"])),
                          colour = c(unlist(input_df[,"colour"]), unlist(input_df[,"colour"])))

    ##actual OmicCircos plotting function
    pcfoc_func <- function(input_df, which_genome, label_df) {
      plot(c(1,800) , c(1,800) , type="n", axes=FALSE, xlab="", ylab="")
      OmicCircos::circos(R = 248,
                         cir = which_genome,
                         type="chr",
                         print.chr.lab = TRUE,
                         W = 5,
                         scale = F)
      OmicCircos::circos(R = 240,
                         cir = which_genome,
                         W = 50,
                         mapping = input_df[,c(1:3,5:7)],
                         type = "link.pg",
                         col = as.vector(input_df$colour))
      OmicCircos::circos(R = 255,
                         cir = which_genome,
                         W = 50,
                         mapping = label_df,
                         type = "label",
                         col = "black",
                         side = "out",
                         cex = .4)
     }

     ##colours per-sample
     nms <- names(table(input_df$sampleID))
     colz <- unlist(lapply(nms, function(f) {
       input_df$colour[match(f, input_df$sampleID)]
     }))

     pdf(pdf_path, width = 9, height = 9)
     par(mar=c(2,2,2,2))
     pcfoc_func(input_df, which_genome, label_df)
     grid::grid.newpage()
     legend(200, 500,
            legend = nms,
            col = colz,
            lty=1, lwd=2)
     dev.off()
  }
}
