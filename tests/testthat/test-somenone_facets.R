testthat::test_that("facets_cna_consensus", {
  fcc_out <- somenone::facets_cna_consensus(
              cncf_match = "test_ex_soma.fit_cncf_jointsegs.tsv",
              pp_match = "test_ex_soma.fit_ploidy_purity.tsv",
              dict_file = "GRCh38_Verily_v1.genome.noChr.dict",
              tag = "testthat_testdata",
              work_dir = system.file(package = "somenone", "extdata"),
              cgc_bed = "cancer_gene_census.bed")
  testthat::expect_s3_class(fcc_out, "data.frame")
})
