library(testthat)
library(Maaslin)


expected_maaslin_results = read.table("maaslin_results.txt",stringsAsFactors=FALSE)
setwd("..") 
maaslin_results <- Maaslin('Maaslin/extdata/maaslin_demo2.tsv','output',strInputConfig='Maaslin/extdata/maaslin_demo2.read.config')
expect_that(expected_maaslin_results$Variable,equals(as.character(maaslin_results$Variable)))
expect_that(expected_maaslin_results$Feature,equals(as.character(maaslin_results$Feature)))
expect_that(expected_maaslin_results$Value,equals(as.character(maaslin_results$Value)))
expect_that(expected_maaslin_results$N,equals(maaslin_results$N))
expect_that(expected_maaslin_results$N.not.0,equals(maaslin_results$N.not.0))
expect_that(expected_maaslin_results$P.value,equals(maaslin_results$P.value))
expect_that(expected_maaslin_results$Q.value,equals(maaslin_results$Q.value))
