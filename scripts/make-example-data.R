## this is a one-off script to make example data for quick starts
## todo, in the long term, replace with real data

library("testthat")
library("SEW")

set.seed(130)
reads_span_n_snps <- 5
n_snps <- 50
L <- sort(sample(1:10000, n_snps))
n_reads <- n_snps * 10
phred_bq <- 8 ## make less accurate
chr <- 10
phasemaster <- matrix(0, nrow = n_snps, ncol = 2)
phasemaster[, 1] <- sample(c(0, 1), n_snps, replace = TRUE)
phasemaster[, 2] <- 1 - phasemaster[, 1]

##
source("~/Google Drive/STITCH/STITCH/R/test-drivers.R")
data_package <- make_acceptance_test_data_package(
    reads_span_n_snps = reads_span_n_snps,
    n_samples = 1,
    n_snps = n_snps,
    n_reads = n_reads,
    seed = 2,
    chr = chr,
    K = 2,
    phasemaster = phasemaster,
    L = L,
    tmpdir = getwd(),
    phred_bq = phred_bq
)

outputdir <- STITCH::make_unique_tempdir()    
SEW(
    chr = data_package$chr,
    bamlist = data_package$bamlist,
    posfile = data_package$posfile,
    outputdir = outputdir,
    phasefile = data_package$phasefile,
    very_verbose = TRUE
)


## now manually tar up package

## manually check using README code (not shown!)

