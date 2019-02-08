#!/usr/bin/env Rscript

## change directory to one up from scripts, no matter how this was called
args <- commandArgs(trailingOnly = FALSE)
for(key in c("--file=", "--f=")) {
    i <- substr(args, 1, nchar(key)) == key
    if (sum(i) == 1) {
        script_dir <- dirname(substr(args[i], nchar(key) + 1, 1000))
        setwd(file.path(script_dir, "../"))
    }
}
Sys.setenv(PATH = paste0(Sys.getenv("PATH"), ":", getwd()))

library("testthat")
library("parallel")
library("STITCH") ## for make_STITCH_cli

## testthat doesn't do what I want outside of package form
## so don't bother wrappping, just fail

cli_function_build <- Sys.getenv("CLI_FUNCTION_BUILD")
if (cli_function_build != "") {
    print(paste0("Using ", cli_function_build))
    dir <- tempdir()
    system(paste0("cp ", cli_function_build, " ", dir, "/"))
    system(paste0("(cd ", dir, " && tar -zxvf ", dir, "/*tar.gz SEW/R/functions.R)"))
    function_file <- file.path(dir, "STITCH/R/functions.R")
} else {
    function_file <- "SEW/R/sew.R"
}

cli_output_file <- "SEW.R"
make_STITCH_cli(
    function_file = "SEW/R/sew.R",
    cli_output_file = cli_output_file,
    integer_vectors = c("unwindIterations", "sample_unwindIterations"),
    other_character_params = c("phasefile", "unwindIterations", "sample_unwindIterations"),
    other_logical_params = c("outputHaplotypeProbabilities", "very_verbose"),
    function_name = "SEW",
    library_name = "SEW"
)
system(paste0("chmod +x ", cli_output_file))

message("test that SEW CLI produces help message")
## behaviour of optparse changed!
## now exits code 0 as one would hope on --help
stdout_file <- tempfile()
stderr_file <- tempfile()
out <- system2(
    cli_output_file,
    args = c(
        "--help"
    ),
    stdout = stdout_file, stderr = stderr_file
)
stderr <- system(paste0("cat ", stderr_file), intern = TRUE)
stdout <- system(paste0("cat ", stdout_file), intern = TRUE)
if (out > 0) {
    message("---stderr---")
    print(stderr)
    message("---stdout---")
    print(stdout)
}
expect_equal(0, out)



##
## make some test data
##
n_snps <- 10
chr <- 10
phasemaster <- matrix(c(rep(0, n_snps), rep(1, n_snps)), ncol = 2)
data_package <- STITCH::make_acceptance_test_data_package(
    reads_span_n_snps = 3,
    n_samples = 1,
    n_snps = n_snps,
    n_reads = n_snps * 10,
    seed = 2,
    chr = chr,
    K = 2,
    phasemaster = phasemaster
)



## this also test character, integer, double and NA variables
message("test that SEW CLI can work")
stdout_file <- tempfile()
stderr_file <- tempfile()
out <- system2(
    cli_output_file,
    args = c(
        paste0("--chr=", data_package$chr),
        paste0("--bamlist=", data_package$bamlist),
        paste0("--posfile=", data_package$posfile),
        paste0("--outputdir=", data_package$outputdir),
        paste0("--phasefile=", data_package$phasefile),
        "--niterations=25",
        "--very_verbose=TRUE"
    ),
    stdout = stdout_file, stderr = stderr_file
)
stderr <- system(paste0("cat ", stderr_file), intern = TRUE)
stdout <- system(paste0("cat ", stdout_file), intern = TRUE)
if (out > 0) {
    message("---stderr---")
    print(stderr)
    message("---stdout---")
    print(stdout)
}
expect_equal(0, out)
