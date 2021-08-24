if ( 1 == 0 ) {

    library("testthat")
    library("SEW")
    dir <- "~/proj/SEW/"
    setwd(paste0(dir, "/SEW/R"))
    a <- dir(pattern = "*R")
    b <- grep("~", a)
    if (length(b) > 0) {
        a <- a[-b]
    }
    o <- sapply(a, source)


}

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



test_that("SEW can phase a sample at high coverage with long reads", {

    outputdir <- STITCH::make_unique_tempdir()    
    SEW(
        chr = data_package$chr,
        bamlist = data_package$bamlist,
        posfile = data_package$posfile,
        outputdir = outputdir,
        phasefile = data_package$phasefile,
        niterations = 25,
        very_verbose = TRUE
    )

    check_sew_output(
        file = file.path(outputdir, paste0("sew.", data_package$chr, ".vcf.gz")),
        data_package
    )

})


test_that("SEW can run without heuristics", {

    outputdir <- STITCH::make_unique_tempdir()    
    SEW(
        chr = data_package$chr,
        bamlist = data_package$bamlist,
        posfile = data_package$posfile,
        outputdir = outputdir,
        phasefile = data_package$phasefile,
        niterations = 25,
        very_verbose = TRUE,
        disable_heuristics = TRUE
    )

    check_sew_output(
        file = file.path(outputdir, paste0("sew.", data_package$chr, ".vcf.gz")),
        data_package
    )

})


test_that("SEW can phase a sample at high coverage with long reads with defined regions", {

    outputdir <- STITCH::make_unique_tempdir()
    regionStart <- 5
    regionEnd <- 7
    SEW(
        chr = data_package$chr,
        bamlist = data_package$bamlist,
        posfile = data_package$posfile,
        outputdir = outputdir,
        phasefile = data_package$phasefile,
        regionStart = regionStart,
        regionEnd = regionEnd,
        buffer = 1
    )

    regionName <- paste0(data_package$chr, ".", regionStart, ".", regionEnd)
    check_sew_output(
        file = file.path(outputdir, paste0("sew.", regionName, ".vcf.gz")),
        data_package,
        which_snps = regionStart:regionEnd
    )

})


test_that("SEW can phase a moderate sample with high bq error with long reads", {

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
    
    ## make
    data_package_local <- STITCH::make_acceptance_test_data_package(
        reads_span_n_snps = reads_span_n_snps,
        n_samples = 1,
        n_snps = n_snps,
        n_reads = n_reads,
        seed = 2,
        chr = chr,
        K = 2,
        phasemaster = phasemaster,
        L = L,
        phred_bq = phred_bq
    )
    ## tmpdir = getwd(),
    
    outputdir <- STITCH::make_unique_tempdir()    
    SEW(
        chr = data_package_local$chr,
        bamlist = data_package_local$bamlist,
        posfile = data_package_local$posfile,
        outputdir = outputdir,
        phasefile = data_package_local$phasefile
    )

    check_sew_output(
        file = file.path(outputdir, paste0("sew.", data_package_local$chr, ".vcf.gz")),
        data_package_local
    )

    ## set to TRUE to build example package
    if (FALSE) {

        setwd("~/Google Drive/SEW")
        unlink("example_package", recursive = TRUE)
        dir.create("example_package")
        
        ## get bam
        system(paste0("mv ", shQuote(data_package_local$bam_files), " example_package"))
        write.table(basename(data_package_local$bam_files), file = "example_package/bamlist.txt", quote = FALSE, col.names = FALSE, row.names = FALSE)
        
        ## get sites, truth
        system(paste0("mv ", shQuote(data_package_local$posfile), " example_package"))
        system(paste0("mv ", shQuote(data_package_local$phasefile), " example_package"))
        
        unlink(data_package_local$outputdir, recursive = TRUE)
        
        ## package
        setwd("example_package")
        system("samtools index *.bam")
        system("tar -czvf ../SEW_example_2019_01_11.tgz *")
        ## change permissions, push to florence
        setwd("../")
        system("chmod 755 SEW_example_2019_01_11.tgz")
        system("rsync -av SEW_example_2019_01_11.tgz florence:~/public_html/ancillary")

    }
    

})
