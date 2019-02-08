test_that("can validate niterations", {

    expect_error(validate_niterations("wer"))
    expect_error(validate_niterations("1"))
    expect_error(validate_niterations(1.1))
    expect_error(validate_niterations(-1))
    expect_error(validate_niterations(0))

    
    expect_null(validate_niterations(1))
    expect_null(validate_niterations(10))
    expect_null(validate_niterations(300))    

    
})

test_that("add required entries output VCF header", {

    command_line <- "SEW(somecommand = 1)"
    sampleNames <- c("samp1", "samp2")
    strandedness <- NULL
    output_vcf_header_file <- tempfile()
    
    write_sew_vcf_header(
        output_vcf_header_file = output_vcf_header_file,
        sampleNames = sampleNames,
        strandedness = strandedness,
        command_line = command_line
    )
    header <- readLines(output_vcf_header_file)

    ## check commands
    expect_equal(
        sum(paste0("##SEW_command=", command_line) == header),
        1
    )
    expect_equal(
        length(grep("##SEW_version", header)) ,
        1
    )
    
})


test_that("C++ version of single_sample_expectation returns the same as the R version", {

    
    sampleReads <- list(
        list(
            1, 0,
            matrix(c(-7, -8), ncol = 1),
            matrix(c(0, 1), ncol = 1)            
        ),
        list(
            0, 0,
            matrix(c(-9), ncol = 1),
            matrix(c(2), ncol = 1)            
        ),
        list(
            2, 0,
            matrix(c(13, 14, 15), ncol = 1),
            matrix(c(0, 1, 2), ncol = 1)            
        )
    )
    N_r <- 3
    nSNPs <- 3
    K <- 2
    eHapsCurrent <- matrix(0, nrow = nSNPs, ncol = K)
    eHapsCurrent[, 1] <- rep(0.1, nSNPs)
    eHapsCurrent[, 2] <- rep(0.9, nSNPs)    
    
    out1 <- single_sample_expectation(
        sampleReads,
        eHapsCurrent,
        nSNPs,
        N_r,
        calculate_updates = TRUE
    )

    out2 <- cpp_single_sample_expectation(
        sampleReads,
        eHapsCurrent,
        nSNPs,
        N_r,
        calculate_updates = TRUE
    )

    expect_equal(out1$log_likelihood, out2$log_likelihood)
    expect_equal(out1$p_reads, out2$p_reads)
    expect_equal(out1$p_reads_given_hap_k, out2$p_reads_given_hap_k)    
    expect_equal(out1$eHapsUpdate_numer, out2$eHapsUpdate_numer)
    expect_equal(out1$eHapsUpdate_denom, out2$eHapsUpdate_denom)
    
})


test_that("both R and C++ single_sample_expectation can use t_break", {

    sampleReads <- list(
        list(
            1, 0,
            matrix(c(-7, -8), ncol = 1),
            matrix(c(0, 1), ncol = 1)            
        ),
        list(
            0, 0,
            matrix(c(-9), ncol = 1),
            matrix(c(2), ncol = 1)            
        ),
        list(
            2, 0,
            matrix(c(13, 14, 15), ncol = 1),
            matrix(c(0, 1, 2), ncol = 1)            
        )
    )
    N_r <- 3
    nSNPs <- 3
    K <- 2
    eHapsCurrent <- matrix(0, nrow = nSNPs, ncol = K)
    eHapsCurrent[, 1] <- rep(0.1, nSNPs)
    eHapsCurrent[, 2] <- rep(0.9, nSNPs)

    t_break <- 1    
    out1 <- single_sample_expectation(
        sampleReads,
        eHapsCurrent,
        nSNPs,
        N_r,
        calculate_updates = TRUE,
        t_break = 1
    )

    out1b <- cpp_single_sample_expectation(
        sampleReads,
        eHapsCurrent,
        nSNPs,
        N_r,
        calculate_updates = TRUE,
        t_break = 1
    )

    expect_equal(out1, out1b)

    ## now - manually do it
    eHapsCurrent[(t_break + 1):nSNPs, 1:2] <- eHapsCurrent[(t_break + 1):nSNPs, 2:1]
    out2 <- single_sample_expectation(
        sampleReads,
        eHapsCurrent,
        nSNPs,
        N_r,
        calculate_updates = TRUE
    )
    
    expect_equal(out1, out2)

})

test_that("C++ version of single_sample_expectation can do eHapsCurrent breaks on the fly", {

    sampleReads <- list(
        list(
            1, 0,
            matrix(c(-7, -8), ncol = 1),
            matrix(c(0, 1), ncol = 1)            
        ),
        list(
            0, 0,
            matrix(c(-9), ncol = 1),
            matrix(c(2), ncol = 1)            
        ),
        list(
            2, 0,
            matrix(c(13, 14, 15), ncol = 1),
            matrix(c(0, 1, 2), ncol = 1)            
        )
    )
    N_r <- 3
    nSNPs <- 3
    K <- 2
    eHapsCurrent <- matrix(0, nrow = nSNPs, ncol = K)
    eHapsCurrent[, 1] <- rep(0.1, nSNPs)
    eHapsCurrent[, 2] <- rep(0.9, nSNPs)    

    t_break <- 1    
    out1 <- single_sample_expectation(
        sampleReads,
        eHapsCurrent,
        nSNPs,
        N_r,
        calculate_updates = FALSE,
        t_break = 1
    )

    out2 <- cpp_single_sample_expectation(
        sampleReads,
        eHapsCurrent,
        nSNPs,
        N_r,
        calculate_updates = FALSE,
        t_break = 1
    )
    
    expect_equal(out1$log_likelihood, out2$log_likelihood)
    
})


test_that("single_sample_expectation can use t_break with t_writing_offset", {

    sampleReads <- list(
        list(
            0, 0,
            matrix(c(-30), ncol = 1),
            matrix(c(0), ncol = 1)            
        ),
        list(
            3, 0,
            matrix(c(-29, -28, -27, -26), ncol = 1),
            matrix(c(1, 2, 3, 4), ncol = 1)            
        ),
        list(
            3, 0,
            matrix(c(31, 32, 33, 34), ncol = 1),
            matrix(c(1, 2, 3, 4), ncol = 1)            
        )
    )
    N_r <- 3
    nSNPs <- 5
    K <- 2
    eHapsCurrent <- matrix(0, nrow = nSNPs, ncol = K)
    eHapsCurrent[, 1] <- rep(0.1, nSNPs)
    eHapsCurrent[, 2] <- rep(0.9, nSNPs)
    eHapsCurrent[4:nSNPs, 1:2] <- eHapsCurrent[4:nSNPs, 2:1]

    t_break <- 3 ## i.e. investigate breaking between SNPs 2 and 3
    which_reads <- 2:3

    t_reading_offset <- 0 
    t_writing_offset <- 1
    
    out1 <- single_sample_expectation(
        sampleReads[which_reads],
        eHapsCurrent,
        nSNPs = nSNPs - t_writing_offset, ## want output to be smaller
        N_r = 2,
        calculate_updates = TRUE,
        t_break = t_break,
        t_reading_offset = t_reading_offset,
        t_writing_offset = t_writing_offset
    )

    eHapsFuture <- make_new_eHaps(
        out1$eHapsUpdate_numer,
        out1$eHapsUpdate_denom
    )

    ## new col1 should only have lower bound
    expect_equal(
        sum(eHapsFuture[, 1] != 0.001),
        0
    )
    ## new col2 should only have upper bound
    expect_equal(
        sum(eHapsFuture[, 2] != 0.999),
        0
    )

    ## check C++ function does the same thing
    out2 <- cpp_single_sample_expectation(
        sampleReads[which_reads],
        eHapsCurrent,
        nSNPs = nSNPs - t_writing_offset, ## want output to be smaller
        N_r = 2,
        calculate_updates = TRUE,
        t_break = t_break,
        t_reading_offset = t_reading_offset,
        t_writing_offset = t_writing_offset
    )
    
    expect_equal(out1, out2)
    

})


test_that("single_sample_expectation can use t_break with t_writing_offset and t_reading_offset", {

    sampleReads <- list(
        list(
            0, 0,
            matrix(c(-30), ncol = 1),
            matrix(c(0), ncol = 1)            
        ),
        list(
            3, 0,
            matrix(c(-28, -27, -26, -25), ncol = 1),
            matrix(c(2, 3, 4, 5), ncol = 1)            
        ),
        list(
            3, 0,
            matrix(c(32, 33, 34, 35), ncol = 1),
            matrix(c(2, 3, 4, 5), ncol = 1)            
        )
    )
    N_r <- 3
    nSNPs <- 6
    K <- 2
    eHapsCurrent <- matrix(0, nrow = nSNPs, ncol = K)
    eHapsCurrent[, 1] <- rep(0.1, nSNPs)
    eHapsCurrent[, 2] <- rep(0.9, nSNPs)
    eHapsCurrent[4:nSNPs, 1:2] <- eHapsCurrent[4:nSNPs, 2:1]

    t_break <- 3 ## i.e. investigate breaking between SNPs 2 and 3
    which_reads <- 2:3

    t_reading_offset <- 2 ## i.e. have smaller eHapsCurrent
    t_writing_offset <- 2 ## this is not really a difference?

    ## have removed t_reading_offset bases
    eHapsCurrent_local <- eHapsCurrent[- (1:t_reading_offset), ]
    out1 <- single_sample_expectation(
        sampleReads[which_reads],
        eHapsCurrent_local,
        nSNPs = nSNPs - t_writing_offset, ## want output to be smaller
        N_r = 2,
        calculate_updates = TRUE,
        t_break = t_break,
        t_reading_offset = t_reading_offset,
        t_writing_offset = t_writing_offset
    )

    eHapsFuture <- make_new_eHaps(
        out1$eHapsUpdate_numer,
        out1$eHapsUpdate_denom
    )

    ## new col1 should only have lower bound
    expect_equal(
        sum(eHapsFuture[, 1] != 0.001),
        0
    )
    ## new col2 should only have upper bound
    expect_equal(
        sum(eHapsFuture[, 2] != 0.999),
        0
    )

    ## check C++ function does the same thing
    out2 <- cpp_single_sample_expectation(
        sampleReads[which_reads],
        eHapsCurrent_local,
        nSNPs = nSNPs - t_writing_offset, ## want output to be smaller
        N_r = 2,
        calculate_updates = TRUE,
        t_break = t_break,
        t_reading_offset = t_reading_offset,
        t_writing_offset = t_writing_offset
    )
    
    expect_equal(out1, out2)

})


test_that("investigate_an_unwinding can work for 2+ iterations", {

    sampleReads <- list(
        list(
            0, 0,
            matrix(c(-30), ncol = 1),
            matrix(c(0), ncol = 1)            
        ),
        list(
            2, 0,
            matrix(c(-28, -27, -26), ncol = 1),
            matrix(c(2, 3, 4), ncol = 1)            
        ),
        list(
            2, 0,
            matrix(c(32, 33, 34), ncol = 1),
            matrix(c(2, 3, 4), ncol = 1)            
        )
    )
    N_r <- 3
    nSNPs <- 5
    K <- 2
    eHapsCurrent <- matrix(0, nrow = nSNPs, ncol = K)
    eHapsCurrent[, 1] <- rep(0.1, nSNPs)
    eHapsCurrent[, 2] <- rep(0.9, nSNPs)
    eHapsCurrent[4:nSNPs, 1:2] <- eHapsCurrent[4:nSNPs, 2:1]

    which_reads <- 2:3
    t_break <- 3

    ## just useful to look at
    reads_span <- get_what_reads_span(
        sampleReads = sampleReads,
        nSNPs = nSNPs,
        unwindIterations = 1
    )
    first_snp <- reads_span[t_break, "SNPstart"]
    last_snp <- reads_span[t_break, "SNPend"]    

    ## set really high to make irrelevant
    out_no <- investigate_an_unwinding(
        sampleReads,
        eHapsCurrent,
        nSNPs,
        t_break = 100,
        which_reads,
        n_unwinding_its = 1,
        first_snp = first_snp,
        last_snp = last_snp
    )
    
    out1 <- investigate_an_unwinding(
        sampleReads,
        eHapsCurrent,
        nSNPs,
        t_break = t_break,
        which_reads,
        n_unwinding_its = 1,
        first_snp = first_snp,
        last_snp = last_snp        
    )

    out2 <- investigate_an_unwinding(
        sampleReads,
        eHapsCurrent,
        nSNPs,
        t_break,
        which_reads,
        n_unwinding_its = 2,
        first_snp = first_snp,
        last_snp = last_snp        
    )

    expect_equal(out1 < out_no, TRUE )
    expect_equal(out2 < out1, TRUE )    

    ## check that the same thing as doing whole thing
    out2NA <- investigate_an_unwinding(
        sampleReads,
        eHapsCurrent,
        nSNPs,
        t_break = t_break,
        which_reads,
        n_unwinding_its = 2,
        first_snp = NA,
        last_snp = NA
    )
    
    expect_equal(out2, out2NA)

})


test_that("can calculate PS", {
    
    sampleReads <- list(
        list(
            0, 0,
            matrix(c(-30), ncol = 1),
            matrix(c(0), ncol = 1)            
        ),
        list(
            3, 0,
            matrix(c(-29, -28, -27, -26), ncol = 1),
            matrix(c(1, 2, 3, 4), ncol = 1)            
        ),
        list(
            3, 0,
            matrix(c(31, 32, 33, 34), ncol = 1),
            matrix(c(1, 2, 3, 4), ncol = 1)            
        )
    )
    L <- c(2, 5, 10, 20, 50)
    inputdir <- tempdir()
    regionName <- "test"
    save(sampleReads, file = STITCH::file_sampleReads(inputdir, 1, regionName))
    phase_set <- calculate_phase_set(inputdir, regionName, L)
    expect_equal(phase_set, c(2, 5, 5, 5, 5))

})


test_that("can calculate more complicated PS", {
    
    sampleReads <- list(
        list(
            3, 0,
            matrix(c(-29, -28, -27, -26), ncol = 1),
            matrix(c(1, 2, 3, 4), ncol = 1)            
        ),
        list(
            3, 0,
            matrix(c(31, 32, 33, 34), ncol = 1),
            matrix(c(1, 2, 3, 4), ncol = 1)            
        ),
        list(
            1, 0,
            matrix(c(31, 32), ncol = 1),
            matrix(c(3, 4), ncol = 1)            
        ),
        list(
            1, 0,
            matrix(c(31, 32), ncol = 1),
            matrix(c(5, 6), ncol = 1)            
        )
    )
    L <- c(2, 5, 10, 20, 50, 100, 200)
    inputdir <- tempdir()
    regionName <- "test"
    save(sampleReads, file = STITCH::file_sampleReads(inputdir, 1, regionName))
    phase_set <- calculate_phase_set(inputdir, regionName, L)
    expect_equal(phase_set, c(2, 5, 5, 5, 5, 100, 100))

})

test_that("can remove buffer from PS", {
    phase_set <- c(2, 2, 2, 2, 2, 10, 10, 10, 10, 10)
    L <- c(2, 3, 4, 5, 6, 10, 11, 12, 13, 14)
    inRegion2 <- c(FALSE, FALSE, FALSE, TRUE, TRUE, TRUE, TRUE, TRUE, FALSE, FALSE)
    expect_equal(
        remove_buffer_from_phase_set(phase_set, L, inRegion2),
        c(5, 5, 10, 10, 10)
    )
})
