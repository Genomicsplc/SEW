check_sew_output <- function(
    file,
    data_package,
    which_snps = NULL
) {
    vcf <- read.table(
        file,
        header = FALSE,
        stringsAsFactors = FALSE
    )
    check_sew_phase(
        vcf,
        data_package$phase,
        which_snps
    )
}


## here, we check phasing is exact
check_sew_phase <- function(vcf, phase, which_snps = NULL) {
    ## GT is first column
    gt <- sapply(strsplit(vcf[, 10], ":"), function(x) x[1])
    gt1 <- substr(gt, 1, 1)
    gt2 <- substr(gt, 3, 3)
    ##
    if (is.null(which_snps)) {
        t1 <- phase[, 1, 1]
        t2 <- phase[, 1, 2]
    } else {
        t1 <- phase[which_snps, 1, 1]
        t2 <- phase[which_snps, 1, 2]
    }
    ## basically, need
    hap1_check <- (sum(gt1 != t1) == 0) | (sum(gt2 != t1) == 0)
    hap2_check <- (sum(gt1 != t2) == 0) | (sum(gt2 != t2) == 0)
    if (!hap1_check | !hap2_check) {
        print(paste0("phasing results are:", gt))
        print(paste0("truth is:", apply(phase[, 1, ], 1, paste, collapse = "|")))
    }
    expect_equal(hap1_check, TRUE)
    expect_equal(hap2_check, TRUE)    
}
