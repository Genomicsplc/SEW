validate_niterations <- function(niterations) {
    if (
    (class(niterations) != "integer") &
    (class(niterations) != "numeric")
    ) {
        stop("reference_iterations must be an integer of at least 1")
    }
    if (round(niterations) != niterations) {
        stop("reference_iterations must be an integer of at least 1")        
    }
    if (niterations < 1) {
        stop("reference_iterations must be an integer of at least 1")                
    }
    return(NULL)
}

validate_unwindIterations <- function(unwindIterations, sample_unwindIterations, usePhaseSet_unwindIterations) {
    if (sum(is.na(unwindIterations)) > 0)
        stop(paste0(
            "unwindIterations must contain no NA but supplied is:",
            paste0(unwindIterations, collapse = ",")
        ))
    if (sum(unwindIterations < 0) > 0)
        stop(paste0("unwindIterations values must all be greater than 0"))
    ##
    if (sum(is.na(sample_unwindIterations)) > 0)
        stop(paste0(
            "sample_unwindIterations must contain no NA but supplied is:",
            paste0(sample_unwindIterations, collapse = ",")
        ))
    if (sum(sample_unwindIterations < 0) > 0)
        stop(paste0("sample_unwindIterations values must all be greater than 0"))
    ##
    if (sum(is.na(usePhaseSet_unwindIterations)) > 0)
        stop(paste0(
            "usePhaseSet_unwindIterations must contain no NA but supplied is:",
            paste0(usePhaseSet_unwindIterations, collapse = ",")
        ))
    if (sum(usePhaseSet_unwindIterations < 0) > 0)
        stop(paste0("usePhaseSet_unwindIterations values must all be greater than 0"))
    ##
    if (length(unwindIterations) != length(sample_unwindIterations))
        stop(paste0(
            "unwindIterations and sample_unwindIterations must have the same length, ",
            "but you have given ", 
            "unwindIterations = c(", paste0(unwindIterations, collapse = ","), ") and ",
            "sample_unwindIterations = c(", paste0(sample_unwindIterations, collapse = ","), ")"
        ))
    if (length(unwindIterations) != length(usePhaseSet_unwindIterations))
        stop(paste0(
            "unwindIterations and usePhasSet_unwindIterations must have the same length, ",
            "but you have given ", 
            "unwindIterations = c(", paste0(unwindIterations, collapse = ","), ") and ",
            "usePhaseSet_unwindIterations = c(", paste0(usePhaseSet_unwindIterations, collapse = ","), ")"
        ))
}


calculate_strandedness <- function(
    p_h_given_O,
    nSNPs,
    L,
    inputdir,
    regionName,
    save_sampleReadsInfo,
    kept_reads
) {
    print_message("Calculate strandedness")
    if (save_sampleReadsInfo == FALSE) {
        return(
            list(
                strandedness = NULL,
                strand_info = NULL
            )
        )
    }
    iBam <- 1
    load(STITCH::file_sampleReads(inputdir, 1, regionName))    
    load(STITCH::file_sampleReadsInfo(inputdir, iBam, regionName))
    ## remove reads previously removed

    sampleReadsInfo <- sampleReadsInfo[kept_reads, ]
    sampleReads <- sampleReads[kept_reads]
    ## 
    strand <- match(as.character(sampleReadsInfo[, -1]), c("-", "+")) - 1
    if (length(strand) != nrow(p_h_given_O)) {
        print_message(paste0("length(strand) = ", length(strand)))
        print_message(paste0("nrow(p_h_given_O) = ", nrow(p_h_given_O))) 
    }
    m <- cbind(
        strand = strand,
        p_h1 = p_h_given_O[, 1],
        p_h2 = p_h_given_O[, 2]
    )
    ## make output matrix
    ## for each SNP - add to 8 counts
    out <- array(0, c(nSNPs, 2, 2, 2))
    dimnames(out)[[1]] <- L
    dimnames(out)[[2]] <- c("hap1", "hap2")
    dimnames(out)[[3]] <- c("-", "+")
    dimnames(out)[[4]] <- c("ref", "alt")
    for(iRead in 1:length(sampleReads)) {
        ## extract read, add to components
        s <- sampleReads[[iRead]]
        b <- s[[3]]
        u <- s[[4]]
        ## info about read
        h <- as.integer(m[iRead, "p_h2"] > 0.5) + 1 ## 1 = hap1, 2 = hap2
        s <- m[iRead, "strand"] + 1
        for(j in 1:length(u)) {
            base <- as.integer(b[j] > 0) + 1 ## 1 = ref, 2 = alt
            t <- u[j] + 1
            ##print(c(t, h, s,base))
            out[t, h, s, base] <- out[t, h, s, base] + 1
        }
    }
    ## calculate metrics?
    strandedness <- array(0, c(nSNPs, 2))
    for(t in 1:nSNPs) {
        ## fraction of haplotype that comes from the forward strand
        strandedness[t, 1] <- sum(out[t, 1, "+", ]) / sum(out[t, 1, , ])
        strandedness[t, 2] <- sum(out[t, 2, "+", ]) / sum(out[t, 2, , ])
        ## binomial p-value within haplotype
        ##p1 <- c(
        ##    binom.test(out[t, 1, , "ref"])$p.value,
        ##    binom.test(out[t, 1, , "alt"])$p.value
        ## )
        ##p2 <- c(
        ##    binom.test(out[t, 2, , "ref"])$p.value,
        ##    binom.test(out[t, 2, , "alt"])$p.value
        ## )
        ## forward strand?
        ##strandedness[t, 1] <- 10 * -log10(min(p1))
        ##strandedness[t, 2] <- 10 * -log10(min(p2))
        ### single metric chi square
        ##o <- out[t, , , ]
        ##e <- calculate_hap_strand_base_expectation(o)
        ##strandedness[t, 1] <- -log10(1 - fisher.test(out[t, , , 1])$p.value)
        ##strandedness[t, 2] <- -log10(1 - fisher.test(out[t, , , 2])$p.value)
        ##R <- t(out[t, 1 , ,])
        ##refRatio <- min(R[1, ]) /  max(R[1, ])
        ##altRatio <- min(R[2, ]) /  max(R[2, ])
    }
    strandedness[strandedness == Inf] <- 0.5
    print_message("Done calculating strandedness")
    return(
        list(
            strandedness = strandedness,
            strand_info = out
        )
    )
}

## manually check a few?

calculate_hap_strand_base_expectation <- function(o) {
    n <- sum(o)
    exp <- lapply(1:3, function(j) {
            if (j==1) return(c(sum(o[1, , ]), sum(o[2, , ])) / n)
            if (j==2) return(c(sum(o[, 1, ]), sum(o[, 2, ])) / n)
            if (j==3) return(c(sum(o[, , 1]), sum(o[, , 2])) / n)
    })
    e <- array(0, c(2,2,2))
    for(j1 in 1:2)
        for(j2 in 1:2)
            for(j3 in 1:2)
                e[j1, j2, j3] <- exp[[1]][j1] * exp[[2]][j2] * exp[[3]][j3]
    e <- e * n
    return(e)
}


calculate_phase_entropy <- function(
    eHapsCurrent,
    method = 1
) {
    if (method == 1) {
        x1 <- eHapsCurrent[, 1]
        x2 <- eHapsCurrent[, 2]
        return(-log10(x1 * x2 * (1 - x1) * (1 - x2)))
    } else if (method == 2) {
        x1 <- eHapsCurrent[, 1]
        x2 <- eHapsCurrent[, 2]
        x3 <- x1 * (1 - x1)
        x4 <- x2 * (1 - x2)
        ## large value bad
        w <- x4 > x3
        w[is.na(w)] <- FALSE
        x3[w] <- x4[w]
        return(-log10(x3))
    }
}



write_sew_vcf_header <- function(
    output_vcf_header_file,
    sampleNames,
    strandedness,
    command_line
) {
    ## set up file names
    ## header
    strandedness_header <- NULL
    if (is.null(strandedness) == FALSE) {
        strandedness_header <- c(
            '##INFO=<ID=SB1,Number=.,Type=Float,Description="Strand balance for haplotype 1 i.e. fraction of reads from forward strand from haplotype 1">\n',
            '##INFO=<ID=SB2,Number=.,Type=Float,Description="Strand balance for haplotype 2 i.e. fraction of reads from forward strand from haplotype 2">\n'
        )
    }
    SEW_version <- sessionInfo()$otherPkgs$SEW$Version
    if (is.null(SEW_version))
        SEW_version <- "interactive"
    header <- paste0(
        '##fileformat=VCFv4.0\n',
        "##SEW_version=", SEW_version, "\n",
        "##SEW_command=", command_line, "\n",
        '##INFO=<ID=PE,Number=.,Type=Float,Description="Phase entropy">\n',
        paste0(strandedness_header, collapse = ""),
        '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n',
        '##FORMAT=<ID=DP,Number=1,Type=Float,Description="Read Depth">\n',
        '##FORMAT=<ID=PS,Number=1,Type=Integer,Description="Phase set identifier">\n',
        '##FORMAT=<ID=D1,Number=1,Type=Float,Description="Dosage for haplotype 1, i.e. probability of alternate allele for haplotype 1">\n',
        '##FORMAT=<ID=D2,Number=1,Type=Float,Description="Dosage for haplotype 2, i.e. probability of alternate allele for haplotype 1">\n',        
        '##FORMAT=<ID=ER1,Number=1,Type=Float,Description="Estimated number of copies of reference alleles on haplotype 1">\n',
        '##FORMAT=<ID=EA1,Number=1,Type=Float,Description="Estimated number of copies of alternate alleles on haplotype 1">\n',
        '##FORMAT=<ID=ER2,Number=1,Type=Float,Description="Estimated number of copies of reference alleles on haplotype 2">\n',
        '##FORMAT=<ID=EA2,Number=1,Type=Float,Description="Estimated number of copies of alternate alleles on haplotype 2">\n'        
    )
    header2 <- paste(
        "#CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO", "FORMAT",
        paste(sampleNames, collapse = "\t", sep="\t"),
        sep="\t"
    )
    ##
    ## add output
    ##
    cat(header, header2, "\n", sep="", file = output_vcf_header_file)
    return(NULL)
}

write_sew_vcf_after_em <- function(
    vcf_output_name,
    outputdir,
    regionName,
    eHapsCurrent,
    eHapsUpdate_numer,
    eHapsUpdate_denom,
    alleleCount,
    phase_entropy,
    pos,
    sampleNames,
    phase_set,
    command_line,
    strandedness = NULL
) {
    ## 
    print_message("Build final vcf")
    output_vcf <- STITCH::get_output_filename(
        output_filename = vcf_output_name,
        outputdir = outputdir,
        regionName = regionName,
        output_format = "bgvcf",
        prefix = "sew"
    )
    print_message(paste0("Output vcf:", output_vcf))
    output_vcf_header_file <- paste0(output_vcf, ".header.gz")
    output_vcf_left <- paste0(output_vcf, ".left.gz")
    ## write header
    write_sew_vcf_header(
        output_vcf_header_file = output_vcf_header_file,
        sampleNames = sampleNames,
        strandedness = strandedness,
        command_line = command_line
    )    
    ## body
    options(scipen=6)
    INFO <- paste0(
        "PE=", roundD(phase_entropy, 3)
    )
    if (is.null(strandedness) == FALSE) {
        INFO <- paste0(
            INFO,
            ";SB1=", roundD(strandedness[, 1], 3),
            ";SB2=", roundD(strandedness[, 2], 3)
        )
    }
    ## sample column
    sample_col <- make_sample_column(
        eHapsCurrent = eHapsCurrent,
        eHapsUpdate_numer = eHapsUpdate_numer,
        eHapsUpdate_denom = eHapsUpdate_denom,
        alleleCount = alleleCount,
        phase_set = phase_set
    )    
    to_out <- matrix(
        paste(
            pos[,1], pos[,2], ".", pos[,3], pos[,4], ".", "PASS", INFO,
            "GT:DP:PS:D1:D2:ER1:EA1:ER2:EA2",
            sample_col,
            sep="\t"
        ), ncol=1
    )
    ## finally, write!
    write.table(
        to_out,
        file = gzfile(output_vcf_left),
        row.names = FALSE,
        col.names = FALSE,
        sep = "",
        quote = FALSE
    )
    ## now merge together
    command <- paste0(
        'bash -c "',
        'paste -d ', shQuote("\t"), ' ',
        '<(gunzip -c ', shQuote(output_vcf_left), ' | cat ) ',
        " | ",
        " cat ", shQuote(output_vcf_header_file), " - | bgzip -c > ",
        shQuote(output_vcf)
       ,'"'
    ) 
    ## probably overkill in this version?
    system(command)
    unlink(output_vcf_header_file)
    unlink(output_vcf_left)
    print_message("Done building final VCF")
}



## calculate phase set based on hets alone
## then fill in the hom sites in an ascending manner
## 
recalculate_phase_set <- function(
    inputdir,
    regionName,
    L,
    eHapsCurrent,
    original_phase_set
) {
    print_message("Recalculate phase set")
    load(STITCH::file_sampleReads(inputdir, 1, regionName))
    is_het_snp <- (round(eHapsCurrent[, 1]) + round(eHapsCurrent[, 2])) == 1
    is_het_snp[is.na(is_het_snp)] <- FALSE
    het_snps <- which(is_het_snp)
    ## define as any two points with >0 connections
    phase_set <- L
    ## do global match here
    ##x <- unlist(sapply(sampleReads, function(x) x[[4]])) + 1
    ##y <- unlist(sapply(1:length(sampleReads), function(i) rep(i, sampleReads[[i]][[1]] + 1)))
    ## now get who for y
    ##base_from_reads <- tapply(X = x, INDEX = y, FUN = I)
    for(iRead in 1:length(sampleReads)) {
        u <- sampleReads[[iRead]][[4]] + 1
        u2 <- u[is_het_snp[u]]
        ##u2 <- intersect(u, het_snps) ## ugh
        if (length(u2) > 0) {
            phase_set[u2] <- min(phase_set[u2], L[u2[1]])
            ##phase_set[u2] <- min(c(min(phase_set[u2]), min(L[u2[1]])))
        }
    }
    ## now - have done the hets, ignored the rest
    ## fill in the blanks with their previous values
    phase_set[-het_snps] <- original_phase_set[-het_snps]
    ## the hets may have caused an increase
    ## if there is a decrease, set it up
    ## i.e. make this non-descending
    for(i in 2:length(phase_set)) {
        if (phase_set[i] < phase_set[i - 1]) {
            phase_set[i] <- phase_set[i - 1]
        }
    }
    print_message("Done recalculating phase set")
    return(phase_set)
}



calculate_phase_set <- function(
    inputdir,
    regionName,
    L
) {
    load(STITCH::file_sampleReads(inputdir, 1, regionName))    
    ## define as any two points with >0 connections
    phase_set <- L
    for(iRead in 1:length(sampleReads)) {
        u <- sampleReads[[iRead]][[4]] + 1
        phase_set[u] <- min(phase_set[u], L[u[1]])
    }
    return(phase_set)
}


remove_buffer_from_phase_set <- function(phase_set, L, inRegion2) {
    new_L <- L[inRegion2]
    new_phase_set <- phase_set[inRegion2]
    x <- match(new_phase_set, new_L)
    if (sum(is.na(x) > 0)) {
        to_redo <- unique(new_phase_set[is.na(x)])
        for(i in 1:length(to_redo)) {
            r <- to_redo[i]
            w <- new_phase_set == r
            a <- new_L[w] ## new positions of this group
            new_phase_set[w] <- min(a)
        }
    }
    return(new_phase_set)
}

make_sample_column <- function(
    eHapsCurrent,
    eHapsUpdate_numer,
    eHapsUpdate_denom,
    alleleCount,
    phase_set,
    n_digits = 4
) {
    ## always try genotype? take most likely marginal value of each haplotype?
    gt <- paste0(round(eHapsCurrent[, 1]), "|", round(eHapsCurrent[, 2]))
    dp <- alleleCount[,2]
    phase_set <- phase_set
    ea <- eHapsUpdate_numer ## col 1 = hap 1, col 2 = hap2    
    er <- eHapsUpdate_denom - eHapsUpdate_numer    
    str <- paste(
        gt,
        roundD(dp, n_digits),
        phase_set,
        roundD(eHapsCurrent[, 1], n_digits),
        roundD(eHapsCurrent[, 2], n_digits),
        roundD(er[, 1], n_digits),
        roundD(ea[, 1], n_digits),
        roundD(er[, 2], n_digits),
        roundD(ea[, 2], n_digits),
        sep = ":"
    )
    return(str)
}

roundD <- function(x, n_digits)
    format(round(x, n_digits), nsmall = n_digits, trim = TRUE)


## run entirety of single sample phasing here
single_sample_phasing_em_complete <- function(
    inputdir,
    regionName,
    nSNPs,
    phase,
    niterations,
    unwindIterations,
    sample_unwindIterations,
    usePhaseSet_unwindIterations,
    n_unwinding_its = 2,
    phase_set,
    very_verbose
) {

    ## set.seed(1)
    load(STITCH::file_sampleReads(inputdir, 1, regionName))
    kept_reads <- array(TRUE, length(sampleReads))
   
    ## initialize
    n_unwound <- NA ## none previously unwound
    K <- 2 ## this will only work for K=2
    eHapsCurrent <- matrix(runif(nSNPs * K), ncol = K)
    eHapsCurrent[, 1] <- 0.1
    eHapsCurrent[, 2] <- 0.9
    eHapsCurrent <- eHapsCurrent / rowSums(eHapsCurrent)
    N_r <- length(sampleReads)
    print_message(paste0("There are ", nSNPs, " SNPs and ", N_r, " reads"))

    ## NULL if no unwindIterations
    reads_span <- get_what_reads_span(
        sampleReads = sampleReads,
        nSNPs = nSNPs,
        unwindIterations = unwindIterations
    )
    print_message("Begin iterations")

    for(iteration in 1:niterations) { 

        ## calculate read probabilities
        out <- cpp_single_sample_expectation(
            sampleReads = sampleReads,
            eHapsCurrent = eHapsCurrent,
            nSNPs = nSNPs,
            N_r = N_r,
            calculate_updates = TRUE
        )
        eHapsUpdate_numer <- out$eHapsUpdate_numer
        eHapsUpdate_denom <- out$eHapsUpdate_denom
        p_reads <- out$p_reads ## for likelihood
        p_reads_given_hap_k <- out$p_reads_given_hap_k ## for read probabilities
        p_h_given_O <- out$p_h_given_O

        ## if any reads have p_reads = 0 (or worse?), remove them
        ## can have NA if reads are insanely long
        ## ultimately we would want these (probably?) for interest of time now exclude
        ## this can happen periodically, usually at initialization with extremely long reads        
        bad_reads <- (p_reads <= 0) | (is.na(p_reads))
        to_remove <- which(bad_reads)
        if (length(to_remove) > 0) {
            print_message(paste0("Removing ", sum(p_reads[to_remove] <= 0, na.rm = TRUE), " reads with P(read) <= 0"))
            print_message(paste0("Removing ", sum(is.na(p_reads)), " reads with P(read) = NA"))
            to_keep <- bad_reads == FALSE
            sampleReads <- sampleReads[to_keep]
            ## here, of those that are still kept, remove as appropriate
            kept_reads[kept_reads][to_remove] <- FALSE
            N_r <- length(sampleReads)            
            p_reads <- p_reads[to_keep]
            p_reads_given_hap_k <- p_reads_given_hap_k[to_keep, ]
            p_h_given_O <- p_h_given_O[to_keep, ]
            ## if any entries are now NA as a result, set to not NA for now. will get over-written later
            to_redo <- which(rowSums(is.na(eHapsUpdate_numer)) > 0 | rowSums(is.na(eHapsUpdate_denom)) > 0)
            eHapsUpdate_numer[to_redo, ] <- runif(length(to_redo) * 2)
            eHapsUpdate_denom[to_redo, ] <- runif(length(to_redo) * 2)
            reads_span <- get_what_reads_span(
                sampleReads = sampleReads,
                nSNPs = nSNPs,
                unwindIterations = unwindIterations
            )
        }

        ## calculate read probabilities
        ## calculate new parameters
        print_message(paste0("iteration = ", iteration, ", ", date()))
        print_message(paste0("Log likelihood = ", -sum(log(p_reads))))
        if (is.na(-sum(log(p_reads)))) {
            save_file <- file.path(inputdir, paste0("debug.", regionName, ".RData"))
            save(N_r, nSNPs, sampleReads, p_reads, eHapsUpdate_denom, eHapsUpdate_numer, p_reads_given_hap_k, p_h_given_O, eHapsCurrent, file = save_file)
            print_message(paste0("Saving some things to:", save_file))
            stop("Error, see above messages")
        }
        eHapsCurrent <- make_new_eHaps(
            eHapsUpdate_numer = eHapsUpdate_numer,
            eHapsUpdate_denom = eHapsUpdate_denom
        )
        
        out <- single_sample_phase_heuristic(
            sampleReads = sampleReads,
            eHapsCurrent = eHapsCurrent,
            nSNPs = nSNPs,
            iteration = iteration,
            p_reads = p_reads,
            reads_span = reads_span,            
            unwindIterations = unwindIterations,
            sample_unwindIterations = sample_unwindIterations,
            usePhaseSet_unwindIterations = usePhaseSet_unwindIterations,
            n_unwinding_its = n_unwinding_its,
            phase_set = phase_set,
            n_unwound = n_unwound,
            very_verbose = very_verbose
        )
        eHapsCurrent <- out$eHapsCurrent
        n_unwound <- out$n_unwound

        ## temporary plotting turned on as necessary
        if ( 1 == 0 ) {
            if (is.na(match(iteration, unwindIterations)) == FALSE) {
                for(t_break in c(21)) {
                    which_reads <- reads_span[t_break, "readStart"]:reads_span[t_break, "readEnd"]
                    which_snps <- reads_span[t_break, "SNPstart"]:reads_span[t_break, "SNPend"]
                    phase_debug_plot(
                        plot_name = paste0("~/tempX.", t_break, ".", iteration, ".pdf"),
                        nSNPs,
                        sampleReads,
                        phase = NULL,
                        eHapsUpdate_numer,
                        eHapsUpdate_denom,
                        which_reads = which_reads,
                        which_snps = which_snps
                    )
                }
            }
        }

        if (is.null(phase) == FALSE) {
            pse <- STITCH::calculate_pse(test = eHapsCurrent, truth = phase[, 1, ])
            print_message(paste0("The PSE is ", pse))
        }
    }

    print_message("Done EM updates")

    return(
        list(
            eHapsCurrent = eHapsCurrent,
            eHapsUpdate_numer = eHapsUpdate_numer,
            eHapsUpdate_denom = eHapsUpdate_denom,
            p_reads = p_reads,
            p_reads_given_hap_k = p_reads_given_hap_k,
            p_h_given_O = p_h_given_O,
            final_likelihood = -sum(log(p_reads)),
            reads_span = reads_span,
            kept_reads = kept_reads
        )
    )
}

single_sample_phase_heuristic <- function(
    sampleReads,
    eHapsCurrent,
    nSNPs,
    iteration,
    p_reads,
    reads_span,
    unwindIterations = NA,
    threshold_unwindIterations = NA,
    usePhaseSet_unwindIterations = NA,
    sample_unwindIterations = NA,
    n_unwinding_its = 2,
    phase_set = NULL,
    n_unwound,
    very_verbose = FALSE
) {

    if (is.na(unwindIterations[1] ))
        return(eHapsCurrent)

    ## for each break, want reads that span it, SNPs that span it
    ## done_unwinding <- is.na(n_unwound) == FALSE & n_unwound == 0
    done_unwinding <- FALSE
    if (is.na(match(iteration, unwindIterations)) == FALSE) {
        if (done_unwinding == FALSE) {

            ## for all possible breaks, calculate change in likelihood
            print_message("Calculate whether unwindings are good")
            reads_span <- calculate_ll_for_breaks(
                sampleReads = sampleReads,
                eHapsCurrent = eHapsCurrent,
                reads_span = reads_span,
                nSNPs = nSNPs,
                p_reads = p_reads,
                n_unwinding_its = n_unwinding_its,
                sample_periodicity = sample_unwindIterations[match(iteration, unwindIterations)]
            )

            print_message("Done calculating whether unwindings are good")
            
            print_message("Begin applying unwindings")
            if (usePhaseSet_unwindIterations[match(iteration, unwindIterations)]) {
                phase_set_to_use <- phase_set
            } else {
                phase_set_to_use <- NULL
            }
            out <- apply_flipping(
                reads_span = reads_span,
                nSNPs = nSNPs,
                eHapsCurrent = eHapsCurrent,
                phase_set = phase_set_to_use,
                very_verbose = very_verbose
            )
            eHapsCurrent <- out$eHapsCurrent
            n_unwound <- out$n_unwound
            
            ## phase_set = phase_set
            
            print_message("Done unwinding")
        } else {
            print_message("No longer unwinding as previous unwinding yielded no unwindings to apply")
        }
    }

    return(
        list(
            eHapsCurrent = eHapsCurrent,
            n_unwound = n_unwound 
        )
    )

}


calculate_ll_for_breaks <- function(
    sampleReads,
    eHapsCurrent,
    reads_span,
    nSNPs,
    p_reads,
    n_unwinding_its = 2,
    sample_periodicity = 1
) {
    for(t_break in 1:(nSNPs - 1)) {
        if ((t_break %% sample_periodicity) == 0) {
            ## calculate new local log likelihood
            if (is.na(reads_span[t_break, "readStart"]) == FALSE) {
                which_reads <- reads_span[t_break, "readStart"]:reads_span[t_break, "readEnd"]
                first_snp <- reads_span[t_break, "SNPstart"]
                last_snp <- reads_span[t_break, "SNPend"]
                ## sum(-log(p_reads[which_reads]))
                reads_span[t_break, "llOri"] <- investigate_an_unwinding(
                    sampleReads = sampleReads,
                    eHapsCurrent = eHapsCurrent,
                    nSNPs = nSNPs,
                    t_break = 0,
                    which_reads = which_reads,
                    n_unwinding_its = n_unwinding_its,
                    first_snp = first_snp,
                    last_snp = last_snp
                )
                reads_span[t_break, "llNew"] <- investigate_an_unwinding(
                    sampleReads = sampleReads,
                    eHapsCurrent = eHapsCurrent,
                    nSNPs = nSNPs,
                    t_break = t_break,
                    which_reads = which_reads,
                    n_unwinding_its = n_unwinding_its,
                    first_snp = first_snp,
                    last_snp = last_snp
                )
            }
        }
    }
    return(reads_span)
}



apply_flipping <- function(
    reads_span,
    nSNPs,
    eHapsCurrent,
    phase_set = NULL,
    very_verbose = FALSE
) {
    phase_setL <- phase_set[-length(phase_set)]
    relative_ll <- reads_span[, "llOri"] - reads_span[, "llNew"]
    ## only consider if 2 orders of magnitude better. this is a very light filter
    to_flip <- relative_ll > 0 ## consider basically any that look better
    to_flip[is.na(to_flip)] <- FALSE
    done <- 0
    remaining_eligible <- to_flip
    use <- c()
    ## all SNPs start off eligible that increase the likelihood
    ## take the SNP with the biggest effect. apply it. remove from eligible
    ## continue until done
    ## only do 1 per region it spans
    while (sum(remaining_eligible) > 0) {
        best <- which.max(relative_ll[remaining_eligible])
        x <- reads_span[remaining_eligible, , drop = FALSE][best, ]
        if (x["SNPend"] == nSNPs) {
            x["SNPend"] <- nSNPs - 1
        }
        t_use <- x["t_pre"]
        use <- c(use, t_use)
        ## filter based on read intersections
        if (is.null(phase_set)) {
            ## remove SNPs in the radius of its radius
            a <- x["SNPstart"]
            b <- x["SNPend"]
            if (a > 1)
                a <- a - 1
            if (b < (nSNPs - 1))
                b <- b + 1
            y <- a:b
            y2 <- min(reads_span[y, "SNPstart"], na.rm=TRUE):max(reads_span[y, "SNPend"], na.rm=TRUE)
            y2 <- y2[y2 != nSNPs]
            remaining_eligible[y2] <- FALSE
        } else {
            remaining_eligible[phase_setL == phase_setL[t_use]] <- FALSE
        }
    }
    ## this could be applied smarter
    ## & length(intersect(1, use)) == 0 & length(intersect(nSNPs, use)) == 0
    n_unwound <- length(use)
    if (length(use) > 0) {
        use <- sort(use)
        start <- c(1, use + 1)
        end <- c(use, nSNPs)
        if (very_verbose) {
            print_message(paste0("Use the following unwindingins:", paste0(use, ":", round(relative_ll[use], 2)), collapse = ", "))
        }
        print_message(paste0(
            "Applying ", length(use), " unwindings of which ",
            sum(relative_ll[use] > log(1e10)),
            " have > 10 orders of magnitude effect"
        ))
        for(i_use in 1:length(start)) {
            flip <- (i_use %% 2) == 1
            if (flip) {
                if (very_verbose) {
                    print_message(paste0("Flipping:", start[i_use], "-", end[i_use]))
                }
                w <- start[i_use]:end[i_use]
                eHapsCurrent[w, 1:2] <- eHapsCurrent[w, 2:1]
            }
        }
        ## add some noise?
        ## K <- 2
        ## eHapsCurrentNoise <- matrix(runif(T * K), ncol = K)
        ## eHapsCurrentNoise <- eHapsCurrentNoise / rowSums(eHapsCurrentNoise)
        ## eHapsCurrent <- 0.8 * eHapsCurrent + 0.2 * eHapsCurrentNoise
    } else {
        print_message("No unwindings were found that improved the likelihood")
    }
    return(
        list(
            eHapsCurrent = eHapsCurrent,
            n_unwound = n_unwound
        )
    )
}




## if first_snp is NA, do updates on full set naively
## if it is not NA, then for the first iteration, write to a smallet set
## then for subsequent iterations, read from the smaller set
investigate_an_unwinding <- function(
    sampleReads,
    eHapsCurrent,
    nSNPs,
    t_break,
    which_reads,
    n_unwinding_its = 2,
    first_snp = NA, ## 0-based
    last_snp = NA ## 0-based
) {
    calculate_updates <- TRUE
    ## default - set to NA
    t_reading_offset <- 0
    t_writing_offset <- 0
    nSNPs_use <- nSNPs
    t_reading_offset <- 0
    ## if this is NA - then write out a smaller amount, when calculating updates
    if (is.na(first_snp) == FALSE) {
        t_writing_offset <- first_snp - 1
        nSNPs_use <- last_snp - first_snp + 1
    }
    for(it in 1:n_unwinding_its) {
        if (it > 1) {
            t_break <- 0
            if (is.na(first_snp) == FALSE)
                t_reading_offset <- first_snp - 1                
        }
        if (it == n_unwinding_its)
            calculate_updates <- FALSE
        sampleReads_local <- sampleReads[which_reads]
        out <- cpp_single_sample_expectation(
            sampleReads = sampleReads_local,
            eHapsCurrent = eHapsCurrent,
            nSNPs = nSNPs_use,
            N_r = length(which_reads),
            calculate_updates = calculate_updates,
            t_break = t_break,
            t_reading_offset = t_reading_offset,
            t_writing_offset = t_writing_offset,
            calculate_read_probabilities = FALSE ## never want these here
        )
        if (calculate_updates)
            eHapsCurrent <- make_new_eHaps(
                eHapsUpdate_numer = out$eHapsUpdate_numer,
                eHapsUpdate_denom = out$eHapsUpdate_denom
            )
    }
    return(out$log_likelihood)
}


make_new_eHaps <- function(
    eHapsUpdate_numer,
    eHapsUpdate_denom
) {
    eHapsFuture <- eHapsUpdate_numer / eHapsUpdate_denom
    eHapsFuture[eHapsFuture < 0.001] <- 0.001
    eHapsFuture[eHapsFuture > (1 - 0.001)] <- (1 - 0.001)
    return(eHapsFuture)
}



get_what_reads_span <- function(
    sampleReads,
    nSNPs,
    unwindIterations
) {

    if (is.na(unwindIterations[1] ))
        return(NULL)

    out <- array(NA, c(nSNPs - 1, 10))
    colnames(out) <- c(
        "t_pre", "t_post", "readStart", "readEnd",
        "readIntersectStart", "readIntersectEnd", 
        "SNPstart", "SNPend", "llOri", "llNew"
    )
    out[, 1] <- 1:(nSNPs-1)
    out[, 2] <- out[, 1] + 1 ## overkill but good for remembering and understanding

    read_origin <- get_where_reads_come_from(sampleReads)

    span_snps <- tapply(read_origin[, 1], read_origin[, 2], I)

    left <- span_snps[match(out[, 1], names(span_snps))]
    right <- span_snps[match(out[, 2], names(span_snps))]
   
    for(t in 1:(nSNPs-1)) {
        ## changed from intersect to union
        ## so any reads that intersect the two sites
        which_reads <- union(left[[t]], right[[t]])
        which_reads_intersect <- intersect(left[[t]], right[[t]])        
        if (length(which_reads) > 0) {
            which_reads_min <- min(which_reads)
            which_reads_max <- max(which_reads)        
            out[t, "readStart"] <- which_reads_min
            out[t, "readEnd"] <- which_reads_max
            if (length(which_reads_intersect) > 0) {
                which_reads_intersect_min <- min(which_reads_intersect)
                which_reads_intersect_max <- max(which_reads_intersect)        
                out[t, "readIntersectStart"] <- which_reads_intersect_min
                out[t, "readIntersectEnd"] <- which_reads_intersect_max
            }
           ## get the SNPs these reads intersect
            which_reads_all <- which_reads_min:which_reads_max
            y <- sort(unique(unlist(lapply(sampleReads[which_reads_all], function(x) x[[4]]))))
            out[t, "SNPstart"] <- min(y) + 1
            out[t, "SNPend"] <- max(y) + 1
        }
    }

    return(out)
}

get_where_reads_come_from <- function(sampleReads) {
    x <- sum(sapply(sampleReads, function(x) {
        u <- x[[4]]
        return(max(u) - min(u) + 1)
    }))
    out2 <- array(0, c(x, 2))
    c <- 1
    for(iRead in 1:length(sampleReads)) {
        sampleRead <- sampleReads[[iRead]]
        ## change to include all SNP that span
        ##J <- sampleRead[[1]]
        u <- sampleRead[[4]]
        u <- matrix(min(u):max(u), ncol = 1)
        J <- nrow(u) - 1
        out2[c + 0:J, 1] <- iRead
        out2[c + 0:J, 2] <- u + 1
        c <- c + J + 1
    }
    colnames(out2) <- c("read", "u")
    return(out2)
}
    


## single iteration, calculate expectation
## t_break - calculate everything assuming there is a 1-based break between
##    eHapsCurrent between SNPs t and t + 1
## t_reading_offset - 0-based offset for reading. when reading eHapsCurrent, it is this much shorter than expected based on u. for example, if normally expected eHapsCurrent to have 6 rows, with a (0-based) u going for 2-3, then t_reading_offset = 2 with an eHapsCurrent having rows 3-6 does the right thing
## t_writing_offset - write results omitting 
single_sample_expectation <- function(
    sampleReads,
    eHapsCurrent,
    nSNPs,
    N_r,
    calculate_updates = TRUE,
    calculate_read_probabilities = TRUE,
    t_break = 0,
    t_reading_offset = 0,
    t_writing_offset = 0
) {
    K <- 2
    log_likelihood <- 0
    nSNPs_output <- 1
    N_r_output <- 1
    if (calculate_updates) {
        nSNPs_output <- nSNPs
    }
    if (calculate_read_probabilities)    
        N_r_output <- N_r
    eHapsUpdate_numer <- array(0, c(nSNPs, K))
    eHapsUpdate_denom <- array(0, c(nSNPs, K))
    p_reads <- array(0, c(1, N_r))        
    p_reads_given_hap_k <- array(0, c(N_r, 2))
    p_h_given_O <- array(0, c(N_r, 2))
    for(r in 1:N_r) {
        sampleRead <- sampleReads[[r]]
        J_r <- sampleRead[[1]]
        bq <- sampleRead[[3]]
        u <- sampleRead[[4]]
        pr <- STITCH::convertScaledBQtoProbs(bq)
        flip_k <- array(0, J_r + 1)
        if (t_break > 0) {
            flip <- (u + 1) > t_break
            flip_k[flip] <- 1            
        }
        p_read_given_hap_k_components <- array(0, c(J_r + 1, K))
        for(k in 1:K) {
            use_k <- k
            for(j in 1:(J_r + 1)) {
                u_j <- u[j] + 1                
                u_j_reading <- u[j] + 1 - t_reading_offset
                if (flip_k[j] == 1)
                    use_k <- 3 - k
                p_read_given_hap_k_components[j, k] <-
                    eHapsCurrent[u_j_reading, use_k] * pr[j, 2] +
                    (1 - eHapsCurrent[u_j_reading, use_k]) * pr[j, 1]
            }
        }
        p_read_given_hap_k <- apply(p_read_given_hap_k_components, 2, prod)
        p_read <- 0.5 * sum(p_read_given_hap_k)
        log_likelihood <- log_likelihood + -sum(log(p_read))
        if (calculate_read_probabilities) {
            p_reads[1, r] <- p_read
            p_reads_given_hap_k[r, ] <- p_read_given_hap_k
        }
        if (calculate_updates) {
            for(j in 1:(J_r + 1)) {
                u_j_reading <- u[j] + 1 - t_reading_offset
                u_j_writing <- u[j] + 1 - t_writing_offset                    
                for(k in 1:K) {
                    use_k <- k
                    if (flip_k[j] == 1)
                        use_k <- 3 - k
                    p_r_no_j <- p_read_given_hap_k[k] / p_read_given_hap_k_components[j,k]
                    p_h_and_g_1 <- p_r_no_j * 0.5 * pr[j, 2] *
                        eHapsCurrent[u_j_reading, use_k] / p_read
                    p_h_and_g_2 <- p_r_no_j * 0.5 * pr[j, 1] *
                        (1 - eHapsCurrent[u_j_reading, use_k]) / p_read
                    if (calculate_read_probabilities)
                        p_h_given_O[r, k] <- p_h_and_g_1 + p_h_and_g_2;
                    eHapsUpdate_numer[u_j_writing, k] <-
                        eHapsUpdate_numer[u_j_writing, k] + p_h_and_g_1
                    eHapsUpdate_denom[u_j_writing, k] <-
                        eHapsUpdate_denom[u_j_writing, k] +
                        p_h_and_g_1 + p_h_and_g_2
                }
            }
        }
    }
    if (!calculate_read_probabilities & calculate_updates) {
        return(
            list(
                eHapsUpdate_numer = eHapsUpdate_numer,
                eHapsUpdate_denom = eHapsUpdate_denom,
                log_likelihood = log_likelihood,
            )
        )
    } else if (calculate_read_probabilities & calculate_updates) {
        return(
            list(
                eHapsUpdate_numer = eHapsUpdate_numer,
                eHapsUpdate_denom = eHapsUpdate_denom,
                log_likelihood = log_likelihood,                
                p_reads = p_reads,
                p_reads_given_hap_k = p_reads_given_hap_k,
                p_h_given_O = p_h_given_O
            )
        )
    } else if (!calculate_updates) {
        return(
            list(
                log_likelihood = log_likelihood
            )
        )
    }
}
    

    



plot_pileup <- function(eHapsUpdate_numer, eHapsUpdate_denom, k, which_snps) {
    par(mar = c(2, 2, 2, 2))
    up <- eHapsUpdate_numer[which_snps, k]
    down <- eHapsUpdate_denom[which_snps, k] - eHapsUpdate_numer[which_snps, k]
    x <- which_snps
    xlim <- range(which_snps)
    ylim <- c(-max(down), max(up))
    plot(
        x = 0, y = 0, xlim = xlim, ylim = ylim,
        xlab = "", ylab = "", axes = TRUE, col = "white"
    )
    ## plot up = alt probabilities
    rect(xleft = x - 0.5, xright = x + 0.5, ybottom = 0, ytop = up, col = "red")
    rect(xleft = x - 0.5, xright = x + 0.5, ybottom = - down, ytop = 0, col = "blue")
}

plot_phase <- function(phase, which_snps) {
    par(mar = c(2, 2, 2, 2))
    x <- 1:length(which_snps)
    xlim <- c(0, length(x))
    ylim <- c(0, 2)
    plot(x = 0, y = 0, xlim = xlim, ylim = ylim, xlab = "", ylab = "", axes = TRUE, col = "white")
    for(k in 1:2) {
        col <- c("blue", "red")[phase[which_snps, 1, k] + 1]
        rect(xleft = x - 0.5, xright = x + 0.5, ybottom = k - 1, ytop = k, col = col)
    }
}



phase_add_points <- function(which_snps, which_reads, sampleReads, max_depth = 100) {
    par(mar = c(2, 2, 2, 2))    
    xlim <- c(0, length(which_snps))
    ylim <- c(0, max_depth) # plot up
    plot(x = 0, y = 0, xlim = xlim, ylim = ylim, xlab = "", ylab = "", axes = FALSE, col = "white")
    axis(side = 1, at = 1:length(which_snps), labels = which_snps)
    axis(side = 2)
    n_reads <- length(which_reads)
    plotted <- array(0, c(max_depth, length(which_snps)))
    for(r in which_reads) {
        ## find where to plot
        sampleRead <- sampleReads[[r]]
        bq <- sampleRead[[3]]
        u <- sampleRead[[4]]
        u <- match(u + 1, which_snps)
        keep <- is.na(u) == FALSE
        bq <- bq[keep]
        u <- u[keep]
        J_r <- length(bq)        
        if (J_r > 0) {
            r_depth <- 1
            while (r_depth < max_depth) {
                if (sum(plotted[r_depth, u]) == 0) {
                    ## can plot here
                    plotted[r_depth, u] <- r
                    lines(x = u, y = rep(r_depth, J_r))                    
                    rect(
                        xleft = u - 0.4, xright = u + 0.4,
                        ybottom = rep(r_depth - 0.4, J_r), ytop = rep(r_depth + 0.4, J_r),
                        col = c("blue", "red")[as.integer(bq > 0) + 1],
                        border = NA
                    )
                    ## add joining line
                    r_depth <- max_depth
                } else {
                    r_depth <- r_depth + 1
                }
            }
        }
    }
    return(plotted)
}



phase_debug_plot <- function(
    plot_name,
    nSNPs,
    sampleReads,
    phase = NULL,
    eHapsUpdate_numer,
    eHapsUpdate_denom,
    which_reads = NA,
    which_snps = NA
) {
    if (is.na(which_reads[1]))
        which_reads <- 1:length(sampleReads)
    if (is.na(which_snps[1])) {
        which_snps <- 1:nSNPs
    }
    pdf(file = plot_name, height = 12, width = 8)#, units = "in", res = 300)
    par(fig = c(0, 1, 0.5, 1), new = FALSE)    
    plotted <- phase_add_points(which_snps, which_reads, sampleReads)
    if (is.null(phase)) {
        par(fig = c(0, 1, 2/6, 3 / 6), new = TRUE)
        plot_phase(phase, which_snps)
    }
    par(fig = c(0, 1, 1/6, 2/6), new = TRUE)
    plot_pileup(eHapsUpdate_numer, eHapsUpdate_denom, k = 1, which_snps)
    par(fig = c(0, 1, 0, 1/6), new = TRUE)
    plot_pileup(eHapsUpdate_numer, eHapsUpdate_denom, k = 2, which_snps)
    dev.off()
}




