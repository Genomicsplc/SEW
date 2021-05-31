#' @title SEW
#' @param chr What chromosome to run. Should match BAM headers
#' @param posfile Where to find file with positions to run. File is tab seperated with no header, one row per SNP, with col 1 = chromosome, col 2 = physical position (sorted from smallest to largest), col 3 = reference base, col 4 = alternate base. Bases are capitalized. Example first row: 1<tab>1000<tab>A<tab>G<tab>
#' @param outputdir What output directory to use
#' @param tempdir What directory to use as temporary directory. If possible, use ramdisk, like /dev/shm/
#' @param bamlist Path to file with bam file locations. File is one row per entry, path to bam files. Bam index files should exist in same directory as for each bam, suffixed either .bam.bai or .bai
#' @param cramlist Path to file with cram file locations. File is one row per entry, path to cram files. cram files are converted to bam files on the fly for parsing into SEW
#' @param reference Path to reference fasta used for making cram files. Only required if cramlist is defined
#' @param genfile Path to gen file with high coverage results. Empty for no genfile. File has a header row with a name for each sample, matching what is found in the bam file. Each subject is then a tab seperated column, with 0 = hom ref, 1 = het, 2 = hom alt and NA indicating missing genotype, with rows corresponding to rows of the posfile. Note therefore this file has one more row than posfile which has no header
#' @param phasefile Path to phase file with truth phasing results. Empty for no phasefile. File has a header row with a name for each sample, matching what is found in the bam file. Each subject is then a tab seperated column, with 0 = ref and 1 = alt, separated by a vertical bar |, e.g. 0|0 or 0|1. Note therefore this file has one more row than posfile which has no header
#' @param outputHaplotypeProbabilities Whether to output haplotype probabilities in files    
#' @param regenerateInput Whether to regenerate input files
#' @param originalRegionName If regenerateInput is FALSE (i.e. using existing data), this is the name of the original region name (chr.regionStart.regionEnd). This is necessary to load past variables
#' @param keepInterimFiles Whether to keep interim parameter estimates
#' @param keepTempDir Whether to keep files in temporary directory
#' @param generateInputOnly Whether to just generate input data then quit
#' @param useSoftClippedBases Whether to use (TRUE) or not use (FALSE) bases in soft clipped portions of reads
#' @param vcf_output_name Override the default VCF output name with this given file name. Please note that this does not change the names of inputs or outputs (e.g. RData, plots), so if outputdir is unchanged and if multiple SEW runs are processing on the same region then they may over-write each others inputs and outputs
#' @param initial_min_hapProb Initial lower bound for probability read comes from haplotype. Double bounded between 0 and 1
#' @param initial_max_hapProb Initial upper bound for probability read comes from haplotype. Double bounded between 0 and 1
#' @param regenerateInputWithDefaultValues If regenerateInput is FALSE and the original input data was made using regionStart, regionEnd and buffer as default values, set this equal to TRUE
#' @param unwindIterations What iterations to search for flips in estimated haplotypes that increase the likelihood
#' @param sample_unwindIterations When performing an unwind iterations, how often to check a SNP for an improvement following unwinding. 1 means check all SNPs
#' @param save_sampleReadsInfo Experimental. Boolean TRUE/FALSE about whether to save additional information about the reads that were extracted
#' @param outputInputInVCFFormat Whether to output the input in vcf format
#' @param downsampleToCov What coverage to downsample individual sites to. This ensures no floating point errors at sites with really high coverage
#' @param downsampleFraction Downsample BAMs by choosing a fraction of reads to retain. Must be value 0<downsampleFraction<1
#' @param chrStart When loading from BAM, some start position, before SNPs occur. Default NA will infer this from either regionStart, regionEnd and buffer, or posfile
#' @param chrEnd When loading from BAM, some end position, after SNPs occur. Default NA will infer this from either regionStart, regionEnd and buffer, or posfile
#' @param regionStart When running imputation, where to start from. The 1-based position x is kept if regionStart <= x <= regionEnd
#' @param regionEnd When running imputation, where to stop
#' @param buffer Buffer of region to perform imputation over. So imputation is run form regionStart-buffer to regionEnd+buffer, and reported for regionStart to regionEnd, including the bases of regionStart and regionEnd
#' @param iSizeUpperLimit Do not use reads with an insert size of more than this value
#' @param bqFilter Minimum BQ for a SNP in a read. Also, the algorithm uses bq<=mq, so if mapping quality is less than this, the read isnt used
#' @param niterations Number of EM iterations.
#' @param plotHapSumDuringIterations Boolean TRUE/FALSE about whether to make a plot that shows the relative number of individuals using each ancestral haplotype in each iteration
#' @param very_verbose Boolean TRUE/FALSE about whether to print additional values which may be helpeful for understanding how the model operates
#' @param keepSampleReadsInRAM Whether to (generally) keep sampleReads in RAM or store them in the temporary directory. SEW will be faster if this is FALSE at the expense of RAM
#' @param use_bx_tag Whether to try and use BX tag in same to indicate that reads come from the same underlying molecule
#' @param bxTagUpperLimit When using BX tag, at what distance between reads to consider reads with the same BX tag to come from different molecules
#' @return Results in properly formatted version
#' @author Robert Davies
#' @export
SEW <- function(
    chr,
    posfile,
    outputdir,
    tempdir = NA,
    bamlist = "",
    cramlist = "",
    reference = "",
    genfile = "",
    phasefile = "",
    outputHaplotypeProbabilities = FALSE,
    regenerateInput = TRUE,
    originalRegionName = NA,
    keepInterimFiles = FALSE,
    keepTempDir = FALSE,
    generateInputOnly = FALSE,
    useSoftClippedBases = FALSE,
    vcf_output_name = NULL,
    initial_min_hapProb = 0.4,
    initial_max_hapProb = 0.6,
    regenerateInputWithDefaultValues = FALSE,
    unwindIterations = c(40, 70, 100, 120, 140, 160, 180, 200),
    sample_unwindIterations = c(5, 5, 1, 1, 1, 1, 1, 1),
    usePhaseSet_unwindIterations = c(rep(FALSE, 8)),
    save_sampleReadsInfo = TRUE,
    outputInputInVCFFormat = FALSE,
    downsampleToCov = 200,
    downsampleFraction = 1,
    chrStart = NA,
    chrEnd = NA,
    regionStart = NA,
    regionEnd = NA,
    buffer = NA,
    iSizeUpperLimit = as.integer(600),
    bqFilter = as.integer(-1),
    niterations = 300,
    plotHapSumDuringIterations = FALSE,
    very_verbose = FALSE,
    keepSampleReadsInRAM = FALSE,
    use_bx_tag = TRUE,
    bxTagUpperLimit = 50000    
 ) {


    ## capture command line
    x <- as.list(environment())
    command_line <- paste0(
        "SEW(",
        paste(names(x), " = ", x, collapse = ", ", sep = ""),
        ")"
    )
    print_message(paste0("Running ", command_line))    
    
    ## do not export
    inputBundleBlockSize <- NA
    nCores <- 1 ## 1 sample, nothing operates on cores
    ##

    ##
    ## specify output description
    ##
    regionName <- chr
    options(scipen = 999) ## Dangit rounding
    if(is.na(regionStart) == FALSE & is.na(regionEnd) == FALSE)
        regionName <- paste0(chr, ".", regionStart,".", regionEnd)



    ##
    ## validate parameters
    ##
    STITCH::validate_nCores(nCores)    
    STITCH::validate_chr(chr)    
    STITCH::validate_posfile(posfile)
    STITCH::validate_outputdir(outputdir)
    STITCH::validate_tempdir(tempdir)
    STITCH::validate_downsampleFraction(downsampleFraction)
    STITCH::validate_regionStart_regionEnd_and_buffer(regionStart, regionEnd, buffer)
    STITCH::validate_bamlist_and_cramlist_for_input_generation(regenerateInput, originalRegionName, bamlist, cramlist, regionStart, regionEnd, buffer, regenerateInputWithDefaultValues)
    STITCH::validate_vcf_output_name(vcf_output_name)
    validate_unwindIterations(unwindIterations, sample_unwindIterations, usePhaseSet_unwindIterations)
    validate_niterations(niterations)
    
    ##
    ##
    ## make necessary directory
    ##
    ##
    out <- STITCH::initialize_directories(
        tempdir = tempdir,
        keepTempDir = keepTempDir,
        outputdir = outputdir
    )
    tempdir <- out$tempdir
    inputdir <- out$inputdir
    
    ##
    ## start
    ##
    print_message("Program start")
    file <- file.path(outputdir, "RData", paste0("start.", regionName ,".RData"))
    if(generateInputOnly == FALSE) # why?
        save(date, file = file)
    

    ##
    ## load the positions, genotypes and phases
    ##
    out <- STITCH::get_and_validate_pos_gen_and_phase(
        posfile = posfile,
        genfile = genfile,
        phasefile = phasefile,
        chr = chr,
        verbose = TRUE
    )
    pos <- out$pos
    gen <- out$gen
    phase <- out$phase
    nSNPs <- out$nSNPs
    L <- out$L

    if (genfile == "" & phasefile != "") {
        ## not sure this is a good idea?        
        gen <- matrix(phase[, , 1] + phase[, , 2], ncol = 1)
    }

    ##
    ## shrink (if regionStart and regionEnd are NA)
    ##
    out <- STITCH::shrink_region(
        regionStart = regionStart,
        regionEnd = regionEnd,
        buffer = buffer,
        L = L,
        pos = pos,
        gen = gen,
        phase = phase
    )
    pos <- out$pos
    gen <- out$gen
    phase <- out$phase
    L <- out$L
    nSNPs <- out$nSNPs
    inRegionL <- out$inRegionL
    start_and_end_minus_buffer <- out$start_and_end_minus_buffer


    ##
    ## determine chrStart and chrEnd
    ##
    out <- STITCH::initialize_chrStart_and_chrEnd(
        chrStart = chrStart,
        chrEnd = chrEnd,
        L = L,
        iSizeUpperLimit = iSizeUpperLimit
    )
    chrStart <- out$chrStart
    chrEnd <- out$chrEnd

    ##
    ##
    ## get sample name (one sample name!)
    ##
    ##
    out <- STITCH::get_sample_names(
        bamlist = bamlist,
        cramlist = cramlist,
        save = FALSE
    )
    N <- out$N
    sampleNames <- out$sampleNames
    bam_files <- out$bam_files
    cram_files <- out$cram_files
    if (N > 1) {
        stop(paste0(
            "More than one bam or cram file supplied. ",
            "SEW currently runs only on one sample at a time"
        ))
    }

    ##
    ## get matches to high coverage gen and phase
    ##
    out <- STITCH::match_gen_and_phase_to_samples(
        sampleNames = sampleNames,
        gen = gen,
        phase = phase
    )
    highCovInLow <- out$highCovInLow
    samples_with_phase <- out$samples_with_phase



    ##
    ## if inputBundleBlockSize is not NA
    ## get bundling matrix for input files
    ##
    bundling_info <- STITCH::get_bundling_position_information(
        N = N,
        nCores = nCores,
        blockSize = inputBundleBlockSize
    )

    ##
    ## either generate the data, or load it from before
    ##
    STITCH::generate_or_refactor_input(regenerateInput = regenerateInput, bundling_info = bundling_info, L = L, pos = pos, nSNPs = nSNPs, bam_files = bam_files, cram_files = cram_files, reference = reference, iSizeUpperLimit = iSizeUpperLimit, bqFilter = bqFilter, chr = chr, outputdir = outputdir, N = N, downsampleToCov = downsampleToCov, sampleNames = sampleNames, inputdir = inputdir, useSoftClippedBases = useSoftClippedBases, regionName = regionName, tempdir = tempdir, chrStart = chrStart, chrEnd = chrEnd, generateInputOnly = generateInputOnly, nCores = nCores, save_sampleReadsInfo = save_sampleReadsInfo, use_bx_tag = use_bx_tag, bxTagUpperLimit = bxTagUpperLimit)    

    ##
    ## if necessary, shrink BAMs, but only if regenerateInput = FALSE
    ##
    STITCH::shrinkReads(N = N, nCores = nCores, originalRegionName = originalRegionName, regionName = regionName, bundling_info = bundling_info, tempdir = tempdir, inputdir = inputdir, inRegionL = inRegionL, regenerateInput = regenerateInput, inputBundleBlockSize = inputBundleBlockSize)
    

    ##
    ## downsample to percent
    ##
    if(downsampleFraction < 1) {
        STITCH::downsampleToFraction(N=N,nCores=nCores,downsampleFraction=downsampleFraction,regionName=regionName,tempdir=tempdir, bundling_info = bundling_info)
    }

    ##
    ## at this point (could be earlier), keep sampleReads in RAM if so desired
    ##
    if (keepSampleReadsInRAM) {
        allSampleReads <- STITCH::load_all_sampleReads_into_memory(
            N = N,
            nCores = nCores,
            tempdir = tempdir,
            regionName = regionName,
            bundling_info = bundling_info
        )
    } else {
        allSampleReads <- NULL
    }


    ##
    ## build alleleCount
    ##
    alleleCount <- STITCH::buildAlleleCount(
        nSNPs = nSNPs,
        N = N,
        nCores = nCores,
        regionName = regionName,
        tempdir = tempdir,
        bundling_info = bundling_info,
        allSampleReads = allSampleReads
    )

    ##
    ## do EM phasing
    ##
    phase_set <- calculate_phase_set(inputdir, regionName, L)    
    out <- single_sample_phasing_em_complete(
        inputdir = inputdir,
        regionName = regionName,
        nSNPs = nSNPs,
        phase = phase,
        niterations = niterations,
        unwindIterations = unwindIterations,
        sample_unwindIterations = sample_unwindIterations,
        usePhaseSet_unwindIterations = usePhaseSet_unwindIterations,
        phase_set = phase_set,
        very_verbose = very_verbose
    )
    eHapsCurrent <- out$eHapsCurrent
    eHapsUpdate_numer <- out$eHapsUpdate_numer
    eHapsUpdate_denom <- out$eHapsUpdate_denom
    p_reads <- out$p_reads
    p_reads_given_hap_k <- out$p_reads_given_hap_k
    p_h_given_O <- out$p_h_given_O
    reads_span <- out$reads_span
    kept_reads <- out$kept_reads

    
    ##
    ## calculate metrics
    ##
    phase_set <- recalculate_phase_set(
        inputdir = inputdir,
        regionName = regionName,
        L = L,
        eHapsCurrent = eHapsCurrent,
        original_phase_set = phase_set
    )    
    phase_entropy <- calculate_phase_entropy(eHapsCurrent)
    out <- calculate_strandedness(
        p_h_given_O = p_h_given_O,
        nSNPs = nSNPs,
        L = L,
        inputdir = inputdir,
        regionName = regionName,
        save_sampleReadsInfo = save_sampleReadsInfo,
        kept_reads = kept_reads
    )
    strandedness <- out$strandedness
    strand_info <- out$strand_info



    ##
    ## remove buffer from things
    ##
    if(is.na(regionStart)==FALSE & is.na(regionEnd)==FALSE) {
        save(
            alleleCount, pos, gen, phase, L,
            eHapsCurrent, eHapsUpdate_numer, eHapsUpdate_denom,
            p_reads, p_h_given_O, p_reads_given_hap_k, phase_entropy,
            strandedness, strand_info, reads_span,
            nCores, N, nSNPs, chr, niterations, unwindIterations,
            file = file.path(
                outputdir, "RData", paste0("EM.all.", regionName, ".withBuffer.RData")
            )
        )
        ## for STITCH>=1.6.0, need eHapsCurrent_tc
        eHapsCurrent_tc <- array(NA, c(2, nSNPs, 1))
        eHapsCurrent_tc[, , 1] <- t(eHapsCurrent)
        ##
        out <- STITCH::remove_buffer_from_variables(L = L,  regionStart = regionStart, regionEnd = regionEnd, pos = pos, gen = gen, phase = phase, alleleCount =  alleleCount, highCovInLow = highCovInLow, eHapsCurrent_tc = eHapsCurrent_tc, eHapsUpdate_numer = eHapsUpdate_numer, eHapsUpdate_denom = eHapsUpdate_denom, strandedness = strandedness, phase_entropy = phase_entropy, gridWindowSize = NA)
        phase_set <- remove_buffer_from_phase_set(
            phase_set = phase_set,
            L = L,
            inRegion2 = out$inRegion2
        )
        pos <- out$pos
        gen <- out$gen
        phase <- out$phase
        alleleCount <- out$alleleCount
        L <- out$L
        nSNPs <- out$nSNPs
        ## ugh, workaround for now
        if (length(dim(out$eHapsCurrent_tc)) == 2) {
            ## release 1.6.2 behaviour
            eHapsCurrent <- t(out$eHapsCurrent_tc)
        } else {
            ## after release 1.6.2 behaviour
            eHapsCurrent <- t(out$eHapsCurrent_tc[, , 1, drop = TRUE])
        }
        eHapsUpdate_numer <- out$eHapsUpdate_numer
        eHapsUpdate_denom <- out$eHapsUpdate_denom
        strandedness <- out$strandedness
        phase_entropy <- out$phase_entropy
    }


    ##
    ## save R objects to disk, for faster access
    ##
    print_message("Save RData objects to disk")
    save(
        alleleCount, pos, gen, phase, L,
        eHapsCurrent, eHapsUpdate_numer, eHapsUpdate_denom,
        p_reads, p_h_given_O, p_reads_given_hap_k, phase_entropy,
        strandedness, strand_info, reads_span,
        nCores, N, nSNPs, chr, niterations, unwindIterations,
        file = file.path(
            outputdir, "RData", paste0("EM.all.", regionName, ".RData")
        )
    )
    
    ##
    ## write out VCF here
    ##
    write_sew_vcf_after_em(
        vcf_output_name = vcf_output_name,
        outputdir = outputdir,
        regionName = regionName,
        eHapsCurrent = eHapsCurrent,
        eHapsUpdate_numer = eHapsUpdate_numer,
        eHapsUpdate_denom = eHapsUpdate_denom,
        alleleCount = alleleCount,
        phase_entropy = phase_entropy,
        pos = pos,
        sampleNames = sampleNames,
        phase_set = phase_set,
        command_line = command_line,
        strandedness = strandedness
    )

    print_message("Program done")
    save(date, file = file.path(outputdir, "RData", paste0("end.",regionName,".RData")))
    
    return(NULL)

}
