#!/usr/bin/env Rscript

if (!suppressPackageStartupMessages(require("optparse")))
    install.packages("optparse", repos="http://cran.rstudio.com/")

option_list <- list(
    make_option(
        "--chr",
        type = "character",
        help = "What chromosome to run. Should match BAM headers"
    ), 
    make_option(
        "--posfile",
        type = "character",
        help = "Where to find file with positions to run. File is tab seperated with no header, one row per SNP, with col 1 = chromosome, col 2 = physical position (sorted from smallest to largest), col 3 = reference base, col 4 = alternate base. Bases are capitalized. Example first row: 1<tab>1000<tab>A<tab>G<tab>"
    ), 
    make_option(
        "--outputdir",
        type = "character",
        help = "What output directory to use"
    ), 
    make_option(
        "--tempdir",
        type = "character",
        help = "What directory to use as temporary directory. If possible, use ramdisk, like /dev/shm/ [default NA] ",
        default = NA
    ), 
    make_option(
        "--bamlist",
        type = "character",
        help = "Path to file with bam file locations. File is one row per entry, path to bam files. Bam index files should exist in same directory as for each bam, suffixed either .bam.bai or .bai [default \"\"] ",
        default = ""
    ), 
    make_option(
        "--cramlist",
        type = "character",
        help = "Path to file with cram file locations. File is one row per entry, path to cram files. cram files are converted to bam files on the fly for parsing into SEW [default \"\"] ",
        default = ""
    ), 
    make_option(
        "--reference",
        type = "character",
        help = "Path to reference fasta used for making cram files. Only required if cramlist is defined [default \"\"] ",
        default = ""
    ), 
    make_option(
        "--genfile",
        type = "character",
        help = "Path to gen file with high coverage results. Empty for no genfile. File has a header row with a name for each sample, matching what is found in the bam file. Each subject is then a tab seperated column, with 0 = hom ref, 1 = het, 2 = hom alt and NA indicating missing genotype, with rows corresponding to rows of the posfile. Note therefore this file has one more row than posfile which has no header [default \"\"] ",
        default = ""
    ), 
    make_option(
        "--phasefile",
        type = "character",
        help = "Path to phase file with truth phasing results. Empty for no phasefile. File has a header row with a name for each sample, matching what is found in the bam file. Each subject is then a tab seperated column, with 0 = ref and 1 = alt, separated by a vertical bar |, e.g. 0|0 or 0|1. Note therefore this file has one more row than posfile which has no header [default \"\"] ",
        default = ""
    ), 
    make_option(
        "--outputHaplotypeProbabilities",
        type = "logical",
        help = "Whether to output haplotype probabilities in files    [default FALSE] ",
        default = FALSE
    ), 
    make_option(
        "--regenerateInput",
        type = "logical",
        help = "Whether to regenerate input files [default TRUE] ",
        default = TRUE
    ), 
    make_option(
        "--originalRegionName",
        type = "character",
        help = "If regenerateInput is FALSE (i.e. using existing data), this is the name of the original region name (chr.regionStart.regionEnd). This is necessary to load past variables [default NA] ",
        default = NA
    ), 
    make_option(
        "--keepInterimFiles",
        type = "logical",
        help = "Whether to keep interim parameter estimates [default FALSE] ",
        default = FALSE
    ), 
    make_option(
        "--keepTempDir",
        type = "logical",
        help = "Whether to keep files in temporary directory [default FALSE] ",
        default = FALSE
    ), 
    make_option(
        "--generateInputOnly",
        type = "logical",
        help = "Whether to just generate input data then quit [default FALSE] ",
        default = FALSE
    ), 
    make_option(
        "--useSoftClippedBases",
        type = "logical",
        help = "Whether to use (TRUE) or not use (FALSE) bases in soft clipped portions of reads [default FALSE] ",
        default = FALSE
    ), 
    make_option(
        "--vcf_output_name",
        type = "character",
        help = "Override the default VCF output name with this given file name. Please note that this does not change the names of inputs or outputs (e.g. RData, plots), so if outputdir is unchanged and if multiple SEW runs are processing on the same region then they may over-write each others inputs and outputs [default NULL] ",
        default = NULL
    ), 
    make_option(
        "--initial_min_hapProb",
        type = "double",
        help = "Initial lower bound for probability read comes from haplotype. Double bounded between 0 and 1 [default 0.4] ",
        default = 0.4
    ), 
    make_option(
        "--initial_max_hapProb",
        type = "double",
        help = "Initial upper bound for probability read comes from haplotype. Double bounded between 0 and 1 [default 0.6] ",
        default = 0.6
    ), 
    make_option(
        "--regenerateInputWithDefaultValues",
        type = "logical",
        help = "If regenerateInput is FALSE and the original input data was made using regionStart, regionEnd and buffer as default values, set this equal to TRUE [default FALSE] ",
        default = FALSE
    ), 
    make_option(
        "--unwindIterations",
        type = "character",
        help = "What iterations to search for flips in estimated haplotypes that increase the likelihood [default c(40, 70, 100, 120, 140, 160, 180, 200)] ",
        default = "c(40, 70, 100, 120, 140, 160, 180, 200)"
    ), 
    make_option(
        "--sample_unwindIterations",
        type = "character",
        help = "When performing an unwind iterations, how often to check a SNP for an improvement following unwinding. 1 means check all SNPs [default c(5, 5, 1, 1, 1, 1, 1, 1)] ",
        default = "c(5, 5, 1, 1, 1, 1, 1, 1)"
    ), 
    make_option(
        "--save_sampleReadsInfo",
        type = "logical",
        help = "Experimental. Boolean TRUE/FALSE about whether to save additional information about the reads that were extracted [default TRUE] ",
        default = TRUE
    ), 
    make_option(
        "--outputInputInVCFFormat",
        type = "logical",
        help = "Whether to output the input in vcf format [default FALSE] ",
        default = FALSE
    ), 
    make_option(
        "--downsampleToCov",
        type = "double",
        help = "What coverage to downsample individual sites to. This ensures no floating point errors at sites with really high coverage [default 200] ",
        default = 200
    ), 
    make_option(
        "--downsampleFraction",
        type = "double",
        help = "Downsample BAMs by choosing a fraction of reads to retain. Must be value 0<downsampleFraction<1 [default 1] ",
        default = 1
    ), 
    make_option(
        "--chrStart",
        type = "integer",
        help = "When loading from BAM, some start position, before SNPs occur. Default NA will infer this from either regionStart, regionEnd and buffer, or posfile [default NA] ",
        default = NA
    ), 
    make_option(
        "--chrEnd",
        type = "integer",
        help = "When loading from BAM, some end position, after SNPs occur. Default NA will infer this from either regionStart, regionEnd and buffer, or posfile [default NA] ",
        default = NA
    ), 
    make_option(
        "--regionStart",
        type = "integer",
        help = "When running imputation, where to start from. The 1-based position x is kept if regionStart <= x <= regionEnd [default NA] ",
        default = NA
    ), 
    make_option(
        "--regionEnd",
        type = "integer",
        help = "When running imputation, where to stop [default NA] ",
        default = NA
    ), 
    make_option(
        "--buffer",
        type = "integer",
        help = "Buffer of region to perform imputation over. So imputation is run form regionStart-buffer to regionEnd+buffer, and reported for regionStart to regionEnd, including the bases of regionStart and regionEnd [default NA] ",
        default = NA
    ), 
    make_option(
        "--iSizeUpperLimit",
        type = "double",
        help = "Do not use reads with an insert size of more than this value [default as.integer(600)] ",
        default = as.integer(600)
    ), 
    make_option(
        "--bqFilter",
        type = "double",
        help = "Minimum BQ for a SNP in a read. Also, the algorithm uses bq<=mq, so if mapping quality is less than this, the read isnt used [default as.integer(-1)] ",
        default = as.integer(-1)
    ), 
    make_option(
        "--niterations",
        type = "integer",
        help = "Number of EM iterations. [default 300] ",
        default = 300
    ), 
    make_option(
        "--plotHapSumDuringIterations",
        type = "logical",
        help = "Boolean TRUE/FALSE about whether to make a plot that shows the relative number of individuals using each ancestral haplotype in each iteration [default FALSE] ",
        default = FALSE
    ), 
    make_option(
        "--very_verbose",
        type = "logical",
        help = "Boolean TRUE/FALSE about whether to print additional values which may be helpeful for understanding how the model operates [default FALSE] ",
        default = FALSE
    ), 
    make_option(
        "--keepSampleReadsInRAM",
        type = "logical",
        help = "Whether to (generally) keep sampleReads in RAM or store them in the temporary directory. SEW will be faster if this is FALSE at the expense of RAM [default FALSE] ",
        default = FALSE
    ), 
    make_option(
        "--use_bx_tag",
        type = "logical",
        help = "Whether to try and use BX tag in same to indicate that reads come from the same underlying molecule [default TRUE] ",
        default = TRUE
    ), 
    make_option(
        "--bxTagUpperLimit",
        type = "integer",
        help = "When using BX tag, at what distance between reads to consider reads with the same BX tag to come from different molecules [default 50000    ] ",
        default = 50000    
    )
)
opt <- suppressWarnings(parse_args(OptionParser(option_list = option_list)))
suppressPackageStartupMessages(library(SEW))
Sys.setenv(PATH = paste0(Sys.getenv("PATH"), ":", getwd()))
SEW(
    chr = opt$chr,
    posfile = opt$posfile,
    outputdir = opt$outputdir,
    tempdir = opt$tempdir,
    bamlist = opt$bamlist,
    cramlist = opt$cramlist,
    reference = opt$reference,
    genfile = opt$genfile,
    phasefile = opt$phasefile,
    outputHaplotypeProbabilities = opt$outputHaplotypeProbabilities,
    regenerateInput = opt$regenerateInput,
    originalRegionName = opt$originalRegionName,
    keepInterimFiles = opt$keepInterimFiles,
    keepTempDir = opt$keepTempDir,
    generateInputOnly = opt$generateInputOnly,
    useSoftClippedBases = opt$useSoftClippedBases,
    vcf_output_name = opt$vcf_output_name,
    initial_min_hapProb = opt$initial_min_hapProb,
    initial_max_hapProb = opt$initial_max_hapProb,
    regenerateInputWithDefaultValues = opt$regenerateInputWithDefaultValues,
    unwindIterations = eval(parse(text=opt$unwindIterations)),
    sample_unwindIterations = eval(parse(text=opt$sample_unwindIterations)),
    save_sampleReadsInfo = opt$save_sampleReadsInfo,
    outputInputInVCFFormat = opt$outputInputInVCFFormat,
    downsampleToCov = opt$downsampleToCov,
    downsampleFraction = opt$downsampleFraction,
    chrStart = opt$chrStart,
    chrEnd = opt$chrEnd,
    regionStart = opt$regionStart,
    regionEnd = opt$regionEnd,
    buffer = opt$buffer,
    iSizeUpperLimit = opt$iSizeUpperLimit,
    bqFilter = opt$bqFilter,
    niterations = opt$niterations,
    plotHapSumDuringIterations = opt$plotHapSumDuringIterations,
    very_verbose = opt$very_verbose,
    keepSampleReadsInRAM = opt$keepSampleReadsInRAM,
    use_bx_tag = opt$use_bx_tag,
    bxTagUpperLimit = opt$bxTagUpperLimit
)
