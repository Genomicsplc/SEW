### Required

* **chr** What chromosome to run. Should match BAM header
* **posfile**  Where to find file with positions to run. File is tab separated with no header, one row per SNP, with col 1 = chromosome, col 2 = physical position (sorted from smallest to largest), col 3 = reference base, col 4 = alternate base. Bases are capitalized. STITCH only handles bi-allelic SNPs
* **outputdir**  What output directory to use / where output files go
* **bamlist** (one of bamlist or cramlist required) Path to file with BAM file location. File is one row per entry, path to BAM files. BAM index files should exist in same directory as for each BAM, suffixed either .bam.bai or .bai

### Optional

* **niterations** Number of EM iterations.
* **unwindIterations** What iterations to search for flips in estimated haplotypes that increase the likelihood. See main README or paper supplement for more details
* **phasefile** Path to phase file with truth phasing results. Empty for no phasefile. File has a header row with a name for each sample, matching what is found in the bam file. Each subject is then a tab seperated column, with 0 = ref and 1 = alt, separated by a vertical bar |, e.g. 0|0 or 0|1. Note therefore this file has one more row than posfile which has no header
* **bqFilter**  Minimum BQ for a SNP in a read. Also, the algorithm uses bq<=mq, so if mapping quality is less than this, the read isnt used
* *downsampleToCov**  What coverage to downsample individual sites to. This ensures no floating point errors at sites with really high coverage
* **regionStart** When running imputation, where to start from
* **regionEnd** When running imputation, where to stop
* **buffer** Buffer of region to perform imputation over. Imputation is run from bases including regionStart - buffer to regionEnd + buffer, including the bases, with 1-based positions. After imputation, the VCF is shrunk to only include positions from regionStart to regionEnd, inclusive

### Output

* VCF named <outputdir>sew.<chr>.<regionStart>.<regionEnd>.vcf.gz, or if no regionStart and regionEnd is given <outputdir>sew.<chr>.vcf.gz
* **vcf_output_name** overrides default name to give VCF output