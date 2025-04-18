# Allele-Specific Copy Number Analysis of Tumors

## Description

This repository provides the ASCAT R package (v3.2.0) that can be used to infer tumour purity, ploidy and allele-specific copy number profiles.

ASCAT is described in detail in: [Allele-specific copy number analysis of tumors. Van Loo P *et al*. *PNAS* (2010)](https://pubmed.ncbi.nlm.nih.gov/20837533).

This repository also contains the code underlying additional publication:
[Allele-specific multi-sample copy number segmentation. Ross EM, Haase K, Van Loo P & Markowetz F. *Bioinformatics* (2020)](https://pubmed.ncbi.nlm.nih.gov/32449758).

## Installation (v3.2.0 version)
Bioconductor package dependencies: [GenomicRanges](https://bioconductor.org/packages/release/bioc/html/GenomicRanges.html) & [IRanges](https://bioconductor.org/packages/release/bioc/html/IRanges.html) (`BiocManager::install(c('GenomicRanges','IRanges'))` with a recent R/BiocManager version).

Processing high-throughput sequencing data: [alleleCounter](https://github.com/cancerit/alleleCount) (C version)

Installing ASCAT using R: `devtools::install_github('VanLoo-lab/ascat/ASCAT')`

## Changes since v2.5.3
### Major changes:
- Default penalty for both ASPCF (`ascat.aspcf`) and ASmultiPCF (`ascat.asmultipcf`) is now **70** (was 25). It is suitable for SNP arrays, as well as WES and WGS data.
- LogR correction can now be used to correct for both GC content (standard requirement) and replication timing (optional). Also, the correction method has been updated (it now uses autosomes to compute correlations with covariates and applies a linear model with *splines* on all chromosomes). Please note that `ascat.correctLogR` should be used from now on (`ascat.GCcorrect` is still there for backward compatibility but is just a wrapper to `ascat.correctLogR`).
- Color scheme has been changed for CNA profiles so it is now colorblind-friendly:
	- Rounded profiles: red is the major allele and blue is the minor allele.
	- Unrounded profiles: purple is the total CN and green is the minor allele.
	- ASPCF plots: red dots are raw values and blue dots are segmented values.
- Because ASCAT leverages genomic information from heterozygous SNPs, the nonPAR region in chromosome X for males is challenging as there are no such SNP, as opposed to PAR1 and PAR2 regions being present on chrX and chrY. We improved CNA calling in chrX by considering specificities between nonPAR and PAR1/PAR2. To this end, `ascat.loadData` has a new argument, `genomeVersion` (either `'hg19'` or `'hg38'`), that enables locating the nonPAR region on chrX. If provided, such information will be considered in the different ASCAT functions. We recommend always providing this information so CNA calling on chrX for males will be more accurate. Since PAR1 and PAR2 are present in both chrX and chrY, a 1+1 status in males refers to 1 copy in X and 1 copy in Y, but 1+0 could either be 1 copy in X (and no copy in Y) or 1 copy of Y (and no copy in X). Also, please note that most platforms have a limited resolution for PAR1 and PAR2 so results should carefully be interpreted in respect to available resolution.
- Ploidy value displayed in CNA profiles no longer comes from the grid search and is now the final tumour ploidy (matching with `ascat.output$ploidy`).

### Minor changes:
- `ascat.plotRawData` and `ascat.plotSegmentedData` have an extra argument, `logr.y_values`, to change Y scale for the logR track. Default is: `c(-2,2)`, whereas previous plots were: `c(-1,1)`.
- '*Aberrant cell fraction*' now refers to '*purity*'. For backward compatibility, `ascat.output$aberrantcellfraction` still exists but we encourage using `ascat.output$purity` instead.
- The density of logR and BAF tracks are now shown in dark red.

### New features in v3:
- New set of instructions, as part of the main `ascat.prepareHTS` function, to derive logR and BAF from high-throughput sequencing (HTS) data. Briefly, [alleleCounter](https://github.com/cancerit/alleleCount) is used to get allele counts at specific loci on a pair of tumour/normal BAM (or CRAM) files, or from a tumor only BAM (or CRAM). This information is then converted into logR and (mirrored) BAF values, based on a similar method than in the [Battenberg package](https://github.com/Wedge-lab/battenberg). Although this method allows running ASCAT on different HTS data:
  - **WES**: we recommend providing a BED file covering sequenced regions of the genome.
  - **WGS**: we recommend running [Battenberg](https://github.com/Wedge-lab/battenberg) for accurate clonal and subclonal allele-specific copy-number alteration calling. However, ASCAT can still be used to get a fast purity/ploidy fit (~30 minutes with 12 CPUs from BAMs to CNA profiles). To this end, we provide a set of files that can be used (see *[ReferenceFiles/WGS](ReferenceFiles/WGS)*).
  - **Targeted sequencing (TS)**: a bespoke function, `ascat.prepareTargetedSeq` has been implemented. Such a function must be run on a batch of normals (no tumours) and will identify high-quality SNPs to investigate. Then, `ascat.prepareHTS` can be used on selected SNPs to process tumour/normal pairs. Because of sparse datapoints, we recommend using `penalty=25` when running `ascat.aspcf`. `ascat.prepareTargetedSeq` was further fine-tuned in the v3.1.1 release and now uses a probabilistic method to infer genotypes (hom/het/noisy) based on counts (instead of fixed VAF in v3.1.0).
  - **For HTS data (WGS, WES and TS), gamma must be set to 1 in `ascat.runASCAT`.**
- A new function to collect metrics of interest has been added: `ascat.metrics`.
- Boundaries can be defined for purity and ploidy (min & max) when running `ascat.runAscat` (arguments: `min_purity`/`max_purity` and `min_ploidy`/`max_ploidy`).
- New function, `ascat.plotAdjustedAscatProfile`, that plots an ASCAT profile with respect to chromosome length (instead of the number of heterozygous SNPs).
- For sequencing data processed with `ascat.prepareHTS`, ASCAT now reports raw (=unmirrored) BAF so true BAF values can be used in downstream analyses (*e.g.* rephase). Please note that such information is available in a couple of new files called `*_rawBAF.txt` (based on `tumourBAF_file` and `normalBAF_file`).

## Testing
We provide some scripts and input data in the *[ExampleData](ExampleData)* folder.

## Reference files
All reference files are hosted on [Zenodo](https://zenodo.org/records/14008443).
- LogR correction files (`ascat.correctLogR`) for standard platforms (Affymetrix SNP 6.0, Affymetrix 250k STY, Illumina 660k and Illumina OmniExpress) can be found in the *[ReferenceFiles/SNParrays](ReferenceFiles/SNParrays)* folder. For other platforms, please use our scripts (in *[LogRcorrection](LogRcorrection)*) to generate such correction files.
- For WGS, we provide logR correction files as well as loci and allele files in *[ReferenceFiles/WGS](ReferenceFiles/WGS)*.
- For WES and TS, we provide logR correction files as well as loci and allele files in: *[ReferenceFiles/WES](ReferenceFiles/WES)*. Please note that reference files for WES and TS contain way more SNPs than the ones for WGS. This is because they will be downsampled so we need to provide an exhaustive list of SNPs to begin with. Do not use such a list for processing WGS data and do not use reference files for WGS to process WES/TS data.

## Supported arrays without matched germline
*Custom10k*, *IlluminaASA*, *IlluminaGSAv3*, *Illumina109k*, *IlluminaCytoSNP*, *IlluminaCytoSNP850k*, *Illumina610k*, *Illumina660k*, *Illumina700k*, *Illumina1M*, *Illumina2.5M*, *IlluminaOmni5*, *IlluminaOmniExpressExome*, *IlluminaGDACyto-8*, *Affy10k*, *Affy100k*, *Affy250k_sty*, *Affy250k_nsp*, *AffyOncoScan*, *AffyCytoScanHD*, *AffySNP6*, *HumanCNV370quad*, *HumanCore12*, *HumanCoreExome24* and *HumanOmniExpress12*.

Because arrays have a defined set of SNP probes, with a fairly constant rate of heterozygous probes across individuals, useful metrics in `ascat.predictGermlineGenotypes` can be inferred from some cases (with no or very few CN changes). However, sequencing data is subjected to massive variations because of design, coverage and/or artefacts. Therefore, we are not able to provide pre-defined metrics for unmatched sequencing data.

We now provide a preset for WGS data under specific conditions: hg38 assembly and >50x coverage. Although this preset might be used in some other conditions, please note that is has not been extensively benchmarked under other conditions so we do not provide any guarantee outside of the scope. Such a preset is called *WGS_hg38_50X*.

## Misc
For more information about ASCAT and other projects of our group, please visit our [website](https://www.crick.ac.uk/research/a-z-researchers/researchers-v-y/peter-van-loo/software/).

# Changes to let ASCAT run on long-read data:
- When running `ascat.prepareHTS()`, the following two parameters must be set (options for alleleCounter): `min_base_qual=10` and `additional_allelecounter_flags="-f 0"`.
- Also, `loci_binsize` must be set to a higher value (default=1, no binning). This activates the binning process and a value of 500 (maybe higher, depending on the average length of your reads) works well in our experience. This reduces autocorrelation in BAF/LogR for long-read sequencing.
- If there was no PCR step in the library preparation, GC correction is not needed.
- Increasing the penalty value in `ascat.aspcf` can help in reducing the noise.
- Example of the new `ascat.prepareHTS` function for long-reads sequencing:

```
ascat.prepareHTS(
  tumourseqfile = tumour_BAM,
  normalseqfile = normal_BAM,
  tumourname = name_tumour,
  normalname = name_normal,
  allelecounter_exe = allelecounter,
  skip_allele_counting_normal = FALSE,
  skip_allele_counting_tumour = FALSE,
  alleles.prefix = G1000_alleles_hg38_chr,
  loci.prefix = G1000_loci_hg38_chr,
  gender = gender,
  genomeVersion = "hg38",
  nthreads = 12,
  tumourLogR_file = "Tumor_LogR.txt",
  tumourBAF_file = "Tumor_BAF.txt",
  loci_binsize = 500,
  min_base_qual= 10,
  additional_allelecounter_flags="-f 0")
```
