# Allele-Specific Copy Number Analysis of Tumors

## Description

This repository provides the ASCAT R package that can be used to infer tumour purity, ploidy and allele-specific copy number profiles.

ASCAT is described in detail in: [Allele-specific copy number analysis of tumors. Van Loo P *et al*. *PNAS* (2010)](http://www.ncbi.nlm.nih.gov/pubmed/20837533).

This repository also contains the code underlying additional publication:
[Allele-specific multi-sample copy number segmentation. Ross EM, Haase K, Van Loo P & Markowetz F. *Bioinformatics* (2020)](https://pubmed.ncbi.nlm.nih.gov/32449758).

## Installation
This pre-v3 version: `devtools::install_github('tlesluyes/ascat/ASCAT',ref='v3.0')`

Standard ASCAT: `devtools::install_github('VanLoo-lab/ascat/ASCAT')`

## Changes since v2.5.2
### Major changes:
- Default penalty for both ASPCF (`ascat.aspcf`) and ASmultiPCF (`ascat.asmultipcf`) is now **70** (was 25). It is suitable for SNP arrays, as well as WES and WGS data.
- LogR correction can now be used to correct for both GC content (standard requirement) and replication timing (optional). Also, the correction method has been updated (it now uses autosomes to compute correlations with covariates and applies a linear model with *splines* on all chromosomes). Please note that `ascat.correctLogR` should be used from now on (`ascat.GCcorrect` is still there for backward compatibility but is just a wrapper to `ascat.correctLogR`).
- Color scheme has been changed for CNA profiles so it is now colorblind-friendly:
	- Rounded profiles: red is the major allele and blue is the minor allele.
	- Unrounded profiles: purple is the total CN and green is the minor allele.
- Because ASCAT leverages genomic information from heterozygous SNPs, the nonPAR region in chromosome X for males is challenging as there are no such SNP, as opposed to PAR1 and PAR2 regions being present on chrX and chrY. We improved CNA calling in chrX by considering specificities between nonPAR and PAR1/PAR2. To this end, `ascat.loadData` has a new argument, `genomeVersion` (either `'hg19'` or `'hg38'`), that enables locating the nonPAR region on chrX. If provided, such information will be considered in the different ASCAT functions. We recommend always providing this information so CNA calling on chrX for males will be more accurate. Since PAR1 and PAR2 are present in both chrX and chrY, a 1+1 status in males refers to 1 copy in X and 1 copy in Y, but 1+0 could either be 1 copy in X (and no copy in Y) or 1 copy of Y (and no copy in X). Also, please note that most platforms have a limited resolution for PAR1 and PAR2 so results should carefully be interpreted in respect to available resolution.
- Ploidy value displayed in CNA profiles no longers comes from the grid search and is now the final tumour ploidy (matching with `ascat.output$ploidy`).

### Minor changes:
- `ascat.plotRawData` and `ascat.plotSegmentedData` have an extra argument, `logr.y_values`, to change Y scale for the logR track. Default is: `c(-2,2)`, whereas previous plots were: `c(-1,1)`.
- '*Aberrant cell fraction*' now refers to '*purity*'. For backward compatibility, `ascat.output$aberrantcellfraction` still exists but we encourage using `ascat.output$purity` instead.

### New features in v3.0:
- New set of instructions, as part of the main `ascat.prepareHTS` function, to derive logR and BAF from high-throughput sequencing (HTS) data (WES, WGS & targeted sequencing). Briefly, [alleleCounter](https://github.com/cancerit/alleleCount) is used to get allele counts at specific loci on a pair of tumour/normal (either BAM or CRAM files). This information is then converted into logR and BAF values, based on a similar method than in the [Battenberg package](https://github.com/Wedge-lab/battenberg). Although this method allows running ASCAT on different HTS data:
  - WES: we recommend providing a BED file covering sequenced regions of the genome.
  - WGS: we recommend running [Battenberg](https://github.com/Wedge-lab/battenberg) for accurate clonal and subclonal allele-specific copy-number alteration calling. However, ASCAT can still be used to get a fast purity/ploidy fit (~30 minutes with 12 CPUs from BAMs to CNA profiles). To this end, we provide a set of files that can be used (see *[ReferenceFiles/WGS](ReferenceFiles/WGS)*).
  - Targeted sequencing: data must be preprocessed using the `ascat.preprocessTargSeq` (to be soon implemented) function to extract loci of interest. Then, such curated loci list may be used as part of the `ascat.prepareHTS` function.
  - **For HTS data, gamma must be set to 1 in `ascat.runASCAT`.**
- A new function to collect metrics of interest has been added: `ascat.metrics`.
- Boundaries can be defined for ploidy (min & max) when running `ascat.runAscat` (arguments: `min_ploidy` and `max_ploidy`).
- New function, `ascat.plotAdjustedAscatProfile`, that plots an ASCAT profile in respect to chromosome length (instead of number of heterozygous SNPs).

## Testing
We provide some scripts and input data in the *[ExampleData](ExampleData)* folder.

## Reference files
- LogR correction files (`ascat.correctLogR`) for standard platforms (Affymetrix SNP 6.0, Affymetrix 250k STY, Illumina 660k and Illumina OmniExpress) can be found in the *[ReferenceFiles/SNParrays](ReferenceFiles/SNParrays)* folder. For other platforms, please use our scripts (in *[LogRcorrection](LogRcorrection)*) to generate such correction files.
- For WGS, we provide logR correction files as well as loci and allele files in *[ReferenceFiles/WGS](ReferenceFiles/WGS)*.
- For WES, we recommend using all biallelic SNPs from the 1,000 Genomes Project. One could reproduce the process described in *[ReferenceFiles/WGS](ReferenceFiles/WGS)*, relaxing allele frequency to increase resolution.

## Supported arrays without matched germline
*Custom10k*, *IlluminaASA*, *IlluminaGSAv3*, *Illumina109k*, *IlluminaCytoSNP*, *IlluminaCytoSNP850k*, *Illumina610k*, *Illumina660k*, *Illumina700k*, *Illumina1M*, *Illumina2.5M*, *IlluminaOmni5*, *Affy10k*, *Affy100k*, *Affy250k_sty*, *Affy250k_nsp*, *AffyOncoScan*, *AffyCytoScanHD*, *AffySNP6*, *HumanCNV370quad*, *HumanCore12*, *HumanCoreExome24*, *HumanOmniExpress12* and *IlluminaOmniExpressExome*.

Because arrays have a defined set of SNP probes, with a fairly constant rate of heterozygous probes accross individuals, useful metrics in `ascat.predictGermlineGenotypes` can be inferred from some cases (with no or very few CN changes). However, sequencing data is subjected to massive variations because of design, coverage and/or artefacts. Therefore, we are not able to provide pre-defined metrics for unmatched sequencing data.

## Misc
For more information about ASCAT and other projects of our group, please visit our [website](https://www.crick.ac.uk/research/a-z-researchers/researchers-v-y/peter-van-loo/software/).
