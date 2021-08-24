# Allele-Specific Copy Number Analysis of Tumors

## Description

This repository provides the ASCAT R package that can be used to infer tumour purity, ploidy and
allele-specific copy number profiles.

ASCAT is described in detail in: [Allele-specific copy number analysis of tumors. Van Loo P *et al*. *PNAS* (2010)](http://www.ncbi.nlm.nih.gov/pubmed/20837533).

This repository also contains the code underlying additional publication:
[Allele-specific multi-sample copy number segmentation. Ross EM, Haase K, Van Loo P & Markowetz F. *Bioinformatics* (2020)](https://pubmed.ncbi.nlm.nih.gov/32449758).

## Installation
`devtools::install_github('VanLoo-lab/ascat/ASCAT')`

## Changes since v2.5.2
### Major changes:
- Default penalty for both ASPCF (*ascat.aspcf*) and ASmultiPCF (*ascat.asmultipcf*) is now 70 (was 25).
- LogR correction (*ascat.GCcorrect*) can now be used to correct for both GC content (standard requirement) and replication timing (optional). Also, the correction method has been updated (it now uses a linear model with *splines*).
- Color scheme has been changed for CNA profiles:
	- Rounded profiles: TBD1 is the major allele and TBD2 is the minor allele.
	- Unrounded profiles: TBD3 is the total CN and TBD4 is the minor allele.

### Minor changes:
- *ascat.plotRawData* and *ascat.plotSegmentedData* have an extra parameter (*logr.y_values*) to change Y scale for the logR track. Default is: c(-2,2), whereas previous plot were: c(-1,1).

### New features in v3.0:
- New set of instructions, as part of the main *ascat.prepareHTS* function, to derive logR and BAF from high-throughput sequencing (HTS) data (WES, WGS & targeted sequencing). Briefly, [alleleCounter](https://github.com/cancerit/alleleCount) is used to get allele counts at specific loci on a pair of tumour/normal (either BAM or CRAM files). This information is then converted into logR and BAF values, based on a similar method than in the [Battenberg package](https://github.com/Wedge-lab/battenberg). Although this method allows running ASCAT on different HTS data:
  - WES: we recommand providing a BED file covering sequenced regions of the genome.
  - WGS: we recommend running [Battenberg](https://github.com/Wedge-lab/battenberg) for accurate clonal and subclonal allele-specific copy-number alteration calling. However, ASCAT can still be used to get a fast purity/ploidy fit (~30 minutes with 12 CPUs from BAMs to CNA profiles). To this end, we provide a set of files that can be used (see *ReferenceFiles/WGS*).
  - Targeted sequencing: data must be preprocessed using the *ascat.preprocessTargSeq* function to extract loci of interest. Then, such curated loci list may be used as part of the *ascat.prepareHTS* function.
  - **For HTS data, gamma must be set to 1 in *ascat.runASCAT*.**
- A new function to collect metrics of interest has been added: *ascat.metrics*.

## Testing
We provide some scripts and input data in the *ExampleData* folder.

## Reference files for logR correction
We provide scripts to generate correction files for any platform in the *LogRcorrection* folder.

## Supported arrays without matched germline
*Custom10k*, *IlluminaASA*, *IlluminaGSAv3*, *Illumina109k*, *IlluminaCytoSNP*, *IlluminaCytoSNP850k*, *Illumina610k*, *Illumina660k*, *Illumina700k*, *Illumina1M*, *Illumina2.5M*, *IlluminaOmni5*, *Affy10k*, *Affy100k*, *Affy250k_sty*, *Affy250k_nsp*, *AffyOncoScan*, *AffyCytoScanHD*, *AffySNP6*, *HumanCNV370quad*, *HumanCore12*, *HumanCoreExome24*, *HumanOmniExpress12* and *IlluminaOmniExpressExome*.

## Misc
For more information about ASCAT and other projects of our group, please visit our [website](https://www.crick.ac.uk/research/a-z-researchers/researchers-v-y/peter-van-loo/software/).
