# Allele-Specific Copy Number Analysis of Tumors

## Description

This repository provides the ASCAT R package that can be used to infer tumour purity, ploidy and
allele-specific copy number profiles.

ASCAT is described in detail in: [Allele-specific copy number analysis of tumors. Van Loo P *et al*. *PNAS* (2010)](http://www.ncbi.nlm.nih.gov/pubmed/20837533).

This repository also contains the code underlying additional publication:
[Allele-specific multi-sample copy number segmentation. Ross EM, Haase K, Van Loo P & Markowetz F. *Bioinformatics* (2020)](https://pubmed.ncbi.nlm.nih.gov/32449758).

## Installation
`devtools::install_github('Crick-CancerGenomics/ascat/ASCAT')`

## Major changes since v2.5.2
- Default penalty for ASPCF is now 70 (was 25).
- LogR correction (*ascat.GCcorrect*) can now be used to correct for both GC content and replication timing. Also, the correction method has been updated (now uses *splines*).
- New set of instructions, as part of the main *ascat.prepareWES* function, to prepare high-throughput sequencing (HTS) data (WES & targeted sequencing). See example in he *ExampleData* folder. Briefly, [alleleCounter](https://github.com/cancerit/alleleCount) is used to get allele counts at 1000G SNP positions on a pair of BAM/CRAM files. This information is then converted into logR and BAF values. This method is derived from the [Battenberg package](https://github.com/Wedge-lab/battenberg). Although this method allows running ASCAT on WGS data, we recommand running Battenberg for WGS and ASCAT for WES and targeted sequencing.
- A QC function (*ascat.QC*) has been added.

### Testing
We provide some scripts and input data in the *ExampleData* folder.

### GC correction script
We provide a method that creates a GC correction file in the *gcProcessing* folder.

### Supported arrays without matched germline
*Custom10k*, *IlluminaASA*, *IlluminaGSAv3*, *Illumina109k*, *IlluminaCytoSNP*, *IlluminaCytoSNP850k*, *Illumina610k*, *Illumina660k*, *Illumina700k*, *Illumina1M*, *Illumina2.5M*, *IlluminaOmni5*, *Affy10k*, *Affy100k*, *Affy250k_sty*, *Affy250k_nsp*, *AffyOncoScan*, *AffyCytoScanHD*, *AffySNP6*, *HumanCNV370quad*, *HumanCore12*, *HumanCoreExome24*, *HumanOmniExpress12* and *IlluminaOmniExpressExome*.

### Misc
For more information about ASCAT and other projects of our group, please visit our [website](https://www.crick.ac.uk/research/a-z-researchers/researchers-v-y/peter-van-loo/software/).
