# Allele-Specific Copy Number Analysis of Tumors

## Description

This repository provides the ASCAT R package that can be used to infer tumour purity, ploidy and allele-specific copy number profiles.

ASCAT is described in detail in: [Allele-specific copy number analysis of tumors. Van Loo P *et al*. *PNAS* (2010)](http://www.ncbi.nlm.nih.gov/pubmed/20837533).

This repository also contains the code underlying additional publication:
[Allele-specific multi-sample copy number segmentation. Ross EM, Haase K, Van Loo P & Markowetz F. *Bioinformatics* (2020)](https://pubmed.ncbi.nlm.nih.gov/32449758).

## Future v3 release (beta-test version)
We are pleased to announce that an early version of the future v3 release is available here: [https://github.com/VanLoo-lab/ascat/tree/v3.0](https://github.com/VanLoo-lab/ascat/tree/v3.0).

This beta version comes with major improvements and new features. Therefore, we recommend users to give it a try, as early adopters, before it becomes the latest (and default) version of ASCAT.

## Installation
`devtools::install_github('VanLoo-lab/ascat/ASCAT')`

### Testing
We provide some scripts and input data in the *ExampleData* folder.

### Output
After running the *ascat.runAscat* function:
`ascat.output$segments`

### GC correction script
We provide a method that creates a GC correction file in the *gcProcessing* folder.

### Supported arrays without matched germline
*Custom10k*, *IlluminaASA*, *IlluminaGSAv3*, *Illumina109k*, *IlluminaCytoSNP*, *IlluminaCytoSNP850k*, *Illumina610k*, *Illumina660k*, *Illumina700k*, *Illumina1M*, *Illumina2.5M*, *IlluminaOmni5*, *Affy10k*, *Affy100k*, *Affy250k_sty*, *Affy250k_nsp*, *AffyOncoScan*, *AffyCytoScanHD*, *AffySNP6*, *HumanCNV370quad*, *HumanCore12*, *HumanCoreExome24*, *HumanOmniExpress12* and *IlluminaOmniExpressExome*.

### Misc
For more information about ASCAT and other projects of our group, please visit our [website](https://www.crick.ac.uk/research/a-z-researchers/researchers-v-y/peter-van-loo/software/).
