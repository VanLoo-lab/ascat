## Standard ASCAT run

```R
library(ASCAT)
ascat.bc = ascat.loadData(Tumor_LogR_file = "Tumor_LogR.txt", Tumor_BAF_file = "Tumor_BAF.txt", Germline_LogR_file = "Germline_LogR.txt", Germline_BAF_file = "Germline_BAF.txt", gender = rep('XX',100), genomeVersion = "hg19") # isTargetedSeq=T for targeted sequencing data
ascat.plotRawData(ascat.bc, img.prefix = "Before_correction_")
ascat.bc = ascat.correctLogR(ascat.bc, GCcontentfile = "GC_example.txt", replictimingfile = "RT_example.txt")
ascat.plotRawData(ascat.bc, img.prefix = "After_correction_")
ascat.bc = ascat.aspcf(ascat.bc) # penalty=25 for targeted sequencing data
ascat.plotSegmentedData(ascat.bc)
ascat.output = ascat.runAscat(ascat.bc, write_segments = T) # gamma=1 for HTS data
QC = ascat.metrics(ascat.bc,ascat.output)
save(ascat.bc, ascat.output, QC, file = 'ASCAT_objects.Rdata')
```

## Minimal ASCAT run (adapted from the standard run: no logR correction, no intermediate plot, no QC)

```R
library(ASCAT)
ascat.bc = ascat.loadData(Tumor_LogR_file = "Tumor_LogR.txt", Tumor_BAF_file = "Tumor_BAF.txt", Germline_LogR_file = "Germline_LogR.txt", Germline_BAF_file = "Germline_BAF.txt", gender = rep('XX',100), genomeVersion = "hg19")
ascat.bc = ascat.aspcf(ascat.bc, out.dir = NA) # penalty=25 for targeted sequencing data
ascat.output = ascat.runAscat(ascat.bc) # gamma=1 for HTS data
save(ascat.bc,ascat.output,file='ASCAT_objects.Rdata')
```


## ASCAT run without matched normal data (`platform` needs to be adapted, see `?ascat.predictGermlineGenotypes`)

```R
library(ASCAT)
ascat.bc = ascat.loadData(Tumor_LogR_file = "Tumor_LogR.txt", Tumor_BAF_file = "Tumor_BAF.txt", gender = rep('XX',100), genomeVersion = "hg19")
ascat.plotRawData(ascat.bc, img.prefix = "Before_correction_")
ascat.bc = ascat.correctLogR(ascat.bc, GCcontentfile = "GC_example.txt", replictimingfile = "RT_example.txt")
ascat.plotRawData(ascat.bc, img.prefix = "After_correction_")
gg = ascat.predictGermlineGenotypes(ascat.bc, platform = "AffySNP6")
ascat.bc = ascat.aspcf(ascat.bc, ascat.gg=gg)
ascat.plotSegmentedData(ascat.bc)
ascat.output = ascat.runAscat(ascat.bc, write_segments = T)
QC = ascat.metrics(ascat.bc,ascat.output)
save(ascat.bc, ascat.output, QC, file = 'ASCAT_objects.Rdata')
```

## ASCAT run with multi-sample segmentation (when shared breakpoints are expected)

```R
library(ASCAT)
ascat.bc = ascat.loadData(Tumor_LogR_file = "Tumor_LogR.txt", Tumor_BAF_file = "Tumor_BAF.txt", Germline_LogR_file = "Germline_LogR.txt", Germline_BAF_file = "Germline_BAF.txt", gender = rep('XX',100), genomeVersion = "hg19")
ascat.plotRawData(ascat.bc, img.prefix = "Before_correction_")
ascat.bc = ascat.correctLogR(ascat.bc, GCcontentfile = "GC_example.txt", replictimingfile = "RT_example.txt")
ascat.plotRawData(ascat.bc, img.prefix = "After_correction_")
ascat.bc = ascat.asmultipcf(ascat.bc)
ascat.plotSegmentedData(ascat.bc)
ascat.output = ascat.runAscat(ascat.bc, write_segments = T) # gamma=1 for HTS data
QC = ascat.metrics(ascat.bc,ascat.output)
save(ascat.bc, ascat.output, QC, file = 'ASCAT_objects.Rdata')
```

## Getting CNA profiles from CEL files using ASCAT

[ASCAT_fromCELfiles.R](ASCAT_fromCELfiles.R)

## Extracting logR and BAF from HTS data and running ASCAT

```R
library(ASCAT)

ascat.prepareHTS(
  tumourseqfile = "Tumour.bam",
  normalseqfile = "Normal.bam",
  tumourname = "Tumour_name",
  normalname = "Normal_name",
  allelecounter_exe = "/PATH/TO/allelecounter",
  alleles.prefix = "G1000_alleles_hg19_chr",
  loci.prefix = "G1000_loci_hg19_chr",
  gender = "XX",
  genomeVersion = "hg19",
  nthreads = 8,
  tumourLogR_file = "Tumor_LogR.txt",
  tumourBAF_file = "Tumor_BAF.txt",
  normalLogR_file = "Germline_LogR.txt",
  normalBAF_file = "Germline_BAF.txt")

ascat.bc = ascat.loadData(Tumor_LogR_file = "Tumor_LogR.txt", Tumor_BAF_file = "Tumor_BAF.txt", Germline_LogR_file = "Germline_LogR.txt", Germline_BAF_file = "Germline_BAF.txt", gender = 'XX', genomeVersion = "hg19")
ascat.plotRawData(ascat.bc, img.prefix = "Before_correction_")
ascat.bc = ascat.correctLogR(ascat.bc, GCcontentfile = "GC_file.txt", replictimingfile = "RT_file.txt")
ascat.plotRawData(ascat.bc, img.prefix = "After_correction_")
ascat.bc = ascat.aspcf(ascat.bc)
ascat.plotSegmentedData(ascat.bc)
ascat.output = ascat.runAscat(ascat.bc, gamma=1, write_segments = T)
QC = ascat.metrics(ascat.bc,ascat.output)
save(ascat.bc, ascat.output, QC, file = 'ASCAT_objects.Rdata')
```

## Extracting logR and BAF from Tumor-only HTS (WGS) data and running ASCAT

```R
library(ASCAT)

ascat.prepareHTS(
  tumourseqfile = "Tumour.bam",
  tumourname = "Tumour_name",
  allelecounter_exe = "/PATH/TO/allelecounter",
  alleles.prefix = "G1000_alleles_hg38_chr",
  loci.prefix = "G1000_loci_hg38_chr",
  gender = "XX",
  genomeVersion = "hg38",
  nthreads = 8,
  tumourLogR_file = "Tumor_LogR.txt",
  tumourBAF_file = "Tumor_BAF.txt")

ascat.bc = ascat.loadData(Tumor_LogR_file = "Tumor_LogR.txt", Tumor_BAF_file = "Tumor_BAF.txt", gender = 'XX', genomeVersion = "hg38")
ascat.plotRawData(ascat.bc, img.prefix = "Before_correction_")
ascat.bc = ascat.correctLogR(ascat.bc, GCcontentfile = "GC_file.txt", replictimingfile = "RT_file.txt")
ascat.plotRawData(ascat.bc, img.prefix = "After_correction_")
gg = ascat.predictGermlineGenotypes(ascat.bc, platform = "HTS_WGS_1000G")
ascat.bc = ascat.aspcf(ascat.bc, ascat.gg=gg)
ascat.plotSegmentedData(ascat.bc)
ascat.output = ascat.runAscat(ascat.bc, gamma=1, write_segments = T)
QC = ascat.metrics(ascat.bc,ascat.output)
save(ascat.bc, ascat.output, QC, file = 'ASCAT_objects.Rdata')
```


## Processing targeted sequencing data
```R
library(ASCAT)

ascat.prepareTargetedSeq(
  Worksheet = "myWorksheet.tsv", # A tab-separated file with specific information. Check format using ?ascat.prepareTargetedSeq
  alleles.prefix = "G1000_alleles_hg19_chr",
  BED_file = "my_targeted_design.bed",
  allelecounter_exe = "/PATH/TO/allelecounter",
  genomeVersion = "hg19",
  nthreads = 8)

ascat.prepareHTS(
  tumourseqfile = "Tumour.bam",
  normalseqfile = "Normal.bam",
  tumourname = "Tumour_name",
  normalname = "Normal_name",
  allelecounter_exe = "/PATH/TO/allelecounter",
  alleles.prefix = "./alleleData/Cleaned/alleleData_chr",
  loci.prefix = "./alleleData/Cleaned/loci_chr",
  gender = "XX",
  genomeVersion = "hg19",
  nthreads = 8,
  tumourLogR_file = "Tumor_LogR.txt",
  tumourBAF_file = "Tumor_BAF.txt",
  normalLogR_file = "Germline_LogR.txt",
  normalBAF_file = "Germline_BAF.txt")
  
ascat.bc = ascat.loadData(Tumor_LogR_file = "Tumor_LogR.txt", Tumor_BAF_file = "Tumor_BAF.txt", Germline_LogR_file = "Germline_LogR.txt", Germline_BAF_file = "Germline_BAF.txt", gender = 'XX', genomeVersion = "hg19", isTargetedSeq=T)
ascat.plotRawData(ascat.bc, img.prefix = "Before_correction_")
ascat.bc = ascat.correctLogR(ascat.bc, GCcontentfile = "GC_file.txt", replictimingfile = "RT_file.txt")
ascat.plotRawData(ascat.bc, img.prefix = "After_correction_")
ascat.bc = ascat.aspcf(ascat.bc, penalty=25)
ascat.plotSegmentedData(ascat.bc)
ascat.output = ascat.runAscat(ascat.bc, gamma=1, write_segments = T)
QC = ascat.metrics(ascat.bc,ascat.output)
save(ascat.bc, ascat.output, QC, file = 'ASCAT_objects.Rdata')
```
